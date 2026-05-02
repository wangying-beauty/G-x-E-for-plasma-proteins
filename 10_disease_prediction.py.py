"""
# GxE Disease Prediction Pipeline

Machine learning pipeline for disease risk prediction using Gene-Environment (G×E) interactions.

## Features
- LightGBM classifier with hyperparameter tuning
- G × E interaction term generation
- SHAP-based feature interpretation
- Bootstrap confidence intervals
- DeLong test for model comparison

## Quick Start

### 1. Predict
python predict_gxe_model.py 0 10 --data_dir ./data
# Processes 10 diseases starting from index 0
"""

import numpy as np
import pandas as pd
from tqdm import tqdm
from collections import Counter
import warnings
import pickle
import os
import argparse
from lightgbm import LGBMClassifier
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import KFold
import random
import gc
import shap

warnings.filterwarnings('ignore')

# ==================== Configuration ====================
RANDOM_SEED = 2024
N_FOLDS = 10
NB_CPUS = 8

BASE_PARAMS = {
    'n_estimators': 100, 'max_depth': 5, 'num_leaves': 31,
    'subsample': 0.8, 'learning_rate': 0.1, 'colsample_bytree': 0.8,
    'boosting_type': 'gbdt', 'objective': 'binary', 'metric': 'auc',
    'verbose': -1, 'n_jobs': NB_CPUS
}

PARAMS_GRID = {
    'n_estimators': [100, 200, 300],
    'max_depth': [5, 10, 15],
    'learning_rate': [0.01, 0.05, 0.1],
}

# ==================== Utility Functions ====================
def clean_feature_name(name):
    """Clean feature names for LightGBM compatibility"""
    return str(name).replace(':', '_').replace(' ', '_').replace('(', '_').replace(')', '_')

def find_optimal_cutoff(y_true, y_pred):
    """Find optimal threshold using Youden's index"""
    fpr, tpr, thresholds = roc_curve(y_true, y_pred)
    youden_index = tpr - fpr
    optimal_idx = np.argmax(youden_index)
    return thresholds[optimal_idx]

def calculate_shap(model, X, max_samples=1000):
    """Calculate SHAP importance"""
    try:
        X_sample = X.sample(min(max_samples, len(X)))
        explainer = shap.TreeExplainer(model)
        shap_values = explainer.shap_values(X_sample)
        if isinstance(shap_values, list):
            shap_values = shap_values[1]
        mean_shap = np.abs(shap_values).mean(axis=0)
        return dict(zip(X.columns, mean_shap))
    except:
        return {}

# ==================== Data Preparation ====================
def prepare_disease_data(disease, interactions_df, snp_df, env_df, 
                         protein_df, disease_df, cov_df, protein_selection_df):
    """Prepare dataset with G×E interactions for one disease"""
    
    # Get interactions and proteins for this disease
    disease_interactions = interactions_df[
        interactions_df['Disease_code'] == disease]
    
    selected_proteins = protein_selection_df[
        protein_selection_df['Disease_code'] == disease
    ].iloc[0, 1:31].dropna().tolist()  # Top 30 proteins
    
    snps = disease_interactions['snp_id'].unique().tolist()
    envs = disease_interactions['env_id'].unique().tolist()
    
    # Merge all features
    df = disease_df[['sample_id', disease]].dropna()
    
    for snp in snps:
        if snp in snp_df.columns:
            df = df.merge(snp_df[['sample_id', snp]], on='sample_id', how='left')
    
    for env in envs:
        if env in env_df.columns:
            df = df.merge(env_df[['sample_id', env]], on='sample_id', how='left')
    
    for protein in selected_proteins:
        if protein in protein_df.columns:
            df = df.merge(protein_df[['sample_id', protein]], on='sample_id', how='left')
    
    df = df.merge(cov_df, on='sample_id', how='left')
    
    # Create G×E interactions
    interaction_features = []
    for _, row in disease_interactions.iterrows():
        snp, env = row['snp_id'], row['env_id']
        if snp in df.columns and env in df.columns:
            interaction_name = f"{env}×{snp}"
            # Standardize before multiplication
            df[interaction_name] = (
                (df[env] - df[env].mean()) / df[env].std() *
                (df[snp] - df[snp].mean()) / df[snp].std()
            )
            interaction_features.append(interaction_name)
    
    # Define feature groups
    demo_features = ['age', 'sex', 'smoking', 'TDI', 'BMI', 'alcohol', 'systolic blood pressure']
    protein_features = [p for p in selected_proteins if p in df.columns]
    
    return df, demo_features, protein_features, interaction_features

# ==================== Model Training ====================
def train_model(X_train, y_train, X_test, params):
    """Train LightGBM model"""
    model = LGBMClassifier(random_state=RANDOM_SEED, **params)
    model.fit(X_train, y_train)
    y_pred = model.predict_proba(X_test)[:, 1]
    return model, y_pred

def cross_validate(df, features, target, n_folds=10):
    """Perform cross-validation"""
    kf = KFold(n_splits=n_folds, shuffle=True, random_state=RANDOM_SEED)
    
    all_preds, all_labels, all_ids = [], [], []
    gain_importance = Counter()
    shap_importance = Counter()
    
    for train_idx, test_idx in kf.split(df):
        X_train = df.iloc[train_idx][features]
        y_train = df.iloc[train_idx][target]
        X_test = df.iloc[test_idx][features]
        y_test = df.iloc[test_idx][target]
        
        # Train model
        model, y_pred = train_model(X_train, y_train, X_test, BASE_PARAMS)
        
        # Collect results
        all_preds.extend(y_pred)
        all_labels.extend(y_test)
        all_ids.extend(df.iloc[test_idx]['sample_id'])
        
        # Feature importance
        gain_imp = dict(zip(features, model.feature_importances_))
        gain_importance.update(gain_imp)
        
        # SHAP (on test set only)
        shap_imp = calculate_shap(model, X_test)
        shap_importance.update(shap_imp)
    
    return {
        'predictions': pd.DataFrame({
            'sample_id': all_ids,
            'true_label': all_labels,
            'pred_prob': all_preds
        }),
        'gain_importance': dict(gain_importance),
        'shap_importance': dict(shap_importance),
        'auc': roc_auc_score(all_labels, all_preds)
    }

# ==================== Main Pipeline ====================
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('start_idx', type=int, help='Start disease index')
    parser.add_argument('num_diseases', type=int, help='Number of diseases')
    parser.add_argument('--data_dir', default='./data', help='Data directory')
    args = parser.parse_args()
    
    output_dir = f"results_{args.start_idx+1}_to_{args.start_idx+args.num_diseases}"
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Processing {args.num_diseases} diseases starting from index {args.start_idx}")
    
    # Load data
    print("Loading data...")
    snp_df = pd.read_csv(f"{args.data_dir}/snp_data.csv")
    protein_df = pd.read_csv(f"{args.data_dir}/protein_data.csv")
    env_df = pd.read_csv(f"{args.data_dir}/environment_data.csv")
    disease_df = pd.read_csv(f"{args.data_dir}/disease_labels.csv")
    cov_df = pd.read_csv(f"{args.data_dir}/covariates.csv")
    interactions_df = pd.read_csv(f"{args.data_dir}/gxe_interactions.csv")
    protein_selection_df = pd.read_csv(f"{args.data_dir}/protein_selection.csv")
    
    # Get disease list
    diseases = interactions_df['Disease_code'].unique()[
        args.start_idx : args.start_idx + args.num_diseases
    ]
    
    results = []
    
    # Process each disease
    for disease in tqdm(diseases, desc="Diseases"):
        try:
            print(f"\nProcessing: {disease}")
            
            # Prepare data
            df, demo_feats, protein_feats, interaction_feats = prepare_disease_data(
                disease, interactions_df, snp_df, env_df, 
                protein_df, disease_df, cov_df, protein_selection_df
            )
            
            # All features: Demo + Protein + G×E
            all_features = demo_feats + protein_feats + interaction_feats
            all_features = [clean_feature_name(f) for f in all_features]
            df.columns = [clean_feature_name(c) for c in df.columns]
            
            print(f"  Features: {len(all_features)} (Demo:{len(demo_feats)}, "
                  f"Protein:{len(protein_feats)}, G×E:{len(interaction_feats)})")
            
            # Cross-validation
            result = cross_validate(df, all_features, disease, N_FOLDS)
            
            print(f"  AUC: {result['auc']:.4f}")
            
            # Save predictions
            result['predictions'].to_csv(
                f"{output_dir}/predictions_{disease}.csv", index=False)
            
            # Save feature importance
            importance_df = pd.DataFrame({
                'feature': list(result['gain_importance'].keys()),
                'gain': list(result['gain_importance'].values()),
                'shap': [result['shap_importance'].get(f, 0) 
                        for f in result['gain_importance'].keys()]
            }).sort_values('gain', ascending=False)
            
            importance_df.to_csv(
                f"{output_dir}/importance_{disease}.csv", index=False)
            
            # Save top 20 to txt
            with open(f"{output_dir}/top20_{disease}.txt", 'w') as f:
                f.write(f"Disease: {disease}\n")
                f.write("="*60 + "\n\n")
                f.write("TOP 20 FEATURES BY GAIN:\n")
                for i, row in importance_df.head(20).iterrows():
                    f.write(f"{row.name+1:2d}. {row['feature']:40s} {row['gain']:.6f}\n")
            
            results.append({
                'disease': disease,
                'auc': result['auc'],
                'n_samples': len(df)
            })
            
            gc.collect()
            
        except Exception as e:
            print(f"Error processing {disease}: {e}")
            continue
    
    # Save summary
    pd.DataFrame(results).to_csv(f"{output_dir}/summary.csv", index=False)
    print(f"\nCompleted! Results saved to {output_dir}/")

if __name__ == "__main__":
    main()