import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import numpy as np
import scanpy as sc
import pandas as pd
from sklearn.model_selection import KFold, StratifiedKFold
from sklearn.preprocessing import MinMaxScaler, RobustScaler
from sklearn.metrics import *
from imblearn.over_sampling import SMOTE, RandomOverSampler
import argparse
import scipy
import warnings
import random
import os
import copy


warnings.filterwarnings('ignore')

seed=42
torch.manual_seed(seed)
np.random.seed(seed)
torch.cuda.manual_seed(seed)
torch.cuda.manual_seed_all(seed)
random.seed(seed)
np.random.seed(seed)
os.environ['PYTHONHASHSEED'] = str(seed)


class SimpleClassifier(nn.Module):
    def __init__(self, input_size):
        super(SimpleClassifier, self).__init__()
        
        self.classifier = nn.Sequential(
            nn.Linear(input_size, args.hdim),
            nn.ReLU(),
            nn.Dropout(args.dropout),
            nn.Linear(args.hdim, 1),
            nn.Sigmoid()
        )

    def forward(self, x):
        return self.classifier(x)

        
def preprocess_bulk_data(bulk_data):
    return np.log1p(bulk_data)

def preprocess_single_cell_data(adata):
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata.X


def train_model(model, train_loader, criterion, optimizer, device):
    model.train()
    running_loss = 0.0
    all_preds = []
    all_labels = []
    
    for inputs, labels in train_loader:
        inputs = inputs.to(device)
        labels = labels.to(device)
        
        optimizer.zero_grad()
        outputs = model(inputs)
        loss = criterion(outputs, labels)
        loss.backward()
        optimizer.step()
        
        running_loss += loss.item()
        
        # Convert predictions to binary (0 or 1)
        preds = (outputs >= 0.5).float()
        all_preds.extend(preds.cpu().numpy())
        all_labels.extend(labels.cpu().numpy())
    
    # Calculate metrics
    return running_loss / len(train_loader)


def test_model(model, test_data, test_labels, device):
    model.eval()
    with torch.no_grad():
        test_tensor = torch.FloatTensor(test_data).to(device)
        outputs = model(test_tensor)
        all_outputs = outputs.cpu().numpy().flatten()
        
        auc = roc_auc_score(test_labels, all_outputs)
        ap = average_precision_score(test_labels, all_outputs)
        fpr, tpr, thresholds = roc_curve(test_labels, all_outputs)
        
        optimal_idx = np.argmax(tpr - fpr)
        best_threshold = thresholds[optimal_idx]
        
        binary_predictions = (all_outputs >= best_threshold).astype(int)
        F1 = f1_score(test_labels, binary_predictions)
        
        # report_dict = classification_report(test_labels, binary_predictions, output_dict=True)
        # report_df = pd.DataFrame(report_dict).T
        metrics = {
            'auc': auc,
            'ap':ap,
            'best_threshold': best_threshold,
            'F1': F1,
            'predictions': binary_predictions,
            'raw_outputs': all_outputs
        }
        
        return metrics


def get_scdeal(current_type, hvg=''):
    
    if current_type == 'scStyleBulk':
        source_data_path = f'/NAS_Storage3/ddochi132/data/scDEAL/LOBICO/HVG/scStyleBulk_{args.data_name}_GDSC_HVG.csv'
        bulk_labels = sc.read_h5ad(f'/NAS_Storage3/ddochi132/data/datasets/GDSC/Multiple/Bulk_comgenes/Bulk_{args.data_name}.h5ad').obs['sensitive']
        
    elif current_type == 'Bulk':
        source_data_path = f'/NAS_Storage3/ddochi132/data/scDEAL/LOBICO/HVG/Bulk_{args.data_name}_GDSC_HVG.csv'
        bulk_labels = sc.read_h5ad(f'/NAS_Storage3/ddochi132/data/datasets/GDSC/Multiple/Bulk_comgenes/Bulk_{args.data_name}.h5ad').obs['sensitive']
        
    bulk_data = pd.read_csv(source_data_path, index_col=0)
    
    single_cell_data = sc.read_h5ad(f'/NAS_Storage3/ddochi132/data/scDEAL/LOBICO/HVG/raw_{args.data_name}_QC_HVG.h5ad')
    single_cell_labels = sc.read_h5ad(f'/NAS_Storage3/ddochi132/data/datasets/GDSC/Multiple/singlecell/{args.data_name}_maxgenes_QC.h5ad').obs['sensitive']
    
    return bulk_data, bulk_labels, single_cell_data, single_cell_labels


def get_scad(current_type, hvg=''):
    
    if current_type == 'scStyleBulk':
        source_data_path = f'/NAS_Storage3/ddochi132/data/SCAD/scStyleBulk/target_raw/HVG/scStyleBulk_{args.data_name}_GDSC_HVG.csv'
        bulk_labels = sc.read_h5ad(f'/NAS_Storage3/ddochi132/data/datasets/GDSC/Multiple/Bulk_comgenes/Bulk_{args.data_name}.h5ad').obs['sensitive']
        
    elif current_type == 'Bulk':
        source_data_path = f'/NAS_Storage3/ddochi132/data/SCAD/scStyleBulk/source_GDSC/HVG/{args.data_name}_GDSC_HVG.csv'
        bulk_labels = sc.read_h5ad(f'/NAS_Storage3/ddochi132/data/datasets/GDSC/Multiple/Bulk_comgenes/Bulk_{args.data_name}.h5ad').obs['sensitive']
        
    bulk_data = pd.read_csv(source_data_path, index_col=0)
    
    single_cell_data = sc.read_h5ad(f'/NAS_Storage3/ddochi132/data/SCAD/scStyleBulk/target_raw/HVG/{args.data_name}_HVG.h5ad')
    single_cell_labels = sc.read_h5ad(f'/NAS_Storage3/ddochi132/data/datasets/GDSC/Multiple/singlecell/{args.data_name}_maxgenes_QC.h5ad').obs['sensitive']
    
    return bulk_data, bulk_labels, single_cell_data, single_cell_labels

   
def get_scRobust_emb(current_type, hvg=''):
    
    if current_type == 'scStyleBulk':
        source_data_path = f'/NAS_Storage3/ddochi132/data/scRobust_WH/embedding_HVG/GDSC_Multiple/scStyleBulk_{args.data_name}_HVG2000.npy'
        bulk_labels = sc.read_h5ad(f'/NAS_Storage3/ddochi132/data/datasets/GDSC/Multiple/Bulk_comgenes/Bulk_{args.data_name}.h5ad').obs['sensitive']
        
    elif current_type == 'Bulk':
        source_data_path = f'/NAS_Storage3/ddochi132/data/scRobust_WH/embedding_HVG/GDSC_Multiple/Bulk_{args.data_name}_HVG2000.npy'
        bulk_labels = sc.read_h5ad(f'/NAS_Storage3/ddochi132/data/datasets/GDSC/Multiple/Bulk_comgenes/Bulk_{args.data_name}.h5ad').obs['sensitive']
        
    bulk_data = np.load(source_data_path)
    
    single_cell_data = np.load(f'/NAS_Storage3/ddochi132/data/scRobust_WH/embedding_HVG/GDSC_Multiple/sc_{args.data_name}_HVG2000.npy')
    single_cell_labels = sc.read_h5ad(f'/NAS_Storage3/ddochi132/data/datasets/GDSC/Multiple/singlecell/{args.data_name}_maxgenes_QC.h5ad').obs['sensitive']
    
    return bulk_data, bulk_labels, single_cell_data, single_cell_labels    


def get_scGPT_emb(current_type, hvg=''):
    
    if current_type == 'scStyleBulk':
        source_data_path = f'/NAS_Storage3/ddochi132/data/scGPT/embedding/GDSC_Multiple/scStyleBulk_{args.data_name}_HVG.npy'
        bulk_labels = sc.read_h5ad(f'/NAS_Storage3/ddochi132/data/datasets/GDSC/Multiple/Bulk_comgenes/Bulk_{args.data_name}.h5ad').obs['sensitive']
        
    elif current_type == 'Bulk':
        source_data_path = f'/NAS_Storage3/ddochi132/data/scGPT/embedding/GDSC_Multiple/Bulk_{args.data_name}_HVG.npy'
        bulk_labels = sc.read_h5ad(f'/NAS_Storage3/ddochi132/data/datasets/GDSC/Multiple/Bulk_comgenes/Bulk_{args.data_name}.h5ad').obs['sensitive']
        
    bulk_data = np.load(source_data_path)
    
    single_cell_data = np.load(f'/NAS_Storage3/ddochi132/data/scGPT/embedding/GDSC_Multiple/raw_{args.data_name}_HVG.npy')
    single_cell_labels = sc.read_h5ad(f'/NAS_Storage3/ddochi132/data/datasets/GDSC/Multiple/singlecell/{args.data_name}_maxgenes_QC.h5ad').obs['sensitive']
    
    return bulk_data, bulk_labels, single_cell_data, single_cell_labels  


def get_scFoundation_emb(current_type, hvg=''):
    
    if current_type == 'scStyleBulk':
        source_data_path = f'/NAS_Storage3/ddochi132/data/scFoundation/embedding_comgenes/GDSC_Multiple/scStyleBulk_{args.data_name}_50M-0.1B-res_singlecell_cell_embedding_t4_resolution.npy'
        bulk_labels = sc.read_h5ad(f'/NAS_Storage3/ddochi132/data/datasets/GDSC/Multiple/Bulk_comgenes/Bulk_{args.data_name}.h5ad').obs['sensitive']
        
    elif current_type == 'Bulk':
        source_data_path = f'/NAS_Storage3/ddochi132/data/scFoundation/embedding_comgenes/GDSC_Multiple/Bulk_{args.data_name}_50M-0.1B-res_bulk_cell_embedding_f1_resolution.npy'
        bulk_labels = sc.read_h5ad(f'/NAS_Storage3/ddochi132/data/datasets/GDSC/Multiple/Bulk_comgenes/Bulk_{args.data_name}.h5ad').obs['sensitive']
        
    bulk_data = np.load(source_data_path)
    
    single_cell_data = np.load(f'/NAS_Storage3/ddochi132/data/scFoundation/embedding_comgenes/GDSC_Multiple/sc_{args.data_name}_50M-0.1B-res_singlecell_cell_embedding_t4_resolution.npy')
    single_cell_labels = sc.read_h5ad(f'/NAS_Storage3/ddochi132/data/datasets/GDSC/Multiple/singlecell/{args.data_name}_maxgenes_QC.h5ad').obs['sensitive']
    
    return bulk_data, bulk_labels, single_cell_data, single_cell_labels  


def main(args):
    device = torch.device("cuda:11" if torch.cuda.is_available() else "cpu")
    current_type = args.type
    combination = f'{current_type}_{args.data_name}_w{args.wait}_b{args.batch}_h{args.hdim}_{args.lr}_d{args.dropout}_{args.bn}_{args.sampling}'
    
    
    if args.foundation_model == 'gexpr':
        bulk_data, bulk_labels, single_cell_data, single_cell_labels = get_scdeal(current_type)
        
        comgenes = bulk_data.columns.intersection(single_cell_data.var_names)
        bulk_data = bulk_data[comgenes] 
        single_cell_data = single_cell_data[:, comgenes] 

        bulk_data_processed = preprocess_bulk_data(bulk_data)
        sc_data_processed_dense = preprocess_single_cell_data(single_cell_data)
    
    elif args.foundation_model == 'scRobust':
        bulk_data_processed, bulk_labels, sc_data_processed_dense, single_cell_labels = get_scRobust_emb(current_type)
    elif args.foundation_model == 'scGPT':    
        bulk_data_processed, bulk_labels, sc_data_processed_dense, single_cell_labels = get_scGPT_emb(current_type)
    elif args.foundation_model == 'scFoundation':  
        bulk_data_processed, bulk_labels, sc_data_processed_dense, single_cell_labels = get_scFoundation_emb(current_type)
        
        
    # Initialize 5-fold cross validation
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    fold_results_f1 = []
    fold_results_auc = []
    fold_results_ap = []
    best_val_f1_lst = []

    
    for fold, (train_idx, val_idx) in enumerate(skf.split(bulk_data_processed, bulk_labels)):
        
        # Split data and labels for this fold
        if args.foundation_model=='gexpr':
            X_train = bulk_data_processed.iloc[train_idx]
            y_train = bulk_labels[train_idx]
            X_val = bulk_data_processed.iloc[val_idx]
            y_val = bulk_labels[val_idx]
            
        else:
            X_train = bulk_data_processed[train_idx]
            y_train = bulk_labels[train_idx]
            X_val = bulk_data_processed[val_idx]
            y_val = bulk_labels[val_idx]

        if args.foundation_model == 'gexpr':
            # Standardize training data
            scaler = MinMaxScaler()
            X_train = scaler.fit_transform(X_train)
            X_val = scaler.transform(X_val)
            sc_data_scaled = scaler.transform(sc_data_processed_dense)
        else:
            sc_data_scaled = sc_data_processed_dense
            
            
        if (len(np.unique(y_train)) >1) and (len(np.unique(y_val)) > 1) and args.sampling == 'smote':
            smote = SMOTE(random_state=42)
            X_train_resampled, y_train_resampled = smote.fit_resample(X_train, y_train)
        elif (len(np.unique(y_train)) >1) and (len(np.unique(y_val)) > 1) and args.sampling == 'up':
            smote = RandomOverSampler(random_state=42)
            X_train_resampled, y_train_resampled = smote.fit_resample(X_train, y_train)
        else:
            X_train_resampled, y_train_resampled = X_train, y_train
        
        # Create data loaders
        train_dataset = TensorDataset(
            torch.FloatTensor(X_train_resampled),
            torch.FloatTensor(y_train_resampled).reshape(-1, 1)
        )
        train_loader = DataLoader(train_dataset, batch_size=args.batch, shuffle=True)
        
        # Initialize model, criterion, and optimizer
        model = SimpleClassifier(input_size=X_train.shape[1]).to(device)
        criterion = nn.BCELoss()
        optimizer = optim.Adam(model.parameters(), lr=args.lr)
        
        best_val_f1 = 0
        best_model = None
        wait = 0
        
        # Training loop
        n_epochs = 100
        for epoch in range(n_epochs):
            # Train
            train_loss = train_model(model, train_loader, criterion, optimizer, device)
            
            # Validation
            val_results = test_model(model, X_val, y_val, device)
            # val_metrics = calculate_metrics(y_val, val_preds)
            
            if val_results['F1'] >= best_val_f1:
                wait = 0
                best_val_f1 = val_results['F1']
                best_model = copy.deepcopy(model.state_dict())
            
            else:
                wait += 1
                if wait >= args.wait:
                    # print("Early Stopping: ", epoch)
                    break
        
        best_val_f1_lst.append(best_val_f1)
        # Load best model for this fold
        
        model = SimpleClassifier(input_size=X_train.shape[1]).to(device)
        model.load_state_dict(best_model)
        
        test_results = test_model(model, sc_data_scaled, single_cell_labels, device)
        
        fold_results_f1.append(test_results['F1'])
        fold_results_auc.append(test_results['auc'])
        fold_results_ap.append(test_results['ap'])
        
    # Average predictions and metrics across folds
    final_f1 = round(np.mean(fold_results_f1, axis=0), 4)
    final_auc = round(np.mean(fold_results_auc, axis=0), 4)
    final_ap = round(np.mean(fold_results_ap, axis=0), 4)
    
    # print('Val F1 : ', round(np.mean(best_val_f1_lst, axis=0),4))
    # print('Val F1 folds: ', [round(f1, 4) for f1 in best_val_f1_lst])
    
    print(combination, ' VAL F1 :', round(np.mean(best_val_f1_lst, axis=0),4), '| TEST F1: ', final_f1, '| TEST AUC: ', final_auc, '| TEST AP: ',final_ap)
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Train Autoencoders and Classifier')

    parser.add_argument('--data_name', type=str, default='')
    parser.add_argument('--type', type=str, default='scStyleBulk')
    parser.add_argument('--drug', type=str, default='GEFITINIB',help='Target drug')
    parser.add_argument('--lr', type=float, default=1e-5)
    parser.add_argument('--reduce', type=str, default='F')
    parser.add_argument('--n_neighbors', type=int, default=20)
    parser.add_argument('--min_dist', type=float, default=0.2)
    
    
    parser.add_argument('--hdim', type=int, default=128)
    parser.add_argument('--wait', type=int, default=5)
    parser.add_argument('--batch', type=int, default=128)
    parser.add_argument('--dropout', type=float, default=0.5)
    parser.add_argument('--bn', type=str, default='False')
    parser.add_argument('--sampling', type=str, default='up')
    parser.add_argument('--foundation_model', type=str, default='gexpr')
    
    # Argument parsing
    args = parser.parse_args()
    main(args)