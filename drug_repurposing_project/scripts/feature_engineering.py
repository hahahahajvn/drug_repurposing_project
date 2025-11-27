import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
import ast

class FeatureEngineer:
    def __init__(self):
        self.scaler = StandardScaler()
        self.feature_names = []
    
    def create_features(self, drugs_df, diseases_df, associations_df):
        """创建药物-疾病对的特征矩阵"""
        print("=== 第二步：特征工程 ===")
        
        features = []
        labels = []
        pairs_info = []
        
        # 只使用已知关联的数据进行训练
        train_associations = associations_df[associations_df['association'] != -1]
        
        for _, assoc in train_associations.iterrows():
            drug_id = assoc['drug_id']
            disease_id = assoc['disease_id']
            
            # 获取药物特征
            drug_row = drugs_df[drugs_df['drug_id'] == drug_id]
            if len(drug_row) == 0:
                continue
            drug_features = drug_row.iloc[0]
            
            # 获取疾病特征
            disease_row = diseases_df[diseases_df['disease_id'] == disease_id]
            if len(disease_row) == 0:
                continue
            disease_features = disease_row.iloc[0]
            
            # 创建特征向量
            feature_vector = self.create_feature_vector(drug_features, disease_features)
            features.append(feature_vector)
            labels.append(assoc['association'])
            pairs_info.append((drug_id, disease_id, drug_features['name'], disease_features['name']))
        
        # 创建特征名称
        self.create_feature_names()
        
        features_array = np.array(features)
        labels_array = np.array(labels)
        
        print(f"创建了 {len(features_array)} 个样本，每个样本 {features_array.shape[1]} 个特征")
        print(f"正样本: {sum(labels_array)}, 负样本: {sum(1 - labels_array)}")
        
        return features_array, labels_array, pairs_info
    
    def create_feature_vector(self, drug, disease):
        """为单个药物-疾病对创建特征向量"""
        features = []
        
        # 1. 基础药物特征
        features.extend([
            drug['molecular_weight'],
            drug['logP'],
            drug['h_bond_donors'],
            drug['h_bond_acceptors'],
            len(drug['targets'].split(',')) if pd.notna(drug['targets']) else 0
        ])
        
        # 2. 基础疾病特征
        features.extend([
            disease['gene_count'],
            disease['avg_gene_length'],
            len(disease['related_genes'].split(',')) if pd.notna(disease['related_genes']) else 0
        ])
        
        # 3. 交互特征（药物靶点与疾病基因的相似性）
        drug_targets = set(drug['targets'].split(',')) if pd.notna(drug['targets']) else set()
        disease_genes = set(disease['related_genes'].split(',')) if pd.notna(disease['related_genes']) else set()
        
        # 计算各种相似性指标
        intersection = drug_targets & disease_genes
        union = drug_targets | disease_genes
        
        # Jaccard相似性
        jaccard = len(intersection) / len(union) if len(union) > 0 else 0
        features.append(jaccard)
        
        # 重叠系数
        overlap = len(intersection) / min(len(drug_targets), len(disease_genes)) if min(len(drug_targets), len(disease_genes)) > 0 else 0
        features.append(overlap)
        
        # 简单计数特征
        features.append(len(intersection))
        features.append(len(drug_targets))
        features.append(len(disease_genes))
        
        return features
    
    def create_feature_names(self):
        """创建特征名称列表"""
        self.feature_names = [
            # 药物特征
            'drug_mol_weight', 'drug_logP', 'drug_h_donors', 'drug_h_acceptors', 'drug_target_count',
            # 疾病特征
            'disease_gene_count', 'disease_avg_gene_length', 'disease_related_gene_count',
            # 交互特征
            'jaccard_similarity', 'overlap_coefficient', 'shared_genes_count', 
            'drug_targets_count', 'disease_genes_count'
        ]
    
    def prepare_prediction_features(self, drugs_df, diseases_df):
        """为所有可能的药物-疾病对准备特征（用于预测）"""
        all_features = []
        all_pairs = []
        
        for _, drug in drugs_df.iterrows():
            for _, disease in diseases_df.iterrows():
                feature_vector = self.create_feature_vector(drug, disease)
                all_features.append(feature_vector)
                all_pairs.append((drug['drug_id'], disease['disease_id'], drug['name'], disease['name']))
        
        return np.array(all_features), all_pairs

def main():
    """特征工程测试"""
    # 加载数据
    drugs_df = pd.read_csv('../data/raw/sample_drugs.csv')
    diseases_df = pd.read_csv('../data/raw/sample_diseases.csv')
    associations_df = pd.read_csv('../data/raw/sample_associations.csv')
    
    # 创建特征
    engineer = FeatureEngineer()
    features, labels, pairs_info = engineer.create_features(drugs_df, diseases_df, associations_df)
    
    print(f"特征矩阵形状: {features.shape}")
    print(f"特征名称: {engineer.feature_names}")
    
    return features, labels, pairs_info, engineer

if __name__ == "__main__":
    features, labels, pairs_info, engineer = main()