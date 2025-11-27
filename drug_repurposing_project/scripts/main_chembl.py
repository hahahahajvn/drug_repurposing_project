import pandas as pd
import numpy as np
import os
import sys

# 添加当前目录到Python路径
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from data_collection_chembl import ChemBLDataCollector
from feature_engineering import FeatureEngineer
from model_training import DrugRepurposingModel

def run_chembl_pipeline():
    """运行基于ChEMBL的完整管道"""
    print("=== 基于ChEMBL的生物医学药物重定位项目 ===\n")
    
    # 第一步：从ChEMBL收集数据
    print("步骤 1/3: 从ChEMBL收集数据")
    collector = ChemBLDataCollector()
    compounds_df, activities_df, targets_df, associations_df = collector.collect_all_data()
    
    # 第二步：特征工程（需要适配ChEMBL数据结构）
    print("\n步骤 2/3: 特征工程")
    engineer = FeatureEngineer()
    
    # 注意：这里需要根据ChEMBL数据结构调整特征工程
    # 暂时使用简化版本
    features, labels, pairs_info = engineer.create_features(
        compounds_df.rename(columns={'chembl_id': 'drug_id', 'pref_name': 'name'}),
        targets_df.rename(columns={'target_chembl_id': 'disease_id', 'pref_name': 'name'}),
        associations_df.rename(columns={'drug_chembl_id': 'drug_id', 'target_chembl_id': 'disease_id'})
    )
    
    # 保存特征数据
    np.save('../data/processed/features.npy', features)
    np.save('../data/processed/labels.npy', labels)
    np.save('../data/processed/feature_names.npy', np.array(engineer.feature_names))
    
    # 第三步：模型训练
    print("\n步骤 3/3: 模型训练")
    model = DrugRepurposingModel()
    best_params = model.train(features, labels, engineer.feature_names)
    
    # 进行预测
    print("\n=== 进行预测 ===")
    prediction_features, prediction_pairs = engineer.prepare_prediction_features(
        compounds_df.rename(columns={'chembl_id': 'drug_id', 'pref_name': 'name'}),
        targets_df.rename(columns={'target_chembl_id': 'disease_id', 'pref_name': 'name'})
    )
    
    predictions = model.predict(prediction_features, prediction_pairs)
    predictions.to_csv('../results/chembl_predictions.csv', index=False)
    
    # 显示最有希望的预测结果
    print("\n=== 最有希望的药物重定位预测 ===")
    promising_predictions = predictions[
        (predictions['probability'] > 0.7) & 
        (predictions['probability'] < 0.95)
    ].head(10)
    
    print(promising_predictions[['drug_name', 'disease_name', 'probability']])
    
    # 保存模型
    model.save_model('../models/chembl_model.pkl')
    
    print("\n=== 项目完成！ ===")
    print("基于ChEMBL的结果已保存到:")
    print("- 预测结果: results/chembl_predictions.csv")
    print("- 模型文件: models/chembl_model.pkl")
    
    return model, predictions

if __name__ == "__main__":
    model, predictions = run_chembl_pipeline()