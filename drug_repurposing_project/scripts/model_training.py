import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, cross_val_score, GridSearchCV
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score, roc_curve
from sklearn.preprocessing import StandardScaler
import joblib
import os

class DrugRepurposingModel:
    def __init__(self):
        self.model = None
        self.scaler = StandardScaler()
        self.feature_names = []
    
    def train(self, features, labels, feature_names):
        """训练模型"""
        print("=== 第三步：模型训练 ===")
        
        self.feature_names = feature_names
        
        # 数据标准化
        features_scaled = self.scaler.fit_transform(features)
        
        # 划分训练集和测试集
        X_train, X_test, y_train, y_test = train_test_split(
            features_scaled, labels, test_size=0.3, random_state=42, stratify=labels
        )
        
        print(f"训练集大小: {X_train.shape[0]}")
        print(f"测试集大小: {X_test.shape[0]}")
        
        # 定义参数网格
        param_grid = {
            'n_estimators': [50, 100],
            'max_depth': [5, 10, None],
            'min_samples_split': [2, 5],
            'min_samples_leaf': [1, 2]
        }
        
        # 使用网格搜索寻找最佳参数
        rf = RandomForestClassifier(random_state=42)
        grid_search = GridSearchCV(rf, param_grid, cv=5, scoring='roc_auc', n_jobs=-1, verbose=1)
        
        print("开始网格搜索...")
        grid_search.fit(X_train, y_train)
        
        self.model = grid_search.best_estimator_
        
        # 评估模型
        self.evaluate_model(X_test, y_test)
        
        # 可视化结果
        self.plot_results(X_test, y_test)
        
        return grid_search.best_params_
    
    def evaluate_model(self, X_test, y_test):
        """评估模型性能"""
        y_pred = self.model.predict(X_test)
        y_pred_proba = self.model.predict_proba(X_test)[:, 1]
        
        print("\n=== 模型性能评估 ===")
        print(f"最佳参数: {self.model.get_params()}")
        print(f"测试集AUC: {roc_auc_score(y_test, y_pred_proba):.3f}")
        
        # 交叉验证
        cv_scores = cross_val_score(self.model, X_test, y_test, cv=5, scoring='roc_auc')
        print(f"5折交叉验证AUC: {cv_scores.mean():.3f} (+/- {cv_scores.std() * 2:.3f})")
        
        print("\n分类报告:")
        print(classification_report(y_test, y_pred))
    
    def plot_results(self, X_test, y_test):
        """绘制结果图表"""
        y_pred_proba = self.model.predict_proba(X_test)[:, 1]
        
        # 创建图表
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # 1. ROC曲线
        fpr, tpr, _ = roc_curve(y_test, y_pred_proba)
        auc_score = roc_auc_score(y_test, y_pred_proba)
        
        axes[0, 0].plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC曲线 (AUC = {auc_score:.3f})')
        axes[0, 0].plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
        axes[0, 0].set_xlim([0.0, 1.0])
        axes[0, 0].set_ylim([0.0, 1.05])
        axes[0, 0].set_xlabel('假正率')
        axes[0, 0].set_ylabel('真正率')
        axes[0, 0].set_title('ROC曲线')
        axes[0, 0].legend(loc="lower right")
        axes[0, 0].grid(True)
        
        # 2. 特征重要性
        importance = self.model.feature_importances_
        indices = np.argsort(importance)[::-1]
        
        axes[0, 1].barh(range(len(importance)), importance[indices])
        axes[0, 1].set_yticks(range(len(importance)))
        axes[0, 1].set_yticklabels([self.feature_names[i] for i in indices])
        axes[0, 1].set_xlabel('重要性')
        axes[0, 1].set_title('特征重要性')
        
        # 3. 混淆矩阵
        y_pred = self.model.predict(X_test)
        cm = confusion_matrix(y_test, y_pred)
        sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=axes[1, 0])
        axes[1, 0].set_xlabel('预测标签')
        axes[1, 0].set_ylabel('真实标签')
        axes[1, 0].set_title('混淆矩阵')
        
        # 4. 概率分布
        axes[1, 1].hist([y_pred_proba[y_test == 0], y_pred_proba[y_test == 1]], 
                       bins=20, alpha=0.7, label=['负样本', '正样本'])
        axes[1, 1].set_xlabel('预测概率')
        axes[1, 1].set_ylabel('频数')
        axes[1, 1].set_title('预测概率分布')
        axes[1, 1].legend()
        
        plt.tight_layout()
        plt.savefig('../results/model_performance.png', dpi=300, bbox_inches='tight')
        plt.show()
    
    def predict(self, features, pairs_info):
        """进行预测"""
        if self.model is None:
            raise ValueError("请先训练模型！")
        
        features_scaled = self.scaler.transform(features)
        probabilities = self.model.predict_proba(features_scaled)[:, 1]
        
        # 创建结果DataFrame
        results = []
        for i, (drug_id, disease_id, drug_name, disease_name) in enumerate(pairs_info):
            results.append({
                'drug_id': drug_id,
                'drug_name': drug_name,
                'disease_id': disease_id,
                'disease_name': disease_name,
                'probability': probabilities[i],
                'prediction': 1 if probabilities[i] > 0.5 else 0
            })
        
        results_df = pd.DataFrame(results)
        results_df = results_df.sort_values('probability', ascending=False)
        
        return results_df
    
    def save_model(self, filepath):
        """保存模型"""
        if self.model is None:
            raise ValueError("没有模型可保存！")
        
        model_data = {
            'model': self.model,
            'scaler': self.scaler,
            'feature_names': self.feature_names
        }
        
        joblib.dump(model_data, filepath)
        print(f"模型已保存到: {filepath}")
    
    def load_model(self, filepath):
        """加载模型"""
        model_data = joblib.load(filepath)
        self.model = model_data['model']
        self.scaler = model_data['scaler']
        self.feature_names = model_data['feature_names']
        print(f"模型已从 {filepath} 加载")

def main():
    """模型训练主函数"""
    # 加载特征数据
    features = np.load('../data/processed/features.npy')
    labels = np.load('../data/processed/labels.npy')
    feature_names = np.load('../data/processed/feature_names.npy')
    
    # 训练模型
    model = DrugRepurposingModel()
    best_params = model.train(features, labels, feature_names)
    
    # 保存模型
    model.save_model('../models/drug_repurposing_model.pkl')
    
    return model, best_params

if __name__ == "__main__":
    model, best_params = main()