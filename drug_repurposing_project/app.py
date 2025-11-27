from flask import Flask, render_template, request, jsonify
import pandas as pd
import numpy as np
import joblib
import os
import sys

# 添加项目路径
sys.path.append('/home/yourusername/mysite/drug_repurposing_project/scripts')

app = Flask(__name__)

# 加载训练好的模型
try:
    model = joblib.load('models/drug_repurposing_model.pkl')
    print("模型加载成功")
except:
    print("模型加载失败，使用模拟模式")
    model = None

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/predict', methods=['POST'])
def predict():
    try:
        data = request.json
        drug_name = data.get('drug_name', '')
        disease_name = data.get('disease_name', '')
        
        # 这里添加预测逻辑
        if model:
            # 实际预测代码
            probability = 0.75  # 模拟数据
        else:
            probability = np.random.uniform(0.1, 0.9)
        
        return jsonify({
            'success': True,
            'drug_name': drug_name,
            'disease_name': disease_name,
            'probability': round(probability, 3),
            'prediction': '有效' if probability > 0.5 else '无效'
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})

@app.route('/api/compounds')
def get_compounds():
    """获取化合物列表API"""
    try:
        compounds_df = pd.read_csv('data/raw/chembl_compounds.csv')
        compounds = compounds_df[['chembl_id', 'pref_name']].to_dict('records')
        return jsonify(compounds)
    except:
        # 返回模拟数据
        return jsonify([
            {'chembl_id': 'CHEMBL1000', 'pref_name': 'METFORMIN'},
            {'chembl_id': 'CHEMBL2000', 'pref_name': 'ASPIRIN'}
        ])

if __name__ == '__main__':
    app.run(debug=True)