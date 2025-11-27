# 药物重定位预测平台

一个基于机器学习和ChEMBL数据库的药物重定位预测系统，通过分析药物-疾病关联数据来发现现有药物的新适应症。

## 项目概述

本项目整合了生物医学数据收集、特征工程、机器学习建模和Web应用开发，提供了一个完整的药物重定位预测解决方案。系统使用ChEMBL数据库的生物活性数据，通过随机森林模型预测药物与疾病之间的潜在关联。

## 主要功能

- 数据收集：从ChEMBL API自动获取化合物、靶点和活性数据
- 特征工程：提取药物和疾病的多种特征，构建交互特征
- 机器学习：使用随机森林算法训练预测模型
- Web界面：提供用户友好的预测界面和结果展示
- RESTful API：支持程序化访问预测功能

## 项目结构
```
drug_repurposing_project/
├── app.py # Flask Web应用主文件
├── index.html # 前端界面
├── requirements.txt # Python依赖包
├── config.py # 配置文件
├── data_collection_chembl.py # ChEMBL数据收集器
├── enhanced_chembl_collector.py # 增强版数据收集器
├── feature_engineering.py # 特征工程模块
├── model_training.py # 模型训练模块
├── main_chembl.py # 主执行管道
├── test_api.py # API测试脚本
├── debug.log # 运行日志
├── data/ # 数据目录
│ ├── raw/ # 原始数据
│ │ ├── chembl_compounds.csv
│ │ ├── chembl_activities.csv
│ │ ├── chembl_targets.csv
│ │ └── chembl_associations.csv
│ └── processed/ # 处理后的数据
├── models/ # 训练好的模型
├── results/ # 结果文件
└── logs/ # 日志文件
```
text

## 安装和配置

### 环境要求

- Python 3.8+
- 依赖包详见 `requirements.txt`

### 安装步骤

1. 克隆项目到本地：
```bash
git clone <repository-url>
cd drug_repurposing_project
安装依赖：

bash
pip install -r requirements.txt
创建必要的目录：

bash
mkdir -p data/raw data/processed models results logs
使用方法
数据收集
运行数据收集脚本从ChEMBL获取数据：

bash
python enhanced_chembl_collector.py
或者使用基础版本：

bash
python data_collection_chembl.py
模型训练
运行完整的训练管道：

bash
python main_chembl.py
或者分步执行：

bash
python feature_engineering.py
python model_training.py
启动Web应用
bash
python app.py
访问 http://localhost:5000 使用Web界面。

API接口
预测接口
POST /predict

请求体：

json
{
  "drug_name": "METFORMIN",
  "disease_name": "DIABETES"
}
响应：

json
{
  "success": true,
  "drug_name": "METFORMIN",
  "disease_name": "DIABETES", 
  "probability": 0.875,
  "prediction": "有效"
}
化合物列表接口
GET /api/compounds

响应：

json
[
  {
    "chembl_id": "CHEMBL1000",
    "pref_name": "METFORMIN"
  },
  {
    "chembl_id": "CHEMBL2000",
    "pref_name": "ASPIRIN"
  }
]
核心模块说明
数据收集模块 (enhanced_chembl_collector.py)
支持从ChEMBL API获取化合物、靶点和活性数据

包含完整的错误处理和重试机制

提供备用数据生成功能

特征工程模块 (feature_engineering.py)
提取的特征包括：

药物特征：分子量、logP、氢键供体/受体数等

疾病特征：相关基因数量、平均基因长度等

交互特征：Jaccard相似性、重叠系数等

模型训练模块 (model_training.py)
使用随机森林分类器

支持网格搜索超参数调优

提供完整的模型评估和可视化

Web应用模块 (app.py)
基于Flask的RESTful API

响应式前端界面

实时预测功能

配置文件
项目的主要配置在 config.py 中，包括：

ChEMBL API端点配置

数据路径设置

样本药物和疾病列表

查询参数配置

数据文件说明
chembl_compounds.csv：化合物基本信息

chembl_activities.csv：生物活性数据

chembl_targets.csv：疾病相关靶点信息

chembl_associations.csv：药物-疾病关联数据

故障排除
如果遇到API连接问题，可以运行测试脚本：

bash
python test_api.py
查看详细日志信息：

bash
tail -f logs/debug.log

本项目仅供学习和研究使用。

text
