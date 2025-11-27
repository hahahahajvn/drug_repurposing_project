# 配置文件 - 专门针对ChEMBL数据库
CHEMBL_CONFIG = {
    'base_url': 'https://www.ebi.ac.uk/chembl/api/data',
    'timeout': 30,
    'max_retries': 3,
    'request_delay': 0.5  # 请求间隔(秒)，避免过快请求
}

# 数据路径配置
DATA_PATHS = {
    'raw_data': '../data/raw/',
    'processed_data': '../data/processed/',
    'models': '../models/',
    'results': '../results/'
}

# 要查询的药物和疾病列表（可以根据需要扩展）
SAMPLE_DRUGS = [
    'METFORMIN', 'ASPIRIN', 'SIMVASTATIN', 'PROPRANOLOL', 'SILDENAFIL',
    'ATORVASTATIN', 'LISINOPRIL', 'METOPROLOL', 'WARFARIN', 'IBUPROFEN'
]

SAMPLE_DISEASES = [
    'DIABETES', 'ALZHEIMER', 'CANCER', 'HYPERTENSION', 'ASTHMA',
    'DEPRESSION', 'ARTHRITIS', 'MIGRAINE', 'EPILEPSY', 'PARKINSON'
]

# ChEMBL API端点配置
API_ENDPOINTS = {
    'compounds': '/compounds.json',
    'activities': '/activities.json',
    'assays': '/assays.json',
    'targets': '/targets.json',
    'mechanisms': '/mechanisms.json',
    'drugs': '/drugs.json',
    'molecules': '/molecules.json'
}

# 查询参数配置
QUERY_PARAMS = {
    'limit': 100,  # 每次查询返回的最大记录数
    'offset': 0,   # 查询偏移量
    'standard_type': 'IC50'  # 标准活动类型
}

# 特征工程配置
FEATURE_CONFIG = {
    # 药物特征
    'drug_features': [
        'molecular_weight', 'alogp', 'hbd', 'hba', 'psa', 'rtb'
    ],
    # 活动数据特征
    'activity_features': [
        'standard_value', 'standard_type', 'standard_units', 'pchembl_value'
    ]
}