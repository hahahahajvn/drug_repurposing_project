import pandas as pd
import numpy as np
import requests
import time
import os
import json
from tqdm import tqdm
from config import CHEMBL_CONFIG, DATA_PATHS, SAMPLE_DRUGS, SAMPLE_DISEASES, API_ENDPOINTS, QUERY_PARAMS

class ChemBLDataCollector:
    def __init__(self):
        self.base_url = CHEMBL_CONFIG['base_url']
        self.session = requests.Session()
        self.setup_directories()
    
    def setup_directories(self):
        """创建必要的文件夹"""
        for path in DATA_PATHS.values():
            os.makedirs(path, exist_ok=True)
    
    def make_api_request(self, endpoint, params=None):
        """通用的API请求函数"""
        url = f"{self.base_url}{endpoint}"
        
        try:
            response = self.session.get(url, params=params, timeout=CHEMBL_CONFIG['timeout'])
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            print(f"API请求失败: {e}")
            return None
        finally:
            # 添加延迟避免过快请求
            time.sleep(CHEMBL_CONFIG['request_delay'])
    
    def search_compounds_by_name(self, compound_name):
        """根据化合物名称搜索"""
        params = {
            'molecule_synonyms__molecule_synonym__iexact': compound_name,
            'limit': QUERY_PARAMS['limit']
        }
        return self.make_api_request(API_ENDPOINTS['compounds'], params)
    
    def get_compound_by_chembl_id(self, chembl_id):
        """根据ChEMBL ID获取化合物详细信息"""
        endpoint = f"/compound/{chembl_id}.json"
        return self.make_api_request(endpoint)
    
    def get_activities_for_compound(self, chembl_id):
        """获取化合物的生物活性数据"""
        params = {
            'molecule_chembl_id': chembl_id,
            'limit': QUERY_PARAMS['limit'],
            'standard_type': QUERY_PARAMS['standard_type']
        }
        return self.make_api_request(API_ENDPOINTS['activities'], params)
    
    def search_targets_by_name(self, target_name):
        """根据靶点名称搜索"""
        params = {
            'pref_name__icontains': target_name,
            'limit': QUERY_PARAMS['limit']
        }
        return self.make_api_request(API_ENDPOINTS['targets'], params)
    
    def get_drug_mechanisms(self, chembl_id):
        """获取药物的作用机制"""
        params = {
            'molecule_chembl_id': chembl_id,
            'limit': QUERY_PARAMS['limit']
        }
        return self.make_api_request(API_ENDPOINTS['mechanisms'], params)
    
    def collect_sample_compounds(self):
        """收集样本化合物的数据"""
        print("=== 从ChEMBL收集化合物数据 ===")
        
        compounds_data = []
        
        for compound_name in tqdm(SAMPLE_DRUGS, desc="收集化合物数据"):
            result = self.search_compounds_by_name(compound_name)
            
            if result and 'compounds' in result:
                for compound in result['compounds']:
                    compound_info = {
                        'chembl_id': compound.get('molecule_chembl_id', ''),
                        'pref_name': compound.get('pref_name', ''),
                        'molecular_weight': compound.get('molecular_weight', None),
                        'alogp': compound.get('alogp', None),
                        'hbd': compound.get('hbd', None),  # 氢键供体
                        'hba': compound.get('hba', None),  # 氢键受体
                        'psa': compound.get('psa', None),  # 极性表面积
                        'rtb': compound.get('rtb', None),   # 可旋转键
                        'smiles': compound.get('molecule_structures', {}).get('canonical_smiles', ''),
                        'inchi': compound.get('molecule_structures', {}).get('standard_inchi', ''),
                        'drug_indication': compound.get('drug_indication', '')
                    }
                    compounds_data.append(compound_info)
            
            # 避免过快请求
            time.sleep(0.5)
        
        return pd.DataFrame(compounds_data)
    
    def collect_activity_data(self, compounds_df):
        """收集化合物的活性数据"""
        print("=== 收集生物活性数据 ===")
        
        activities_data = []
        
        for chembl_id in tqdm(compounds_df['chembl_id'].unique(), desc="收集活性数据"):
            result = self.get_activities_for_compound(chembl_id)
            
            if result and 'activities' in result:
                for activity in result['activities']:
                    activity_info = {
                        'chembl_id': chembl_id,
                        'assay_chembl_id': activity.get('assay_chembl_id', ''),
                        'target_chembl_id': activity.get('target_chembl_id', ''),
                        'target_pref_name': activity.get('target_pref_name', ''),
                        'standard_type': activity.get('standard_type', ''),
                        'standard_value': activity.get('standard_value', None),
                        'standard_units': activity.get('standard_units', ''),
                        'pchembl_value': activity.get('pchembl_value', None),
                        'relation': activity.get('relation', ''),
                        'assay_type': activity.get('assay_type', '')
                    }
                    activities_data.append(activity_info)
        
        return pd.DataFrame(activities_data)
    
    def collect_disease_related_targets(self):
        """收集与疾病相关的靶点"""
        print("=== 收集疾病相关靶点数据 ===")
        
        targets_data = []
        
        for disease_keyword in tqdm(SAMPLE_DISEASES, desc="收集疾病靶点"):
            result = self.search_targets_by_name(disease_keyword)
            
            if result and 'targets' in result:
                for target in result['targets']:
                    target_info = {
                        'target_chembl_id': target.get('target_chembl_id', ''),
                        'pref_name': target.get('pref_name', ''),
                        'target_type': target.get('target_type', ''),
                        'organism': target.get('organism', ''),
                        'disease_keyword': disease_keyword
                    }
                    targets_data.append(target_info)
        
        return pd.DataFrame(targets_data)
    
    def create_sample_associations(self, compounds_df, activities_df, targets_df):
        """创建样本药物-疾病关联数据"""
        print("=== 创建药物-疾病关联数据 ===")
        
        associations = []
        
        # 基于活性数据创建已知关联
        known_associations = activities_df.dropna(subset=['pchembl_value']).copy()
        known_associations = known_associations[known_associations['pchembl_value'] > 5]  # 选择高活性数据
        
        # 创建正样本（已知有效的关联）
        for _, activity in known_associations.iterrows():
            associations.append({
                'drug_chembl_id': activity['chembl_id'],
                'target_chembl_id': activity['target_chembl_id'],
                'disease_keyword': activity['target_pref_name'],
                'association': 1,
                'evidence': 'activity_data',
                'pchembl_value': activity['pchembl_value']
            })
        
        # 创建负样本（随机配对，假设无关联）
        all_drugs = compounds_df['chembl_id'].unique()
        all_targets = targets_df['target_chembl_id'].unique()
        
        # 随机选择一些药物-靶点对作为负样本
        np.random.seed(42)  # 保证可重复性
        num_negative = min(50, len(all_drugs) * len(all_targets))
        
        for _ in range(num_negative):
            drug = np.random.choice(all_drugs)
            target = np.random.choice(all_targets)
            
            # 确保这不是已知的正样本
            if not ((known_associations['chembl_id'] == drug) & 
                   (known_associations['target_chembl_id'] == target)).any():
                associations.append({
                    'drug_chembl_id': drug,
                    'target_chembl_id': target,
                    'disease_keyword': 'random_pair',
                    'association': 0,
                    'evidence': 'random_negative',
                    'pchembl_value': None
                })
        
        return pd.DataFrame(associations)
    
    def save_data(self, compounds_df, activities_df, targets_df, associations_df):
        """保存收集的数据"""
        compounds_df.to_csv(os.path.join(DATA_PATHS['raw_data'], 'chembl_compounds.csv'), index=False)
        activities_df.to_csv(os.path.join(DATA_PATHS['raw_data'], 'chembl_activities.csv'), index=False)
        targets_df.to_csv(os.path.join(DATA_PATHS['raw_data'], 'chembl_targets.csv'), index=False)
        associations_df.to_csv(os.path.join(DATA_PATHS['raw_data'], 'chembl_associations.csv'), index=False)
        
        print("数据保存完成！")
    
    def collect_all_data(self):
        """收集所有数据的主函数"""
        # 收集化合物数据
        compounds_df = self.collect_sample_compounds()
        print(f"收集到 {len(compounds_df)} 个化合物记录")
        
        # 收集活性数据
        activities_df = self.collect_activity_data(compounds_df)
        print(f"收集到 {len(activities_df)} 个活性记录")
        
        # 收集靶点数据
        targets_df = self.collect_disease_related_targets()
        print(f"收集到 {len(targets_df)} 个靶点记录")
        
        # 创建关联数据
        associations_df = self.create_sample_associations(compounds_df, activities_df, targets_df)
        print(f"创建了 {len(associations_df)} 个关联记录")
        
        # 保存数据
        self.save_data(compounds_df, activities_df, targets_df, associations_df)
        
        return compounds_df, activities_df, targets_df, associations_df

def main():
    """主函数"""
    collector = ChemBLDataCollector()
    
    # 收集所有数据
    compounds_df, activities_df, targets_df, associations_df = collector.collect_all_data()
    
    print("\n=== 数据汇总 ===")
    print(f"化合物数据: {compounds_df.shape}")
    print(f"活性数据: {activities_df.shape}")
    print(f"靶点数据: {targets_df.shape}")
    print(f"关联数据: {associations_df.shape}")
    
    # 显示数据样例
    print("\n化合物数据样例:")
    print(compounds_df.head())
    
    print("\n活性数据样例:")
    print(activities_df.head())
    
    return compounds_df, activities_df, targets_df, associations_df

if __name__ == "__main__":
    compounds_df, activities_df, targets_df, associations_df = main()