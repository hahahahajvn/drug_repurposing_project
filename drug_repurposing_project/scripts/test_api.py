# 快速测试脚本 test_api.py
import requests

def test_chembl_api():
    """测试ChEMBL API连接"""
    base_url = "https://www.ebi.ac.uk/chembl/api/data"
    
    # 测试不同的端点
    endpoints = [
        '/compound.json',
        '/compounds.json', 
        '/molecule.json',
        '/activity.json',
        '/activities.json',
        '/target.json',
        '/targets.json'
    ]
    
    for endpoint in endpoints:
        url = base_url + endpoint
        try:
            response = requests.get(url, params={'limit': 1}, timeout=10)
            print(f"{endpoint}: {response.status_code}")
            if response.status_code == 200:
                data = response.json()
                print(f"  成功! 找到 {len(data)} 条记录")
        except Exception as e:
            print(f"{endpoint}: 错误 - {e}")

if __name__ == "__main__":
    test_chembl_api()