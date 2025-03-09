"""
常用光学材料折射率数据库
包含各种常见材料在可见光范围内的折射率数据，支持波长相关的折射率
"""

import numpy as np
from scipy.interpolate import interp1d
import json
import os

# 基础材料数据（这些材料始终存在，不会被覆盖）
BASE_MATERIALS = {
    'Air': {
        'description': '空气',
        'data': 'constant',
        'values': 1.0
    }
}

# 初始化材料数据库
MATERIALS = BASE_MATERIALS.copy()

def load_materials_from_file(filename="custom_materials.json"):
    """从文件加载材料数据库"""
    global MATERIALS
    try:
        if os.path.exists(filename):
            with open(filename, 'r', encoding='utf-8') as f:
                materials_data = json.load(f)
                # 保持基础材料不变
                MATERIALS = BASE_MATERIALS.copy()
                # 添加自定义材料
                MATERIALS.update(materials_data)
    except Exception as e:
        print(f"加载材料数据库时出错: {e}")

def save_materials_to_file(filename="custom_materials.json"):
    """保存材料数据库到文件（不保存基础材料）"""
    try:
        # 将材料数据转换为可序列化的格式，排除基础材料
        materials_data = {}
        for name, material in MATERIALS.items():
            if name not in BASE_MATERIALS:  # 不保存基础材料
                if material['data'] == 'table':
                    materials_data[name] = {
                        'description': material['description'],
                        'data': material['data'],
                        'values': {
                            'wavelength': material['values']['wavelength'],
                            'n': material['values']['n'],
                            'k': material['values']['k']
                        }
                    }
        
        # 保存到文件
        with open(filename, 'w', encoding='utf-8') as f:
            json.dump(materials_data, f, ensure_ascii=False, indent=4)
            
    except Exception as e:
        print(f"保存材料数据库时出错: {e}")

# 其他函数保持不变
def get_material_names():
    """返回所有材料名称列表"""
    return list(MATERIALS.keys())

def get_material_refractive_index(material_name, wavelength=550):
    """
    根据材料名称和波长返回折射率
    
    参数:
    material_name: 材料名称
    wavelength: 波长(nm)，默认550nm
    
    返回:
    复折射率 (n + kj)
    """
    if material_name in MATERIALS:
        material = MATERIALS[material_name]
        
        if material['data'] == 'constant':
            if material['values'] is None:
                return None
            return material['values']
        
        elif material['data'] == 'table':
            values = material['values']
            if values is None:
                return None
            wavelengths = values['wavelength']
            n_values = values['n']
            k_values = values['k']
            
            try:
                # 使用插值计算指定波长的折射率
                n_interp = interp1d(wavelengths, n_values, kind='cubic', bounds_error=False, fill_value='extrapolate')
                k_interp = interp1d(wavelengths, k_values, kind='cubic', bounds_error=False, fill_value='extrapolate')
                
                n = float(n_interp(wavelength))
                k = float(k_interp(wavelength))
                
                # 对n和k进行四舍五入，保留4位小数，避免出现过长的小数
                n = round(n, 4)
                k = round(k, 4)
                
                return complex(n, k)
            except Exception as e:
                print(f"插值计算折射率时出错: {e}")
                return None
        
        elif material['data'] == 'formula':
            # 未来可以添加使用公式计算折射率的功能
            pass
        
    return None

def get_material_n_k_data(material_name):
    """
    获取材料的完整n,k数据表
    
    参数:
    material_name: 材料名称
    
    返回:
    字典，包含波长、n值和k值数组
    """
    try:
        if material_name in MATERIALS:
            material = MATERIALS[material_name]
            
            if material['data'] == 'constant':
                # 对于常数折射率，创建一个简单的数据表
                value = material['values']
                if value is None:
                    return None
                    
                if isinstance(value, complex):
                    n, k = value.real, value.imag
                else:
                    n, k = value, 0
                    
                return {
                    'wavelength': [400, 800],
                    'n': [n, n],
                    'k': [k, k]
                }
                
            elif material['data'] == 'table':
                if material['values'] is None:
                    return None
                return material['values']
                
        return None
    except Exception as e:
        print(f"获取材料n,k数据时出错: {e}")
        return None

def get_material_description(material_name):
    """根据材料名称返回描述"""
    try:
        if material_name in MATERIALS:
            return MATERIALS[material_name]['description']
        return ""
    except Exception as e:
        print(f"获取材料描述时出错: {e}")
        return ""

def get_material_info(material_name):
    """根据材料名称返回完整信息"""
    try:
        if material_name in MATERIALS:
            return MATERIALS[material_name]
        return None
    except Exception as e:
        print(f"获取材料信息时出错: {e}")
        return None

def get_refractive_index_at_wavelength(material_name, wavelength):
    """
    获取指定材料在指定波长的折射率
    
    参数:
    material_name: 材料名称
    wavelength: 波长(nm)
    
    返回:
    复折射率 (n + kj)
    """
    try:
        result = get_material_refractive_index(material_name, wavelength)
        if result is not None and isinstance(result, complex):
            # 确保复数的实部和虚部都被四舍五入
            real_part = round(result.real, 4)
            imag_part = round(result.imag, 4)
            return complex(real_part, imag_part)
        return result
    except Exception as e:
        print(f"获取指定波长的折射率时出错: {e}")
        return None

# 在模块初始化时加载材料数据
load_materials_from_file() 