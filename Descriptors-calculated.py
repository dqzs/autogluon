#!/usr/bin/env python
# coding: utf-8

# In[1]:


import rdkit
import pandas as pd
import numpy as np
from collections import Counter

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from mordred import Calculator, descriptors
from warnings import simplefilter
simplefilter(action='ignore', category=FutureWarning)


# In[4]:


experimental_data = pd.read_excel('Data/荧光分子数据总.xlsx')
experimental_data


# In[6]:


import os
folder_path = 'Data/sdf/'
files = []
for filename in os.listdir(folder_path):
    if os.path.isfile(os.path.join(folder_path, filename)):
        files.append(filename)
print(files)


# In[7]:


mols = []
fail = []
for i in files:
    suppl = Chem.SDMolSupplier(folder_path + str(i),strictParsing=False)
    if suppl[0] is not None:
        mols.append(suppl[0])
    else:
        fail.append(i)
mols
print("有%d个读取失败的,文件名称%s"%(len(fail),fail))


# In[8]:


print("有%d个读取失败的,文件名称%s"%(len(fail),fail))
for i in fail:
    files.remove(i)
print("删除后有%d个分子"%(len(files)))


# In[9]:


# 比较Rdkit与Mordred分子描述符 
calc = Calculator(descriptors, ignore_3D=True)
Mordred_description = []
Rdkit_description = [x[0] for x in Descriptors._descList]
for i in calc.descriptors:
    Mordred_description.append(i.__str__())
for i in Mordred_description:
    if i in Rdkit_description:
        Rdkit_description.remove(i)
        
Molecular_descriptor = []

descriptor_calculator = MoleculeDescriptors.MolecularDescriptorCalculator(Rdkit_description)
j =0
for i in mols:
    Calculator_descript = pd.DataFrame(calc.pandas([i]))
    rdkit_descriptors = pd.DataFrame([descriptor_calculator.CalcDescriptors(i)],columns=Rdkit_description)
    Calculator_descript = Calculator_descript.join(rdkit_descriptors)
    Molecular_descriptor.append(Calculator_descript)
    j+=1
    print(j)
    
a = Molecular_descriptor[0]
for i in Molecular_descriptor[1:]:
    a = a.append(i)
a = a.reset_index(drop=True)
# 删除计算失败的值
a = a.drop(labels=a.dtypes[a.dtypes == "object"].index,axis=1)
a


# In[10]:


# 将分子描述符与预测值结合
NAME = [x.replace("-3d","").replace(".sdf","") for x in files]
a.insert(0,"NAME",NAME)


# In[11]:


a.to_excel("a.xlsx")


# In[12]:


k = 0
for i,j in zip(Test_data.iloc[0,:],Test_data.iloc[1,:]):
    if i!=j:
        k+=1
        print(i,j)
k    


# In[12]:


experimental_data


# In[14]:


import math

# lef_series_index = [x.replace("\xa0","") for x in experimental_data["3dSD"].values]
y_series_index = []
for x in experimental_data["NAME"].values:
    if type(x) == str:
        y_series_index.append(x.replace("\xa0",""))
    else:
        y_series_index.append(x)
        
for i in range(len(y_series_index)):
    if type(y_series_index[i])!=str:
        y_series_index[i] = y_series_index[i-1]
y_series_value = [float(x) for x in experimental_data["EM"].values] 
y_series = pd.Series(y_series_value)
y_series.index = y_series_index
y_series = pd.DataFrame(y_series)
y_series


# In[18]:


final_data = pd.DataFrame()
for i in a["NAME"]:
    temp_data = pd.DataFrame()
    for j in range(y_series[y_series.index == i].shape[0]):
        temp_data = temp_data.append(a[a["NAME"] ==i])
    temp_data.insert(1,"SEF",y_series[y_series.index ==i].values)
    final_data = final_data.append(temp_data)
        
final_data   


# In[15]:


final_data = pd.merge(a, y_series, left_on="NAME", right_index=True)
final_data = final_data.reset_index(drop=True)
final_data = final_data.rename(columns={0: "EM"})
final_data


# In[16]:


final_data = pd.merge(a, y_series, left_on="NAME", right_index=True, how='left') 
final_data = final_data.reset_index(drop=True) 
final_data = final_data.rename(columns={0: "EM"}) 
final_data


# In[17]:


# 获取最后一列的列名
last_column = final_data.columns[-1]

# 将最后一列移到第二列位置
columns = list(final_data.columns)
columns.insert(1, last_column)
columns = columns[:-1]

# 重新排列列的顺序
final_data = final_data[columns]


# In[18]:


final_data


# In[19]:


final_data.to_excel("final_data_总.xlsx")

