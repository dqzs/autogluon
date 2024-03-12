#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
import math

from collections import Counter


# In[3]:


raw_data = pd.read_excel('DATA/flu_data.xlsx',index_col=0)
raw_data


# In[4]:


# 删除0值和空值
consolidated_data  = raw_data.dropna(axis=1,how='any')
data0=(consolidated_data==0).sum(axis = 0)
number0 = data0[data0>=raw_data.shape[0]/2]
number0.index
consolidated_data = consolidated_data.drop(labels=number0.index,axis=1)
print(consolidated_data.shape)
consolidated_data.head()


# In[5]:


# 删除低方差的值
var_data = consolidated_data.iloc[:,2:].var()
del_col = var_data[var_data<0.1].index
consolidated_data = consolidated_data.drop(labels=del_col,axis=1)
print(consolidated_data.shape)
consolidated_data.head()


# In[6]:


# 删除重复值多的列
Duplicated_series=pd.Series(np.zeros(len(consolidated_data.columns)))
Duplicated_series.index = consolidated_data.columns.values
for i in consolidated_data.columns:
    Duplicated_series[i] = len(Counter(consolidated_data[i]))
    
Duplicated_series = Duplicated_series.sort_values(ascending=True) 
# Duplicated_series60 = Duplicated_series[Duplicated_series>60]    
Duplicated_series15 = Duplicated_series[Duplicated_series<=raw_data.shape[0]/10]
consolidated_data = consolidated_data.drop(labels=Duplicated_series15.index,axis=1)
print(consolidated_data.shape)
consolidated_data.head()


# In[7]:


correlation = consolidated_data.iloc[:,2:].corr('spearman')
correlation


# In[8]:


df_bool = (abs(correlation) > 0.85)
correlation[df_bool]


# In[9]:


DN_correlation = consolidated_data.iloc[:,1:].corr('spearman')["EM"]
print(DN_correlation)


# In[10]:


bb = []
col_index = correlation.index
for i in range(0,len(col_index)):
    for j in range(i+1,len(col_index)):
        bb.append([col_index[i],col_index[j]])

for i in bb:
    if not math.isnan(correlation[df_bool].loc[i[0],i[1]]):
        print(i)
    


# In[18]:


k = 0
del_list = []
for i in bb:
    if not math.isnan(correlation[df_bool].loc[i[0],i[1]]) and ('EM'not in i):
        k+=1
        if abs(DN_correlation[i[0]])>abs(DN_correlation[i[1]]):
            if i[1] not in del_list:
                del_list.append(i[1])
        else:
            if i[0] not in del_list:
                del_list.append(i[0])
            
print(del_list)
k   


# In[19]:


Test_data = consolidated_data.drop(labels=del_list,axis=1)
Test_data


# In[21]:


Test_data.to_excel("DATA/预处理后-PCE-Fluorescent-probes_data.xlsx")


# In[ ]:




