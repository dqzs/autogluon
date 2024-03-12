#!/usr/bin/env python
# coding: utf-8

# In[2]:


import autogluon as ag
import pandas as pd
import numpy as np
import os,urllib
from sklearn.metrics import r2_score
from sklearn.model_selection import train_test_split


# In[3]:


from autogluon.tabular import TabularPredictor as task


# In[4]:


# 加载数据集
data= pd.read_excel("DATA/probes_data.xlsx",index_col=0)
random_data = data.sample(frac =1).reset_index(drop=True)


# In[5]:


data


# In[6]:


data_new =random_data.iloc[:,1:]
data_new =data_new.astype(float)


# In[ ]:





# In[7]:


# 第二步：分割数据为训练集、验证集和测试集
train_data, temp_data = train_test_split(data_new, test_size=0.3, random_state=42) 
valid_data, test_data = train_test_split(temp_data, test_size=0.5, random_state=42)


# In[8]:


label_column = train_data.columns[0]
label_column


# In[ ]:





# In[ ]:


from autogluon.tabular import TabularPredictor

# 定义feature_prune_kwargs参数
feature_prune_kwargs = {
    
    'force_prune': True,
    'n_repeats': 3,
    'n_trials_per_model': 10,
    'prune_threshold': 0.001,
    'feature_importance_type': 'autogluon'
}

# 第三步：配置并启动AutoGluon训练
predictor = TabularPredictor(label=label_column,eval_metric="r2",
                             problem_type="regression"
                            ).fit(data_new,
                                  use_bag_holdout=True,
                                  #time_limit = 10000,
                                  #feature_prune_kwargs=feature_prune_kwargs,
                                  feature_prune_kwargs={'prune_threshold': 0.01}, 
                                  auto_stack=True, 
                                  hyperparameters='multimodal',
                                  presets='medium_quality',
                                 )


# In[ ]:


# 输出最优模型
predictor.get_model_best()


# In[ ]:


predictor.path


# In[9]:


# 查看特征重要性
feature_importance = predictor.feature_importance(data_new)


# In[ ]:





# In[18]:


#有的数据会出问题，如果出问题，用下面的数据
#1：
#predictor.evaluate(test_data.astype(float))['r2']
#2：
y_pred = predictor.predict(test_data)
r2_score(test_data[label_column], y_pred)


# In[19]:


predictor.evaluate_predictions(test_data[label_column].astype(float), y_pred.astype(float))


# In[32]:


# 根据特征重要性的排名重新训练数据，并绘制得分图
scores = []
for i in range(1, 21):
    # 选择前i个特征
    features = feature_importance.index[:i].tolist()
    #数据集进行链接
    df1 = pd.DataFrame(data_new.iloc[:,:1])
    df2=pd.DataFrame(data[features])
    data_1= pd.concat([df1, df2],axis=1)
    data_1
    size = int(0.8*data_1.shape[0])
    train_df_1 = data_1[:size]
    test_df_1 = data_1[size:]
    train_data_1 = task.Dataset(data=train_df_1)
    test_data_1 = task.Dataset(data=test_df_1)
    # 训练模型
    predictor = TabularPredictor(label=label_column,eval_metric="r2").fit(train_data_1,
                                                                          hyperparameters='multimodal',
                                                                          presets='good_quality',)
    
    # 计算得分
    score = predictor.evaluate(test_data_1.astype(float))["r2"]
    
    #y_pred = predictor.predict(test_data)
    #score = r2_score(test_data[label_column], y_pred)
    scores.append(score)
    
    # 显示进度
    print(f'Training with top {i} features, score = {score:.4f}')


# In[16]:


# 绘制得分图
import matplotlib.pyplot as plt
plt.plot(range(1, len(feature_importance) + 1), scores)
plt.xlabel('Number of features')
plt.ylabel('Score')
plt.show()


# In[30]:


from pylab import xticks,yticks,np
# 绘制得分图
import matplotlib.pyplot as plt
# 查看前20个特征的得分
plt.figure(figsize=(11,4))
plt.title("Features-R2")
#  选择的特征数量
plt.xlabel("Number of features selected")
# 交叉验证得分
plt.ylabel("R2_Score")

# 修改横轴坐标刻度
xticks(np.linspace(1,31,31,endpoint=True))
# yticks(np.linspace(0.7,1,7,endpoint=True))

# 画出各个特征的得分zzzzz
plt.plot(range(1, len(scores[:30])+1),scores[:30])
plt.grid()
plt.show()


# In[21]:


max_value = max(scores)
index = scores.index(max_value)
index


# ## 根据上面的 图，得知前个特征最好

# In[31]:


# 选择前i个特征
features = feature_importance.index[:4].tolist()
#数据集进行链接
df_select_1 = pd.DataFrame(data.iloc[:,:2])
df_select_2=pd.DataFrame(data[features])
feature_select= pd.concat([df_select_1, df_select_2],axis=1)


# In[32]:


feature_importance.index[:4].tolist()


# In[33]:


features


# In[34]:


feature_select


# In[35]:


feature_select.to_excel("DATA/data-feature-.xlsx")


# In[ ]:




