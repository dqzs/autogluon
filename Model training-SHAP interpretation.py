#!/usr/bin/env python
# coding: utf-8

# In[9]:


import autogluon as ag
import pandas as pd
import numpy as np
import os,urllib
import matplotlib.pyplot as plt
from autogluon.tabular import TabularDataset, TabularPredictor
from sklearn.metrics import r2_score
from sklearn.model_selection import train_test_split


# In[124]:


# 加载数据集
data= pd.read_excel("Cp-data-10.xlsx",index_col=0,)
random_data = data.sample(frac =1).reset_index(drop=True)


# In[125]:


random_data_de_name = random_data.iloc[:,1:]


# In[126]:


random_data_de_name


# In[10]:


# 第二步：分割数据为训练集、验证集和测试集
train_data, temp_data = train_test_split(random_data_de_name, test_size=0.3, random_state=42) 
valid_data, test_data = train_test_split(temp_data, test_size=0.5, random_state=42)


# In[12]:


label = 'Cp'


# In[ ]:


# 第三步：配置并启动AutoGluon训练
predictor = TabularPredictor(label=label,eval_metric="r2").fit(train_data, tuning_data=valid_data, presets='best_quality',use_bag_holdout=True)


# In[19]:


print('Evaluation on test data:')
scores = predictor.evaluate(test_data)
print(scores)


# ### 评估模型

# In[21]:


# 第四步：评估模型
y_train_pred = predictor.predict(train_data)
y_test_pred = predictor.predict(test_data)
y_valid_pred = predictor.predict(valid_data)


# In[22]:


print('r2 score on train data: ', r2_score(train_data[label], y_train_pred))
print('r2 score on test data: ', r2_score(test_data[label], y_test_pred))
print('r2 score on valid data: ', r2_score(valid_data[label], y_valid_pred))


# In[23]:


predictor.leaderboard()


# In[24]:


predictor.feature_importance(valid_data)


# In[25]:


# 第六步：图形化模型解释
importance_df = predictor.feature_importance(valid_data)
plt.barh(importance_df.index, importance_df['importance'])
plt.xlabel('Feature Importance')
plt.show() 


# ### 进阶版

# In[26]:


predictor.fit_summary()


# In[27]:


predictor.leaderboard(test_data, extra_metrics=['mse', 'rmse', 'mae', 'median_absolute_error',  'mape', 'spearmanr', 'pearsonr'], silent=True)


# ### 保留最佳模型

# In[28]:


# 删除其余模型（减少内存开销）
predictor.delete_models(models_to_keep='best')
# 输出最优模型
predictor.get_model_best()


# In[44]:


predictor.fit_summary(show_plot=True)


# In[143]:


predictor.leaderboard(extra_info=True, silent=True)


# ### 合并预测值列表

# In[ ]:


data_10_nolab = data_select_10_new.drop(labels = [label],axis = 1)
print(data_10_nolab)


# In[34]:


data_nolab = random_data_de_name.drop(labels=[label],axis = 1)
print(data_nolab)


# In[35]:


y_predict = predictor.predict(data_nolab)


# In[39]:


y_predict


# In[40]:


data_10_merge =random_data.copy()


# In[41]:


data_10_merge.insert(2,'Predict',y_predict)


# In[42]:


data_10_merge


# ### 准确率图

# In[47]:


#创建 DataFrame 存储预测结果和真实值
train_df = pd.DataFrame({'train_predictions': y_train_pred, 'y_train': train_data[label]})
valid_df = pd.DataFrame({'valid_predictions': y_valid_pred, 'y_valid': valid_data[label]})
test_df = pd.DataFrame({'test_predictions': y_test_pred, 'y_test': test_data[label]})


# In[48]:


plt.figure(figsize=(10, 6))
plt.subplot(131)
plt.scatter(train_df['train_predictions'], train_df['y_train'], label="Train samples", c="#d95f02") 
plt.xlabel("Predicted value") 
plt.ylabel("True value") 
plt.grid(True, linestyle = '--') 
plt.legend() 
plt.plot([0, 500], [0, 500], c="k", zorder=0) 
plt.xlim([0, 500]) 
plt.ylim([0, 500])
plt.title(f'Train R2: {r2_score(train_df["y_train"], train_df["train_predictions"]):.2f}')

plt.subplot(132)
plt.scatter(valid_df['valid_predictions'], valid_df['y_valid'], label="Validation samples", c="#1b9e77") 
plt.xlabel("Predicted value") 
plt.ylabel("True value") 
plt.grid(True, linestyle = '--') 
plt.legend() 
plt.plot([0, 500], [0, 500], c="k", zorder=0) 
plt.xlim([0, 500]) 
plt.ylim([0, 500])
plt.title(f'Validation R2: {r2_score(valid_df["y_valid"], valid_df["valid_predictions"]):.2f}')

plt.subplot(133)
plt.scatter(test_df['test_predictions'], test_df['y_test'], label="Test samples", c="#7570b3") 
plt.xlabel("Predicted value") 
plt.ylabel("True value") 
plt.grid(True, linestyle = '--') 
plt.legend() 
plt.plot([0, 500], [0, 500], c="k", zorder=0) 
plt.xlim([0, 500]) 
plt.ylim([0, 500])
plt.title(f'Test R2: {r2_score(test_df["y_test"], test_df["test_predictions"]):.2f}')

plt.tight_layout()
plt.show()


# In[138]:


#创建四个子图用于展示不同数据集的结果
fig, axs = plt.subplots(ncols=4, figsize=(15,4))

#训练集绘图
axs[0].scatter(train_df['train_predictions'], train_df['y_train'], label="Train samples", c="#d95f02") 
axs[0].set_xlabel("Predicted value") 
axs[0].set_ylabel("True value") 
axs[0].grid(True, linestyle = '--') 
axs[0].legend() 
axs[0].plot([0, 500], [0, 500], c="k", zorder=0) 
axs[0].set_xlim([0, 500]) 
axs[0].set_ylim([0, 500]) 
axs[0].set_title(f'Train R2: {r2_score(train_df["y_train"], train_df["train_predictions"]):.2f}')

#验证集绘图
axs[1].scatter(valid_df['valid_predictions'], valid_df['y_valid'], label="Validation samples", c="#1b9e77") 
axs[1].set_xlabel("Predicted value") 
axs[1].set_ylabel("True value") 
axs[1].grid(True, linestyle = '--') 
axs[1].legend() 
axs[1].plot([0, 500], [0, 500], c="k", zorder=0) 
axs[1].set_xlim([0, 500]) 
axs[1].set_ylim([0, 500]) 
axs[1].set_title(f'Validation R2: {r2_score(valid_df["y_valid"], valid_df["valid_predictions"]):.2f}')

#测试集绘图
axs[2].scatter(test_df['test_predictions'], test_df['y_test'], label="Test samples", c="#7570b3") 
axs[2].set_xlabel("Predicted value") 
axs[2].set_ylabel("True value") 
axs[2].grid(True, linestyle = '--') 
axs[2].legend() 
axs[2].plot([0, 500], [0, 500], c="k", zorder=0) 
axs[2].set_xlim([0, 500]) 
axs[2].set_ylim([0, 500]) 
axs[2].set_title(f'Test R2: {r2_score(test_df["y_test"], test_df["test_predictions"]):.2f}')



# 组合所有数据集的预测结果
all_predictions = np.concatenate((y_train_pred, y_valid_pred, y_test_pred)) 
all_y = np.concatenate((train_data[label], valid_data[label], test_data[label]))

# 绘制测试值
axs[3].scatter(y_test_pred, test_data[label], label='Test samples', c="#FF0000")
# 绘制训练值
axs[3].scatter(y_train_pred, train_data[label], label='Train samples', c="#00FF00")
# 绘制真实值
axs[3].scatter(y_valid_pred, valid_data[label], label='Validation samples', c="#0000FF")

axs[3].set_xlabel("Predicted value") 
axs[3].set_ylabel("True value") 
axs[3].grid(True, linestyle='--') 
axs[3].legend() 
axs[3].plot([0, 500], [0, 500], c="k", zorder=0) 
axs[3].set_xlim([0, 500]) 
axs[3].set_ylim([0, 500]) 
# 组合所有数据集的预测结果
all_predictions = np.concatenate((y_train_pred, y_valid_pred, y_test_pred)) 
all_y = np.concatenate((train_data[label], valid_data[label], test_data[label]))
axs[3].set_title(f'All R2: {r2_score(all_y, all_predictions):.2f}')
#添加标题和调整子图间距
fig.suptitle('Regression Results', fontsize=16) 
fig.tight_layout(rect=[0, 0.03, 1, 0.95])
#展示图形
plt.show()


# ### 模型解释

# In[64]:


import shap
import time
shap.initjs()

import warnings
warnings.filterwarnings('ignore')


# In[165]:


feature_names = train_data.iloc[:,1:].columns
train_data = train_data.copy()

val_data = test_data.iloc[:,1:].copy()

display(train_data.head())
display(val_data)


# In[86]:


class AutogluonWrapper:
    def __init__(self, predictor, feature_names):
        self.ag_model = predictor
        self.feature_names = feature_names
    
    def predict(self, X):
        if isinstance(X, pd.Series):
            X = X.values.reshape(1,-1)
        if not isinstance(X, pd.DataFrame):
            X = pd.DataFrame(X, columns=self.feature_names)
        return self.ag_model.predict(X)


# In[87]:


ag_wrapper = AutogluonWrapper(predictor, feature_names)


# In[92]:


feature_names


# In[93]:


# 创建SHAP解释器
explainer = shap.KernelExplainer(ag_wrapper.predict,shap.kmeans(train_data.iloc[:,1:], 10))


# In[ ]:





# In[94]:


NSHAP_SAMPLES = 100  # how many samples to use to approximate each Shapely value, larger values will be slower
N_VAL = 100  # how many datapoints from validation data should we interpret predictions for, larger values will be slower


# In[98]:


train


# In[162]:


#ROW_INDEX = 0  # index of an example datapoint
single_datapoint = train.iloc[:100,:]
single_prediction = ag_wrapper.predict(single_datapoint)

shap_values_single = explainer.shap_values(single_datapoint, nsamples=NSHAP_SAMPLES)
shap.force_plot(explainer.expected_value, shap_values_single, train.iloc[:100,:],show=False)
plt.savefig('shap/force_plot.jpg')


# In[166]:


#我们还可以绘制在许多预测中聚合的内核SHAP解释，比如在验证数据的第一个N_VAL数据点。
shap_values = explainer.shap_values(val_data.iloc[0:N_VAL,:], nsamples=NSHAP_SAMPLES)
shap.force_plot(explainer.expected_value, shap_values, val_data.iloc[0:N_VAL,:],show=False)
plt.savefig('shap/force_plot_1.jpg')


# In[167]:


shap_values = explainer.shap_values(val_data.iloc[0:N_VAL,:], nsamples=NSHAP_SAMPLES)


# In[168]:


import matplotlib.pyplot as plt
import shap

fig, ax = plt.subplots(figsize=(10, 6))
shap.summary_plot(shap_values, val_data.iloc[0:N_VAL,:], color_bar=False, show=False)
# 创建一个ScalarMappable对象
mappable = plt.cm.ScalarMappable(cmap=plt.cm.RdBu)
mappable.set_array([shap_values.min(), shap_values.max()])
# 添加Colorbar
cb = plt.colorbar(mappable, ax=ax)
cb.set_label('SHAP Value', fontsize=12)
cb.ax.tick_params(labelsize=10)

plt.show()
plt.savefig('shap/shap_value.jpg')


# In[169]:


shap.summary_plot(shap_values, val_data, plot_type='bar',color_bar=True,color='red',show=False)


# In[170]:


# 创建数据框
shap_df = pd.DataFrame(shap_values, columns=val_data.columns.tolist())

# 对每列特征按SHAP值绝对值进行排序
shap_df_abs_mean = shap_df.abs().mean().sort_values(ascending=True)

# 按重要性排序后的特征名称和SHAP值绝对值
feature_names = shap_df_abs_mean.index.tolist()
abs_means = shap_df_abs_mean.values.tolist()

# 生成柱状图
shap.summary_plot(shap_values, val_data, plot_type='bar', color_bar=True, color='red', show=False)

# 添加数据值
for i, val in enumerate(abs_means):
    plt.text(val + 0.03, i, round(val, 2))
plt.show()
plt.savefig('shap/shap_bar.jpg')


# In[107]:


import shap
shap.dependence_plot("ATS0pe", shap_values, val_data.iloc[0:N_VAL,:])


# In[109]:


#val_data[label] = y_valid  # add labels to validation DataFrame
predictor.feature_importance(test_data)


# ### waterfall

# In[110]:


# not support waterfall() but waterfall_legacy()
# shap.plots.waterfall(kernel_shap_values[0])
shap.plots._waterfall.waterfall_legacy(
    explainer.expected_value, 
    shap_values[0], 
    feature_names=val_data.columns)


# ### dependence_plot

# In[111]:


shap.dependence_plot("ATS0pe", shap_values, val_data[:N_VAL])


# ### decision_plot

# In[112]:


shap.decision_plot(explainer.expected_value, shap_values, val_data[:N_VAL],show= False)
plt.savefig('shap/decision_plot/jpg')


# In[113]:


predictor.get_model_best()


# In[114]:


predictor.fit_weighted_ensemble()


# In[115]:


predictor.leaderboard(test_data,extra_info = False)


# In[127]:


random_data_de_name


# In[140]:


#相关矩阵的热力图
import seaborn as sns
#d = random_data_de_name.loc[:,["Cp","ATS0pe","NumValenceElectrons","nRot","mZagreb1","nHBAcc","NssCH2","GGI4","Kappa2","Temperature","SlogP_VSA2"]].corr()

plt.subplots(figsize = (10,10))
sns.heatmap(random_data_de_name.corr(),annot = True,vmax = 1,square = True,cmap = "viridis")
plt.show()


# In[142]:


import seaborn as sns 
d = random_data_de_name.loc[:,["Cp","ATS0pe","NumValenceElectrons","nRot","mZagreb1","nHBAcc","NssCH2","GGI4","Kappa2","Temperature","SlogP_VSA2"]].corr()

plt.subplots(figsize = (10,10)) 
sns.heatmap(random_data_de_name.corr(),annot = True,vmax = 1,square = True,cmap = "viridis", mask = np.triu(np.ones_like(d, dtype=bool))) 
plt.show()


# ### eda

# In[195]:


import autogluon.eda.analysis as eda
import autogluon.eda.visualization as viz
import autogluon.eda.auto as auto

auto.analyze(
    train_data=train_data, label="Cp",
    anlz_facets=[
        # Puts output into the root namespace
        eda.interaction.Correlation(),
        # Puts output into the focus namespace
        eda.Namespace(namespace='focus', children=[
            eda.interaction.Correlation(focus_field='Fare', focus_field_threshold=0.3),
        ])
    ],
    viz_facets=[
        # Renders correlations from the root namespace
        viz.interaction.CorrelationVisualization(),
        # Renders correlations from the focus namespace
        viz.interaction.CorrelationVisualization(namespace='focus'),
    ]
)

plt.savefig('shap/hot_plot.jpg')
