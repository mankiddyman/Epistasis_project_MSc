import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from data_wrangling import *
import seaborn as sns

# code that generates violin and pair plots for single mutant data

#%%
# loading single mutant dataframe
data = pd.read_excel('../data/SM_params.xlsx')
# cleaning dataframe, can be applied to large dataframes as well
# removing wt row 
data.drop(data[data['Unnamed: 0'] == 'WT'].index, inplace = True)
# removing rss column 
df = data.drop(columns=['rss'])
# renaming first column to Node
df.columns = df.columns.str.replace('Unnamed: 0', 'Node')
# explain 

#%%
new_df = df
# categorising mutants by their names into Output, Regulator and Sensor 
new_df.loc[new_df["Node"].str.contains("Output"), "Node"] = "Output"
new_df.loc[new_df["Node"].str.contains("Regulator"), "Node"] = "Regulator"
new_df.loc[new_df["Node"].str.contains("Sensor"), "Node"] = "Sensor"

# Violin Plots
param_names= list(new_df.columns.values)
param_names.remove('Node')
for i in param_names:
    plt.figure()
    vp = sns.violinplot(data=new_df, x=i, y="Node", hue="Node")
    plt.title("Violin Plot for parameter " + i)
#%%
# modifying dataframe for pairplots 
# creating separate dataframes for each node 

# Output Pair Plot

# output dataframe
# removing regulator and sensor data
pp_df = df
output_data = pp_df
output_data = output_data[output_data["Node"].str.contains("Regulator")==False]    
output_data = output_data[output_data["Node"].str.contains("Sensor")==False] 
output_df = output_data.drop(columns=['Node'])
output_pairplot = output_df
sns.pairplot(output_pairplot)

# Regulator Pair Plot

# regulator dataframe
# removing output and sensor data
regulator_data = pp_df
regulator_data = regulator_data[regulator_data["Node"].str.contains("Output")==False]    
regulator_data = regulator_data[regulator_data["Node"].str.contains("Sensor")==False] 
regulator_df = regulator_data.drop(columns=['Node'])
regulator_pairplot = regulator_df
sns.pairplot(regulator_pairplot)

# Sensor Pair Plot

# sensor dataframe
# removing output and regulator data
sensor_data = pp_df
sensor_data = sensor_data[sensor_data["Node"].str.contains("Output")==False]  
sensor_data = sensor_data[sensor_data["Node"].str.contains("Regulator")==False]  
sensor_df = sensor_data.drop(columns=['Node'])
sensor_pairplot = sensor_df
sns.pairplot(sensor_pairplot)