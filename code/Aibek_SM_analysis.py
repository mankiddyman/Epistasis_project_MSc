import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from data_wrangling import *
import seaborn as sns

# code that generates violin and pair plots for single mutant data
# automate process for different/bigger dataframes
#%%
# loade single mutant dataframe
data = pd.read_excel('../data/SM_params.xlsx')
# clean dataframe, can be applied to large dataframes as well
# remove wt row 
data.drop(data[data['Unnamed: 0'] == 'WT'].index, inplace = True)
# remove rss column 
df = data.drop(columns=['rss'])
# rename first column to Node
df.columns = df.columns.str.replace('Unnamed: 0', 'Node')

#%%
new_df = df
# categorise mutants by their names into Output, Regulator and Sensor 
new_df.loc[new_df["Node"].str.contains("Output"), "Node"] = "Output"
new_df.loc[new_df["Node"].str.contains("Regulator"), "Node"] = "Regulator"
new_df.loc[new_df["Node"].str.contains("Sensor"), "Node"] = "Sensor"

# Violin Plots
final_df = new_df
param_names= list(new_df.columns.values)
param_names.remove('Node')
new_param_names = param_names
for j in new_param_names:
    if j.find("A") != -1:
        new_param_names[new_param_names.index(j)]='A'
    elif j.find("B") != -1:
        new_param_names[new_param_names.index(j)]='B'
    elif j.find("C") != -1:
        new_param_names[new_param_names.index(j)]='C'
    elif j.find("N") != -1:
        new_param_names[new_param_names.index(j)]='N'
    elif j.find("F") != -1:
        new_param_names[new_param_names.index(j)]='F'
        
# remove repetitions in new_param_names      
new_param_names = list(dict.fromkeys(new_param_names))

# create new dataframe from lists
# columns of new dataframe: Node, A, B, C, N, F       
# final_df.columns[final_df.columns.str.contains("A")] = "A"

# change columns to lists
# truncate lists 
# automate process for different/bigger dataframes
F_o_list = final_df['F_o'].tolist()

# Sensor
sensor_vp_df = final_df
sensor_vp_df = sensor_vp_df[sensor_vp_df["Node"].str.contains("Output")==False] 
sensor_vp_df = sensor_vp_df[sensor_vp_df["Node"].str.contains("Regulator")==False] 
A_s_list = sensor_vp_df['A_s'].tolist()
B_s_list = sensor_vp_df['B_s'].tolist()
C_s_list = sensor_vp_df['C_s'].tolist()
N_s_list = sensor_vp_df['N_s'].tolist()
f_o_sensor_list = F_o_list[20:30]

# Regulator
regulator_vp_df = final_df
regulator_vp_df = regulator_vp_df[regulator_vp_df["Node"].str.contains("Output")==False] 
regulator_vp_df = regulator_vp_df[regulator_vp_df["Node"].str.contains("Sensor")==False] 
A_r_list = regulator_vp_df['A_r'].tolist()
B_r_list = regulator_vp_df['B_r'].tolist()
C_r_list = regulator_vp_df['C_r'].tolist()
N_r_list = regulator_vp_df['N_r'].tolist()
f_o_regulator_list = F_o_list[10:20]

# Output
output_vp_df = final_df
output_vp_df = output_vp_df[output_vp_df["Node"].str.contains("Regulator")==False] 
output_vp_df = output_vp_df[output_vp_df["Node"].str.contains("Sensor")==False] 
A_o_list = output_vp_df['A_o'].tolist()
B_o_list = output_vp_df['B_o'].tolist()
C_o_list = output_vp_df['C_o'].tolist()
N_o_list = output_vp_df['N_o'].tolist()
f_o_output_list = F_o_list[0:10]

# Output half
output_half_vp_df = final_df
output_half_vp_df = output_half_vp_df[output_half_vp_df["Node"].str.contains("Regulator")==False]
output_half_vp_df = output_half_vp_df[output_half_vp_df["Node"].str.contains("Sensor")==False] 
A_h_list = output_half_vp_df['A_h'].tolist()
B_h_list = output_half_vp_df['B_h'].tolist()
C_h_list = output_half_vp_df['C_h'].tolist()
N_h_list = output_half_vp_df['N_o'].tolist()
f_o_half_list = F_o_list[0:10]

# create columns for the new dataframe

# column 'A'
A_o_list.extend(A_h_list)
A_o_list.extend(A_r_list)
A_o_list.extend(A_s_list)
A_list = A_o_list

# column 'B'
B_o_list.extend(B_h_list)
B_o_list.extend(B_r_list)
B_o_list.extend(B_s_list)
B_list = B_o_list

# column 'C'
C_o_list.extend(C_h_list)
C_o_list.extend(C_r_list)
C_o_list.extend(C_s_list)
C_list = C_o_list

# column 'N'
N_o_list.extend(N_h_list)
N_o_list.extend(N_r_list)
N_o_list.extend(N_s_list)
N_list = N_o_list

# column 'F'
f_o_output_list.extend(f_o_half_list)
f_o_output_list.extend(f_o_regulator_list)
f_o_output_list.extend(f_o_sensor_list)
F_list = f_o_output_list 

# column 'node'
node_list = final_df['Node'].tolist()
for item in node_list:
    if item == "Output":
        node_list.insert(10, "Output Half")
    else:
        break

# check length of columns
a = len(A_list) 
b = len(B_list)
c = len(C_list)
n = len(N_list)
f = len(F_list)
node = len(node_list)
# create new dataframe from lists
# rename columns in dataframe
param_df = pd.DataFrame(list(zip(node_list, A_list, B_list, C_list, N_list, F_list)),columns =['Node', 'A', 'B', 'C', 'N', 'F'])

#%%
# combined violin plots
for i in new_param_names:
    plt.figure()
    vp = sns.violinplot(data=param_df, x=i, y="Node", hue="Node")
    plt.title("Violin Plot for parameter " + i)

#%%
# separate violin plots
for i in param_names:
    plt.figure()
    vp = sns.violinplot(data=new_df, x=i, y="Node", hue="Node")
    plt.title("Violin Plot for parameter " + i)
    
#%%

single_node_list = node_list
for i in range(40):
    single_node_list.insert(10, "Output")
for j in range(40):
    single_node_list.insert(50, "Output Half")
for k in range(40):
    single_node_list.insert(100, "Regulator")
for l in range(40):
    single_node_list.insert(150, "Sensor")
# Output
output_list = output_vp_df['A_o'].tolist()
output_list.extend(output_vp_df['B_o'].tolist())
output_list.extend(output_vp_df['C_o'].tolist())
output_list.extend(output_vp_df['N_o'].tolist())
output_list.extend(F_o_list[0:10])
# Output Half
output_list.extend(output_half_vp_df['A_h'].tolist())
output_list.extend(output_half_vp_df['B_h'].tolist())
output_list.extend(output_half_vp_df['C_h'].tolist())
output_list.extend(output_half_vp_df['N_o'].tolist())
output_list.extend(F_o_list[0:10])
# Regulator
output_list.extend(regulator_vp_df['A_r'].tolist())
output_list.extend(regulator_vp_df['B_r'].tolist())
output_list.extend(regulator_vp_df['C_r'].tolist())
output_list.extend(regulator_vp_df['N_r'].tolist())
output_list.extend(F_o_list[10:20])
# Sensor
output_list.extend(sensor_vp_df['A_s'].tolist())
output_list.extend(sensor_vp_df['B_s'].tolist())
output_list.extend(sensor_vp_df['C_s'].tolist())
output_list.extend(sensor_vp_df['N_s'].tolist())
output_list.extend(F_o_list[20:30])

# range every 10 items 
# if index = 10,20,.. plot in range
# while 
# index range
A_o = ['A']*10
B_o = ['B']*10
C_o = ['C']*10
N_o = ['N']*10
F_o = ['F']*10
A_h = ['A']*10
B_h = ['B']*10 
C_h = ['C']*10
N_h = ['N']*10
F_h = ['F']*10
A_r = ['A']*10
B_r = ['B']*10
C_r = ['C']*10
N_r = ['N']*10
F_r = ['F']*10
A_s = ['A']*10
B_s = ['B']*10
C_s = ['C']*10
N_s = ['N']*10
F_s = ['F']*10
param_list = [B_o,C_o,N_o,F_o,A_h,B_h,C_h,N_h,F_h,A_r,B_r,C_r,N_r,F_r,A_s,B_s,C_s,N_s,F_s]
for i in param_list:
    A_o.extend(i)

# log
log_output = list(np.log(output_list))    
single_df = pd.DataFrame(list(zip(single_node_list, A_o, log_output)),columns =['Node', 'Parameter', 'Value'])
sns.set_style("dark")
sns.set_context("paper", rc={"lines.linewidth": 1.1})
plt.figure()
vp = sns.violinplot(data=single_df, x="Node", y="Value", hue="Parameter", scale="count")
sns.move_legend(vp, "lower right")
vp.set_ylim([-20, 15])  
[vp.axvline(x+.5,color='k') for x in vp.get_xticks()]
plt.title("Violin Plots", fontweight='bold')
 
#%%
# modify dataframe for pairplots 
# create separate dataframes for each node 

# Output Pair Plot

# output dataframe
# remove output half, regulator, sensor data
pp_df = param_df
output_data = pp_df
output_data = output_data[output_data["Node"].str.contains("Regulator")==False]    
output_data = output_data[output_data["Node"].str.contains("Sensor")==False] 
output_data = output_data[output_data["Node"].str.contains("Output Half")==False] 
output_df = output_data.drop(columns=['Node'])
output_pairplot = output_df
sns.pairplot(output_pairplot)

# Regulator Pair Plot

# regulator dataframe
# remove output, output half, sensor data
regulator_data = pp_df
regulator_data = regulator_data[regulator_data["Node"].str.contains("Output")==False]  
regulator_data = regulator_data[regulator_data["Node"].str.contains("Output Half")==False]   
regulator_data = regulator_data[regulator_data["Node"].str.contains("Sensor")==False] 
regulator_df = regulator_data.drop(columns=['Node'])
regulator_pairplot = regulator_df
sns.pairplot(regulator_pairplot)

# Sensor Pair Plot

# sensor dataframe
# remove output, output half, regulator data
sensor_data = pp_df
sensor_data = sensor_data[sensor_data["Node"].str.contains("Output")==False]  
sensor_data = sensor_data[sensor_data["Node"].str.contains("Output Half")==False]
sensor_data = sensor_data[sensor_data["Node"].str.contains("Regulator")==False]  
sensor_df = sensor_data.drop(columns=['Node'])
sensor_pairplot = sensor_df
sns.pairplot(sensor_pairplot)

# Output Half Pair Plot

# output half dataframe
# remove output, regulator, sensor
output_half_data = pp_df
output_half_data.drop(output_half_data[output_half_data['Node'] == 'Output'].index, inplace = True) 
output_half_data = output_half_data[output_half_data["Node"].str.contains("Regulator")==False] 
output_half_data = output_half_data[output_half_data["Node"].str.contains("Sensor")==False]  
output_half_df = output_half_data.drop(columns=['Node'])
output_half_pairplot = output_half_df
sns.pairplot(output_half_pairplot)