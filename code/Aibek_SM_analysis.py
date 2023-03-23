import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

"""
# code that generates violin and pair plots for single mutant data
# automate process for different/bigger dataframes
"""

#%%
# Hill 

def hill_vp(data):
    data.drop(data[data['Unnamed: 0'] == 'WT'].index, inplace = True)
    df = data.drop(columns=['rss'])
    df = data.drop(columns=['time_elapsed_s'])
    df.columns = df.columns.str.replace('Unnamed: 0', 'Node')
    new_df = df
    new_df.loc[new_df["Node"].str.contains("Output"), "Node"] = "Output"
    new_df.loc[new_df["Node"].str.contains("Regulator"), "Node"] = "Regulator"
    new_df.loc[new_df["Node"].str.contains("Sensor"), "Node"] = "Sensor"
    dff = new_df
    param_names= list(dff.columns.values)
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
    new_param_names = list(dict.fromkeys(new_param_names))
    F_o_list = dff['F_o'].tolist()
    # Sensor
    sensor_vp_df = dff
    sensor_vp_df = sensor_vp_df[sensor_vp_df["Node"].str.contains("Output")==False] 
    sensor_vp_df = sensor_vp_df[sensor_vp_df["Node"].str.contains("Regulator")==False] 
    A_s_list = sensor_vp_df['A_s'].tolist()
    B_s_list = sensor_vp_df['B_s'].tolist()
    C_s_list = sensor_vp_df['C_s'].tolist()
    N_s_list = sensor_vp_df['N_s'].tolist()
    f_o_sensor_list = F_o_list[20:30]
    # Regulator
    regulator_vp_df = dff
    regulator_vp_df = regulator_vp_df[regulator_vp_df["Node"].str.contains("Output")==False] 
    regulator_vp_df = regulator_vp_df[regulator_vp_df["Node"].str.contains("Sensor")==False] 
    A_r_list = regulator_vp_df['A_r'].tolist()
    B_r_list = regulator_vp_df['B_r'].tolist()
    C_r_list = regulator_vp_df['C_r'].tolist()
    N_r_list = regulator_vp_df['N_r'].tolist()
    f_o_regulator_list = F_o_list[10:20]
    # Output
    output_vp_df = dff
    output_vp_df = output_vp_df[output_vp_df["Node"].str.contains("Regulator")==False] 
    output_vp_df = output_vp_df[output_vp_df["Node"].str.contains("Sensor")==False] 
    A_o_list = output_vp_df['A_o'].tolist()
    B_o_list = output_vp_df['B_o'].tolist()
    C_o_list = output_vp_df['C_o'].tolist()
    N_o_list = output_vp_df['N_o'].tolist()
    f_o_output_list = F_o_list[0:10]
    # create columns for the new dataframe
    # column 'A'
    A_o_list.extend(A_r_list)
    A_o_list.extend(A_s_list)
    A_list = A_o_list
    # column 'B'
    B_o_list.extend(B_r_list)
    B_o_list.extend(B_s_list)
    B_list = B_o_list
    # column 'C'
    C_o_list.extend(C_r_list)
    C_o_list.extend(C_s_list)
    C_list = C_o_list
    # column 'N'
    N_o_list.extend(N_r_list)
    N_o_list.extend(N_s_list)
    N_list = N_o_list
    # column 'F'
    f_o_output_list.extend(f_o_regulator_list)
    f_o_output_list.extend(f_o_sensor_list)
    F_list = f_o_output_list 
    # column 'node'
    node_list = dff['Node'].tolist()
    param_df = pd.DataFrame(list(zip(node_list, A_list, B_list, C_list, N_list, F_list)),columns =['Node', 'A', 'B', 'C', 'N', 'F'])
    single_node_list = node_list
    for i in range(40):
        single_node_list.insert(10, "Output")
    for j in range(40):
        single_node_list.insert(50, "Regulator")
    for k in range(40):
        single_node_list.insert(100, "Sensor")
    # Output
    output_list = output_vp_df['A_o'].tolist()
    output_list.extend(output_vp_df['B_o'].tolist())
    output_list.extend(output_vp_df['C_o'].tolist())
    output_list.extend(output_vp_df['N_o'].tolist())
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
    param_list = ['A']*10 + ['B']*10 + ['C']*10 + ['N']*10 + ['F']*10 + ['A']*10 + ['B']*10 + ['C']*10 + ['N']*10 + ['F']*10 + ['A']*10 + ['B']*10 + ['C']*10 + ['N']*10 + ['F']*10 
    # log
    log_output =map(math.log10,output_list)   
    single_df = pd.DataFrame(list(zip(single_node_list, param_list, log_output)),columns =['Node', 'Parameter', 'Value'])
    sns.set_style("dark")
    sns.set_context("paper", rc={"lines.linewidth": 1.1})
    plt.figure()
    vp = sns.violinplot(data=single_df, x="Node", y="Value", hue="Parameter", scale="count",bw=0.55)
    sns.move_legend(vp, "lower right")
    vp.set_ylim([-20, 20])  
    plt.xlabel('Node',fontsize=11)
    plt.ylabel('Parameter Value (log10 scale)',fontsize=11)
    #vp.set(xlabel='Node',ylabel='Parameter Value (log10 scale)')
    [vp.axvline(x+.5,color='k') for x in vp.get_xticks()]
    plt.title("Violin Plots for Hill Model", fontweight='bold',fontsize=12)
#%%
# Thermodynamic

def td_vp(data):
    data.drop(data[data['Unnamed: 0'] == 'WT'].index, inplace = True)
    # remove rss column 
    df = data.drop(columns=['time_elapsed_s'])
    dff = df.drop(columns=['rss'])
    # rename first column to Node
    dff.columns = dff.columns.str.replace('Unnamed: 0', 'Node')
    new_df = dff
    # categorise mutants by their names into Output, Regulator and Sensor 
    new_df.loc[new_df["Node"].str.contains("Output"), "Node"] = "Output"
    new_df.loc[new_df["Node"].str.contains("Regulator"), "Node"] = "Regulator"
    new_df.loc[new_df["Node"].str.contains("Sensor"), "Node"] = "Sensor"
    final_df = new_df
    param_names= list(new_df.columns.values)
    param_names.remove('Node')
    new_param_names = param_names
    for j in new_param_names:
        if j.find("A") != -1:
            new_param_names[new_param_names.index(j)]='A'
        elif j.find("C") != -1:
            new_param_names[new_param_names.index(j)]='C'
        elif j.find("K") != -1:
            new_param_names[new_param_names.index(j)]='K'
        elif j.find("P") != -1 and j.find("P_b") != -1:
            new_param_names[new_param_names.index(j)]='P_b'
        elif j.find("P") != -1 and j.find("P_b") == -1:
            new_param_names[new_param_names.index(j)]='P'
        elif j.find("F") != -1:
            new_param_names[new_param_names.index(j)]='F'
    new_param_names = list(dict.fromkeys(new_param_names))
    F_o_list = final_df['F_o'].tolist()
    # Sensor
    sensor_vp_df = final_df
    sensor_vp_df = sensor_vp_df[sensor_vp_df["Node"].str.contains("Output")==False] 
    sensor_vp_df = sensor_vp_df[sensor_vp_df["Node"].str.contains("Regulator")==False] 
    Pb_s_list = sensor_vp_df['P_b'].tolist()
    P_s_list = sensor_vp_df['P_u'].tolist()
    K_s_list = sensor_vp_df['K_12'].tolist()
    C_s_list = sensor_vp_df['C_pa'].tolist()
    A_s_list = sensor_vp_df['A_s'].tolist()
    f_o_sensor_list = F_o_list[20:30]
    # Regulator
    regulator_vp_df = final_df
    regulator_vp_df = regulator_vp_df[regulator_vp_df["Node"].str.contains("Output")==False] 
    regulator_vp_df = regulator_vp_df[regulator_vp_df["Node"].str.contains("Sensor")==False] 
    Pb_r_list = sensor_vp_df['P_b'].tolist()
    P_r_list = sensor_vp_df['P_r'].tolist()
    K_r_list = sensor_vp_df['K_t'].tolist()
    C_r_list = sensor_vp_df['C_pt'].tolist()
    A_r_list = sensor_vp_df['A_r'].tolist()
    f_o_regulator_list = F_o_list[10:20]
    # Output
    output_vp_df = final_df
    output_vp_df = output_vp_df[output_vp_df["Node"].str.contains("Regulator")==False] 
    output_vp_df = output_vp_df[output_vp_df["Node"].str.contains("Sensor")==False] 
    Pb_o_list = sensor_vp_df['P_b'].tolist()
    P_o_list = sensor_vp_df['P_o'].tolist()
    K_o_list = sensor_vp_df['K_l'].tolist()
    C_o_list = sensor_vp_df['C_pl'].tolist()
    A_o_list = sensor_vp_df['A_o'].tolist()
    f_o_output_list = F_o_list[0:10]
    # column 'A'
    A_o_list.extend(A_r_list)
    A_o_list.extend(A_s_list)
    A_list = A_o_list
    # column 'P'
    P_o_list.extend(P_r_list)
    P_o_list.extend(P_s_list)
    P_list = P_o_list
    # column 'C'
    C_o_list.extend(C_r_list)
    C_o_list.extend(C_s_list)
    C_list = C_o_list
    # column 'K'
    K_o_list.extend(K_r_list)
    K_o_list.extend(K_s_list)
    K_list = K_o_list
    # column 'F'
    f_o_output_list.extend(f_o_regulator_list)
    f_o_output_list.extend(f_o_sensor_list)
    F_list = f_o_output_list 
    # column 'P_b'
    Pb_o_list.extend(Pb_r_list)
    Pb_o_list.extend(Pb_s_list)
    Pb_list = Pb_o_list
    # column 'node'
    node_list = final_df['Node'].tolist()
    # new dataframe
    param_df = pd.DataFrame(list(zip(node_list, P_list, Pb_list, C_list, K_list, A_list, F_list)),columns =['Node', 'P','P_b','C','K','A','F'])
    single_node_list = node_list
    for i in range(50):
        single_node_list.insert(10, "Output")
    for j in range(50):
        single_node_list.insert(60, "Regulator")
    for k in range(50):
        single_node_list.insert(120, "Sensor")
    # Output
    output_list = output_vp_df['P_o'].tolist()
    output_list.extend(output_vp_df['P_b'].tolist())
    output_list.extend(output_vp_df['C_pl'].tolist())
    output_list.extend(output_vp_df['K_l'].tolist())
    output_list.extend(output_vp_df['A_o'].tolist())
    output_list.extend(F_o_list[0:10])
    # Regulator
    output_list.extend(regulator_vp_df['P_r'].tolist())
    output_list.extend(regulator_vp_df['P_b'].tolist())
    output_list.extend(regulator_vp_df['C_pt'].tolist())
    output_list.extend(regulator_vp_df['K_t'].tolist())
    output_list.extend(regulator_vp_df['A_r'].tolist())
    output_list.extend(F_o_list[10:20])
    # Sensor
    output_list.extend(sensor_vp_df['P_u'].tolist())
    output_list.extend(sensor_vp_df['P_b'].tolist())
    output_list.extend(sensor_vp_df['C_pa'].tolist())
    output_list.extend(sensor_vp_df['K_12'].tolist())
    output_list.extend(sensor_vp_df['A_s'].tolist())
    output_list.extend(F_o_list[20:30])
    P_o = ['P']*10+['P_b']*10+['C']*10+['K']*10+['A']*10+['F']*10+['P']*10+['P_b']*10+['C']*10+['K']*10+['A']*10+['F']*10+['P']*10+['P_b']*10+['C']*10+['K']*10+['A']*10+['F']*10
    log_output =map(math.log10,output_list)     
    single_df = pd.DataFrame(list(zip(single_node_list, P_o, log_output)),columns =['Node', 'Parameter', 'Value'])
    sns.set_style("dark")
    sns.set_context("paper", rc={"lines.linewidth": 1.1})
    plt.figure()
    single_dff = single_df.drop([80,82,85,89,101,106,107])
    vp = sns.violinplot(data=single_dff, x="Node", y="Value", hue="Parameter", scale="count", bw=0.35)
    sns.move_legend(vp, "lower right")
    vp.set_ylim([-20, 20])  
    [vp.axvline(x+.5,color='k') for x in vp.get_xticks()]
    plt.xlabel('Node',fontsize=11)
    plt.ylabel('Parameter Value (log10 scale)',fontsize=11)
    plt.title("Violin Plots for Thermodynamic Model", fontweight='bold',fontsize=12)
#%%
# violin plots
# data_hm = pd.read_excel('../results/model_hill.modelSM_params_all.xlsx')
# data_td = pd.read_excel('../results/model_thermodynamic.modelSM_params_sensor.xlsx')
# hm_vp = hill_vp(data_hm)
# td_vp = td_vp(data_td)