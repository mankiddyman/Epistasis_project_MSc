import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

#%% read in the output, regulator and sensor mutant data to a single dataframe
#MT_sensor_df, MT_regulator_df and MT_output_df
path_main = "../data/initial_data/mutants_separated"

#create lists of column names for the dataframes that we will put the data into
data_list=['inducer']
for i in range(10):
    data_list.append('Output_'+str(i+1))
    data_list.append('Regulator_'+str(i+1))
    data_list.append('Sensor_'+str(i+1))

#create empty dataframes for the data including repeats using the column names created
MT_df_reps = pd.DataFrame(columns = data_list)

#fill the dataframes with the mutant data for each of output, regulator and sensor
for i in range(10):
    path_out = path_main + "/Output" + str(i+1) + ".dat"
    path_reg = path_main + "/Regulator" + str(i+1) + ".dat"
    path_sens = path_main + "/Sensor" + str(i+1) + ".dat"
    
    sens_df_col = "Sensor_"+str(i+1)
    reg_df_col = "Regulator_"+str(i+1)
    out_df_col = "Output_"+str(i+1)
    
    MT_df_reps["inducer"]= (pd.read_table(path_sens, delim_whitespace= True, header=None).iloc[:,0])
    MT_df_reps[sens_df_col] = (pd.read_table(path_sens, delim_whitespace= True, header=None).iloc[:,1])
    MT_df_reps[reg_df_col] = (pd.read_table(path_reg, delim_whitespace= True, header=None).iloc[:,1])
    MT_df_reps[out_df_col] = (pd.read_table(path_out, delim_whitespace= True, header=None).iloc[:,1])

#create empty dataframes for the data averaging for repeats using the column names created
MT_df = pd.DataFrame(columns = data_list)
for i in range(16):
    row = MT_df_reps.iloc[[0+i,16+i,32+i],:].sum()/3
    MT_df = pd.concat([MT_df, row.to_frame().T], ignore_index=True)

#%% 
#%% send dataframe to data file
MT_df....
#%%  