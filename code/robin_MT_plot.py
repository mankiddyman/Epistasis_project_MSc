import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

#read in the output, regulator and sensor data into their own dataframes 
path_main = "../data/initial_data/mutants_separated"

Outputs_list=[]
Regulator_list = []
Sensor_list=[]
for i in range(10):
    Outputs_list.append('Output_'+str(i+1))
    Regulator_list.append('Regulator'+str(i+1))
    Sensor_list.append('Sensor_'+str(i+1))

MT_sensor_df = pd.DataFrame(columns = Sensor_list)
MT_regulator_df = pd.DataFrame(columns = Regulator_list)
MT_output_df = pd.DataFrame(columns = Outputs_list)

for i in range(10):
    path_out = path_main + "/Output" + str(i+1) + ".dat"
    path_reg = path_main + "/Regulator" + str(i+1) + ".dat"
    path_sens = path_main + "/Sensor" + str(i+1) + ".dat"
    
    sens_df_col = "Sensor_"+str(i+1)
    reg_df_col = "Regulator_"+str(i+1)
    out_df_col = "Output_"+str(i+1)
    
    MT_sensor_df[sens_df_col] = pd.read_table(path_sens, delim_whitespace= True, header=None).iloc[:,1]
    MT_regulator_df[reg_df_col] = pd.read_table(path_reg, delim_whitespace= True, header=None).iloc[:,1]
    MT_output_df[out_df_col] = pd.read_table(path_out, delim_whitespace= True, header=None).iloc[:,1]

