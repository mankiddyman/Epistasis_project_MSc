import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

#read in the output, regulator and sensor data into their own dataframes 
path_main = "../data/initial_data/mutants_separated"

#create lists of column names for the dataframes that we will put the data into
Outputs_list=[]
Regulator_list = []
Sensor_list=[]
for i in range(10):
    Outputs_list.append('Output_'+str(i+1))
    Regulator_list.append('Regulator'+str(i+1))
    Sensor_list.append('Sensor_'+str(i+1))

#create empty dataframes using the column nmaes created
MT_sensor_df = pd.DataFrame(columns = Sensor_list)
MT_regulator_df = pd.DataFrame(columns = Regulator_list)
MT_output_df = pd.DataFrame(columns = Outputs_list)
MT_inducer_df = pd.DataFrame(columns = ["Inducer_conc"])

#inducer concentration is the same for each mutant so only one column is needed
MT_inducer_df["Inducer_conc"] = pd.read_table(path_main +"/Output1.dat",delim_whitespace= True, header=None).iloc[:,0]

#fill the dataframes with the mutant data for each of output, regulator and sensor
for i in range(10):
    path_out = path_main + "/Output" + str(i+1) + ".dat"
    path_reg = path_main + "/Regulator" + str(i+1) + ".dat"
    path_ind = path_main + "/Sensor" + str(i+1) + ".dat"
    
    sens_df_col = "Sensor_"+str(i+1)
    reg_df_col = "Regulator_"+str(i+1)
    out_df_col = "Output_"+str(i+1)
    
    MT_sensor_df[sens_df_col] = pd.read_table(path_sens, delim_whitespace= True, header=None).iloc[:,1]
    MT_regulator_df[reg_df_col] = pd.read_table(path_reg, delim_whitespace= True, header=None).iloc[:,1]
    MT_output_df[out_df_col] = pd.read_table(path_out, delim_whitespace= True, header=None).iloc[:,1]


#define scatter plotting function with log scales
def WT_Plotter(ax,y, params):
    for mut in y:

    out = ax.scatter(wt_inducer, wt_df[y], **params, marker = 'o')
    xScale = ax.set_xscale('log')
    yScale = ax.set_yscale('log')
    xlab = ax.set_xlabel("Inducer concentration (%)")
    ylab = ax.set_ylabel(y)
    #xAxis = ax.xticks([1E-6, 1E-5, 1E-3, 1E-1], [0, 1E-5, 1E-3, 1E-1])
    return out, xScale, yScale, #x