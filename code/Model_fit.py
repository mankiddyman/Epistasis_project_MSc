import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize
from Models import *

#use model chosen from Models.py to obtain fit of model to data

wt_df = pd.read_csv('../data/WT_single.csv')
wt_I=wt_df.S
wt_Sensor=wt_df.Sensor
wt_Regulator=wt_df.Regulator
wt_Output=wt_df.Output
wt_Stripe=wt_df.Stripe

# minimizing sum of sqrd deviations of log model and log data
data = wt_df
params = [2,1,1,1,1,1,1,1,1,1]
def min_fun(params,data=data):
    log_sen = np.log10(data.Sensor)
    log_reg = np.log10(data.Regulator)
    log_out = np.log10(data.Output)
    ind = data.S
    # params is a list of list where each sublist contains parameters for one equation
    Sensor_est, Regulator_est, Output_est = model_leaky(ind,params[0],params[1],params[2],params[3], params[4],params[5],params[6],params[7],params[8],params[9],params[10],params[11],option='full')
    log_sen_est = np.log10(Sensor_est)
    log_reg_est = np.log10(Regulator_est)
    log_out_est = np.log10(Output_est)
    result = 0
    for i in range(0,len(data.Sensor)):
        result += np.power((log_sen[i] - log_sen_est[i]), 2)
        result += np.power((log_reg[i] - log_reg_est[i]), 2)
        result += np.power((log_out[i] - log_out_est[i]), 2)
    return result




#now gonna do the scipy.optimize.minimize
#%%
min_result=minimize(min_fun,x0=[311,17220,1081,0.9,2000,6000,0.0421,0.308,2000,4054,5.25E-5,3.17],method='CG')

#plotting the predictions now


#generating estimates
Sensor_est_array,Regulator_est_array,Output_est_array=init_model(data.S,min_result.x[0],min_result.x[1],min_result.x[2],min_result.x[3],min_result.x[4],min_result.x[5],min_result.x[6],min_result.x[7],min_result.x[8],min_result.x[9],option="Full")

wt_Signal=wt_I
def WT_fit_plot(ax, y, params):
    return ax.plot(wt_Signal, y, **params)


#define scatter plotting function with log scales
def WT_Plotter(ax,y, params):
    out = ax.scatter(wt_Signal, wt_df[y], **params, marker = 'o')
    xScale = ax.set_xscale('log')
    yScale = ax.set_yscale('log')
    xlab = ax.set_xlabel("Inducer concentration (%)")
    ylab = ax.set_ylabel(y)
    #xAxis = ax.xticks([1E-6, 1E-5, 1E-3, 1E-1], [0, 1E-5, 1E-3, 1E-1])
    return out, xScale, yScale, #xAxis


fig, (sensor, regulator, Stripe) = plt.subplots(1,3, constrained_layout=True)

WT_Plotter(sensor, "Sensor", {'color':'red'})
WT_fit_plot(sensor, Sensor_est_array, {'color':'red'})
WT_Plotter(regulator, "Regulator", {'color': 'blue'})
WT_fit_plot(regulator, Regulator_est_array, {'color':'blue'})
WT_Plotter(Stripe,"Stripe", {'color': 'green'})
WT_fit_plot(Stripe, Output_est_array, {'color':'green'})
plt.suptitle("Wild type data plots")

plt.show()

# %%
