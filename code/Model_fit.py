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
#best guess from first fit = [311,17220,1081,0.9,2000,6000,0.0421,0.308,2000,4054,5.25E-5,3.17]
#%%
coeffs = {'A_s': 1,'B_s': 1,'C_s':1,'N_s':0, 'A_r':1,'B_r':1,'C_r':1,'N_r':0,'A_h':1, 'B_h':1, 'C_h':1, 'A_o':1,'B_o':1,'C_o': 1,'N_o': 0}
import random
start_coeffs = np.random.rand(15) * (1e4-1e-4) + 0.5
print(start_coeffs)
# function to input paramaters as a list
def coeff_input(param_list):
    if len(param_list) != len(coeffs):
        print('input paramaters list is not the same length as coefficents: COEFFIECIENTS UNCHANGED')
    else:
        for i,key in enumerate(list(coeffs.keys())):
            coeffs[str(key)] = param_list[i]
    return coeffs
#%%

def min_fun(coefficients,data=data):
    log_sen = np.log10(data.Sensor)
    log_reg = np.log10(data.Regulator)
    log_out = np.log10(data.Output)
    log_stripe = np.log10(data.Stripe)   
    ind = data.S
    
    Sensor_est, Regulator_est, Output_est, Stripe_est = model_4_pred(ind,*coefficients)
    log_sen_est = np.log10(Sensor_est)
    log_reg_est = np.log10(Regulator_est)
    log_out_est = np.log10(Output_est)
    log_stripe_est = np.log10(Stripe_est)

    result = np.power((log_sen - log_sen_est), 2)
    result += np.power((log_reg - log_reg_est), 2)
    result += np.power((log_out - log_out_est), 2)
    result += np.power((log_stripe - log_stripe_est), 2)
    return np.sum(result)
#now gonna do the scipy.optimize.minimize
#%%
#start_coeffs =stats.uniform(0.001, 100).rvs(15)
#print(start_coeffs)
# start_coeffs=[311,17220,1081,0.9 #As,Bs,Cs,Ns
# , 6.80803615e+02,1.88728919e+04,3.91470236e-02,3.74001669e-01, #Ar,Br,Cr,Nr 
# +1.02584227e+03,9.47437625e+04,3.96678625e-03, # Ah, Bh, Ch
# 1,8.58578541e-01,1.05191555e+01,7.15818491e+01,8.30769275e-01] #Fo,Ao,Bo,Co,No
start_coeffs=[1]*16
start_coeffs=[7.31790329e+02, 1.63890722e+04, 1.21558902e+03, 1.23445789e+00,2.19715875e+03, 5.09396044e+04, 8.76632583e-03, 1.37272209e+00, 3.97046404e+03, 2.81037530e+04, 5.99908397e-04, 8.61568305e-01, 7.03425130e-01, 7.57153375e+00, 1.25692066e+00, 3.39280741e+00]
print("starting parameters")
print(start_coeffs)
# adding bounds such that all estimated parameters are +ve
bnds = ((0, None), (0, None), (0, None), (0, None), (0, None), (0, None), (0, None),(0,None), (0, None), (0, None), (0, None), (0, None), (0, None), (0, None), (0, None), (0, None))
min_result=minimize(min_fun,x0= start_coeffs ,method='Nelder-Mead',bounds=bnds,options={"maxiter":1e5,"disp":True})

#plotting the predictions now

#generating estimates
Sensor_est_array,Regulator_est_array,Output_est_array, Stripe_est_array = model_4_pred(data.S,*min_result.x)

Sensor_est_array_initial,Regulator_est_array_initial,Output_est_array_initial, Stripe_est_array_initial = model_4_pred(data.S,*start_coeffs)


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
    return out, xScale, yScale, 


fig, ((Sensor, Regulator), (Output, Stripe)) = plt.subplots(2,2, constrained_layout=True)

WT_Plotter(Sensor, "Sensor", {'color':'red'})
WT_fit_plot(Sensor, Sensor_est_array, {'color':'red'})
WT_fit_plot(Sensor, Sensor_est_array_initial, {'color':'black'})

WT_Plotter(Regulator, "Regulator", {'color': 'blue'})
WT_fit_plot(Regulator, Regulator_est_array, {'color':'blue'})
WT_fit_plot(Regulator, Regulator_est_array_initial, {'color':'black'})

WT_Plotter(Output, "Output", {'color': 'purple'})
WT_fit_plot(Output, Output_est_array, {'color':'purple'})
WT_fit_plot(Output, Output_est_array_initial, {'color':'black'})


WT_Plotter(Stripe,"Stripe", {'color': 'green'})
WT_fit_plot(Stripe, Stripe_est_array, {'color':'green'})
WT_fit_plot(Stripe, Stripe_est_array_initial, {'color':'black'})


plt.suptitle("Wild type data plots")

plt.show()
print(min_result)
print("converged parameters \n",min_result.x)
# %%
