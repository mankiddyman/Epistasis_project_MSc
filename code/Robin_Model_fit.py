import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize
from Models import *
from data_wrangling import *

#use model chosen from Models.py to obtain fit of model to data

df = pd.read_csv('../data/WT_single.csv')

#substitute mutation data into the WT dataframe for a particular set of mutantions
#for single mutations to start:
def get_data_SM(mutation):
    df_MT = df
    data = meta_dict["SM"]
    data = data.loc[data['Mutant_ID'] == mutation]
    #WT data missing measurements for inducer = 0.2 so drop last columnn
    data = data[:-1]

    a = re.search("[0-9]",mutation).start()
    mut_col = f"{mutant[:a]}"
    mutation_means = f"{mutant[:a]}_mean"
    df_MT[[mut_col, "Stripe"]] = data[[mutation_means, "Stripe_mean"]].values
    return df_MT

#example of above function to get SM dat for Sensor1 mutation:
S1_df = get_data_SM("Sensor1")

# minimizing sum of sqrd deviations of log model and log data
data = S1_df

#paramater estimates from WT:
WT_params = [6.08397103e+02, 1.52504577e+04, 1.66805905e+03, 1.19893355e+00, 6.87964693e+02, 2.34976114e+04, 6.23671728e-02, 3.91730917e-01, 5.90606548e+02, 3.52871257e+04, 5.29890033e-04, 8.29829791e-01, 4.28817019e+00, 3.13322189e+00, 1.80901848e+00]
#%%
coeffs = {'A_s': 1,'B_s': 1,'C_s':1,'N_s':0, 'A_r':1,'B_r':1,'C_r':1,'N_r':0,'A_h':1, 'B_h':1, 'C_h':1, 'A_o':1,'B_o':1,'C_o': 1,'N_o': 0}
# function to input paramaters as a list
def coeff_input(param_list):
    if len(param_list) != len(coeffs):
        print('input paramaters list is not the same length as coefficents: COEFFIECIENTS UNCHANGED')
    else:
        for i,key in enumerate(list(coeffs.keys())):
            coeffs[str(key)] = param_list[i]
    return coeffs
#%%
#change dictionary to list of coeffs
def coef_dict_to_list(coef_dict):
    return list(coef_dict.values())

init_coeffs = coeff_input(WT_params)

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
print("starting paramaters:", WT_params)
bnds = ((0, None), (0, None), (0, None), (0, None), (0, None), (0, None), (0, None), (0, None), (0, None), (0, None), (0, None), (0, None), (0, None), (0, None), (0, None))
min_result=minimize(min_fun,x0= WT_params ,method='Nelder-Mead',bounds=bnds)
#plotting the predictions now

#generating estimates
Sensor_est_array,Regulator_est_array,Output_est_array, Stripe_est_array = model_4_pred(data.S,*min_result.x)

Sensor_est_array_initial,Regulator_est_array_initial,Output_est_array_initial, Stripe_est_array_initial = model_4_pred(data.S,*WT_params)


Signal=data.S
def WT_fit_plot(ax, y, params):
    return ax.plot(Signal, y, **params)


#define scatter plotting function with log scales
def WT_Plotter(ax,y, params):
    out = ax.scatter(Signal, data[y], **params, marker = 'o')
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

title = ["SM data type data plots for mutation", "Sensor 1"]
plt.suptitle(title)

plt.show()
print("final parameter estimates:", min_result.x)
# %%
