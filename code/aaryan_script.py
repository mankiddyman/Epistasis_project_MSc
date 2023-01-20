import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

#beginning our work by recreating the non linear least squares fitting of wild-type circuits

#S collumn is signal 
#Sensor collumn is the concentration of sensor in the network with topology: inducer -> sensor (GFP output)
#Regulator collumn is the concentration of regulator in the following network:inducer ->S -|R (GFP output)
#Output is conc of output in :
#inducer -> S -| Output (GFP)
#stripe is the full network measuring the conc of output
wt_df = pd.read_csv('../data/WT_single.csv')
wt_Signal=wt_df.S
wt_Sensor=wt_df.Sensor
wt_Regulator=wt_df.Regulator
wt_Output=wt_df.Output
wt_Stripe=wt_df.Stripe

#%%
fig, axs= plt.subplots(2,2, figsize=(8,8),constrained_layout=True)

axs[0,0].set_xscale("log")
axs[0,0].set_yscale("log")
axs[0,1].set_xscale("log")
axs[0,1].set_yscale("log")
axs[1,0].set_xscale("log")
axs[1,0].set_yscale("log")
axs[1,1].set_xscale("log")
axs[1,1].set_yscale("log")
axs[0,0].set_xlabel("Signal Concentration")
axs[0,1].set_xlabel("Signal Concentration")
axs[1,0].set_xlabel("Signal Concentration")
axs[1,1].set_xlabel("Signal Concentration")

axs[0,0].plot(wt_Signal,wt_Sensor, linestyle="",marker="o")
axs[0,0].set_title(r'inducer -> sensor (GFP output)')
axs[0,0].set_ylabel("Sensor concentration")

axs[0,1].plot(wt_Signal,wt_Regulator, linestyle="",marker="o")
axs[0,1].set_title(r'inducer ->S -|R (GFP output)')
axs[0,1].set_ylabel("Regulator concentration")

axs[1,0].plot(wt_Signal,wt_Output, linestyle="",marker="o")
axs[1,0].set_title(r'inducer -> S -| Output (GFP)')
axs[1,0].set_ylabel("Output concentration")

axs[1,1].plot(wt_Signal,wt_Stripe, linestyle="",marker="o")
axs[1,1].set_title(r'Full circuit with stripe')
axs[1,1].set_ylabel("Output concentration")

#plots make sense
#should investigate the relationships between inducer and sensor (calculated using regression of inducer )
#function to do subsequent non linear regression is :
#%%
#we will now create model functions for each of the steady state values of Sensor, Regulator and Output.
#Henceforth referred to as S, R and O

 
def S_func_full_network(I_conc,A_s,B_s,C_s,N_s):
    val=(A_s+B_s*(C_s*I_conc)**N_s)/(1+(C_s*I_conc)**N_s)
    return val
def R_func_full_network(S_conc,A_r,C_r,N_r):
    return ((A_r)/(1+(C_r*S_conc)**N_r))

def O_func_full_network(S_and_R_conc,A_o,C_o,N_o):
    return ((A_o)/(1+(C_o*(S_and_R_conc))**N_o))

def O_func_half_network(S_conc,A_o,C_o,N_o):
    return ((A_o)/(1+(C_o*S_conc)**N_o))
#for full network
popt_S, pcov = curve_fit(S_func_full_network, wt_Signal,wt_Sensor)

S_conc_estimate=S_func_full_network(np.array(wt_Signal),*popt_S)

popt_R, pcov = curve_fit(R_func_full_network, S_conc_estimate,wt_Regulator,p0=np.array([200,20,1.5]))

R_conc_estimate=R_func_full_network(np.array(wt_Sensor),*popt_R)
#for half network of output conc

popt_O, pcov= curve_fit(O_func_half_network,S_conc_estimate,wt_Output,p0=np.array([20,0.1,1]))

O_conc_half_network_estimate=O_func_half_network(np.array(wt_Sensor),*popt_O)

sum_conc_estimate=S_conc_estimate+R_conc_estimate
#%%
popt_O_full, pcov= curve_fit(O_func_full_network,sum_conc_estimate,wt_Stripe,p0=np.array([60,2,1.5]))

sum_conc_stripe=wt_Sensor+R_conc_estimate

O_conc_full_network_estimate=O_func_full_network(np.array(sum_conc_estimate),*popt_O_full)


fig, axs= plt.subplots(2,2, figsize=(8,8),constrained_layout=True)

axs[0,0].set_xscale("log")
axs[0,0].set_yscale("log")
axs[0,1].set_xscale("log")
axs[0,1].set_yscale("log")
axs[1,0].set_xscale("log")
axs[1,0].set_yscale("log")
axs[1,1].set_xscale("log")
axs[1,1].set_yscale("log")
axs[0,0].set_xlabel("Signal Concentration")
axs[0,1].set_xlabel("Signal Concentration")
axs[1,0].set_xlabel("Signal Concentration")
axs[1,1].set_xlabel("Signal Concentration")


#we can use the full network eqn for this because sensor only responds to inducer
popt, pcov = curve_fit(S_func_full_network, wt_Signal,wt_Sensor)
axs[0,0].plot(wt_Signal,wt_Sensor, linestyle="",marker="o")
axs[0,0].plot(wt_Signal,S_func_full_network(np.array(wt_Signal),*popt))
axs[0,0].set_title(r'inducer -> sensor (GFP output)')
axs[0,0].set_ylabel("Sensor concentration")

axs[0,1].plot(wt_Signal,wt_Regulator, linestyle="",marker="o")
axs[0,1].plot(wt_Signal,R_conc_estimate)
axs[0,1].set_title(r'inducer ->S -|R (GFP output)')
axs[0,1].set_ylabel("Regulator concentration")

axs[1,0].plot(wt_Signal,wt_Output, linestyle="",marker="o")
axs[1,0].plot(wt_Signal,O_conc_half_network_estimate)
axs[1,0].set_title(r'inducer -> S -| Output (GFP)')
axs[1,0].set_ylabel("Output concentration")

#output funciton requires the concentrations of regulator and sensor
axs[1,1].plot(wt_Signal,wt_Stripe, linestyle="",marker="o")
axs[1,1].plot(wt_Signal,O_conc_full_network_estimate)
axs[1,1].set_title(r'Full circuit with stripe')
axs[1,1].set_ylabel("Output concentration")

# %%
