import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

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
