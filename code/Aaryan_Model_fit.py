#%%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize
from Models import *
from data_wrangling import *
from itertools import chain
from itertools import repeat
from PyPDF2 import PdfMerger
import inspect
from Model_fitting_functions import *
 #change dictionary to list of coeffs
def coef_dict_to_list(coef_dict):
    return list(coef_dict.values())
#the single mutant to be studied
   # function to input paramaters as a list

def get_data_SM(mutation:str):
        df_MT = df
        data = meta_dict["SM"]
        data = data.loc[data['Mutant_ID'] == mutation]
        #WT data missing measurements for inducer = 0.2 so drop last columnn
        data = data[:-1]

        a = re.search("[0-9]",mutation).start()
        mut_col = f"{mutation[:a]}"
        mutation_means = f"{mutation[:a]}_mean"
        df_MT[[mut_col, "Stripe"]] = data[[mutation_means, "Stripe_mean"]].values
        df_MT[["Mutant_ID"]]=mutation
        if mutation.startswith("Output"):
            df_MT["Sensor"]=meta_dict["WT"].Sensor

        return df_MT
#data=pd.read_csv('../data/WT_single.csv')

    #now gonna do the scipy.optimize.minimize
def WT_fit_plot(ax, y,params,label:str):
        return ax.plot(I_conc, y, **params,label=f"{label}")


    #define scatter plotting function with log scales
def WT_Plotter(ax,y, params):
    out = ax.scatter(I_conc, data_[y], **params, marker = 'o')
    xScale = ax.set_xscale('log')
    yScale = ax.set_yscale('log')

    return out, xScale, yScale 
def dict_to_list(params_dict,return_keys=False):
    if return_keys==True:
       a=[list(i.keys()) for i in list(params_dict.values())]
    elif return_keys==False:
        a=[list(i.values()) for i in list(params_dict.values())]

    return list(chain.from_iterable(a))
        

all_files=os.listdir(os.path.join(path_main,"mutants_separated"))
dat_files=[x for x in all_files if '.dat' in x]

SM_names=[x[:-4] for x in dat_files]
#SM_names=SM_names[0:10] #remove later

#paramater estimates from WT:
#WT_params = [6.59963463e+02, 1.63471394e+04, 1.25925588e+03, 1.16043953e+00, 1.99831036e+03, 2.04000859e+11, 2.77180774e+06, 8.37522575e-01, 5.47787795e-06, 6.71081447e+04, 1.41294256e-03, 5.41433758e+07,2.12643877e+00, 2.72060461e+00, 1.25044310e+00]

#WT_params = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]   


model_type_in=model_hill
params_dict={"sen_params":{"As":618.047086,"Bs":16278.856600,"Cs":1300.653790,"Ns":1.096541},"reg_params":{"Ar":1916.175610,"Br":18874.240800,"Cr":0.009030,"Nr":0.820433},"out_h_params":{"Ah":683.835638,"Bh":32464.380200,"Ch":0.000473},"out_params":{"Fo":2.821352,"Ao":0.632148,"Bo":0.972768 ,"Co":2.640174,"No":1.919339}}
SM_mutant_of_interest="Regulator1"
n_iter = float=1e5
#%%
#thermodynamic
params_dict={"sen_params":{"a_s":1,"Kp_s":1,"P":1,"Ki_s":1,"C_pi_s":1},"reg_params":{"a _r":1,"Kp_r":1,"Ks_r":1},"out_params":{"a_o":1,"Kp_o":1,"K_lacI_o":1,"C_po_lacI_o":1}}

#model_hill_shaky testing
################################
params_dict={"sen_params":{"As":1,"Bs":1,"Cs":1,"Ns":1},"reg_params":{"Ar":1,"Br":1,"Cr":1,"Nr":1},"out_h_params":{"Ah":1,"Bh":1,"Ch":1},"out_params":{"Fo":1,"Ao":1,"Bo":1,"Co":1,"No":1}}

WT_params=param_dictToList(params_dict)
#WT_params=list(np.array(WT_params)+100)
#WT_params=[1]*13
#WT_params=list(np.array(WT_params)+100)
bnds = tuple(repeat((0,None), len(WT_params)))
model_type = model_hill_shaky
n_iter = float=1e5

WT_params=dict_to_list(params_dict)
WT_params=list(np.array(WT_params)+100)
WT_params=[1]*13
WT_params=list(np.array(WT_params)+100)
bnds = ((0, None), (0, None), (0, None), (0, None), (0, None), (0, None), (0, None),(0, None), (0, None), (0, None),(1, None), (0,1), (0,1))

#fitting

#thermodynamic model
#rewrite this
params_dict={"sen_params":{"a_s":1,"Kp_s":1,"P":1,"Ki_s":1,"C_pi_s":1},"reg_params":{"a_r":1,"Kp_r":1,"Cr":1,"Nr":1},"out_h_params":{"Ah":1,"Bh":1,"Ch":1}}


params_dict={"sen_params":{"As":1,"Bs":1,"Cs":1,"Ns":1},"reg_params":{"Ar":1,"Br":1,"Cr":1,"Nr":1},"out_h_params":{"Ah":1,"Bh":1,"Ch":1},"free_params":{"Fo":1},"out_params":{"Ao":1,"Bo":1,"Co":1,"No":1}}




#now generating initial guesses

model_type_in=thermodynamic_model
#%%

#functionality to add ,
#  automatically find out number of arguements, automatically set bounds on 


#bound generation for model_hill
model_type=model_hill
params_dict={"sen_params":{"As":1,"Bs":1,"Cs":1,"Ns":1},"reg_params":{"Ar":1,"Br":1,"Cr":1,"Nr":1},"out_h_params":{"Ah":1,"Bh":1,"Ch":1},"out_params":{"Ao":1,"Bo":1,"Co":1,"No":1},"free_params":{"Fo":1},}
model_list=dict_to_list(params_dict,False)
model_dict_keys=dict_to_list(params_dict,True)
custom_settings=[[2,2,2],[4,4,4],["Ao","Bo","Co"]]
generate_bounds(params_dict,node="Regulator",custom_settings=custom_settings)