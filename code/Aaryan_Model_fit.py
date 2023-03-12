#%%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize
from Models import *
from Models import model_hill
from data_wrangling import *
from itertools import chain
from itertools import repeat
from PyPDF2 import PdfMerger
import inspect
from Model_fitting_functions import *
 #change dictionary to list of coeffs
#%%

#steady hill model_1
params_dict_hill={"sen_params":{"A_s":2,"B_s":1,"C_s":1,"N_s":1},"reg_params":{"A_r":1,"B_r":1,"C_r":1,"N_r":1},"out_h_params":{"A_h":1,"B_h":1,"C_h":1},"out_params":{"A_o":1,"B_o":1,"C_o":1,"N_o":1},"free_params":{"F_o":1}}

params_list_hill=dict_to_list(params_dict_hill)

I_conc=meta_dict["WT"].S
hill=model_hill(params_list_hill,I_conc)
func=model_hill.model

# params_list_hill=get_WT_params(model_type=func,start_guess=params_list_hill,n_iter=1e10,method="Nelder-Mead",params_dict=params_dict_hill,custom_settings=[],tol=.01)
# print(params_list_hill)
# params_list_hill=get_WT_params(model_type=func,start_guess=params_list_hill,n_iter=1e10,method="Powell",params_dict=params_dict_hill)
# params_list_hill=get_WT_params(model_type=func,start_guess=params_list_hill,n_iter=1e10,method="Powell",params_dict=params_dict_hill)
# print(params_list_hill)

# params_list_hill=get_WT_params(model_type=func,start_guess=params_list_hill,n_iter=1e10,method="Nelder-Mead",params_dict=params_dict_hill)
# #first 8 parameters are good fits, will reset last 8 to 1
# for i in range(0,len(params_list_hill)):
#     if i>=8:
#         params_list_hill[i]=1
# print(params_list_hill)
# params_list_hill=get_WT_params(model_type=func,start_guess=params_list_hill,n_iter=1e10,method="Powell",params_dict=params_dict_hill)


params_list_hill=[618.05, 16278.86, 1300.65, 1.23445789e+00,2.19715875e+03, 5.09396044e+04, 8.76632583e-03, 1.37272209e+00, 3.97046404e+03, 2.81037530e+04, 5.99908397e-04, 8.61568305e-01, 7.03425130e-01, 7.57153375e+00, 1.25692066e+00, 3.39280741e+00]
converged_params_list_hill=get_WT_params(model_type=func,start_guess=params_list_hill,n_iter=1e5,method="Nelder-Mead",params_dict=params_dict_hill,custom_settings=[],tol=0.0001)
converged_params_dict_hill=list_to_dict(old_dict=params_dict_hill,new_values=converged_params_list_hill)

#first get wt params

#then apply bounds to estimate for all single mutants, exporting converged params to excel

model_fitting_SM(model_type=func,n_iter=1e5,params_dict=converged_params_dict_hill)
#%%
#steady_hill model 2 with changed fluorescence and output equations
hill_2=model_hill([1]*13,meta_dict["WT"].S)
params_dict_hill_2=hill_2.example_dict_model_2
params_list_hill_2=dict_to_list(params_dict_hill_2)
func=model_hill.model_2
converged_params_list_hill_2=get_WT_params(model_type=func,start_guess=params_list_hill_2,n_iter=1e5,method="Nelder-Mead",params_dict=params_dict_hill_2,custom_settings=[],tol=10)
converged_params_list_hill_2=get_WT_params(model_type=func,start_guess=converged_params_list_hill_2,n_iter=1e5,method="TNC",params_dict=params_dict_hill_2,custom_settings=[],tol=1)
converged_params_list_hill_2=[618.05, 16278.86, 1300.65, 1.23445789e+00,2.19715875e+03, 5.09396044e+04, 8.76632583e-03, 1.37272209e+00, 3.97046404e+03, 2.81037530e+04, 7.57153375e+00, 1.25692066e+00, 3.39280741e+00]
converged_params_list_hill_2=get_WT_params(model_type=func,start_guess=converged_params_list_hill_2,n_iter=1e5,method="Nelder-Mead",params_dict=params_dict_hill_2,custom_settings=[],tol=1)
converged_params_list_hill_2=get_WT_params(model_type=func,start_guess=converged_params_list_hill_2,n_iter=1e5,method="Nelder-Mead",params_dict=params_dict_hill_2,custom_settings=[],tol=1)
converged_params_dict_hill_2=list_to_dict(old_dict=params_dict_hill_2,new_values=converged_params_list_hill_2)
model_fitting_SM(model_type=func,n_iter=1e5,params_dict=converged_params_dict_hill_2)
#%%

#thermodynamic model

therm_model=model_thermodynamic(params_list=[1]*14,I_conc=meta_dict["WT"].S)
params_therm_dict=therm_model.example_dict
params_therm_list=dict_to_list(params_therm_dict)
func=therm_model.model
params_therm_dict={"sen_params":{"P_b":1,"P_u":1,"K_12":1,"C_pa":1,"A_s":1},"reg_params":{"P_r":1,"C_pt":1,"K_t":1,"A_r":1},"out_h_params":{},"out_params":{"P_o":1,"C_pl":1, "K_l":1,"A_o":1},"free_params":{},"fixed_params":{"F_o":1}}
params_therm_list=dict_to_list(params_therm_dict)
params_therm_list=get_WT_params(model_type=func,start_guess=params_therm_list,n_iter=1e5,method="Nelder-Mead",params_dict=params_therm_dict,custom_settings=[[1,0,0,0,0],[None,1,1,None,None],["C_pa","C_pt","C_pl","P_p","F_o"]],tol=1)
params_therm_list=get_WT_params(model_type=func,start_guess=params_therm_list,n_iter=1e5,method="Nelder-Mead",params_dict=params_therm_dict,custom_settings=[[1,0,0,0,0],[None,1,1,None,None],["C_pa","C_pt","C_pl","P_p","F_o"]],tol=1)
params_therm_list=get_WT_params(model_type=func,start_guess=params_therm_list,n_iter=1e5,method="TNC",params_dict=params_therm_dict,custom_settings=[[1,0,0,0,0],[None,1,1,None,None],["C_pa","C_pt","C_pl","P_p","F_o"]],tol=1)
#I like the sensor so i will reinitialise parameters of sensor
a=params_therm_list[0:3]
b=params_therm_list[10:12]
params_therm_list=dict_to_list(params_therm_dict)
params_therm_list[0:3]=a
params_therm_list[10:12]=b
params_therm_dict=list_to_dict(old_dict=params_therm_dict,new_values=params_therm_list)
params_therm_list=get_WT_params(model_type=func,start_guess=params_therm_list,n_iter=1e5,method="Nelder-Mead",params_dict=params_therm_dict,custom_settings=[[1,0,0,params_therm_dict['fixed_params']['F_o'],params_therm_dict['fixed_params']['P_p']],[None,1,1,params_therm_dict['fixed_params']['F_o'],params_therm_dict['fixed_params']['P_p']],["C_pa","C_pt","C_pl","P_p","F_o"]],node="Regulator",tol=1)
params_therm_list=[1e-5]*12
params_therm_list=get_WT_params(model_type=func,start_guess=params_therm_list,n_iter=1e5,method="Nelder-Mead",params_dict=params_therm_dict,custom_settings=[[1,0,0],[None,1,1],["C_pa","C_pt","C_pl"]],tol=0.0011)

#competition model
#%%
model_comp=CompDeg(params_list=[1]*15,I_conc=meta_dict["WT"].S)
func=model_comp.model_2
params_dict_comp=model_comp.example_dict_model_2
params_list_comp=dict_to_list(params_dict_comp)
params_list_comp=[4.20504769e+02, 1.47356534e+04, 1.67916264e+03, 1.27219156e+00, 23, 2.36888092e+04, 1.03576041e-02 ,9.10072254e-01, 1.67389421e+00, 8.95345e+01, 2.65769862e+00 ,1.37995266e+00 ,23612.8375e+00, 5, 100]


params_list_comp=get_WT_params(model_type=func,start_guess=params_list_comp,n_iter=1e5,method="Nelder-Mead",params_dict=params_dict_comp,custom_settings=[[1e-5,1e-5,1e-5],[None,None,None],["Deg","A_s","C_o"]],tol=1)

#%%
#%%
#getting wt params
func=model_hill_shaky.model
a=model_hill_shaky([1]*16,meta_dict["WT"].S)
I_conc=meta_dict["WT"].S
params_list_hill_shaky=dict_to_list(a.example_dict)
params_dict_hill_shaky=a.example_dict
params_list_hill_shaky=get_WT_params(model_type=func,start_guess=params_list_hill_shaky,n_iter=1e5,method="TNC",params_dict=params_dict_hill_shaky,custom_settings=[],tol=0.001)
#applying bounds to fit for sm mutants
        

all_files=os.listdir(os.path.join(path_main,"mutants_separated"))
dat_files=[x for x in all_files if '.dat' in x]

SM_names=[x[:-4] for x in dat_files]
#SM_names=SM_names[0:10] #remove later

#paramater estimates from WT:
#WT_params = [6.59963463e+02, 1.63471394e+04, 1.25925588e+03, 1.16043953e+00, 1.99831036e+03,    2.04000859e+11, 2.77180774e+06, 8.37522575e-01, 5.47787795e-06, 6.71081447e+04, 1.41294256e-03, 5.41433758e+07,2.12643877e+00, 2.72060461e+00, 1.25044310e+00]

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