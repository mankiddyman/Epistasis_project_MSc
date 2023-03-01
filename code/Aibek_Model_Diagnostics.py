from Models import *
from Epistasis_calc_functions import *
from data_wrangling import *
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy import stats

"""
a program that:
1) uses Epistasis_calc_functions to calculate mean epistasis for the models defined in Models.py
2) plots histograms of experimental and model epistasis
3) performs statistical tests 
"""
#----------------------------------------------------------------------------------------------------------
#%%
model_list = ['model_leaky', 'model_hill', 'thermodynamic_model', 'model_hill_shaky']
data = pd.read_excel('../data/SM_params.xlsx') # change to model sm 

leaky = model_leaky #1
hill = model_hill #2
td = thermodynamic_model #3
hill_shaky = model_hill_shaky #4
 
# epistasis for Hill Model
# eps_hill_df = 
# eps_ind_hill_df = 

# epistasis for Leaky Model
# eps_leaky_df =
# eps_ind_leaky_df =

# epistasis for Hill Shaky Model
# eps_shaky_df =
# eps_ind_shaky_df =

# epistasis for Thermodynamic Model
# eps_td_df =
# eps_ind_td_df = 

# test
#%%
# rss histogram per model
# experimental epistasis, model epistasis
exp_eps = pd.read_excel('../data/Source_Data.xlsx', sheet_name='Figure 2')
exp_eps.drop(exp_eps[exp_eps['Unnamed: 0'] == 'genotype'].index, inplace = True)
exp_eps.columns = exp_eps.columns.str.replace('Unnamed: 9', 'Epistasis')
ela_df = exp_eps[['Epistasis']]
ela_eps = ela_df["Epistasis"]
#ela_eps.hist()
#fig, axes = plt.subplots(1, 2)
#ela_eps.hist(ax=axes[0])
#ela_eps.hist(ax=axes[1])
#axes[0].set_title("Experimental")
#axes[1].set_title("Leaky Model")
eps_hill_df = get_Eps(model_hill) # important
eps_hill= eps_hill_df['Ep']

#%%
print('The following models are available')
for i in model_list:
    print(model_list.index(i)+1, end=')')
    print(" ", i)
model = input('Select a model: ')
# handle incorrect input
while model != '1' and model != '2' and model != '3':
    print('WRONG INPUT! PLEASE SELECT A MODEL')
    print('The following models are available')
    for i in model_list:
        print(model_list.index(i)+1, end=')')
        print(" ", i)
    model = input('Select a model: ')
# handle correct input
while model == '1' or model == '2' or model == '3':
    if model == '1':
        print('You have selected the Leaky model')
        cont = input('Do you want to select a new model? y/n')
        if cont == 'n':
            print('Residual sum of squares')
            data["rss"].hist()
            plt.title("Histogram of residual sum of squares")
            print('Diagnostic 1: epistasis comparison')
            fig, axes = plt.subplots(1, 2)
            ela_eps.hist(ax=axes[0])
            ela_eps.hist(ax=axes[1])
            axes[0].set_title("Experimental")
            axes[1].set_title("Leaky Model")
            break
    elif model == '2':
        print('You have selected the Hill model')
        cont = input('Do you want to select a new model? y/n')
        if cont == 'n':
            print('Diagnostic 1: residual sum of squares')
            print('Diagnostic 2: epistasis comparison')
            fig, axes = plt.subplots(1, 2)
            ela_eps.hist(ax=axes[0])
            eps_hill.hist(ax=axes[1])
            axes[0].set_title("Experimental")
            axes[1].set_title("Hill Model")
            fig, axes = plt.subplots(1,1)
            plt.hist(ela_eps, alpha=0.5, label='Experimental')
            plt.hist(eps_hill, alpha=0.5, label='Hill Model')
            plt.legend(loc='upper right')
            plt.show()
    elif model == '3':
        print('You have selected the Thermodynamic model')
        cont = input('Do you want to select a new model? y/n')
        if cont == 'n':
            print('Residual sum of squares')
            data["rss"].hist()
            plt.title("Histogram of residual sum of squares")
            print('Diagnostic 1: epistasis comparison')
            # fig, axes = plt.subplots(1, 2)
            # ela_eps.hist(ax=axes[0])
            # ela_eps.hist(ax=axes[1])
            # axes[0].set_title("Experimental")
            # axes[1].set_title("Thermodynamic Model")
            # fig, axes = plt.subplots(1,1)
            # plt.hist(ela_eps, alpha=0.5, label='Experimental')
            # plt.hist(eps_td, alpha=0.5, label='Thermodynamic Model')
            # plt.legend(loc='upper right')
            # plt.show()
    model = input('Select a model: ')
    if model != '1' and model != '2' and model != '3':
        print('WRONG INPUT! PLEASE SELECT A MODEL')
        print('The following models are available')
        for i in model_list:
            print(model_list.index(i)+1, end=')')
            print(" ", i)
        model = input('Select a model: ')