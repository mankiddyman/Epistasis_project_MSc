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
#td = thermodynamic_model #3
hill_shaky = model_hill_shaky #4

rss_data = pd.read_excel('../data/model_hill_SM_params.xlsx')
rss_df = rss_data[["rss","time_elapsed_s"]]
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
#ela_eps.hist()
#fig, axes = plt.subplots(1, 2)
#ela_eps.hist(ax=axes[0])
#ela_eps.hist(ax=axes[1])
#axes[0].set_title("Experimental")
#axes[1].set_title("Leaky Model")
# changed because class 
eps_hill_df = get_Eps(model_hill.model) # important
eps_hill= eps_hill_df['Ep']

eps_exp_df = get_Eps('observed')
eps_exp = eps_exp_df['Ep']

eps_exp_list = list(eps_exp)
eps_exp_test = np.array(eps_exp_list)
eps_hill_list = list(eps_hill)
eps_hill_test = np.array(eps_hill_list)
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
            #print('Residual sum of squares')
            #data["rss"].hist()
            #plt.title("Histogram of residual sum of squares")
            #print('Diagnostic 1: epistasis comparison')
            #fig, axes = plt.subplots(1, 2)
            #ela_eps.hist(ax=axes[0])
            #ela_eps.hist(ax=axes[1])
            #axes[0].set_title("Experimental")
            #axes[1].set_title("Leaky Model")
            break
    elif model == '2':
        print('You have selected the Hill model')
        cont = input('Do you want to select a new model? y/n')
        if cont == 'n':
            print('Diagnostic 1: residual sum of squares')
            print('Diagnostic 2: epistasis comparison')
            fig, axes = plt.subplots(1, 2)
            eps_exp.hist(ax=axes[0])
            eps_hill.hist(ax=axes[1])
            axes[0].set_title("Experimental")
            axes[1].set_title("Hill Model")
            fig, axes = plt.subplots(1,1)
            plt.hist(eps_exp, alpha=0.5, label='Experimental')
            plt.hist(eps_hill, alpha=0.5, label='Hill Model')
            plt.legend(loc='upper right')
            #fig, axes = plt.subplots(1,1) 
            #rss_df.plot(x="time_elapsed_s",y="rss")
            plt.show()
            var_ela = np.var(eps_exp)
            var_hill = np.var(eps_hill)
            ratio_var = np.divide(var_ela, var_hill)
            # 1.42<4
            alpha = 0.05
            test = stats.ttest_ind(a=eps_exp_test, b=eps_hill_test, equal_var=True)
            p_value = getattr(test, 'pvalue')
            if (p_value<alpha):
                print("p value = " +str(p_value)+  " < 0.05")
                print("We reject null hypothesis and conclude that the experimental and model epistasis means are not equal")
            else:
                print("p value " +str(p_value)+ "> 0.05")
                print("we fail to reject null hypothesis and conclude that the experimental and model epistasis means are equal")
            break
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