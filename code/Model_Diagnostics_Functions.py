import numpy as np
from scipy import stats
from Models import *
from Epistasis_calc_functions import *
import matplotlib.pyplot as plt

"""
This program contains a set of function used in Aibek_Model_Diagnostics.py
get_Data() calculates and returns epistasis for a chosen model
get_pvalue() conducts a two sample t test for two models and returns a p value with a conclusion for a hypothesis test
get_histogram() returns diagnostic plots 
"""
#%%
def get_Data(model):
    # observed epistasis
    if model == 'observed':
        eps_exp_df = get_Eps('observed')
        eps = eps_exp_df['Ep']
    # model hill epistasis    
    elif model == 'model_hill':
        eps_hill_df = get_Eps(model_hill.model)
        eps = eps_hill_df['Ep']
    # elif model == 'other model'
    return eps 

#%%
def get_pvalue(model1, model2):
    eps1 = get_Data(model1) # move outside function 
    eps2 = get_Data(model2)
    eps1_list = list(eps1)
    eps2_list = list(eps2)
    # observed epistasis
    eps1_test = np.array(eps1_list)
    # model epistasis
    eps2_test = np.array(eps2_list)
    # check equal variance assumption
    var_eps1 = np.var(eps1)
    var_eps2 = np.var(eps2)
    ratio_var = np.divide(var_eps1, var_eps2)
    # conduct two sample t test 
    alpha = 0.05
    # change model name
    if ratio_var < 4:
        test = stats.ttest_ind(a=eps1_test, b=eps2_test, equal_var=True)
        p_value = getattr(test, 'pvalue')
        if model2 == 'model_hill':
            model_name = 'Hill model'
        elif model2 == 'model_hill_shaky':
            model_name = 'Hill Shaky model'
        elif model2 == 'model_thermodynamic':
            model_name = 'Thermodynamic model'
        elif model2 == 'compDeg':
            model_name == 'Competition Model'
        if p_value<alpha:
            message = "p value = " +str(p_value)+  " < 0.05"
            message = message + " .We reject the null hypothesis. We have sufficient evidence to conclude that the experimental and " + model_name+ " epistasis means are not equal."
        else:
            message = "p value = " +str(p_value)+  " > 0.05"
            message = message + " .We fail to reject the null hypothesis. We do not have sufficient evidence to conclude that the experimental and " + model_name+ " epistasis means are not equal."
    return message

#%%
def get_histograms(model1, model2): # x, y
    model_name = model2
    # model_name = model.removesuffix('.model')
    if model_name == 'model_hill':
        label_name = 'Hill Model'
    elif model_name == 'model_hill_shaky':
        label_name = 'Hill Shaky Model'
    elif model_name == 'model_thermodynamic':
        label_name = 'Thermodynamic Model'
    elif model_name == 'compDeg':
        label_name = 'Competition Model'
    # separate plots
    eps1 = get_Data(model1) 
    eps2 = get_Data(model2)
    fig1, axes = plt.subplots(1, 2)
    eps1.hist(ax=axes[0])
    eps2.hist(ax=axes[1])
    axes[0].set_title("Experimental")
    axes[1].set_title("Hill Model") 
    # combined plots
    fig2, axes = plt.subplots(1,1)
    plt.hist(eps1, alpha=0.5, label='Experimental')
    plt.hist(eps2, alpha=0.5, label=label_name)
    plt.legend(loc='upper right')
    return fig1, fig2

#%% 
# variables
#eps_obs = get_Data('observed')
#eps_hill = get_Data('model_hill')
# observed, model_hill
plot1 = get_histograms('observed', 'model_hill')
# observed, model_hill_shaky
pvalue1 = get_pvalue('observed', 'model_hill')