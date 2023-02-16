import math
import pandas as pd
import numpy as np
from data_wrangling import meta_dict
from Models import *
import re
from itertools import permutations
import matplotlib.pyplot as plt
from scipy import stats

#Calculate epistasis for all double and triple mutants from observed and expected flourcsence at low, medium and high inducer concs

df_S = meta_dict['SM']
df_DM = meta_dict['DM']
df_TM = meta_dict['TM']
df_M = (pd.concat([df_DM, df_TM], ignore_index = True)).drop_duplicates()

df_fits = pd.read_excel('../data/SM_params.xlsx').rename(columns={'Unnamed: 0': 'mutant'})

I_conc = {'low':0.0, 'medium': 0.0002, 'high':0.2}
I_conc_np = np.array(list(I_conc.values()))
I_conc_np[1] = 0.000195 #medium Inducer concentration is rounded in df_S, want more accurate value (0.000195) when fitting using a model
WT_params = df_fits.iloc[0].values.flatten().tolist()[1:-1]

#stripe output at low, medium and high incucer concentration for WT
g_WT = np.array(df_DM['obs_fluo_mean'][df_DM['genotype']=='WT'])
g_WT_sd = np.array(df_DM['obs_SD'][df_DM['genotype']=='WT'])

#calculates the GFP for a single mutant relative to WT
def g(mutant:str):
    g_single = np.empty(len(I_conc))
    g_single_std = np.empty(len(I_conc))
    for i, conc in enumerate(list(I_conc.values())):
        g_single[i] = df_S['Stripe_mean'][(df_S['Inducer'] == conc) & (df_S['Mutant_ID'] == mutant)]
        g_single_std[i] = df_S['Stripe_stdev'][(df_S['Inducer'] == conc) & (df_S['Mutant_ID'] == mutant)]
    g_single = np.divide(g_single, g_WT)
    return g_single , g_single_std

#expected flourecence for a set set of parameters and a chosen model at concentrations in I_conc dictionary relative to WT
def g_hat(mutant:str, model):
    params = df_fits[df_fits['mutant']== mutant].values.flatten().tolist()[1:-1]
    g_hat = np.divide(model(I_conc_np, *params)[-1],g_WT)
    return g_hat

#expected log fold flourecence for mutants assuming log additivity of single mutants
def G_hat_logadd(mutants:list):
    G_hat_la = g(mutants[0])[0]
    G_hat_la_std = np.power(g(mutants[0])[1], 2)
    for mut in mutants[1:]:
        G_hat_la = np.multiply(g(mut)[0],G_hat_la)
        G_hat_la_std += np.power(g(mut)[1], 2)
    G_hat_la = np.log10(G_hat_la)
    x = np.power(G_hat_la_std, 1/2)
    G_hat_la_std = np.multiply(G_hat_la, x)
    return G_hat_la, G_hat_la_std

#log additive expected GFP for double or triple mutants for a given model
def G_hat(model, mutants:list):
    #copy df_fits and replace relevant parameters with mutated ones in df_fits row for first mutant in "mutants" list
    df_fits1 = df_fits.copy(deep = False)
    for mut in mutants[1:]:
        param_id = f"_{mut[0].lower()}"
        columns_OI = df_fits1.filter(like=param_id).columns #parameters of interest to be replaced with mutant parameters
        df_fits1.loc[df_fits1['mutant'] == mutants[0],columns_OI] = list(df_fits1[df_fits1['mutant'] == mut].filter(like=param_id).values)
    #get parameter values for mutant combination
    G_hat = np.log10(g_hat(mutants[0], model))
    return G_hat

#get list of possible mutant id names - e.g ['Sensor1','Regulator1'] --> ['S1_R1' , 'R1_S1']
def get_mut_ids(mutants:list):
    muts = []
    for mut in mutants:
        mut_index = re.search("[0-9]",mut).start()
        mutant1_id = f"{mut[:1]}{mut[mut_index:]}"
        muts.extend([mutant1_id])
    #permutations of mutants to search against in df_M
    mut_perms = list(permutations(muts))
    for i,ids in enumerate(mut_perms):
        string = ''
        for j in range(len(ids)):
            string += ids[j] +'_'
        mut_perms[i] = string[:-1]
    return mut_perms

#the opposite of get_mut_ids that turns id into list of names e.g. S1_O1 --> ['Sensor1','Regulator1']
def get_mut_names(mut_id:str):
    mutant_names = []
    ids = mut_id.split('_')
    for mut in ids:
        if mut.startswith('S'):
            mut = f"Sensor{mut[1:]}"
            mutant_names.extend([mut])
        if mut.startswith('R'):
            mut = f"Regulator{mut[1:]}"
            mutant_names.extend([mut])
        if mut.startswith('O'):
            mut = f"Output{mut[1:]}"
            mutant_names.extend([mut])
    return mutant_names

#observed GFP for double or triple mutants
#takes list of possible combinations of mutation ids as parameter e.g. ['O1_R1', 'R1_O1']
def G_obs(mutants:list):
    mut_perms = get_mut_ids(mutants)
    df_mutant = df_M.loc[df_M['genotype'] == mut_perms[0]]
    for mut in mut_perms[1:]:
        dfa = df_M.loc[df_M['genotype'] == mut]
        df_mutant = pd.concat([df_mutant, dfa])    
    #get mean and std for low, med, high inducer concs
    #NB - assumes low inducer concs come above medium, which come before high in mutant data
    MT_mean = np.array(df_mutant['obs_fluo_mean'])
    MT_sd = np.array(df_mutant['obs_SD'])
    G_obs = np.log10(np.divide(MT_mean, g_WT))
    return G_obs, MT_sd

def Epistasis(mutants:list, model = 'logadd'):
    G = G_obs(mutants)
    G_obs_mean = G[0]
    G_obs_std = G[1]
    if model == 'logadd':
        Ghat_logadd = G_hat_logadd(mutants)
        Ghat_logadd_mean = Ghat_logadd[0]
        Ghat_logadd_std = Ghat_logadd[1]
        Epsilon =  G[0] - Ghat_logadd[0]
        p_val = stats.ttest_ind_from_stats(mean1 = Ghat_logadd_mean,std1 = Ghat_logadd_std, nobs1 = 3, mean2 = G_obs_mean, std2 = G_obs_std, nobs2 = 3)[1]  
    else:
        Ghat = G_hat(model, mutants)
        Epsilon = G[0] - Ghat
        p_val = 0
    return Epsilon, p_val

#calculate epistasis for all double mutants
def get_Eps(model='logadd'):
    Ep_low = []
    Ep_low_pVal = []
    Ep_medium = []
    Ep_med_pVal = []
    Ep_high = []
    Ep_high_pVal = []
    print("printing every 100th mutant processed to show progress...")
    for i, mut_id in enumerate(df_M['genotype'][(df_M['inducer level'] == 'low')][1:]):
        if i % 100 == 0:
            print(mut_id)
        mut_names = get_mut_names(mut_id)
        #error: stops at Sensor7 - problem searching df_S for "0.0002"
        Ep_low += [Epistasis(mut_names, model)[0][0]]
        Ep_low_pVal += [Epistasis(mut_names, model)[0][0]]
        Ep_medium += [Epistasis(mut_names, model)[0][1]]
        Ep_med_pVal += [Epistasis(mut_names, model)[1][1]]
        Ep_high += [Epistasis(mut_names, model)[0][2]]
        Ep_high_pVal += [Epistasis(mut_names, model)[1][2]]

    Eps = [math.nan]*3 +Ep_low + Ep_medium + Ep_high
    Eps_pVal = [math.nan]*3 +Ep_low_pVal + Ep_med_pVal + Ep_high_pVal
    return Eps, Eps_pVal

Epsistases = get_Eps()
df_M["mean_epistasis"] = Epsistases[0]
df_M["pVal_epistasis"] = Epsistases[1]

df_Eps = df_M.loc[(df_M['inducer level'] == 'low') & (df_M['genotype category'] != 'single')].copy()

Eps_ml = np.subtract(df_M['mean_epistasis'][df_M['inducer level'] == 'medium'].to_numpy() , df_M['mean_epistasis'][df_M['inducer level'] == 'low'].to_numpy())
Eps_hm = np.subtract(df_M['mean_epistasis'][df_M['inducer level'] == 'high'].to_numpy() , df_M['mean_epistasis'][df_M['inducer level'] == 'medium'].to_numpy())

df_Eps[['med-low', 'high-med']] = [Eps_ml[1:], Eps_hm[1:]]
df_Eps['med-low']= list(Eps_ml[1:])

plt.hist(df_M['mean epistasis'][df_M['genotype category']== 'triple'])

np.mean(df_M['mean epistasis'][df_M['genotype category']== 'pairwise'])
np.power(np.mean(df_M['mean epistasis'][df_M['genotype category']== 'triple']),10)