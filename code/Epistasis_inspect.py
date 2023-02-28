import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Epistasis_calc import df_Eps

df_M = pd.read_excel('../data/df_M.xlsx')
df_Eps = pd.read_excel('../data/df_Eps.xlsx')

#plot significant epistasis
cols = np.where(df_M['Sig_Epistasis'] == True, 'darkgray', np.where(df_M['genotype category'] == 'pairwise', 'palegreen', 'lightskyblue'))
alphas = np.where(df_M['Sig_Epistasis'] == True, 0.5, 0.8)
fig, ax = plt.subplots()
ax.scatter(df_M['pVal_epistasis'], df_M['mean_epistasis'], c = cols,alpha = alphas, zorder = 2.5, edgecolor = 'black' )
ax.axvline(0.05, c = 'darkgray', ls = '--')
plt.show

#plot inducer dependant epistasis
cols = np.where(df_Eps['Sig_Epistasis'] == True, 'darkgray', np.where(df_Eps['genotype category'] == 'pairwise', 'palegreen', 'lightskyblue'))
alphas = np.where(df_Eps['Sig_Epistasis'] == True, 0.5, 0.8)
n_notSig = np.count_nonzero(cols == 'darkgray')
n_pairwise = np.count_nonzero(cols == 'darkgray')
n_triplet = np.count_nonzero(cols == 'darkgray')
prop_upperLeft = np.divide(np.count_nonzero((df_Eps['med-low'] >0) & (df_Eps['high-med'] <0)), len(df_Eps['med-low']))
prop_lowerRight = np.divide(np.count_nonzero((df_Eps['med-low'] <0) & (df_Eps['high-med'] >0)), len(df_Eps['med-low']))
prop_upperRight = np.divide(np.count_nonzero((df_Eps['med-low'] >0) & (df_Eps['high-med'] >0)), len(df_Eps['med-low']))
prop_lowerLeft = np.divide(np.count_nonzero((df_Eps['med-low'] <0) & (df_Eps['high-med'] <0)), len(df_Eps['med-low']))

#%%
fig, ax = plt.subplots()
scatter = ax.scatter(df_Eps['med-low'], df_Eps['high-med'], c = cols, alpha = alphas, zorder = 2.5, edgecolor = 'black')
legend1 = ax.legend(*scatter.legend_elements(), loc="lower left")
ax.add_artist(legend1)
ax.axhline(y=0, c = 'darkgray', ls = '--')
ax.axvline(x=0, c = 'darkgray', ls = '--')
ax.set_xlim([-1,1])
ax.set_ylim([-1,1])
ax.set_xlabel('$\u03B5_{medium}$ - $\u03B5_{low}$')
ax.set_ylabel('$\u03B5_{high}$ - $\u03B5_{medium}$')
ax.set_title('inducer dependant Epistasis for '+ 'logadd model')

#proportins mutants in each quadrant
ax.text(-0.75, -0.75, 'n = '+ str(np.round(prop_lowerLeft*100, 0))+'%', verticalalignment='center', horizontalalignment='center')
ax.text(0.75, -0.75, 'n = '+ str(np.round(prop_lowerRight*100, 0))+'%', verticalalignment='center', horizontalalignment='center')
ax.text(-0.75, 0.75, 'n = '+ str(np.round(prop_upperLeft*100, 0))+'%', verticalalignment='center', horizontalalignment='center')
ax.text(0.75, 0.75, str(np.round(prop_upperRight*100, 0))+'%', verticalalignment='center', horizontalalignment='center')
plt.show()