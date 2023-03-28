import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from Chi_Figure3_func import Sort_mutants, convertDF

#%%
def epsVP_df():
    O_obs, R_obs, S_obs = Sort_mutants(model = 'hill', strat = 'all')
    N_O_obs = convertDF(O_obs,'Output', 'observed')
    N_R_obs = convertDF(R_obs,'Regulator', 'observed')
    N_S_obs = convertDF(S_obs,'Sensor', 'observed')
    Out, Reg, Sen = Sort_mutants('hill', strat = 'all')
    New_Out = convertDF(Out,'Output','hill')
    New_Reg = convertDF(Reg,'Regulator','hill')
    New_Sen = convertDF(Sen,'Sensor','hill')
    
    hill_df = pd.concat([N_O_obs,N_R_obs,N_S_obs,New_Out,New_Reg,New_Sen])


    OutT, RegT, SenT = Sort_mutants('thermodynamic', 'sensor')
    New_OutT = convertDF(OutT,'Output','thermodynamic')
    New_RegT = convertDF(RegT,'Regulator','thermodynamic')
    New_SenT = convertDF(SenT,'Sensor','thermodynamic')
    
    therm_df = pd.concat([N_O_obs,N_R_obs,N_S_obs,New_OutT,New_RegT,New_SenT])

    return hill_df, therm_df
#%%
Hill, therm = epsVP_df()
#%%
# Hill Model
hill_out = Hill
hill_out = hill_out[hill_out["Genotype"].str.contains("R")==False] 
hill_out = hill_out[hill_out["Genotype"].str.contains("S")==False] 
hill_reg = Hill
hill_reg = hill_reg[hill_reg["Genotype"].str.contains("O")==False] 
hill_reg = hill_reg[hill_reg["Genotype"].str.contains("S")==False] 
hill_sen = Hill
hill_sen = hill_sen[hill_sen["Genotype"].str.contains("O")==False] 
hill_sen = hill_sen[hill_sen["Genotype"].str.contains("R")==False] 
hill_sen_p = hill_sen
hill_sen_t = hill_sen
hill_sen_p = hill_sen_p[hill_sen_p["Category"].str.contains("triplet")==False] 
hill_sen_t = hill_sen_t[hill_sen_t["Category"].str.contains("pairwise")==False] 
hill_out_p = hill_out
hill_out_t = hill_out
hill_out_p = hill_out_p[hill_out_p["Category"].str.contains("triplet")==False] 
hill_out_t = hill_out_t[hill_out_t["Category"].str.contains("pairwise")==False] 
hill_reg_p = hill_reg
hill_reg_t = hill_reg
hill_reg_p = hill_reg_p[hill_reg_p["Category"].str.contains("triplet")==False] 
hill_reg_t = hill_reg_t[hill_reg_t["Category"].str.contains("pairwise")==False]
hill_out_p.loc[hill_out_p["model"].str.contains("hill"), "model"] = "hill"
hill_sen_p.loc[hill_sen_p["model"].str.contains("hill"), "model"] = "hill"
hill_reg_p.loc[hill_reg_p["model"].str.contains("hill"), "model"] = "hill" 
hill_out_t.loc[hill_out_t["model"].str.contains("hill"), "model"] = "hill"
hill_sen_t.loc[hill_sen_t["model"].str.contains("hill"), "model"] = "hill"
hill_reg_t.loc[hill_reg_t["model"].str.contains("hill"), "model"] = "hill"
# Double
fig = plt.figure()
plt.subplot(3,1,1)
plt.title("Hill Model Pairwise Epistasis", fontweight='bold')
hp1 = sns.violinplot(data=hill_out_p, x="Genotype", y="Epistasis", hue="model", split=True)
hp1.set(ylabel=None)
sns.move_legend(hp1, "upper left", bbox_to_anchor=(1, 1))
plt.subplot(3,1,2)
hp2 = sns.violinplot(data=hill_reg_p, x="Genotype", y="Epistasis", hue="model", split=True)
hp2.legend([],[], frameon=False)
plt.subplot(3,1,3)
hp3 = sns.violinplot(data=hill_sen_p, x="Genotype", y="Epistasis", hue="model", split=True)
hp3.set(ylabel=None)
hp3.legend([],[], frameon=False)
# Triple
fig = plt.figure()
plt.subplot(3,1,1)
hp1 = sns.violinplot(data=hill_out_t, x="Genotype", y="Epistasis", hue="model", split=True)
hp1.set(ylabel=None)
plt.title("Hill Model Triple Epistasis", fontweight='bold')
sns.move_legend(hp1, "upper left", bbox_to_anchor=(1, 1))
plt.subplot(3,1,2)
hp2 = sns.violinplot(data=hill_reg_t, x="Genotype", y="Epistasis", hue="model", split=True)
hp2.legend([],[], frameon=False)
plt.subplot(3,1,3)
hp3 = sns.violinplot(data=hill_sen_t, x="Genotype", y="Epistasis", hue="model", split=True)
hp3.set(ylabel=None)
hp3.legend([],[], frameon=False)

#%%
# Thermodynamic Model
td_out = therm
td_out = td_out[td_out["Genotype"].str.contains("R")==False] 
td_out = td_out[td_out["Genotype"].str.contains("S")==False] 
td_reg = therm
td_reg = td_reg[td_reg["Genotype"].str.contains("O")==False] 
td_reg = td_reg[td_reg["Genotype"].str.contains("S")==False] 
td_sen = therm
td_sen = td_sen[td_sen["Genotype"].str.contains("O")==False] 
td_sen = td_sen[td_sen["Genotype"].str.contains("R")==False] 
td_sen_p = td_sen
td_sen_t = td_sen
td_sen_p = td_sen_p[td_sen_p["Category"].str.contains("triplet")==False] 
td_sen_t = td_sen_t[td_sen_t["Category"].str.contains("pairwise")==False] 
td_out_p = td_out
td_out_t = td_out
td_out_p = td_out_p[td_out_p["Category"].str.contains("triplet")==False] 
td_out_t = td_out_t[td_out_t["Category"].str.contains("pairwise")==False] 
td_reg_p = td_reg
td_reg_t = td_reg
td_reg_p = td_reg_p[td_reg_p["Category"].str.contains("triplet")==False] 
td_reg_t = td_reg_t[td_reg_t["Category"].str.contains("pairwise")==False]
td_out_p.loc[td_out_p["model"].str.contains("therm"), "model"] = "thermodynamic"
td_sen_p.loc[td_sen_p["model"].str.contains("therm"), "model"] = "thermodynamic"
td_reg_p.loc[td_reg_p["model"].str.contains("therm"), "model"] = "thermodynamic" 
td_out_t.loc[td_out_t["model"].str.contains("therm"), "model"] = "thermodynamic"
td_sen_t.loc[td_sen_t["model"].str.contains("therm"), "model"] = "thermodynamic"
td_reg_t.loc[td_reg_t["model"].str.contains("therm"), "model"] = "thermodynamic"
# Double
fig = plt.figure()
plt.subplot(3,1,1)
plt.title("Thermodynamic Model Pairwise Epistasis", fontweight='bold')
tp1 = sns.violinplot(data=td_out_p, x="Genotype", y="Epistasis", hue="model", split=True)
tp1.set(ylabel=None)
sns.move_legend(tp1, "upper left", bbox_to_anchor=(1, 1))
plt.subplot(3,1,2)
tp2 = sns.violinplot(data=td_reg_p, x="Genotype", y="Epistasis", hue="model", split=True)
tp2.legend([],[], frameon=False)
plt.subplot(3,1,3)
tp3 = sns.violinplot(data=td_sen_p, x="Genotype", y="Epistasis", hue="model", split=True)
tp3.set(ylabel=None)
tp3.legend([],[], frameon=False)
# Triple
fig = plt.figure()
plt.subplot(3,1,1)
plt.title("Thermodynamic Model Triple Epistasis", fontweight='bold')
tp1 = sns.violinplot(data=td_out_t, x="Genotype", y="Epistasis", hue="model", split=True)
tp1.set(ylabel=None)
sns.move_legend(tp1, "upper left", bbox_to_anchor=(1, 1))
plt.subplot(3,1,2)
tp2 = sns.violinplot(data=td_reg_t, x="Genotype", y="Epistasis", hue="model", split=True)
tp2.legend([],[], frameon=False)
plt.subplot(3,1,3)
tp3 = sns.violinplot(data=td_sen_t, x="Genotype", y="Epistasis", hue="model", split=True)
tp3.set(ylabel=None)
tp3.legend([],[], frameon=False)