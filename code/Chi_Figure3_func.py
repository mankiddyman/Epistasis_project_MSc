import numpy as np
from scipy import stats
from Models import *
from Epistasis_calc_functions import *
# from Model_Diagnostics_Functions import *
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import cm
import statistics as s
import matplotlib.patches as mpatches

"Functions to plot figure 3 from Ruben's paper. Use Sort_mutants() + convertDF() is you need a dataframe"
"with Epistasis values and genotype info. Otheriwse, use Epistasis_hist() to get plots. May need to change fetching"
"path and file name for Episatsis files but can also use get_Eps instead if working."

#%%
def Sort_mutants(model): #observed, model_hill, and thermodynamic
    # if model == str('observed'):
    #     Eps_df = pd.read_excel('../results/Eps_observed.xlsx') #takes from results.
    # elif model == str('model_hill'):
    #     Eps_df = pd.read_excel('../results/hill_epistasis.xlsx') #saves time if already generated, just modift the df.
    # elif model == str('thermodynamic'):
    #     Eps_df = pd.read_excel('../results/thermodynamic_epistasis.xlsx')
    # else:
    #     print('the inputted model is invalid, please selec from observed, model_hill or thermodynamic')

    if model == str('observed'):
        Eps_df = get_Eps('observed') 
    elif model == str('model_hill'):
        Eps_df = get_Eps(model=model_hill) 
    elif model == str('thermodynamic'):
        Eps_df = get_Eps(model=model_thermodynamic)
    else:
        print('the inputted model is invalid, please selec from observed, model_hill or thermodynamic')
    #output
        #filter out O10
    o10 = Eps_df[Eps_df['genotype'].str.contains(str('O10'))] 
        #filter out 01 by removing duplicates from o10, merge two dataframes and assign boolean column for duplicated then filter out.
    o1 = Eps_df[Eps_df['genotype'].str.contains(str('O1'))] 
    merge = pd.concat([o1,o10], axis=0)
    merge['dups'] = merge.duplicated(keep=False)
    o1 = merge[(merge['dups']==False)]


    number = [2,3,4,5,6,7,8,9] #Filter out rest of mutant genotypes
    Output_df = pd.DataFrame()
    Output_df['O1'] = o1['Ep']
    Output_df.reset_index(drop=True, inplace=True) #gets rid of index from original table
    

    for num in number:
        genotype = 'O'+ str(num)
        temp = Eps_df[Eps_df['genotype'].str.contains(genotype)] #could simply with .query for genotypes
        temp.reset_index(drop=True, inplace=True)
        Output_df[genotype] = temp['Ep']
    o10.reset_index(drop=True, inplace=True)
    Output_df['O10'] = o10['Ep']

    #Sensor
    s10 = Eps_df[Eps_df['genotype'].str.contains(str('S10'))] 
    s1 = Eps_df[Eps_df['genotype'].str.contains(str('S1'))] 
    mergeS = pd.concat([s1,s10], axis=0) #merge dataframes
    mergeS['dups'] = mergeS.duplicated(keep=False)
    s1 = mergeS[(mergeS['dups']==False)]

    Sensor_df = pd.DataFrame()
    Sensor_df['S1'] = s1['Ep']
    Sensor_df.reset_index(drop=True, inplace=True)
    for num in number:
        genotype = 'S'+ str(num)
        temp = Eps_df[Eps_df['genotype'].str.contains(genotype)]
        temp.reset_index(drop=True, inplace=True)
        Sensor_df[genotype] = temp['Ep']
    s10.reset_index(drop=True, inplace=True)
    Sensor_df['S10'] = s10['Ep']

    #Regulator
    r10 = Eps_df[Eps_df['genotype'].str.contains(str('R10'))] 
    r1 = Eps_df[Eps_df['genotype'].str.contains(str('R1'))] 
    mergeR = pd.concat([r1,r10], axis=0)
    mergeR['dups'] = mergeR.duplicated(keep=False)
    r1 = mergeR[(mergeR['dups']==False)]

    #Filter out rest of mutant genotypes
    Regulator_df = pd.DataFrame()
    Regulator_df['R1'] = r1['Ep']
    Regulator_df.reset_index(drop=True, inplace=True)
    for num in number:
        genotype = 'R'+ str(num)
        temp = Eps_df[Eps_df['genotype'].str.contains(genotype)]
        temp.reset_index(drop=True, inplace=True)
        Regulator_df[genotype] = temp['Ep']
    r10.reset_index(drop=True, inplace=True)
    Regulator_df['R10'] = r10['Ep']

    return Output_df, Regulator_df, Sensor_df

#%%
def convertDF(dataframe,node): #need to be dimensions of 360, node is "Output","Regulator", or "Sensor"
#final dataframe
    if node == 'Output':
        node = 'O'
    elif node == 'Regulator':
        node = 'R'
    elif node == 'Sensor':
        node = 'S'
    else:
        print("incorrect node, please enter either Output, Regulator or Sensor")
        node = input("enter a node: ")
    col_list = []
    for column in dataframe:
        eps = dataframe[column]
        eps_list = list(eps)
        for i in eps_list:
            col_list.append(i)
    Final_df = pd.DataFrame(columns=['Epistasis','Genotype','Category','LMH'])
    Final_df['Epistasis'] = col_list
    
    gen_list = []
    for i in range(1,11):
        name = node + str(i)
        gen_list += [name]*360
    Final_df['Genotype'] = gen_list

    Category_list = []
    Category_list = (['pairwise'] * 60 + ['triplet'] * 300)*10 #pairwise = 0, triplet = 1
    Final_df['Category'] = Category_list

    LMH_list = []
    LMH_list = (['low']*20 + ['medium']*20 + ['high']*20 + ['low']*100 + ['medium']*100 + ['high']*100)*10
    Final_df['LMH'] = LMH_list

    return Final_df

# %%
def Epistasis_hist(model): #plots epistasis values for all genotypes + median +stdev (need to plot mean instead)

    # model = 'observed'

    Out, Reg, Sen = Sort_mutants(model)
    New_Out = convertDF(Out,'Output')
    New_Reg = convertDF(Reg,'Regulator')
    New_Sen = convertDF(Sen,'Sensor')
    
    figure, axis = plt.subplots(3,1)

    #output
    group = 'Genotype'
    column = 'Epistasis'
    Cato = 'Category'

    grouped = New_Out.groupby(group)

    names, vals, cats, xs, mean, std = [], [], [], [], [], []

    for i, (name, subdf) in enumerate(grouped):
        names.append(name)
        vals.append(subdf[column].tolist())
        cats.append(subdf[Cato].tolist())
        xs.append(np.random.normal(i+1, 0.04, subdf.shape[0]))

    names = list(Out.keys()) #proper order of genotypes
    for x, val, cat in zip(xs, vals, cats):
        pair = val[0:60]
        trip = val[60:]
        x1 = x[0:60]
        x2 = x[60:]
        # temp = s.mean(val)
        # mean.append([temp])
        # temp = s.stdev(val)
        # std.append([temp])
        std = s.stdev(val)
            # axis[0].scatter(rand_jitter(x1), pair, c= 'yellow' , alpha=1/5, marker='o')
            # axis[0].scatter(rand_jitter(x2), trip, c= 'paleturquoise', alpha=1/5, marker='o')
        sns.regplot(x = x1, y = pair, fit_reg = False, x_jitter = 0.2, 
                    color= 'grey', marker='.', ax=axis[0], scatter_kws={'alpha':0.7})
        sns.regplot(x = x2, y = trip, fit_reg = False, x_jitter = 0.2, color= 'mediumseagreen',
                    marker='.', ax=axis[0], scatter_kws={'alpha':0.1})
        #axis[0].errorbar(x1, pair, yerr=std, fmt='none', lw = 3)
    medianp = dict(color='black',linewidth=1.5)
    axis[0].boxplot(vals, labels=names, showfliers=False, patch_artist=False, showbox=False, medianprops=medianp)

    xmin, xmax, ymin, ymax = axis[0].axis()
    lims = 'y-max:' + str(round(ymax,2)) + ' y-min:' + str(round(ymin,2))
    axis[0].annotate(lims, (8.8,ymax*0.8), fontsize = 5)

    # axis[0].errorbar(np.arange(1,10), mean, std, fmt='ok', lw = 3)

    axis[0].set(xlabel='Output Genotypes', title = 'Variation of Epistasis '+ '| model: '+ model)
    axis[0].tick_params(axis='x', which='major', labelsize=8)
    axis[0].hlines(y= 0, xmin= 0.5, xmax= 10.5, linestyles='dashed', linewidth = 1, alpha = 0.45, colors='k')

    pairwise = mpatches.Patch(color= 'grey', label='Pairwise')
    Triplet = mpatches.Patch(color= 'mediumseagreen', label='Triplet (O)')
    TripletR = mpatches.Patch(color= 'royalblue', label='Triplet (R)')
    TripletS = mpatches.Patch(color= 'lightcoral', label='Triplet (S)')
    axis[0].legend(handles=[pairwise,Triplet,TripletR, TripletS], bbox_to_anchor=(1, 1), title = "Legend")

     #Regulator
    group = 'Genotype'
    column = 'Epistasis'
    Cato = 'Category'

    grouped = New_Reg.groupby(group)

    names, vals, cats, xs, mean, std = [], [], [], [], [], []

    for i, (name, subdf) in enumerate(grouped):
        names.append(name)
        vals.append(subdf[column].tolist())
        cats.append(subdf[Cato].tolist())
        xs.append(np.random.normal(i+1, 0.04, subdf.shape[0]))

    names = list(Reg.keys()) #proper order of genotypes
    for x, val, cat in zip(xs, vals, cats):
        pair = val[0:60]
        trip = val[60:]
        x1 = x[0:60]
        x2 = x[60:]
        # temp = s.mean(val)
        # mean.append([temp])
        # temp = s.stdev(val)
        # std.append([temp])
        std = s.stdev(val)
            # axis[0].scatter(rand_jitter(x1), pair, c= 'yellow' , alpha=1/5, marker='o')
            # axis[0].scatter(rand_jitter(x2), trip, c= 'paleturquoise', alpha=1/5, marker='o')
        sns.regplot(x = x1, y = pair, fit_reg = False, x_jitter = 0.2, 
                    color= 'grey', marker='.', ax=axis[1], scatter_kws={'alpha':0.7})
        sns.regplot(x = x2, y = trip, fit_reg = False, x_jitter = 0.2, color= 'royalblue',
                    marker='.', ax=axis[1], scatter_kws={'alpha':0.1})
        #axis[0].errorbar(x1, pair, yerr=std, fmt='none', lw = 3)
    medianp = dict(color='black',linewidth=1.5)
    axis[1].boxplot(vals, labels=names, showfliers=False, patch_artist=False, showbox=False, medianprops=medianp)

    # axis[0].errorbar(np.arange(1,10), mean, std, fmt='ok', lw = 3)

    axis[1].set(ylabel='Epistasis (\u03B5)')
    axis[1].tick_params(axis='x', which='major', labelsize=8)
    axis[1].hlines(y= 0, xmin= 0.5, xmax= 10.5, linestyles='dashed', linewidth = 1, alpha = 0.45, colors='k')

    xmin, xmax, ymin, ymax = axis[1].axis()
    lims = 'y-max:' + str(round(ymax,2)) + ' y-min:' + str(round(ymin,2))
    axis[1].annotate(lims, (8.8,ymax*0.8), fontsize = 5)
    #Sensor
    group = 'Genotype'
    column = 'Epistasis'
    Cato = 'Category'

    grouped = New_Sen.groupby(group)

    names, vals, cats, xs, mean, std = [], [], [], [], [], []

    for i, (name, subdf) in enumerate(grouped):
        names.append(name)
        vals.append(subdf[column].tolist())
        cats.append(subdf[Cato].tolist())
        xs.append(np.random.normal(i+1, 0.04, subdf.shape[0]))

    names = list(Sen.keys()) #proper order of genotypes
    for x, val, cat in zip(xs, vals, cats):
        pair = val[0:60]
        trip = val[60:]
        x1 = x[0:60]
        x2 = x[60:]
        # temp = s.mean(val)
        # mean.append([temp])
        # temp = s.stdev(val)
        # std.append([temp])
        std = s.stdev(val)
            # axis[0].scatter(rand_jitter(x1), pair, c= 'yellow' , alpha=1/5, marker='o')
            # axis[0].scatter(rand_jitter(x2), trip, c= 'paleturquoise', alpha=1/5, marker='o')
        sns.regplot(x = x1, y = pair, fit_reg = False, x_jitter = 0.2, 
                    color= 'grey', marker='.', ax=axis[2], scatter_kws={'alpha':0.7})
        sns.regplot(x = x2, y = trip, fit_reg = False, x_jitter = 0.2, color= 'lightcoral',
                    marker='.', ax=axis[2], scatter_kws={'alpha':0.1})
        #axis[0].errorbar(x1, pair, yerr=std, fmt='none', lw = 3)
    medianp = dict(color='black',linewidth=1.5)
    axis[2].boxplot(vals, labels=names, showfliers=False, patch_artist=False, showbox=False, medianprops=medianp)

    # axis[0].errorbar(np.arange(1,10), mean, std, fmt='ok', lw = 3)

    axis[2].set(xlabel='Genotype')
    axis[2].tick_params(axis='x', which='major', labelsize=8)
    axis[2].hlines(y= 0, xmin= 0.5, xmax= 10.5, linestyles='dashed', linewidth = 1, alpha = 0.45, colors='k')

    xmin, xmax, ymin, ymax = axis[2].axis()
    lims = 'y-max:' + str(round(ymax,2)) + ' y-min:' + str(round(ymin,2))
    axis[2].annotate(lims, (8.8,ymax*0.8), fontsize = 5)
    
    # plt.gcf().set_size_inches(7,5)
    plt.savefig("../results/" + model + "_Epistasis_Histograms.pdf", format="pdf", bbox_inches='tight')

    return 

# %%
######################
# test
# Out_thermo, Reg_thermo, Sen_thermo, = Sort_mutants('thermodynamic')

# Out_thermo_df = convertDF(Out_thermo,'Output')
# Reg_thermo_df = convertDF(Reg_thermo,'Regulator')
# Sen_thermo_df = convertDF(Sen_thermo,'Sensor')

# filepath = '../results/SensorDF.csv'
# Sen_thermo_df.to_csv(filepath, index=False)

#save plots here

Epistasis_hist('observed')

Epistasis_hist('model_hill')

Epistasis_hist('thermodynamic')



#%%