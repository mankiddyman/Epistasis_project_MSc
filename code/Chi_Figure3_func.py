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

Eps_toExcel('observed')
Eps_toExcel(model = model_hill)
Eps_toExcel(model = model_thermodynamic)
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
        Eps_df = pd.read_excel('../results/Eps_observed.xlsx') 
    elif model == str('model_hill'):
        Eps_df = pd.read_excel('../results/Eps_model_hill.xlsx') 
    elif model == str('thermodynamic'):
        Eps_df = pd.read_excel('../results/Eps_model_thermodynamic.xlsx') 
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
    Category_list = ((['pairwise'] * 20 + ['triplet'] * 100)*3)*10 #pairwise = 0, triplet = 1
    Final_df['Category'] = Category_list

    LMH_list = []
    LMH_list = (['low']*120 + ['medium']*120 + ['high']*120)*10
    Final_df['LMH'] = LMH_list

    return Final_df

# %%
def Epistasis_hist(model): #plots epistasis values for all genotypes + median +stdev (need to plot mean instead)

    Out, Reg, Sen = Sort_mutants(model)
    New_Out = convertDF(Out,'Output')
    New_Reg = convertDF(Reg,'Regulator')
    New_Sen = convertDF(Sen,'Sensor')
    
    figure, axis = plt.subplots(3,1)

    #output
    group = 'Genotype'
    column = 'Epistasis'
    Cato = 'Category'

    grouped = New_Out.groupby([group, Cato], sort = False)


    names, vals, xs, means, stds = [], [], [], [], []

    for i, (name, subdf) in enumerate(grouped):
        vals.append(subdf[column].tolist())
        means.append(s.mean(vals[i]))
        stds.append(s.stdev(vals[i]))
        xs.append(np.random.normal(i+1, 0.04, subdf.shape[0]))

    x1, x2, vals_pair, vals_trip ,mean_pair, std_pair, mean_trip, std_trip = [] , [] , [], [], [], [], [], []
    abc = np.arange(1,len(vals)+1,2)
    edf = np.arange(2,len(vals)+1,2) #proper order of genotypes

    names = list(Out.keys())

    for i in range(0,len(vals),2):
        vals_pair.append(vals[i])
        x1.append(xs[i])
        mean_pair.append(means[i])
        std_pair.append(stds[i])
    for i in range(1,len(vals),2):
        vals_trip.append(vals[i])
        x2.append(xs[i])
        mean_trip.append(means[i])
        std_trip.append(stds[i])

    for val_p, val_t, x_1, x_2, mean_p, std_p, mean_t, std_t, a, e in zip(vals_pair, vals_trip, x1, x2, mean_pair, std_pair, mean_trip, std_trip, abc, edf):
        sns.regplot(x = x_1, y = val_p, fit_reg = False, x_jitter = 0.2, 
                        color= 'seagreen', marker='.', ax=axis[0], scatter_kws={'alpha':0.7})
        sns.regplot(x = x_2, y = val_t, fit_reg = False, x_jitter = 0.2, 
                        color= 'mediumseagreen', marker='.', ax=axis[0], scatter_kws={'alpha':0.7})
        axis[0].errorbar(a , mean_p, std_p, linestyle='None',
                          fmt='_', c = 'k', elinewidth = 1.5, capsize = 1.5)
        axis[0].errorbar(e , mean_t, std_t, linestyle='None',
                          fmt='_', c = 'k', elinewidth = 1.5, capsize = 1.5)

    x=[]    
    for i in range(0,len(xs),2):
        avg = (s.mean(xs[i]) + (s.mean(xs[i+1])))/2
        x.append(avg)

    plt.sca(axis[0])
    plt.xticks(x,names)

    xmin, xmax, ymin, ymax = axis[0].axis()
    lims = 'y-max:' + str(round(ymax,2)) + ' y-min:' + str(round(ymin,2))
    axis[0].annotate(lims, (max(x2[9]-3),ymax*0.8), fontsize = 5)


    axis[0].set(xlabel='Output Genotypes', title = 'Variation of Epistasis '+ '| model: '+ model)
    axis[0].tick_params(axis='x', which='major', labelsize=8)
    axis[0].axhline(0, linestyle='dashed', linewidth = 1, alpha = 0.45, c='k')
    

    pairwiseO = mpatches.Patch(color= 'seagreen', label='Pairwise (O)')
    pairwiseR = mpatches.Patch(color= 'royalblue', label='Pairwise (R)')
    pairwiseS = mpatches.Patch(color= 'orangered', label='Pairwise (S)')
    TripletO = mpatches.Patch(color= 'mediumseagreen', label='Triplet (O)')
    TripletR = mpatches.Patch(color= 'cornflowerblue', label='Triplet (R)')
    TripletS = mpatches.Patch(color= 'lightcoral', label='Triplet (S)')
    axis[0].legend(handles=[pairwiseO, TripletO, pairwiseR, TripletR, pairwiseS, TripletS], bbox_to_anchor=(1, 1), title = "Legend")
     
    #Regulator
    grouped = New_Reg.groupby([group, Cato], sort = False)

    names, vals, xs, means, stds = [], [], [], [], []

    for i, (name, subdf) in enumerate(grouped):
        vals.append(subdf[column].tolist())
        means.append(s.mean(vals[i]))
        stds.append(s.stdev(vals[i]))
        xs.append(np.random.normal(i+1, 0.04, subdf.shape[0]))

    x1, x2, vals_pair, vals_trip ,mean_pair, std_pair, mean_trip, std_trip = [] , [] , [], [], [], [], [], []
    abc = np.arange(1,len(vals)+1,2)
    edf = np.arange(2,len(vals)+1,2) #proper order of genotypes

    names = list(Reg.keys()) #proper order of genotypes

    for i in range(0,len(vals),2):
        vals_pair.append(vals[i])
        x1.append(xs[i])
        mean_pair.append(means[i])
        std_pair.append(stds[i])
    for i in range(1,len(vals),2):
        vals_trip.append(vals[i])
        x2.append(xs[i])
        mean_trip.append(means[i])
        std_trip.append(stds[i])

    for val_p, val_t, x_1, x_2, mean_p, std_p, mean_t, std_t, a, e in zip(vals_pair, vals_trip, x1, x2, mean_pair, std_pair, mean_trip, std_trip, abc, edf):
        sns.regplot(x = x_1, y = val_p, fit_reg = False, x_jitter = 0.2, 
                        color= 'royalblue', marker='.', ax=axis[1], scatter_kws={'alpha':0.7})
        sns.regplot(x = x_2, y = val_t, fit_reg = False, x_jitter = 0.2, 
                        color= 'cornflowerblue', marker='.', ax=axis[1], scatter_kws={'alpha':0.7})
        axis[1].errorbar(a , mean_p, std_p, linestyle='None',
                          fmt='_', c = 'k', elinewidth = 1.5, capsize = 1.5)
        axis[1].errorbar(e , mean_t, std_t, linestyle='None',
                          fmt='_', c = 'k', elinewidth = 1.5, capsize = 1.5)

    axis[1].set(ylabel='Epistasis (\u03B5)')
    axis[1].tick_params(axis='x', which='major', labelsize=8)

    x=[]    
    for i in range(0,len(xs),2):
        avg = (s.mean(xs[i]) + (s.mean(xs[i+1])))/2
        x.append(avg)

    plt.sca(axis[1])
    plt.xticks(x,names)
        
    

    xmin, xmax, ymin1, ymax1 = axis[1].axis()
    lims = 'y-max:' + str(round(ymax1,2)) + ' y-min:' + str(round(ymin1,2))
    axis[1].annotate(lims, (max(x2[9]-3),ymax1*0.8), fontsize = 5)
    axis[1].axhline(0, linestyle='dashed', linewidth = 1, alpha = 0.45, c='k')

    #Sensor
    grouped = New_Sen.groupby([group, Cato], sort = False)

    names, vals, xs, means, stds = [], [], [], [], []

    for i, (name, subdf) in enumerate(grouped):
        vals.append(subdf[column].tolist())
        means.append(s.mean(vals[i]))
        stds.append(s.stdev(vals[i]))
        xs.append(np.random.normal(i+1, 0.04, subdf.shape[0]))

    x1, x2, vals_pair, vals_trip ,mean_pair, std_pair, mean_trip, std_trip = [] , [] , [], [], [], [], [], []
    abc = np.arange(1,len(vals)+1,2)
    edf = np.arange(2,len(vals)+1,2) #proper order of genotypes

    names = list(Sen.keys()) #proper order of genotypes
    for i in range(0,len(vals),2):
        vals_pair.append(vals[i])
        x1.append(xs[i])
        mean_pair.append(means[i])
        std_pair.append(stds[i])
    for i in range(1,len(vals),2):
        vals_trip.append(vals[i])
        x2.append(xs[i])
        mean_trip.append(means[i])
        std_trip.append(stds[i])

    for val_p, val_t, x_1, x_2, mean_p, std_p, mean_t, std_t, a, e in zip(vals_pair, vals_trip, x1, x2, mean_pair, std_pair, mean_trip, std_trip, abc, edf):
        sns.regplot(x = x_1, y = val_p, fit_reg = False, x_jitter = 0.2, 
                        color= 'orangered', marker='.', ax=axis[2], scatter_kws={'alpha':0.7})
        sns.regplot(x = x_2, y = val_t, fit_reg = False, x_jitter = 0.2, 
                        color= 'lightcoral', marker='.', ax=axis[2], scatter_kws={'alpha':0.7})
        axis[2].errorbar(a , mean_p, std_p, linestyle='None',
                          fmt='_', c = 'k', elinewidth = 1.5, capsize = 1.5)
        axis[2].errorbar(e , mean_t, std_t, linestyle='None',
                          fmt='_', c = 'k', elinewidth = 1.5, capsize = 1.5)

    
    axis[2].set(xlabel='Genotype')
    axis[2].tick_params(axis='x', which='major', labelsize=8)
    x=[]    
    for i in range(0,len(xs),2):
        avg = (s.mean(xs[i]) + (s.mean(xs[i+1])))/2
        x.append(avg)

    plt.sca(axis[2])
    plt.xticks(x,names)

    xmin, xmax, ymin2, ymax2 = axis[2].axis()
    lims = 'y-max:' + str(round(ymax2,2)) + ' y-min:' + str(round(ymin2,2))
    axis[2].annotate(lims, (max(x2[9]-3),ymax2*0.8), fontsize = 5)


    axis[2].set(xlabel='Genotypes')
    axis[2].tick_params(axis='x', which='major', labelsize=8)
    axis[2].axhline(0, linestyle='dashed', linewidth = 1, alpha = 0.45, c='k')

    # # plt.gcf().set_size_inches(7,5)
    plt.savefig("../results/" + model + "_Epistasis_Boxplots_3ab.pdf", format="pdf", bbox_inches='tight')

    # return 

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
def get_epsInfo(eps, df):
        # find row 
        # in that row access genotype category, genotype, inducer level
        # return dictionary
        i = np.isclose(df['Ep'], eps).argmax()
        output = []
        gen_cat = df.iloc[i]['genotype category']
        gen = df.iloc[i]['genotype']
        ind_lvl = df.iloc[i]['inducer level']
        output.append(gen_cat)
        output.append(gen)
        output.append(ind_lvl)
        return output

def compare_df(df1, df2):
        new_df = pd.DataFrame(list(zip(df1, df2)),columns =['old', 'new'])
        new_df['different'] = np.where(abs((new_df['old']/new_df['new'])-1)>0.2, 'yes', '-')
        return new_df
#%%

get_epsInfo(-0.350322053, New_Eps)