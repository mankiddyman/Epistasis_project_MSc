import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as patches

observed = Figures(model = 'model_thermodynamic_sensor',split_pt = False) #defined below
indEps_observed = indComp('model_thermodynamic')

#set what colour you want to display the model in
observed_col = 'red'
model_hill_all_col = 'green'
model_therm_sens_col = 'navy'

#%% use this to define the function "Figures"
#Figures function gives figures of standard deviation and mean epistasis for mutants grouped by inducer conc and node mutated
def Figures(model= 'model_hill_all', split_pt = False, model_col:str = model_hill_all_col):
    #function to get Eps for a given model
    def get_modelEps(model:str = 'observed'):
        columns = [ 'genotype category', 'inducer level', 'genotype','Ep']
        df = pd.read_excel(f"../results/Eps_{model}.xlsx")[columns]
        df_Eps = pd.DataFrame(columns[:-2]+['std', 'node']).set_index(0).T
        #get all instances of each mutant in its own dataframe and group to find the mean and std
        for node in ['Sensor', 'Regulator', 'Output']:
            for i in range(1,11):
                mut = f"{node[0]}{i}"
                search = df.loc[((df['genotype'].str.startswith(mut + '_')) | (df['genotype'].str.endswith('_' + mut))| (df['genotype'].str.contains('_' + mut + '_')))].groupby(['genotype category', 'inducer level']) 
                df_mut_mean = search.mean().rename(columns={"Ep": "mean"})
                df_mut_std = search.std().rename(columns={"Ep": "std"})
                df_mut = pd.concat([df_mut_mean,df_mut_std], axis = 1).reset_index()
                #need to reorder inducer level column for plotting later
                df_mut['order'] = [3,1,2]*(int(len(df_mut)/3))
                df_mut = df_mut.sort_values(by = "order").drop('order', axis = 1)
                #add a label for the node
                df_mut['node'] = [node]*len(df_mut)
                #add mutant df to df with all mutants in
                df_Eps = pd.concat([df_Eps, df_mut])
        #add a column indicating the model
        df_Eps['model'] = [f"{model}"]*len(df_Eps)
        #prefer 'LMH' to 'inducer level' since only one word, same for genotype category
        return df_Eps.rename(columns = {'inducer level': 'LMH', 'genotype category': 'cat'}).reset_index(drop = True)

    df_Eps_obs = get_modelEps() #get obseved eps
    df_Eps_mod = get_modelEps(model)#get model data
    df_Eps = pd.concat([df_Eps_obs, df_Eps_mod])
    #stick node, and LMH columns together so that seaborn can group them properly
    df_Eps["desc"] = df_Eps["node"].map(str) + ", " + df_Eps["LMH"]

    #manually get the x tick positions and names
    x_ticks = []
    x_names = []
    for i in range(df_Eps["desc"].nunique()):
        adj = 0.25
        if i % 3 == 1:
            x_ticks.append(i)
        elif i % 3 == 2:
            x_ticks.append(i-adj)
        else:
            x_ticks.append(i+adj)
    x_names = list(df_Eps["LMH"].unique())*3
    sec_x_names = list(df_Eps["node"].unique())

    def Plotter(ax, data, y):
        palette ={"observed": observed_col, model: model_hill_all_col}
        sns.violinplot(ax = ax, data=data, x='desc', y=y, hue="model", split=True, dodge=False, palette=palette)
        #ax.xaxis.set_ticks(x_ticks)
        ax.set_xticklabels(x_names, fontsize=axis_size)
        ax.set_xlabel('')
        ax.add_patch(patches.Rectangle((-0.5,-100), 9/3, 1000, alpha = 0.1, color = 'r', zorder = 0.1)) 
        ax.add_patch(patches.Rectangle((-0.5+9/3,-100), 9/3, 1000, alpha = 0.1, color = 'g', zorder = 0.1)) 
        ax.add_patch(patches.Rectangle((-0.5+(9*2)/3,-100), 9/3, 1000, alpha = 0.1, color = 'b', zorder = 0.1)) 

    def set_axs(ax1, ax2):
        ax1.axhline(zorder = 0.5, c = 'darkgrey', linestyle = '--')
        #axis Ticks
        ax1.set_xticks([])
        secax = ax2.secondary_xaxis('bottom')
        secax.set_xticks(x_ticks[1::3], sec_x_names, fontsize = axis_size)
        secax.tick_params(pad=pad)
        #ax.set_xticklabels( names )
        for ticklabel, tickcolor in zip(secax.get_xticklabels(), ['r', 'g', 'b']):
            ticklabel.set_color(tickcolor)
        #axis labels
        ax1.set_ylabel(f"mean of $\epsilon$", fontsize = axis_size)
        ax2.set_ylabel(f"standard deviation of $\epsilon$", fontsize = axis_size)
        #legend
        ax1.get_legend().remove()
        ax2.get_legend().remove()

    if split_pt == False:
        title_size = 15
        axis_size = 14
        pad = 20
        Fig, (ax1, ax2) = plt.subplots(2, 1, figsize = (8.27,7))
        Plotter(ax1, data = df_Eps, y= "mean")
        Plotter(ax2, data = df_Eps, y= "std")
        set_axs(ax1, ax2)
        file_path = "../results/"+model+"Ep_compare_pt.jpg"

    else:
        title_size = 30
        axis_size = 25
        pad = 40
        Fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize = (16,14))
        Plotter(ax1, data = df_Eps.loc[df_Eps.cat == "pairwise"], y= "mean")
        Plotter(ax3, data = df_Eps.loc[df_Eps.cat == "pairwise"], y= "std")
        Plotter(ax2, data = df_Eps.loc[df_Eps.cat == "triplet"], y= "mean")
        Plotter(ax4, data = df_Eps.loc[df_Eps.cat == "triplet"], y= "std")
        set_axs(ax1, ax3)
        set_axs(ax1 = ax2, ax2 = ax4)
        ax2.set_ylabel('')
        ax4.set_ylabel('')
        ax3.set_xticklabels(['L', 'M', 'H']*3)
        ax4.set_xticklabels(['L', 'M', 'H']*3)
        file_path = "../results/"+model+"Ep_compare_pt.jpg"
      
    handles, labels = ax2.get_legend_handles_labels()
    Fig.legend(handles, labels, loc="upper left", bbox_to_anchor=(0.9, 0.9), fontsize = axis_size)
    Fig.suptitle(f"Mean and Varience of Calculated versus Obseved Epistasis for \n{model.split('_')[1]} model with sampling strategy '{model.split('_')[2]}'".title(), fontsize = title_size)
    
    Fig.savefig(file_path)
    return Fig 


#plot comparing epistasis at different inducer concs
#first, get a dataframe with Ep value at high, hedium and low I conc on one row
def indComp(model='observed'):
    df = pd.read_excel('../results/Eps_'+str(model)+'.xlsx', index_col=0).set_index('genotype')
    #reorganise data so that epistasis at low, medium and high I concs are easily compared
    Ep_medium = pd.DataFrame(df['Ep'].loc[(df['inducer level']=='medium')]).rename(columns ={'Ep':'Ep_medium'})
    Ep_high = pd.DataFrame(df['Ep'].loc[(df['inducer level']=='high')]).rename(columns={'Ep':'Ep_high'})
    df = df[['Ep', 'genotype category']].loc[(df['inducer level']=='low')]
    df = df.rename(columns={'Ep':'Ep_low'})
    df = df.join(Ep_high.join(Ep_medium))

    #get proportions that are in each quadrant of the plot using the following dataframes
    df_up = df[((df['Ep_low'] < df['Ep_medium']) & (df['Ep_medium'] < df['Ep_high']))]
    df_down = df[(df['Ep_low'] > df['Ep_medium']) & (df['Ep_medium'] > df['Ep_high'])]
    df_peak = df[(df['Ep_low'] < df['Ep_medium']) & (df['Ep_medium'] > df['Ep_high'])]
    df_trough = df[(df['Ep_low'] > df['Ep_medium']) & (df['Ep_medium'] < df['Ep_high'])]
    #totals:
    n,n_pair, n_triplet = len(df),len(df[df['genotype category']=='pairwise']), len(df[df['genotype category']=='triplet'])

    n_up = [len(df_up), len(df_up[df_up['genotype category']=='pairwise']), len(df_up[df_up['genotype category']=='triplet'])]

    n_down = [len(df_down), len(df_down[df_down['genotype category']=='pairwise']), len(df_down[df_down['genotype category']=='triplet'])]

    n_peak = [len(df_peak), len(df_peak[df_peak['genotype category']=='pairwise']), len(df_peak[df_peak['genotype category']=='triplet'])]
    n_trough = [len(df_trough), len(df_trough[df_trough['genotype category']=='pairwise']), len(df_trough[df_trough['genotype category']=='triplet'])]

    # % calculator
    def percent(x,n):
        return np.round(x*100/n,0)

    #now plot a scatter plot of med-low vs high-med
    x = df.Ep_medium - df.Ep_low
    y= df.Ep_high - df.Ep_medium
    cols = np.where(df['genotype category'] == 'pairwise', 'purple', 'orange')

    fig1, ax = plt.subplots()
    scatter = ax.scatter(x, y, c = cols, zorder = 2.5, edgecolor = 'black', alpha = 0.65, label = ['pairwise', 'triple'])
    ax.axhline(y=0, c = 'darkgray', ls = '--')
    ax.axvline(x=0, c = 'darkgray', ls = '--')
    xlim = max(abs(min(x)), abs(max(x)))
    ylim = max(abs(min(y)), abs(max(y)))
    ax.set_xlim([-xlim*1.05,xlim*1.05])
    ax.set_ylim([-ylim*1.05,ylim*1.05])
    ax.set_xlabel('$\u03B5_{medium}$ - $\u03B5_{low}$')
    ax.set_ylabel('$\u03B5_{high}$ - $\u03B5_{medium}$')
    ax.set_title('inducer dependant Epistasis for '+ str(model)+' data')

    #proportins mutants in each quadrant
    def annotate(direction):
        if direction == n_up:
            a, b, c = 1, 1, 0.7 #a & b position %s horizontally & vertically, c adjusts pair/trip %s
        elif direction == n_down:
            a, b, c = -1,-1, 0.9
        elif direction == n_peak:
            a, b, c = 1,-1, 0.9
        else:
            a, b, c = -1,1,0.7 
        x_pos, y_pos = a*xlim,b*ylim
        ax.text(x_pos*0.8,y_pos*0.8 , str(percent(direction[0],n))+'%', verticalalignment='center', horizontalalignment='center', size = 20)
        ax.text(x_pos*0.8, y_pos*c, str(percent(direction[1],n_pair))+'%', verticalalignment='top', horizontalalignment='right', size = 9, c = 'purple')
        ax.text(x_pos*0.8, y_pos*c, str(percent(direction[2],n_triplet))+'%', verticalalignment='top', horizontalalignment='left', size = 9, c = 'orange')
    annotate(n_up)
    annotate(n_down)
    annotate(n_peak)
    annotate(n_trough)
    #and a basic legend
    ax.text(-xlim, -ylim*1.25, 'Total, n = '+str(n), verticalalignment='top', horizontalalignment='left', size = 15, c= 'k')
    ax.text(-xlim, -ylim*1.36, 'Pairwise, n = '+str(n_pair), verticalalignment='top', horizontalalignment='left', size = 15, c= 'purple')
    ax.text(-xlim, -ylim*1.46, 'Triple, n = '+str(n_triplet), verticalalignment='top', horizontalalignment='left', size = 15, c= 'orange')
    plt.show()

    #plot a key for the inducer dependence plot
    x = np.array([1/6, 0.5, 5/6])
    y_peak = np.array([0.3, 0.7, 0.3])
    y_up = np.array([0.2, 0.5, 0.8])
    lim = (-1.05,1.05)

    fig2, ((ax1), (ax2)) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [5, 7]}, figsize = (3, 5))
    def plot(ax, x,y):
        return ax.plot(x,y, marker= 'o', markerfacecolor= 'lightgray',c= 'k', linestyle='-')
    plot(ax2, -1*x,y_peak)
    plot(ax2, x, -1*y_peak)
    plot(ax2, x, y_up)
    plot(ax2, -1*x, -1*np.flip(y_up))
    plot(ax1, x, y_up)
    plot(ax2, [-1/3, 0, 1/3], [0,0,0])
    ax2.axhline(zorder = 0.5, c = 'lightgrey', linestyle = '--')
    ax2.axvline(zorder = 0.5, c = 'lightgrey', linestyle = '--')
    ax2.set(xlim= lim, ylim = lim,  xticks= [0], yticks=[0], xlabel = '$\epsilon_{medium} - \epsilon_{low}$', ylabel = '$\epsilon_{high} - \epsilon_{medium}$')
    ax1.set(xlim = (0,1), ylim = (0,1))
    ax1.set_axis_off()
    ax1.text(0.499, 0.5, '{', verticalalignment='bottom', horizontalalignment='right', size = 35, c= 'k')
    ax1.text(0.37, 0.65, '$\epsilon_{high} - \epsilon_{medium}$', verticalalignment='center', horizontalalignment='right', size = 15, c= 'k')
    ax1.text(0.51, 0.5, '}', verticalalignment='top', horizontalalignment='left', size = 35, c= 'k')
    ax1.text(0.63, 0.32, '$\epsilon_{medium} - \epsilon_{low}$', verticalalignment='center', horizontalalignment='left', size = 15, c= 'k')
    ax2.text(0.4,0, 'inducer \nindependence', size = '7', verticalalignment = 'center')
    fig2.suptitle('Inducer dependence of $\epsilon$')

    #now to get totals and proportions for more specific plots
    #lmh indicates Ep < 0 at low, medium and high I conc
    up_lmh = df_up[(df_up.Ep_low < 0) & ( df_up.Ep_medium < 0) & (df_up.Ep_high < 0)] 
    up_lm = df_up[(df_up.Ep_low < 0) & ( df_up.Ep_medium < 0) & (df_up.Ep_high > 0)] 
    up_l = df_up[(df_up.Ep_low < 0) & ( df_up.Ep_medium > 0) & (df_up.Ep_high > 0)] 
    up_over = df_up[(df_up.Ep_low > 0) & ( df_up.Ep_medium > 0) & (df_up.Ep_high > 0)] 

    peak_lmh = df_peak[(df_peak.Ep_low < 0) & ( df_peak.Ep_medium < 0) & (df_peak.Ep_high < 0)] 
    peak_lm = df_peak[(df_peak.Ep_low < 0) & ( df_peak.Ep_medium < 0) & (df_peak.Ep_high > 0)] 
    peak_l = df_peak[(df_peak.Ep_low < 0) & ( df_peak.Ep_medium > 0) & (df_peak.Ep_high > 0)] 
    peak_over = df_peak[(df_peak.Ep_low > 0) & ( df_peak.Ep_medium > 0) & (df_peak.Ep_high > 0)] 

    trough_lmh = df_trough[(df_trough.Ep_low < 0) & ( df_trough.Ep_medium < 0) & (df_trough.Ep_high < 0)] 
    trough_lm = df_trough[(df_trough.Ep_low < 0) & ( df_trough.Ep_medium < 0) & (df_trough.Ep_high > 0)] 
    trough_l = df_trough[(df_trough.Ep_low < 0) & ( df_trough.Ep_medium > 0) & (df_trough.Ep_high > 0)] 
    trough_over = df_trough[(df_trough.Ep_low > 0) & ( df_trough.Ep_medium > 0) & (df_trough.Ep_high > 0)] 

    down_lmh = df_down[(df_down.Ep_low < 0) & ( df_down.Ep_medium < 0) & (df_down.Ep_high < 0)] 
    down_lm = df_down[(df_down.Ep_low < 0) & ( df_down.Ep_medium < 0) & (df_down.Ep_high > 0)] 
    down_l = df_down[(df_down.Ep_low < 0) & ( df_down.Ep_medium > 0) & (df_down.Ep_high > 0)] 
    down_over = df_down[(df_down.Ep_low > 0) & ( df_down.Ep_medium > 0) & (df_down.Ep_high > 0)] 
    return fig1, fig2

#plot significant epistasis
def Sig_Ep(model='observed'):
    df_M = pd.read_excel('../results/Eps_'+str(model)+'.xlsx', index_col=0)
    cols = np.where(df_M['Sig_Epistasis'] == True, 'darkgray', np.where(df_M['genotype category'] == 'pairwise', 'palegreen', 'lightskyblue'))
    alphas = np.where(df_M['Sig_Epistasis'] == True, 0.5, 0.8)
    fig, ax = plt.subplots()
    ax.scatter(df_M['Ep_pVal'], df_M['Ep'], c = cols,alpha = alphas, zorder = 2.5, edgecolor = 'black' )
    ax.axvline(0.05, c = 'darkgray', ls = '--')
    plt.show
    return fig
    
#%%
