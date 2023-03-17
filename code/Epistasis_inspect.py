import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb
from Chi_Figure3_func import Sort_mutants, convertDF


Figures() #defined below
indComp('model_hill')

#%% use this to define the function "Figures"
#Figures function gives figures of standard deviation and mean epistasis for mutants grouped by inducer conc and node mutated
def Figures(model= 'observed'):
    #first, get dataframe with epistasis values for each mutant
    df = Sort_mutants(model)
    df1 = convertDF(df[0], 'Output')
    df2 = convertDF(df[1], 'Regulator')
    df3 = convertDF(df[2], 'Sensor')
    #group df and calculate mean and varience for each mutation, grouped by inducer conc, node and pairwise/triplet category
    df_Eps = pd.concat([pd.concat([df1, df2], axis = 0), df3], axis = 0).reset_index().drop('index', axis = 1)
    df_Eps['node'] = np.where(df_Eps['Genotype'].str.startswith('S'), 'Sensor', np.where(df_Eps['Genotype'].str.startswith('O'), 'Output', 'Regulator'))
    #reorder categories; want low then med then high on the plot
    node_order = ['Sensor', 'Regulator', 'Output']
    LMH_order = ['low', 'medium', 'high']
    #blank dataframe to fill with reordered values
    df_Eps1 = pd.DataFrame(df_Eps.columns.values).set_index(0).T
    for i in node_order:
        dfi = df_Eps.loc[df_Eps['node'] == i]
        for j in LMH_order:
            dfij = dfi.loc[dfi['LMH'] == j]
            df_Eps1 = pd.concat([df_Eps1, dfij], ignore_index=True)
    df_Eps = df_Eps1

    #get means and standard deviation for each mutant
    df_means = df_Eps.groupby(['LMH','node', 'Genotype' , 'Category'], sort = False).mean().reset_index()
    df_vars = df_Eps.groupby(['LMH', 'node','Genotype', 'Category'], sort = False).std().reset_index()

    #split pairwise and triplet, means and std
    df_means_p = df_means.loc[df_means['Category']== 'pairwise']
    df_means_t = df_means.loc[df_means['Category']== 'triplet']
    df_vars_p = df_vars.loc[df_means['Category']== 'pairwise']
    df_vars_t = df_vars.loc[df_means['Category']== 'triplet']
    to_plots = [df_means_p, df_means_t, df_vars_p, df_vars_t]

    group = ['node','LMH']
    col_choice = ['r', 'g', 'b']
    #group data by node and inducer cons ready to be plotted on a figure
    def Plot_setup(df:pd.DataFrame):
        grouped = df.groupby(group, sort=False)
        names, vals, xs, bar_pos, cols = [], [],[], [], []
        for i, (name, subdf) in enumerate(grouped):
            names.append(name[1])
            vals.append(subdf['Epistasis'].tolist())
            bar_pos += [i+1]
            xs.append(np.random.normal(i+1, 0.04, subdf.shape[0]))
            if i < 3:
                cols += [[col_choice[0]]*subdf.shape[0]]
            elif i < 6:
                cols += [[col_choice[1]]*subdf.shape[0]]
            else:
                cols += [[col_choice[2]]*subdf.shape[0]]
        #posttions for boxplots and x values for scatter plots
        for i in range(3):
            adj = 0.3
            xs[3*i] += adj
            bar_pos[3*i] += adj
            xs[3*i+2] -= adj
            bar_pos[3*i+2] -= adj
        return names, vals, xs, bar_pos, cols

    #defina a plotting function
    def Plotter(ax, names, vals, xs, bar_pos, cols,isMean = False, triple = False):
        #ax.boxplot(vals,labels=None, showfliers = False, positions = bar_pos, meanline = True, whis = 0,showcaps = False,  medianprops= {'c':'black'})
        if triple == True:
            title = 'triple'
            title_col = 'orange'
        else:
            title = 'pairwise'
            title_col = 'purple'
        ax.set_title(title, c = title_col) 
        for x, val, col in zip(xs, vals, cols):
            ax.scatter(x, val, alpha=0.4, c = col)
            ax.errorbar(np.mean(x) , np.mean(val), np.std(val), linestyle='None', fmt='_', c = 'k', elinewidth = 1.5, capsize = 1.5)
        if isMean == True:
            ax.axhline(zorder = 0.5, c = 'lightgrey', linestyle = '--')
            ax.set_ylim(min(df_means.Epistasis)*1.05, max(df_means.Epistasis)*1.05)
        else:
            ax.set_ylim(0, max(df_vars.Epistasis)*1.05)
        ax.set_xticklabels(names,rotation=45, fontsize=8)
        secax = ax.secondary_xaxis('bottom')
        secax.set_xticks(bar_pos[1::3], node_order)
        secax.tick_params(pad=40)
        #ax.set_xticklabels( names )
        for ticklabel, tickcolor in zip(secax.get_xticklabels(), col_choice):
            ticklabel.set_color(tickcolor)
        return ax

    #Define the figures, mean first,
    Fig_obs_Mean, (means_p, means_t) = plt.subplots(1,2)
    mean_P = Plot_setup(df_means_p)
    Plotter(means_p, *mean_P, True)
    mean_T = Plot_setup(df_means_t)
    Plotter(means_t,*mean_T, True, triple = True)
    means_p.set_ylabel('mean $\epsilon$')
    Fig_obs_Mean.suptitle(str(model))
    plt.show()

    #then varience
    Fig_obs_Var, (vars_p, vars_t) = plt.subplots(1,2)
    var_P = Plot_setup(df_vars_p)
    Plotter(vars_p, *var_P)
    var_T = Plot_setup(df_vars_t)
    Plotter(vars_t,*var_T, triple = True)
    vars_p.set_ylabel('standard deviation of $\epsilon$')
    Fig_obs_Var.suptitle(str(model))
    plt.show()
    return Fig_obs_Mean, Fig_obs_Var


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

    x = np.array([1/6, 0.5, 5/6])
    y_peak = np.array([0.3, 0.7, 0.3])
    y_up = np.array([0.2, 0.5, 0.8])
    lim = (-1.05,1.05)

    fig2, ((ax1), (ax2)) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [1, 3]}, figsize = (2, 5))
    def plot(ax, x,y):
        return ax.plot(x,y, marker= 'o', markerfacecolor= 'k',c= 'gray', linestyle='-')
    plot(ax2, -1*x,y_peak)
    plot(ax2, x, -1*y_peak)
    plot(ax2, x, y_up)
    plot(ax2, -1*x, -1*np.flip(y_up))
    plot(ax1, x, y_up)
    ax2.axhline(zorder = 0.5, c = 'lightgrey', linestyle = '--')
    ax2.axvline(zorder = 0.5, c = 'lightgrey', linestyle = '--')
    ax2.set(xlim= lim, ylim = lim,  xticks= [], yticks=[], xlabel = '$\epsilon_{medium} - \epsilon_{low}$', ylabel = '$\epsilon_{high} - \epsilon_{medium}$')
    ax1.set(xlim = (0,1), ylim = (0,1), xticks= [], yticks=[])
    ax1.set_axis_off()


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
    return fig1

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

# %%
