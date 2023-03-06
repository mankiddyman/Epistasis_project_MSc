from Models import *
from Model_Diagnostics_Functions import *

"""
a program that:
1) uses Epistasis_calc_functions to calculate mean epistasis for the models defined in Models.py
2) plots histograms of experimental and model epistasis
3) performs statistical tests 
"""
#----------------------------------------------------------------------------------------------------------
#%%
model_list = ['model_hill', 'model_hill_shaky', 'thermodynamic_model', 'compDeg']
data = pd.read_excel('../data/SM_params.xlsx') # change to model sm 


# hill = model_hill #1
# hill_shaky = model_hill_shaky #2
# td = thermodynamic_model #3
# comp = compDeg #4

# test
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
        print('You have selected the Hill model')
        cont = input('Do you want to select a new model? y/n')
        if cont == 'n':
            print('Diagnostic: epistasis comparison')
            get_histograms('observed', 'model_hill')
            print(pvalue1) 
            break
    elif model == '2':
        print('You have selected the Hill Shaky model')
        cont = input('Do you want to select a new model? y/n')
        if cont == 'n':
            print('Diagnostic: epistasis comparison')
            # get_Data('model_hill_shaky') 
            # get_Data('observed') 
            # get_histograms('observed', 'model_hill_shaky')
            # get_pvalue('observed', 'model_hill_shaky') 
            break
    elif model == '3':
        print('You have selected the Thermodynamic model')
        cont = input('Do you want to select a new model? y/n')
        if cont == 'n':
            print('Diagnostic: epistasis comparison')
            # get_Data('thermodynamic_model') 
            # get_Data('observed') 
            # get_histograms('observed', 'thermodynamic_model')
            # get_pvalue('observed', 'thermodynamic_model') 
            break
    elif model == '4':
        print('You have selected the Competition model')
        cont = input('Do you want to select a new model? y/n')
        if cont == 'n':
            print('Diagnostic: epistasis comparison')
            # get_Data('compDeg') 
            # get_Data('observed') 
            # get_histograms('observed', 'compDeg')
            # get_pvalue('observed', 'compDeg') 
            break
    model = input('Select a model: ')
    if model != '1' and model != '2' and model != '3':
        print('WRONG INPUT! PLEASE SELECT A MODEL')
        print('The following models are available')
        for i in model_list:
            print(model_list.index(i)+1, end=')')
            print(" ", i)
        model = input('Select a model: ')