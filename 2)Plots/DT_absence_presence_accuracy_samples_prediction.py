#%%
# Source: https://hands-on.cloud/implementation-of-support-vector-machine-svm-using-python/
# importing required libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D 
import random
import os

# sklearn imports 
from sklearn.metrics import accuracy_score,confusion_matrix
from sklearn import model_selection
from sklearn.tree import DecisionTreeClassifier
from sklearn.tree import _tree

# =============================================
# ---------- EXPORT in Rules Format -----------
# =============================================

def get_rules(tree, feature_names, class_names):
    tree_ = tree.tree_
    feature_name = [
        feature_names[i] if i != _tree.TREE_UNDEFINED else "undefined!"
        for i in tree_.feature
    ]

    paths = []
    path = []
    
    def recurse(node, path, paths):
        
        if tree_.feature[node] != _tree.TREE_UNDEFINED:
            name = feature_name[node]
            threshold = tree_.threshold[node]
            p1, p2 = list(path), list(path)
            p1 += [f"({name} <= {np.round(threshold, 3)})"]
            recurse(tree_.children_left[node], p1, paths)
            p2 += [f"({name} > {np.round(threshold, 3)})"]
            recurse(tree_.children_right[node], p2, paths)
        else:
            path += [(tree_.value[node], tree_.n_node_samples[node])]
            paths += [path]
            
    recurse(0, path, paths)

    # sort by samples count
    samples_count = [p[-1][1] for p in paths]
    ii = list(np.argsort(samples_count))
    paths = [paths[i] for i in reversed(ii)]
    
    rules = []
    for path in paths:
        rule = "if "
        
        for p in path[:-1]:
            if rule != "if ":
                rule += " and "
            rule += str(p)
        rule += " then "
        if class_names is None:
            rule += "response: "+str(np.round(path[-1][0][0][0],3))
        else:
            classes = path[-1][0][0]
            l = np.argmax(classes)
            rule += f"class: {class_names[l]} (proba: {np.round(100.0*classes[l]/np.sum(classes),2)}%)"
        rule += f" | based on {path[-1][1]:,} samples"
        rules += [rule]
        
    return rules

def DT_Features(rules):

    rules = str(rules)
    # create list of OTUs in rules
    splitted = rules.split(')')
    features = []

    # get all mentions of OTUs
    for i in range(0,len(splitted)-1):
        x = splitted[i]
        x1 = x.split('(')[1]
        x2 = x1.split(' ')[0]
        
        if 'OTU_' in x2:
            features.append(x2)

    # count frequency of OTUs in list
    frequency = {}
    for word in features:

        count = frequency.get(word,0)
        frequency[word] = count + 1

    # open file for writing, "w" is writing
    w = []

    # loop over dictionary keys and values
    for key, val in frequency.items():
        # write every key and value to file
        w.append([key, val])

    w = pd.DataFrame(w)
    w.columns = ['OTU','value']
    # sort on descending order
    # sort = pd.read_csv(name, names=['OTU','value'])
    sort = w.sort_values(by='value', ascending=False)

    return sort

def zonation (values, threshold):
    # ============
    # Give labels
    # ============
    Zonation = [2]*len(values)

    for i,value in enumerate(values):

        if np.isnan(value):
            Zonation[i] = float('nan')
        elif value > threshold:
            Zonation[i] = 1
        else:
            Zonation[i] = 0
    
    return Zonation


variables = ['O2conc','NOxconc','Mnconc','NH4conc']
thresholds = [0.0,0.003,0.005,0.010,0.015]
thresholds = np.array(thresholds)
accuracy = [[float('nan'),float('nan'),float('nan'),float('nan'),float('nan')],
            [float('nan'),float('nan'),float('nan'),float('nan'),float('nan')],
            [float('nan'),float('nan'),float('nan'),float('nan'),float('nan')],
            [float('nan'),float('nan'),float('nan'),float('nan'),float('nan')]]

accuracy_AMOR_pres = np.array(accuracy)
accuracy_AMOR_abs = np.array(accuracy)
accuracy_WNAG_pres = np.array(accuracy)
accuracy_WNAG_abs = np.array(accuracy)
accuracy_MAR_pres = np.array(accuracy)
accuracy_MAR_abs = np.array(accuracy)
accuracy_NWPO_pres = np.array(accuracy)
accuracy_NWPO_abs = np.array(accuracy)
accuracy_SCS_pres = np.array(accuracy)
accuracy_SCS_abs = np.array(accuracy)
accuracy_SPO_pres = np.array(accuracy)
accuracy_SPO_abs = np.array(accuracy)

number_AMOR_pres = np.array(accuracy)
number_AMOR_abs = np.array(accuracy)
number_WNAG_pres = np.array(accuracy)
number_WNAG_abs = np.array(accuracy)
number_MAR_pres = np.array(accuracy)
number_MAR_abs = np.array(accuracy)
number_NWPO_pres = np.array(accuracy)
number_NWPO_abs = np.array(accuracy)
number_SCS_pres = np.array(accuracy)
number_SCS_abs = np.array(accuracy)
number_SPO_pres = np.array(accuracy)
number_SPO_abs = np.array(accuracy)

for i, var in enumerate(variables):

    for j, thres in enumerate(thresholds):
         
        print("========================== \n")       
        print(f"Variable is: {var} \n")
        print(f"Threshold is: {thres} \n")
        print("========================== \n")    
        # ===================
        # Data selection
        # ===================
        data = pd.read_csv(os.getcwd()+"/../1)Preprocessing/4)Significant_gene/Data_prep_ML_FAMILY_datacomparison.csv",
                delimiter=",", dtype=str, index_col=None)
        data.columns = data.columns.str.replace(data.columns[0], 'Core')
        # data['Core'] = data['Core'].astype(float).astype(int)
        # data.reset_index(drop=True, inplace=True)
        data = data.astype(float)
        data['Core'] = data['Core'].astype(int)
        data = data.dropna(subset = [var], how ='any')
        data.reset_index(drop=True, inplace=True) 
        data['zonation'] = zonation(data[var].values,thres)
        data['zonation']


        # ==================================================
        # ---------- REMOVE sample in TRANSITION -----------
        # ==================================================
        predictor_loc = np.where(data.columns == 'zonation')[0][0]

        # remove cores that are not of interest (GS19GC25, GS21GC09, GS20GC20, GS20GC21)
        data = data[data['Core'] != 18]
        data = data[data['Core'] != 19]
        data = data[data['Core'] != 22]
        data = data[data['Core'] != 21]

        data.reset_index(drop=True, inplace=True)

        # select relevant features according to q-value
        # # ------------------------------------------------------
        feature_selection1 = np.array(pd.read_csv(os.getcwd()+'/../1)Preprocessing/4)Significant_gene/Significant_Features_FAMILY_comparison_AMOR_WNAG_MAR_SCS_NWPO_SPO_oxygen.csv',header=None))
        feature_selection2 = np.array(pd.read_csv(os.getcwd()+'/../1)Preprocessing/4)Significant_gene/Significant_Features_FAMILY_comparison_AMOR_WNAG_MAR_SCS_NWPO_SPO_nitrate.csv',header=None))
        feature_selection3 = np.array(pd.read_csv(os.getcwd()+'/../1)Preprocessing/4)Significant_gene/Significant_Features_FAMILY_comparison_AMOR_WNAG_MAR_SCS_NWPO_SPO_manganese.csv',header=None))
        feature_selection4 = np.array(pd.read_csv(os.getcwd()+'/../1)Preprocessing/4)Significant_gene/Significant_Features_FAMILY_comparison_AMOR_WNAG_MAR_SCS_NWPO_SPO_ammonium.csv',header=None))

        feature_selection = np.concatenate([feature_selection1,feature_selection2,feature_selection3,feature_selection4])
        feature_selection = np.unique(feature_selection)

        loc = []
        for f in range(len(feature_selection)):
            
            condition = data.columns == feature_selection[f]
            selection = np.where(condition)
            loc.append(selection[0][0])

        print('Feature selection was all present in dataset: '+ str(len(feature_selection)==len(loc)))

        # feature selection procedure completion
        # --------------------------------------
        feature_selection = data.columns[loc]

        # ============================
        # ---- Look up taxa ----
        #============================

        lvl = 'FAMILY'
        TAX = pd.read_csv(os.getcwd()+"/../1)Preprocessing/3)Merge_OTUtables/TAXtable_merged_family.csv")
        features = [np.array(TAX[TAX.iloc[:,0]==OTU][lvl])[0] for OTU in feature_selection]

        feature_names_order = [x.strip() for x in features]
        len(feature_names_order)

        # ===================================
        # ---------- BALANCE DATA -----------
        # ===================================

        # core GS19GC10 (core == 17) represent the dataset the most (207 sample of 658)
        # here I will separate part of the data
        # take 2/3 of the data
        data.reset_index(drop=True, inplace=True)

        list_GS19GC10 = data[data['Core'] == 17].index
        list_random = []
        for c in range(int(list_GS19GC10.shape[0]/4)):
            list_random.append(random.choice(list_GS19GC10))

        X_seperation = data.iloc[list_random,loc].values.astype(float)
        y_seperation = data.iloc[list_random,predictor_loc]
        y_seperation = np.array(y_seperation)
        core_seperation = np.array(data.iloc[list_random]['Core'])

        # ================================
        # Assigning X and Y variables 
        # ================================
        data_temp = data.drop(list_random)
        data_temp.reset_index(drop=True, inplace=True)

        #Extracting Independent and dependent Variable  
        X = data_temp.iloc[:, loc].values.astype(float)
        y = data_temp.iloc[:, predictor_loc]
        y = np.array(y)
        core = np.array(data_temp['Core'])


        # =================================
        # ---------- SPLIT DATA -----------
        # =================================

        seed = 650                    # Fix random seed for reproducibility
        # Shuffle and split the data into train and a 
        # concatenation of validation and test sets with a ratio of 0.7/0.3
        X_train, X_test, y_train, y_test, core_train, core_test = model_selection.train_test_split(
            X, y, core,
            test_size=0.3, shuffle=True, random_state=seed
        )

        X_test = np.vstack([X_test, X_seperation])
        y_test = np.append(y_test, y_seperation)
        core_test = np.append(core_test,core_seperation)

        sum(core==17)/X.shape[0]
        sum(core_train==17)/X_train.shape[0]
        sum(core_test==17)/X_test.shape[0]


        # store pacific in another variable
        # and remove 26 - 39 from test Data
        pacific = np.arange(27,40)

        pacific_data = np.where(core_test == 26)[0]
        for core in pacific:
            if core in core_test:
                pacific_data = np.concatenate([pacific_data,np.where(core_test == core)[0]])

        X_pacific = X_test[pacific_data,:]
        y_pacific = y_test[pacific_data]

        # store atlantic in another variable (40)
        atlantic_data = np.where((core_test == 40) | (core_test == 41) | (core_test == 42))[0]
        X_atlantic = X_test[atlantic_data,:]
        y_atlantic = y_test[atlantic_data]

        # store north pond in another variable (41 and 42)
        NP_data = np.where((core_test == 43) | (core_test == 44))[0]
        X_NP = X_test[NP_data,:]
        y_NP = y_test[NP_data]

        # store south china sea in another variable (45 up to 55)
        SCS = np.arange(46,56)

        SCS_data = np.where(core_test == 45)[0]
        for core in SCS:
            if core in core_test:
                SCS_data = np.concatenate([SCS_data,np.where(core_test == core)[0]])

        X_SCS = X_test[SCS_data,:]
        y_SCS = y_test[SCS_data]

        # store South Pacific Ocean in another variable (55 up to 67)
        SPO = np.arange(56,68)

        SPO_data = np.where(core_test == 55)[0]
        for core in SPO:
            if core in core_test:
                SPO_data = np.concatenate([SPO_data,np.where(core_test == core)[0]])

        X_SPO = X_test[SPO_data,:]
        y_SPO = y_test[SPO_data]

        # AMOR
        AMOR = np.arange(2,26)
        AMOR_data = np.where(core_test == 1)[0]

        for core in AMOR:
            if core in core_test:
                AMOR_data = np.concatenate([AMOR_data,np.where(core_test == core)[0]])


        X_AMOR = X_test[AMOR_data,:]
        y_AMOR = y_test[AMOR_data]


        iterations = 1

        for it in range(iterations):
            # ===========================================
            # ---------- TRAIN AND TEST MODEL -----------
            # ===========================================

            # ------------------------------------
            # Test choosen model on test dataset
            # ------------------------------------
            print("pre training")
            alpha = 0
            max_dep = 500

            dtree_sklearn = DecisionTreeClassifier(criterion='gini', max_depth=max_dep, ccp_alpha=alpha)
            dtree_sklearn.fit(X_train,y_train)  

            print("============ Train data ============")
            y_pred_train = dtree_sklearn.predict(X_train)

            print("============ All together data ============")
            y_pred_test = dtree_sklearn.predict(X_test)
            print("after training")

        # -------------
        # AMOR
        # -------------
        print("AMOR")
        P = np.where(y_AMOR == 1)
        A = np.where(y_AMOR == 0)

        y_AMOR_pred = dtree_sklearn.predict(X_AMOR)

        if len(A[0]) > 0:
            acc_absent = accuracy_score(y_AMOR[A], y_AMOR_pred[A])
            accuracy_AMOR_abs[i,j] = acc_absent
            number_AMOR_abs[i,j] = len(y_AMOR[A])
        if len(P[0]) > 0:
            acc_present = accuracy_score(y_AMOR[P], y_AMOR_pred[P])
            accuracy_AMOR_pres[i,j] = acc_present
            number_AMOR_pres[i,j] = len(y_AMOR[P])
        
        # -------------
        # Atlantic (WNAG)
        # -------------
        print("WNAG")
        if var == 'O2conc':
            P = np.where(y_atlantic == 1)
            A = np.where(y_atlantic == 0)

            y_atlantic_pred = dtree_sklearn.predict(X_atlantic)

            acc_absent = accuracy_score(y_atlantic[A], y_atlantic_pred[A])
            acc_present = accuracy_score(y_atlantic[P], y_atlantic_pred[P])

            accuracy_WNAG_pres[i,j] = acc_present
            accuracy_WNAG_abs[i,j] = acc_absent
            number_WNAG_pres[i,j] = len(y_atlantic[P])
            number_WNAG_abs[i,j] = len(y_atlantic[A])

        # -------------
        # North Pond (MAR)
        # -------------
        print("MAR")
        if (var == 'O2conc') | (var == 'NOxconc'):
            P = np.where(y_NP == 1)
            A = np.where(y_NP == 0)

            y_NP_pred = dtree_sklearn.predict(X_NP)
            if len(A[0]) > 0:
                acc_absent = accuracy_score(y_NP[A], y_NP_pred[A])
                accuracy_MAR_abs[i,j] = acc_absent
                number_MAR_abs[i,j] = len(y_NP[A])
            if len(P[0]) > 0:
                acc_present = accuracy_score(y_NP[P], y_NP_pred[P])
                accuracy_MAR_pres[i,j] = acc_present
                number_MAR_pres[i,j] = len(y_NP[P])
            
        # -------------
        # Pacific (NWPO)
        # -------------
        print("NWPO")
        if (var == 'O2conc') | (var == 'NOxconc') | (var == 'NH4conc'):
            P = np.where(y_pacific == 1)
            A = np.where(y_pacific == 0)

            y_pacific_pred = dtree_sklearn.predict(X_pacific)

            if len(A[0]) > 0:
                acc_absent = accuracy_score(y_pacific[A], y_pacific_pred[A])
                accuracy_NWPO_abs[i,j] = acc_absent
                number_NWPO_abs[i,j] = len(y_pacific[A])
            if len(P[0]) > 0:
                acc_present = accuracy_score(y_pacific[P], y_pacific_pred[P])
                accuracy_NWPO_pres[i,j] = acc_present
                number_NWPO_pres[i,j] = len(y_pacific[P])
   
        # SCS
        # -------------
        print("SCS")
        P = np.where(y_SCS == 1)
        A = np.where(y_SCS == 0)

        y_SCS_pred = dtree_sklearn.predict(X_SCS)
        if len(A[0]) > 0:
            acc_absent = accuracy_score(y_SCS[A], y_SCS_pred[A])
            accuracy_SCS_abs[i,j] = acc_absent
            number_SCS_abs[i,j] = len(y_SCS[A])
        if len(P[0]) > 0:
            acc_present = accuracy_score(y_SCS[P], y_SCS_pred[P])
            accuracy_SCS_pres[i,j] = acc_present
            number_SCS_pres[i,j] = len(y_SCS[P])

   
        # SPO
        # -------------
        print("SPO")
        P = np.where(y_SPO == 1)
        A = np.where(y_SPO == 0)

        y_SPO_pred = dtree_sklearn.predict(X_SPO)
        if len(A[0]) > 0:
            acc_absent = accuracy_score(y_SPO[A], y_SPO_pred[A])
            accuracy_SPO_abs[i,j] = acc_absent
            number_SPO_abs[i,j] = len(y_SPO[A])
        if len(P[0]) > 0:
            acc_present = accuracy_score(y_SPO[P], y_SPO_pred[P])
            accuracy_SPO_pres[i,j] = acc_present
            number_SPO_pres[i,j] = len(y_SPO[P])
        
# fig.delaxes(ax[2][1])
# plt.tight_layout()
# plt.show()

fig, ax = plt.subplots(2,3,figsize=(20, 10))

locations = ['AMOR','WNAG','MAR','NWPO',"SCS","SPO"]
col = ['orange','red','darkgreen','purple']
variables = [r'$O_2$',r'$NO_x$',r'$Mn$',r'$NH_4$']
r = 0
c = 0

accuracy_abs = [accuracy_AMOR_abs,accuracy_WNAG_abs,accuracy_MAR_abs,accuracy_NWPO_abs,accuracy_SCS_abs,accuracy_SPO_abs]
accuracy_pres = [accuracy_AMOR_pres,accuracy_WNAG_pres,accuracy_MAR_pres,accuracy_NWPO_pres,accuracy_SCS_pres,accuracy_SPO_pres]
number_abs = [number_AMOR_abs,number_WNAG_abs,number_MAR_abs,number_NWPO_abs,number_SCS_abs,number_SPO_abs]
number_pres = [number_AMOR_pres,number_WNAG_pres,number_MAR_pres,number_NWPO_pres,number_SCS_pres,number_SPO_pres]

for l, loc in enumerate(locations):
    ax[r,c].grid(axis='both')
    for v, var in enumerate(variables):
        ax[r,c].plot(thresholds*1000, accuracy_pres[l][v,:], color = col[v],linestyle='-')
        ax[r,c].scatter(thresholds*1000, accuracy_pres[l][v,:], color = col[v])
        ax[r,c].plot(thresholds*1000, accuracy_abs[l][v,:], color = col[v], linestyle='--')
        ax[r,c].scatter(thresholds*1000, accuracy_abs[l][v,:], color = col[v],  label = var)

    ax[r,c].set_title(loc, fontweight='bold')
    ax[r,c].set_ylim(-0.05,1.05)
    ax[r,c].set_xlabel('Threshold (μM)')
    ax[r,c].set_ylabel('Accuracy (%)')
    
    if c < 2:
        c += 1
    else:
        r += 1
        c = 0

 
handles, labels = ax[0,0].get_legend_handles_labels()

linePres = Line2D([0], [0], label='Present', color='darkgrey', linestyle = '-')
lineAbs = Line2D([0], [0], label='Absent', color='darkgrey', linestyle = '--')

handles.extend([linePres,lineAbs])

fig.subplots_adjust(top=0.88, bottom=0.11, left=0.125, right=0.9, hspace=0.245, wspace=0.2)
plt.savefig(os.getcwd()+'/Sensitivity_analysis_accuracy.png')

fig.show()

fig, ax = plt.subplots(2,3,figsize=(20, 10))

locations = ['AMOR','WNAG','MAR','NWPO',"SCS","SPO"]
col = ['orange','red','darkgreen','purple']
variables = [r'$O_2$',r'$NO_x$',r'$Mn$',r'$NH_4$']
r = 0
c = 0

for l, loc in enumerate(locations):
    ax[r,c].grid(axis='both')
    for v, var in enumerate(variables):
        ax[r,c].plot(thresholds*1000, number_pres[l][v,:], color = col[v],linestyle='-')
        ax[r,c].scatter(thresholds*1000, number_pres[l][v,:], color = col[v])
        ax[r,c].plot(thresholds*1000, number_abs[l][v,:], color = col[v], linestyle='--')
        ax[r,c].scatter(thresholds*1000, number_abs[l][v,:], color = col[v],  label = var)

    ax[r,c].set_title(loc, fontweight='bold')
    ax[r,c].set_xlabel('Threshold (μM)')
    ax[r,c].set_ylabel('Number of samples')
    
    if c < 2:
        c += 1
    else:
        r += 1
        c = 0

handles, labels = ax[0,0].get_legend_handles_labels()

linePres = Line2D([0], [0], label='Present', color='darkgrey', linestyle = '-')
lineAbs = Line2D([0], [0], label='Absent', color='darkgrey', linestyle = '--')

handles.extend([linePres,lineAbs])

fig.subplots_adjust(top=0.88, bottom=0.11, left=0.125, right=0.9, hspace=0.245, wspace=0.2)

plt.savefig(os.getcwd()+'/Sensitivity_analysis_samples.png')
fig.show()
# %%
