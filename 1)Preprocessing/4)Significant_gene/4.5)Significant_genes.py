
####
# FIX PROBLEM WITH UNPUBLISHED DATA OUT OF HERE
###

# %%
# Libraries
import seaborn as sns
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import os 
from pathlib import Path
from scipy.interpolate import splrep, BSpline,CubicSpline
from tqdm import tqdm
import composition_stats as compst
import os 

var = 'oxygen'
# var = 'manganese'
# var = 'nitrate'
# var = 'ammonium'

dic_name = {'oxygen':'O2conc', 'nitrate':'NOxconc','manganese':'Mnconc','ammonium':'NH4conc'}
dic_loc = {'oxygen':2, 'nitrate':4,'manganese':10,'ammonium':6,'Depth':1,'Core':0}

# ===================
# Data selection
# ===================

TAX = pd.read_csv(os.getcwd()+"/../3)Merged_OTUtables/TAXtable_merged_family.csv")

Data = pd.read_csv(os.getcwd()+"/../4)Significant_gene/Data_prep_ML_FAMILY_datacomparison_clr.csv", delimiter=",", dtype=str, index_col=None)
Data.columns = Data.columns.str.replace(Data.columns[0], 'Core')
Data['Core'] = Data['Core'].astype(float).astype(int)

Data = Data.dropna(subset = [dic_name[var]], how ='any') 
Data[dic_name[var]] = Data[dic_name[var]].astype(float)

loc_features = []
for i in range(len(Data.columns)):
    if 'OTU_' in Data.columns[i]:
        loc_features.append(i)

feature_names = Data.iloc[:,loc_features].columns

# remove cores that are not of interest (GS19GC25, GS20GC20, GS20GC21)
Data = Data[Data['Core'] != 18]
Data = Data[Data['Core'] != 22]
Data = Data[Data['Core'] != 21]

Data.reset_index(drop=True, inplace=True)

#%%
# --------------------
# ====================
# Compute two-sample t-statistics for differential expression of each gene between ER positive and ER negative groups.
# ====================
# --------------------

# Devide the Data in two groups
# ------------------------------
if (False):
    # The dependent variable
    ER_stat = np.array(Data[dic_name[var]])*1000
    # Define threshold
    threshold_uM = 0
    val = ER_stat == threshold_uM
    # CLR transform chemical data
    ER_stat = compst.clr(np.array(ER_stat+0.001))
    # find CLR transformed threshold
    threshold = ER_stat[val][0]
else:
    ER_stat = np.array(Data[dic_name[var]])*1000
    threshold = 0

# Find locations
samples_ER_neg = np.where(ER_stat <= threshold)[0]
samples_ER_pos = np.where(ER_stat > threshold)[0]

# Data with positive or negative ER status
Data_pos_neg = Data.iloc[:,loc_features].astype(float)
gene_names = feature_names

# remove OTUs not present in any data
zeros_columns = np.where(np.all(Data_pos_neg == 0, axis=0))[0]

if len(zeros_columns) >= 1:
    Data_pos_neg = Data_pos_neg.drop(Data_pos_neg.columns[zeros_columns], axis = 1)

# Devide the Data
Data_neg = Data_pos_neg.iloc[samples_ER_neg,:].astype(float)
Data_pos = Data_pos_neg.iloc[samples_ER_pos,:].astype(float)

# t[i] = (mean(x[i][2]) - mean(x[i][1]))/np.sqrt(s^2[i][1]/n[1] + s^2[i][1]/n[1]) -> as described in "Statistical significance for genomewide studies"
# where i is each gene
t_values = []
for i in range(Data_pos_neg.shape[1]):

    x_neg = np.array(Data_neg.iloc[:,i])
    x_pos = np.array(Data_pos.iloc[:,i])
    var_neg = np.var(x_neg)/len(x_neg)
    var_pos = np.var(x_pos)/len(x_pos)
    t_values.append((np.mean(x_neg) - np.mean(x_pos))/np.sqrt(var_pos+var_neg))

t_values = np.array(t_values)

#%%
# --------------------
# ====================
# Compute theoretical and empirical p-values for each gene (use the permutation procedure from Remark C in the paper).
# ====================
# --------------------


# number of permutations
n_permutations = 100
# create permutated t-test for every permutation
t_jb = np.zeros((Data_pos_neg.shape[1], n_permutations))

for b in tqdm(range(n_permutations)):

    # Permutation of array labels
    permutation_locations = np.random.permutation(Data_pos_neg.shape[0])
    Data_perm = np.array(Data_pos_neg)
    Data_perm = Data_perm[permutation_locations,:] 
    x_neg = Data_perm[:Data_neg.shape[0],:]
    x_pos = Data_perm[Data_neg.shape[0]:,:]
    
    var_neg = np.var(x_neg,axis=0)/x_neg.shape[0]
    var_pos = np.var(x_pos,axis=0)/x_pos.shape[0]
    t_jb[:,b] = (np.mean(x_neg,axis=0) - np.mean(x_pos,axis=0))/(np.sqrt(var_pos + var_neg))


        
# calculate p-value for each gene (theoretical)
p_values_theoretical = []
for i in range(t_jb.shape[0]):

    p_ib = []
    for b in range(n_permutations):
        count = np.sum(abs(t_jb[:,b]) >= abs(t_jb[i,1]))
        p_ib.append(count/(t_jb.shape[0]*n_permutations))
    
    p_values_theoretical.append(np.sum(p_ib))

p_values_theoretical = np.array(p_values_theoretical)


# calculate p-value for each gene (empirical)
p_values = []
for i in range(t_jb.shape[0]):

    p_ib = []
    for b in range(n_permutations):
        count = np.sum(abs(t_jb[:,b]) >= abs(t_values[i]))
        p_ib.append(count/(t_jb.shape[0]*n_permutations))
    
    p_values.append(np.sum(p_ib))

p_values = np.array(p_values)

#%%

# --------------------
# ====================
# Plot the p-value histograms. Do the theoretical and empirical distributions differ?
# ====================
# --------------------

plt.figure()
plt.hist(p_values_theoretical)
plt.xlabel('P value')
plt.ylabel('Number of genes')
plt.title('Theoretical distribution (Null hypothesis)')
plt.show()

# figure 1 of paper
plt.figure()
plt.hist(p_values)
plt.xlabel('P value')
plt.ylabel('Number of genes')
plt.title('Empirical distribution')
plt.show()

# The empirical distribution differs from the theoretical distribution


#%%
# --------------------
# ====================
# Implement the algorithm for estimating q-values from Remark B in the paper.
# ====================
# --------------------

## step 1
sorted1 = p_values.argsort()
p_sort = p_values[sorted1]

## step 2
lamb = np.arange(0.,.96,0.01)

pi = []
counts = np.array([(p_sort > i).sum() for i in lamb])
for i, l in enumerate(lamb):
    pi.append(counts[i]/(len(p_sort)*(1-l)))

pi = np.array(pi)

## step 3 
# the natural cubic spline can be achieved using splrep, BSpline from scipy or poly1d, however the fitting is very poorly.
# therefore, cubic spline will be used for pi0 estimation (cs).
poly = np.poly1d(np.polyfit(lamb,pi, 4))
tck = splrep(lamb,pi, s=3)
# fit cubic spline
cs = CubicSpline(lamb,pi, axis=0, bc_type='natural', extrapolate=True)

# figure 3 of paper
fig, ax = plt.subplots(figsize=(6.5, 4))
ax.plot(lamb, poly(lamb), c= 'blue', label = '3 degree polynomial',alpha=0.5)
ax.plot(lamb, BSpline(*tck)(lamb), c = 'red', label='BSpline',alpha=0.5)
ax.plot(lamb, cs(lamb), c = 'green', label='cubic spline')
ax.scatter(lamb,pi,s=3, c = 'black')
plt.xlabel(r'$\lambda$')
plt.ylabel(r'$\pi_0(\lambda$)')
plt.legend()
plt.show()

## step 4
# taken a bit lower since the last point in the graph goes up

if var == 'manganese':
    pi0 = BSpline(*tck)(0.7) 
    
if var == 'oxygen' or var == 'nitrate' or var == 'ammonium':
    pi0 = BSpline(*tck)(.8)

print(f'$\pi_0$ is {pi0}')

## step 5
# determine the last value of q (q(p[m]))
q_value = pi0 * p_sort
q_value[-1] = min(q_value[-1], 1.0)

## step 6
# overwrite the other q_values with the following algorithm
for i in range(len(p_sort)-2, 0, -1):
    q_value[i] = min(pi0*len(p_sort)*p_sort[i]/(i), q_value[i+1])     



#%%
# --------------------
# ====================
# Reproduce Fig. 2 from the paper for the ER+/- differential expression p-values.
# ====================
# --------------------

fig, ax = plt.subplots(2,2,figsize=(6.5, 4))

ax[0,0].scatter(t_values[sorted1],q_value, s=.1)
ax[0,0].set_xlabel(r't-statistics')
ax[0,0].set_ylabel(r'q-values')

ax[0,1].plot(p_sort,q_value)
ax[0,1].set_xlabel(r'p-values')
ax[0,1].set_ylabel(r'q-values')

q_value_sum = []
for i in range(q_value.shape[0]):
    q_value_sum.append(np.sum(q_value <= q_value[i]))

q_value_sum = np.array(q_value_sum)

ax[1,0].plot(q_value, q_value_sum)
ax[1,0].set_ylabel(r'number of significant genes')
ax[1,0].set_xlabel(r'q-values')
ax[1,0].set_xlim(-0.001,0.1)


expected_FP = q_value*len(q_value)
ax[1,1].plot(q_value_sum,expected_FP)
ax[1,1].set_xlabel(r'number of significant genes')
ax[1,1].set_ylabel(r'number of expected FP')
plt.tight_layout()
plt.show()

#%%
# Find significant gene names
significant_genes_loc = max(np.where(expected_FP ==0)[0])
print(f'There are {q_value_sum[significant_genes_loc]} significant genes with q-value {q_value[significant_genes_loc]}')

genes = gene_names[sorted1]

OTUs_significant = genes[:q_value_sum[significant_genes_loc]]

features = [np.array(TAX[TAX.iloc[:,0]==OTU]['FAMILY'])[0] for OTU in OTUs_significant]

feature_names_order = [x.strip() for x in features]


sorted(feature_names_order)

# %%
np.savetxt(os.getcwd()+"/../4)Significant_gene/Significant_Features_FAMILY_comparison_AMOR_WNAG_MAR_SCS_NWPO_SPO_"+var+"_clr_03062024.csv", OTUs_significant, delimiter=",", fmt="%s")

