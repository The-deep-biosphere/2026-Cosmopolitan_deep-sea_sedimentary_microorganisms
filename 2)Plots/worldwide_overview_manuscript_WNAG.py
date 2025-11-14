# Libraries
import seaborn as sns
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import os

var = 'oxygen'
var2 = 'manganese'
var3 = 'nitrate'

# var = 'manganese'
# var2 ='oxygen'

dic_name = {'oxygen':'O2conc', 'nitrate':'NO3conc','manganese':'Mnconc','NH4':'NH4conc'}
dic_loc = {'oxygen':2, 'nitrate':4,'manganese':10,'NH4':6,'Depth':1,'Core':0}

selected_core = 2
# ===================
# Data selection
# ===================
TAX = pd.read_csv(os.getcwd()+"/../1)Preprocessing/3)Merge_OTUtables/TAXtable_merged_family.csv")
TAX.columns = TAX.columns.str.replace(TAX.columns[0], 'OTU')

data = pd.read_csv(os.getcwd()+"/../1)Preprocessing/4)Significant_gene/Data_prep_ML_FAMILY_datacomparison.csv",
                delimiter=",", dtype=str, index_col=None)
data.columns = data.columns.str.replace(data.columns[0], 'Core')

data = data.dropna(subset = [dic_name[var]], how ='any') 
# data = data.dropna(subset = [dic_name[var2]], how ='any') 
# data = data.dropna(subset = [dic_name[var3]], how ='any') 

loc_features = []
for i in range(len(data.columns)):
    if 'OTU_' in data.columns[i]:
        loc_features.append(i)

feature_names = data.iloc[:,loc_features].columns

# remove cores that are not of interest (GS19GC25, GS21GC09, GS20GC20, GS20GC21)
data = data[data['Core'] != '18.0']
data = data[data['Core'] != '19.0']
data = data[data['Core'] != '22.0']
data = data[data['Core'] != '21.0']

data['Core'] = data['Core'].astype(float).astype(int)
# ==============
# Colour coding
# ==============
data.reset_index(drop=True, inplace=True)
colour_code = np.array(data.Core)
colour_label = {0:'AMOR',1:"NWPO",2:'WNAG',3:'MAR',4:'SCS',5:'SPO'}

for i,c in enumerate(data.Core):
    
    # AMOR
    if c <= 25:
        colour_code[i] = 0

    # NWPO
    if c in np.arange(26,40):
        colour_code[i] = 1
    
    # WNAG
    if c in np.arange(40,43):
        # store atlantic in another variable (40, 41, 42)
        colour_code[i] = 2
    
    # MAR
    if c in np.arange(43,45):
        # store north pond in another variable (43 and 44)
        colour_code[i] = 3

    # SCS
    if c in np.arange(45,55):
        # store south china sea in another variable (45 up to 55)
        colour_code[i] = 4

    # SPO
    if c in np.arange(55,68):
        # store south pacific ocean in another variable (56 up to 68)
        colour_code[i] = 5

core_selection = np.where(colour_code == selected_core)[0]
colour_code = colour_code[core_selection]
data = data.iloc[core_selection,:]
data.reset_index(drop=True, inplace=True)

# find core names based on core ID
core_names = pd.read_csv(os.getcwd()+"/Core_Coordinates.csv")

# ================================
# -------- Sorting data ----------
# ================================

# sort on environemental variable, core, depth, 
# where the first variable has the highest priority 
data.Depth = data.Depth.astype(float)
data[dic_name[var]] = data[dic_name[var]].astype(float)*1000
data[dic_name[var]] = data[dic_name[var]].astype(int)

if var == 'manganese':
    data = data.sort_values([dic_name[var],'Depth','Core'], ascending=False)

else:
    data = data.sort_values([dic_name[var],'Depth','Core'], ascending=[True,False,True])
data.reset_index(drop=True, inplace=True)

X = data.iloc[:, loc_features].astype(float)
Y = data[dic_name[var]]
Y2 = data[dic_name[var2]].astype(float)*1000
Y3 = data[dic_name[var3]].astype(float)*1000

NH4 = data[dic_name['NH4']].astype(float)*1000
core = data['Core'].astype(float).astype(int)
depth = data['Depth'].astype(float)

X_temp = data.iloc[:, loc_features].astype(float)

# select relevant features according to q-value
# ------------------------------------------------------
feature_selection = ['OTU_1010', 'OTU_10037', 'OTU_10', 'OTU_11724', 'OTU_10033', 'OTU_10146', 'OTU_10100', 'OTU_1', 'OTU_1000', 'OTU_10511', 'OTU_10002', 'OTU_10260',
      'OTU_10512', 'OTU_10054', 'OTU_10022', 'OTU_10246','OTU_11074','OTU_10014','OTU_1033','OTU_10379']

prior = []
for i in range(len(feature_selection)):

    condition = feature_names == feature_selection[i]
    selection = np.where(condition)
    prior.append(selection[0][0])

loc = np.unique(prior)#np.concatenate([prior,other])
len(loc)

# feature selection procedure completion
# --------------------------------------
X = X.iloc[:,loc]
X.columns = feature_names[loc]


# change OTU label to order taxonomy
# -------------------------------------
TAX_features = [np.array(TAX[TAX.iloc[:,0]==OTU]['ORDER'])[0] for OTU in feature_names[loc]]
TAX_features_class = [np.array(TAX[TAX.iloc[:,0]==OTU]['CLASS'])[0] for OTU in feature_names[loc]]
TAX_features_family = [np.array(TAX[TAX.iloc[:,0]==OTU]['FAMILY'])[0] for OTU in feature_names[loc]]

# remove spacing
feature_names_order = [x.strip() for x in TAX_features]
feature_names_class = [x.strip() for x in TAX_features_class]
feature_names_family = [x.strip() for x in TAX_features_family]

# attach new feature names
X.columns = feature_names_family# [f'{feature_names_class[OTU]} ({feature_names_order[OTU]})' for OTU in range(len(feature_names_order))]
X.columns = X.columns.str.replace('Unclassified_', '')
X = X.sort_index(axis=1) # alphabatical order of the families

# ===============================================
# Order rows (X index) on Y, high to low oxygen
# ===============================================

X_select = X


rows = np.array(range(Y.shape[0]))+0.5
conc = Y
core_numbers = np.sort(data['Core'].astype(float).astype(int).unique())
# colours2 = sns.color_palette("tab20",np.max(core_numbers.astype(int))+1)
core_ID = [np.array(core_names[core_names['ID'] == i]['core'])[0] for i in core_numbers]

colours = sns.color_palette("tab20",np.max(core_numbers.astype(int))+1)

plt.rcParams['font.size'] = '25'
plt.rcParams["figure.figsize"] = (25,40)
size_dot = 50
f,(axs) = plt.subplots(1,3, gridspec_kw={'width_ratios':[0.2,0.2,1]})

axdepth = axs[0]
ax3 = axs[1]
ax4 = axs[2]

axdepth.grid(axis='both')
ax3.grid(axis='both')

f.subplots_adjust(wspace=0.1)

for i, c in enumerate(core_numbers):
    axdepth.scatter(depth[data.Core == c],rows[data.Core == c]+.5, color = colours[i], s=size_dot+3, label=core_ID[i])


axdepth.set_yticks([])
axdepth.set_ylabel('Samples', fontsize = '30')
axdepth.set_ylim(0,Y.shape[0])
axdepth.set_xticks(np.arange(0, max(depth)+500., 2000.))
axdepth.set_xlim(-5,max(depth)+500)

axdepth.set_xlabel('Depth [cm]', fontsize = '25')
axdepth.tick_params(axis="x", labelsize=15)
axdepth.legend(loc='upper center', bbox_to_anchor=(6, 1.05),
        ncol=5, fancybox=True, shadow=False)

ax3.plot(conc,rows, c = 'c')
ax3.yaxis.tick_right()
ax3.set_yticks([])
ax3.tick_params(axis="x", labelsize=15)
ax3.set_ylim(0,Y.shape[0])
ax3.set_xlabel(r'O$_2$ [µM]', fontsize = '25')
ax3.set_xlim(0,300)
ax3.set_xticks(np.arange(0, 300, 100.0))
line = np.max(np.where(conc == 0))
ax3.text(80, line, '0 µM')
ax3.axhspan(0, line, facecolor='0.5', alpha=0.5)

g1 = sns.heatmap(X_select, xticklabels=True,ax=ax4, cbar=True, cmap='rocket_r') # YlOrBr_r
g1.set_ylabel('')
g1.set_xlabel('')
g1.set_yticks([])
g1.set_ylim(0,np.max(rows))

plt.subplots_adjust(left = .03, right = .97, bottom = .2, top = .95, wspace=.2, hspace=.2)

plt.savefig(os.getcwd()+f'/overview_manuscript_{var}-sites_WNAG.png', dpi=300)

plt.show()


