#%%
import numpy as np
import pandas as pd
# import matplotlib.pyplot as plt
import pylab
import matplotlib as mpl
import os

# ===================
# Data selection
# ===================

TAX = pd.read_csv(os.getcwd()+"/../1)Preprocessing/3)Merge_OTUtables/TAXtable_merged_family.csv")

TAX.columns = TAX.columns.str.replace(TAX.columns[0], 'OTU')

data = pd.read_csv(os.getcwd()+"/../1)Preprocessing/4)Significant_gene/Data_prep_ML_FAMILY_datacomparison.csv",
                delimiter=",", dtype=str, index_col=None)

data.columns = data.columns.str.replace(data.columns[0], 'Core')

# remove cores that are not of interest (GS19GC25, GS21GC09, GS20GC20, GS20GC21)
data = data[data['Core'] != '18.0']
data = data[data['Core'] != '19.0']
data = data[data['Core'] != '22.0']
data = data[data['Core'] != '21.0']

data.reset_index(drop=True, inplace=True)

# find core names based on core ID
core_names = pd.read_csv(os.getcwd()+"/Core_Coordinates.csv")
data.Depth = data.Depth.astype(float)
data.Core = data.Core.astype(float).astype(int)
names = [core_names.core[np.where(np.array(core_names.ID) == i)[0][0]] for i in np.array(data.Core)]


colour_code = np.array(data.Core)
colour_label = {0:'AMOR',1:"NWPO",2:'WNAG',3:'MAR',4:'SCS',5:'SPO'}

for i,c in enumerate(np.array(data.Core)):
    
    if c <= 25:
        colour_code[i] = 0

    if c in np.arange(26,40):
        colour_code[i] = 1
    
    if c in np.arange(40,43):
        # store atlantic in another variable (40)
        colour_code[i] = 2
    
    if c in np.arange(43,45):
        # store north pond in another variable (41 and 42)
        colour_code[i] = 3

    if c in np.arange(45,55):
        # store south china sea in another variable (45 up to 55)
        colour_code[i] = 4

    if c in np.arange(55,68):
        # store south pacific ocean in another variable (56 up to 68)
        colour_code[i] = 5


# --------------------------
# OTU selection Pie Chart
# --------------------------

OTU_Aerobic_AMMOX = ['OTU_1010','OTU_10037',] 
OTU_ANAMMOX = ['OTU_10','OTU_11724']
OTU_Anaerobic_Hetro = ['OTU_10246','OTU_11074','OTU_1033','OTU_10379','OTU_10033', 'OTU_10146','OTU_10100', 'OTU_1', 'OTU_1000','OTU_10054']
OTU_Aerobic_Hetro = ['OTU_10511','OTU_10002', 'OTU_10260','OTU_10512', 'OTU_10014', 'OTU_10022']

ALL = np.concatenate([OTU_Aerobic_AMMOX,OTU_ANAMMOX,OTU_Aerobic_Hetro,OTU_Anaerobic_Hetro])

TAX_features_family = [np.array(TAX[TAX.iloc[:,0]==OTU]['FAMILY'])[0] for OTU in ALL]


def colors (palette,N):
    c = sns.color_palette(palette,N)
    return c

def colorFader(c1,c2,N): #fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1): mix = sample x color is x/total sample colours
    c1=np.array(mpl.colors.to_rgb(c1))
    c2=np.array(mpl.colors.to_rgb(c2))

    col = []
    for i in range(N):
        col.append(mpl.colors.to_hex((1-(i+1)/N)*c1 + (i+1)/N*c2))
    return col


OTUs = ([OTU_Aerobic_AMMOX,OTU_ANAMMOX,OTU_Aerobic_Hetro,OTU_Anaerobic_Hetro,ALL])
name_list = ['Aerobic Ammonium Oxidation','Anaerobic Ammonium Oxidation','Aerobic OM degradation','Anaerobic OM degradation','All']
import seaborn as sns
colours1 = sns.color_palette("RdPu",20)
colours2 = sns.color_palette("tab10",len(ALL)-20)
colours = colours1+colours2


colours = np.concatenate([colorFader('#f7e07f', '#c9a409',len(OTU_Aerobic_AMMOX)),colorFader('#eab7bb', '#ed2e3b',len(OTU_ANAMMOX)),
                          colorFader('#b7eab8', '#218423',len(OTU_Aerobic_Hetro)),colorFader('#d3d6f2', '#162699',len(OTU_Anaerobic_Hetro))])


colour_groups = np.concatenate([['#c9a409']*len(OTU_Aerobic_AMMOX),['#ed2e3b']*len(OTU_ANAMMOX),['#218423']*len(OTU_Aerobic_Hetro),['#162699']*len(OTU_Anaerobic_Hetro)])
hatch_groups = np.concatenate([['/']*len(OTU_Aerobic_AMMOX),['.']*len(OTU_ANAMMOX),['+']*len(OTU_Aerobic_Hetro),['O']*len(OTU_Anaerobic_Hetro)])

# fix colour per taxa
taxa = [np.array(TAX[TAX.iloc[:,0]==OTU]['FAMILY'])[0] for OTU in ALL]
taxa = [x.strip() for x in taxa]
dic_colour = {}

dic_colour_groups = {} 
dic_hatch_groups={}
for i,key in enumerate(taxa):
    dic_colour_groups[key] = colour_groups[i]
    dic_hatch_groups[key] = hatch_groups[i]
    dic_colour[key] = colours[i]

pylab.rcParams["figure.figsize"] = [20, 20]
pylab.rcParams["figure.autolayout"] = True
f,(axs) = pylab.subplots(len(OTUs),1)

for j in np.arange(len(OTUs)):

    loc_features = []
    for i in range(len(OTUs[j])):
        if OTUs[j][i] in data.columns:
            loc = np.where(data.columns == OTUs[j][i])[0][0]
            loc_features.append(loc)

    # get taxanomy family level
    TAX_features_family = [np.array(TAX[TAX.iloc[:,0]==OTU]['FAMILY'])[0] for OTU in OTUs[j]]
    feature_names_family = [x.strip() for x in TAX_features_family]
    axs[j].set_title(name_list[j]+' \n', fontsize = 16, fontweight='bold')
    axs[j].axis('off')

    
    for i in np.unique(colour_code):
        # select location and features from data
        # ----------------------------------------
        temp_data = np.array(data.iloc[np.where(colour_code == i)[0],loc_features].astype(float))
        y = temp_data.sum(axis=0)/temp_data.sum()

        ax = f.add_subplot(len(OTUs),max(colour_code)+1,(j)*(max(colour_code)+1)+i+1)
        
        patches1, texts1 = ax.pie(y, startangle=90, colors=[dic_colour[i] for i in feature_names_family])
            
        if name_list[j] == 'All':
            patches, texts = ax.pie(y, startangle=90, colors=[dic_colour_groups[i] for i in feature_names_family], hatch=[dic_hatch_groups[i] for i in feature_names_family])

        ax.set_title(colour_label[i])


pylab.subplots_adjust(left = .008, right = .608, bottom = .013, top = .839, wspace=0, hspace=.21)

figLegend = pylab.figure(figsize = (7, 20))

# produce a legend for the objects in the other figure
pylab.figlegend(patches1, feature_names_family, loc = 'upper left')

# save the two figures to files
# f.savefig("piechart10.png")
# figLegend.savefig("legend_piechart10.png")

pylab.show()
# %%
