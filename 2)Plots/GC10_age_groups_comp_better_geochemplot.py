#%%
# Source: https://hands-on.cloud/implementation-of-support-vector-machine-svm-using-python/
# importing required libraries
import numpy as np
import pandas as pd
# import matplotlib.pyplot as plt
import matplotlib as mpl
import os


cores = pd.read_csv(os.getcwd()+"/Core_Coordinates.csv",delimiter=",", dtype=str, index_col=None)
TAX = pd.read_csv(os.getcwd()+"/../1)Preprocessing/3)Merge_OTUtables/TAXtable_merged_family.csv")

variables = ['O2_Coxrate','NO3_Coxrate','MnO2_Coxrate','SO4_Coxrate','nitrite_oxidation_rate','ammonium_oxidation_rate','Mn_anammox_rate','Mnox_rate','NO3_Mnox_rate','Anammox_rate']

# possible other cores: array([ 6,  9, 10, 11, 13, 14, 15, 16, 17])
pred_core = 17
pred_core_name = cores.core[cores.ID.astype(int) == pred_core].values[0]

directory = '/home/renee/Renee-PhD/Projects/Sequencing_Data_ML/Writing/Geochemical_Zonation/Preprocessing/'
data = pd.read_csv(os.getcwd()+"/../1)Preprocessing/4)Significant_gene/Data_prep_ML_FAMILY_datacomparison_owncode_NWPO_SPO_SCS_WNAG_MAR_AMOR_remove_duplicates.csv", #Data_prep_ML_FAMILY_datacomparison_owncode_NWPO_SPO_SCS_WNAG_MAR_AMOR
                delimiter=",", dtype=str, index_col=None)

data.columns = data.columns.str.replace(data.columns[0], 'Core')
data['Core'] = data['Core'].astype(float).astype(int)
data = data.astype(float)
data.reset_index(drop=True, inplace=True)

# select relevant features 
OTU_Aerobic_AMMOX = ['OTU_1010','OTU_10037',] 
OTU_ANAMMOX = ['OTU_10','OTU_11724']
OTU_Anaerobic_Hetro = ['OTU_10246','OTU_11074','OTU_1033','OTU_10379','OTU_10033', 'OTU_10146','OTU_10100', 'OTU_1', 'OTU_1000','OTU_10054']
OTU_Aerobic_Hetro = ['OTU_10511','OTU_10002', 'OTU_10260','OTU_10512', 'OTU_10014', 'OTU_10022']
            
feature_selection = np.concatenate([OTU_Aerobic_AMMOX,OTU_ANAMMOX,OTU_Aerobic_Hetro,OTU_Anaerobic_Hetro])

loc = []
for i in range(len(feature_selection)):

    condition = data.columns == feature_selection[i]
    selection = np.where(condition)
    loc.append(selection[0][0])

print('Feature selection was all present in dataset: '+ str(len(feature_selection)==len(loc)))

# feature selection procedure completion
# --------------------------------------
feature_selection = data.columns[loc]

features = [np.array(TAX[TAX.iloc[:,0]==OTU]['FAMILY'])[0] for OTU in feature_selection]
feature_names_order = [x.strip() for x in features]


def ratios(OTU_list, TAXtable, Data, Core_number):
    
    N_OTUs = len(OTU_list)
    OTU_dic = {}
    OTUs = []
    # OTU_ratios = []
    
    # select core 
    data_temp = Data[Data['Core'] == Core_number]
    
    for i in range(N_OTUs):
        OTU_dic[OTU_list[i]] = np.array(TAXtable[TAXtable.iloc[:,0]==OTU_list[i]]['FAMILY'])[0]
        OTUs.append(np.array(data_temp[OTU_list[i]]).astype(float))
        
    return OTUs, OTU_dic, data_temp

OTU_ratios, OTU_dic, data_temp = ratios(feature_selection, TAX, data, pred_core)

width = 1       # the width of the bars
set_size = 30

# Set general font size
import pylab
pylab.rcParams['font.size'] = '25'

fig, ax = pylab.subplots(1,2,figsize=(20, 15), gridspec_kw={'width_ratios':[0.3,0.7]},
                         constrained_layout=True, sharey=True)

ax[0].grid(axis='both')
ax[1].grid(axis='both')

# fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1): mix = sample x color is x/total sample colours
def colorFader(c1,c2,N): 
    c1=np.array(mpl.colors.to_rgb(c1))
    c2=np.array(mpl.colors.to_rgb(c2))

    col = []
    for i in range(N):
        col.append(mpl.colors.to_hex((1-(i+1)/N)*c1 + (i+1)/N*c2))
    return col


colours = np.concatenate([colorFader('#f7e07f', '#c9a409',len(OTU_Aerobic_AMMOX)),colorFader('#eab7bb', '#ed2e3b',len(OTU_ANAMMOX)),
                          colorFader('#b7eab8', '#218423',len(OTU_Aerobic_Hetro)),colorFader('#d3d6f2', '#162699',len(OTU_Anaerobic_Hetro))])

 
lns1 = ax[0].plot(data_temp.NOxconc*1000, data_temp.Depth, c='black', label = r'$NO_3^-$',linewidth=3)
lns2 = ax[0].plot(data_temp.NH4conc*1000, data_temp.Depth, c='darkgrey', label = r'$NH_4^+$',linewidth=3)
ax2 = ax[0].twiny()
lns3 = ax2.plot(data_temp.O2conc*1000, data_temp.Depth, c='blue', label = r'$O_2$',linewidth=3)

# add these three lines to the legend
lns = lns1+lns2+lns3
labs = [l.get_label() for l in lns]
ax[0].legend(lns, labs, loc=4)

ax[0].xaxis.tick_top()
ax[0].xaxis.set_label_position('top') 
ax2.set_xlabel(r'$O_2$ [$\mu M$]', size = set_size, c='blue')
ax2.xaxis.label.set_color('blue')
ax[0].set_xlabel(f"NO$_3^-$ and NH$_4^+$ [$\mu M$]\n", size = set_size)


ax2.xaxis.set_ticks_position('bottom')
ax2.xaxis.set_label_position('bottom') 

from mpl_axes_aligner import align
align.yaxes(ax[0], 0, ax2, 0, 0.5)

OTU_ratios_total = np.array(np.abs(OTU_ratios)).sum(axis=0)
left = np.zeros(len(data_temp.Depth))

ci = 0
for i in range(len(feature_selection)):

    ax[1].barh(data_temp.Depth, np.abs(OTU_ratios[i]), width, color=colours[ci], label=OTU_dic[feature_selection[i]], left=left)
    left += np.abs(OTU_ratios[i])

    ci += 1

ax[1].set_ylabel('Depth (cm)', size = set_size)
ax[1].xaxis.tick_top()
ax[1].xaxis.set_label_position('top') 

ax[0].set_ylim(np.max(data_temp.Depth),0)
ax[1].set_ylim(np.max(data_temp.Depth),0)

pylab.tight_layout()
figLegend = pylab.figure(figsize = (11, 30))

# produce a legend for the objects in the other figure
pylab.figlegend(*ax[1].get_legend_handles_labels(), loc = 'upper left')

# save the two figures to files
# fig.savefig("G10_barplot15.png")
# figLegend.savefig("G10_barplot_legend15.png")


pylab.show()

