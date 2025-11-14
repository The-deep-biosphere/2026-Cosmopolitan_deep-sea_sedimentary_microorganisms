# %%
# import libraries
import pandas as pd
import numpy as np
from skbio.stats.composition import clr
import os

# ==============
# Loading data
# ==============

# reading csv file and drop samples with NAs in oxygen data
# -----------------------------------------------------------
data = pd.read_csv(os.getcwd()+"/../4)Significant_gene/MetaGeoSeqData_with_interpolation_MergedData.csv.csv")
data = data.dropna(subset = ["Oxygen_interpolated","NO3_interpolated","NO2_interpolated","Ammonia_interpolated","Fe_interpolated","Mn_interpolated","Sulfate_interpolated"], how ='all')  # filter out the NAs if all of the listed datasets is NA

# get column locations for creating dataset GEO + SEQ
col_Core = data.columns.get_loc('Core')
col_Depth = data.columns.get_loc('Depth')
col_Oxygen_interpolated = data.columns.get_loc('Oxygen_interpolated')
col_Nitrate_interpolated = data.columns.get_loc('NO3_interpolated')
col_Nitrite_interpolated = data.columns.get_loc('NO2_interpolated')
col_Ammonia_interpolated = data.columns.get_loc('Ammonia_interpolated')
col_Fe_interpolated = data.columns.get_loc('Fe_interpolated')
col_Mn_interpolated = data.columns.get_loc('Mn_interpolated')
col_Sulfate_interpolated = data.columns.get_loc('Sulfate_interpolated')

data = np.array(data)

# Sample names with GEO data available
# ------------------------------------
GEOsample_names = data[:,0]

# =================================================
# Preparation OTUtable before merging with GEOdata
# =================================================

# load OTU table 
OTUtable = pd.read_csv(os.getcwd()+"/../3)Merged_OTUtables/OTUtable_merged_family_clr.csv")
OTUtable.columns = OTUtable.columns.str.replace(OTUtable.columns[0], 'OTU')
sample_names = OTUtable.columns

# Filter out sample column for which the geodata is not available for
# --------------------------------------------------------------------
select = 10
for i in range(len(sample_names)):
    if sample_names[i] in GEOsample_names:
        select = np.append(select, 1)
    else:
        select = np.append(select, 0)
        print("No GEO data available for "+sample_names[i])

# Change column with OTU numbers to 1 to keep this column
select[1] = 1

# onlys select samples in OTUtable where GEOdata is available
select = select[1:] == 1
OTUtable = OTUtable.iloc[:,select]

# Create updated sample names list
sample_names = OTUtable.columns[1:]

OTUtable = np.array(OTUtable)

# Remove OTUs that are not present in any data
# ----------------------------------------------
# ls_remove = []
# for i in range(len(OTUtable[:,0])):
#     if np.count_nonzero(OTUtable[i,1:]) == 0:
#         ls_remove.append(i)
#         print("ATTENTION: row",str(i),", OTU "+ str(OTUtable[i,0]) +" is removed due to lack of reads")

# OTUtable = np.delete(np.array(OTUtable), (ls_remove), axis=0)

# %%
# ====================================
# Calculate the relative abundances
# ====================================

#####################
#### not for CLR ####
#####################
OTUtable[:,1:OTUtable.shape[1]] = OTUtable[:,1:OTUtable.shape[1]].astype(float)

for j in range(1,OTUtable.shape[1]):
    sum_column = OTUtable[:,j].sum()
    OTUtable[:,j] = OTUtable[:,j]/sum_column

# %%
# ====================================
# remove string from OTU number
# ====================================

OTU_number = []

for i in range(len(OTUtable[:,0])):
    OTU_number.append(OTUtable[i,0].split("OTU_")[1])

# =======================================
# Combine the OTU table and the rates
# =======================================

# Create new data frame that will be used for ML
# Each row represent a sample, sample charactaristics and OTU
# The dataframe will become:
# the amount of samples analyzed *times* the amount of OTUs in the OTUtable variable

N_OTU = len(OTUtable[:,0]) # number of OTUs 
combined = np.zeros((len(sample_names),N_OTU+9))

for j in range(len(sample_names)):
    
    # Select all OTUs of one sample
    sample = sample_names[j] 
    print(sample)
    
    for i in range(N_OTU):
        
        # Select row of GEOdata that consist of the info of the sample selected
        select = sample == GEOsample_names
        
        combined[j,0] = data[select,col_Core][0]    # Core ID number
        combined[j,1] = data[select,col_Depth][0]   # Depth [cm]
        # combined[j,2] = OTUtable[i,j+1]                # Relative abundance
        # combined[j,3] = OTUtable[i,0].split("OTU_")[1]                 # OTUtable[i,0]; if OTU_ before number: OTUtable[i,0].split("OTU_")[1]
        combined[j,2] = data[select,col_Oxygen_interpolated][0]  # O2 concentration [uM]
        combined[j,3] = data[select,col_Nitrate_interpolated][0]  # Nitrate concentration [uM]
        combined[j,4] = data[select,col_Nitrite_interpolated][0]  # Nitrite concentration [uM]
        combined[j,5] = data[select,col_Ammonia_interpolated][0]  # Ammonia concentration [uM]
        combined[j,6] = data[select,col_Fe_interpolated][0]  # Fe concentration [uM]
        combined[j,7] = data[select,col_Mn_interpolated][0]  # Mn concentration [uM]
        combined[j,8] = data[select,col_Sulfate_interpolated][0]  # SO4 concentration [mM]
        combined[j,9:N_OTU+9] = OTUtable[:,j+1]                        # Abundance

# %%
# Save file
np.savetxt(os.getcwd()+"/Data_prep_ML_FAMILY_datacomparison_clr.csv", combined, delimiter=",", fmt="%s",
           header="Core,Depth,O2conc,NO3conc,NO2conc,NH4conc,Feconc,Mnconc,SO4conc,"+','.join(OTUtable[:,0]))
