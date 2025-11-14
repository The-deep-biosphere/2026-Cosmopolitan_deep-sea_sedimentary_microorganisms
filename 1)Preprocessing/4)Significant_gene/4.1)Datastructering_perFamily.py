#%%
import pandas as pd
import numpy as np
import os

# Load data
otutable = pd.read_csv(os.getcwd()+"/../3)Merge_OTUtables/OTUtable_merged.csv")
taxtable = pd.read_csv(os.getcwd()+"/../3)Merge_OTUtables/TAXtable_merged.csv")
taxtable = taxtable.iloc[:,:12]

# rename first column to overcome errors with names
otutable = otutable.rename(columns = {otutable.columns[0]:'OTU'})
taxtable = taxtable.rename(columns = {taxtable.columns[0]:'OTU'})
col_order = taxtable.columns.get_loc('FAMILY')

# replace unclassified archaea with right classification
list_arch = pd.read_csv(os.getcwd()+"/../3)Merge_OTUtables/reclassification_archaea/loki_arch_list_OTUs.csv",header=None)
replacement = ['root','Main genome','Archaea','Archaea (superkingdom)','Asgard','Lokiarchaeota','Lokiarchaeia','Unclassified_Lokiarchaeia','Unclassified_Lokiarchaeia','Unclassified_Lokiarchaeia','Unclassified_Lokiarchaeia']

for OTU in list_arch[0]:
    loc = np.where(taxtable['OTU'] == list_arch[0][0])[0][0]
    taxtable.iloc[loc,1:] = replacement

# =======================================
# Renew OTU numbers per order
# =======================================
unique = taxtable['FAMILY'].unique()
# unique = unique[unique != 'Unclassified']                    # remove unclassified to prevent it from collapsing together
names_taxtable = list(taxtable.columns[:col_order+1])

# =======================================
# Renew OTU numbers per order
# =======================================
unique = taxtable['FAMILY'].unique()
# unique = unique[unique != 'Unclassified']                    # remove unclassified to prevent it from collapsing together
names_taxtable = list(taxtable.columns[:col_order+1])
new_taxtable = []        # create a new taxtable with new OTU number and the relating ORDER 

#%%
for tax_order in unique:

    # take new OTU number (first of the existing OTUs with the same order)
    OTU = taxtable[taxtable.iloc[:,col_order] == tax_order]
    OTU_number = list(OTU['OTU'])[0]
    # print(OTU_number)
    
    # store the new OTU number and relating order
    new_taxtable.append(list(OTU.iloc[0,:col_order+1]))
    
    # select all row equal to the order
    for row in np.where([taxtable['FAMILY'] == tax_order])[1]:
        
        # replace OTU number for all the equal orders
        otu = taxtable['OTU'][row]
        # print(otu)
        # replace all OTUs with the same order with the same number
        otutable['OTU'] = otutable['OTU'].replace(otu, OTU_number)

new_taxtable = pd.DataFrame(new_taxtable)
new_otutable = otutable.groupby(['OTU']).sum()

def Header_otutable(otutable):
    names = "OTU"
    for i in range(len(otutable)-1):
        names = str(names+","+otutable[i+1])
    return names

#%%
# Save family level taxa and OTU table
np.savetxt(os.getcwd()+"/../3)Merge_OTUtables/TAXtable_merged_family.csv", new_taxtable, delimiter=",", fmt="%s", header = Header_otutable(names_taxtable))
new_otutable.to_csv(os.getcwd()+"/../3)Merge_OTUtables/OTUtable_merged_family.csv")