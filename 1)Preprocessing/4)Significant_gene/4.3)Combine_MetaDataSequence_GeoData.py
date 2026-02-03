#%%
# importing required libraries
import pandas as pd # version 2.3.2
import numpy as np # version 1.26.4
from scipy.interpolate import interp1d # scipy version 1.10.1; interpolating oxygen profile; scipy version 1.15.0 doesnt have it, then change to np.interp instead
import matplotlib.pyplot as plt # plotting interpolated oxygen profile as check
import os
 
# import metadata of sequencing comparison
MDSeq = pd.read_csv(os.getcwd()+"/../3)Merge_OTUtables/Metadata_Merged_rm_duplicates.csv",
                 delimiter=",", dtype=str)

# import metadata of geochemistry comparison
MDGeo = pd.read_csv(os.getcwd()+"/../4)Significant_gene/All_MICROB_GEO_data.csv",
                 delimiter=",", dtype=str)


#%%
# remove blancs, errors and NAs sequencing Metadata
MDSeq = MDSeq.dropna(subset = ["Depth"])
MDSeq = MDSeq[MDSeq.Depth != "B"]
MDSeq = MDSeq[MDSeq.Depth != "R"]

#%%
# Create columns in MDSeq for the upcoming loop
MDSeq['CruiseID'] = "NA"
MDSeq['CoreID'] = "NA"

for i in range(len(MDSeq.index)):
    
    # -------------------------------------------------------
    # Separate sample info over multiple columns
    # -------------------------------------------------------
    # row are removed so index is not the row number; therefore, the first index is picked
    seq_name = MDSeq.Name[MDSeq.index[i]] 
    
    # split the Sequence Name
    split = seq_name.split("_")
    
    # The first item is the cruise ID
    MDSeq.CruiseID[MDSeq.index[i]] = split[0]
    
    # The second is the core ID
    MDSeq.CoreID[MDSeq.index[i]] = split[1]
    
    # The third part is the depth and should be equal to the depth noted
    if MDSeq.Depth[MDSeq.index[i]] != split[2]:
        print("In sample "+seq_name+" an error has occurred; depths are not equal. Depth is "+MDSeq.Depth[MDSeq.index[i]]+" In line ")
        print(i)
# %%
MDSeq['Oxygen'] = "nan"
MDSeq['Oxygen_interpolated'] = "nan"
MDSeq['NO3'] = "nan"
MDSeq['NO3_interpolated'] = "nan"
MDSeq['Ammonia'] = "nan"
MDSeq['Ammonia_interpolated'] = "nan"
MDSeq['Fe_interpolated'] = "nan"
MDSeq['Mn_interpolated'] = "nan"


for i in np.array(MDSeq['Core'].unique()): 
    CruiseID = np.array(MDSeq[MDSeq['Core'] == str(i)].CruiseID)[1]
    CoreID = np.array(MDSeq[MDSeq['Core'] == str(i)].CoreID)[1]
    Cruise_Core = CruiseID+"_"+CoreID

    # -----------------------------------------------------------------------
    # Create variables for Interpolating 
    # -----------------------------------------------------------------------
    # Make dataframe without NaN's in Oxygen
    MDOx = MDGeo[MDGeo["Oxygen"] != "nan"]
    MDOx = MDOx.dropna(subset = ["Oxygen"])

    # Make dataframe without NaN's in NO3
    MDNO = MDGeo[MDGeo["NO3"] != "nan"]
    MDNO = MDNO.dropna(subset = ["NO3"])
    
    # Make dataframe without NaN's in Ammonia
    MDNH = MDGeo[MDGeo["NH4"] != "nan"]
    MDNH = MDNH.dropna(subset = ["NH4"])
    
    # Make dataframe without NaN's in Fe
    MDFe = MDGeo[MDGeo["Fe"] != "nan"]
    MDFe = MDFe.dropna(subset = ["Fe"])
    
    # Make dataframe without NaN's in Manganese
    MDMn = MDGeo[MDGeo["Mn"] != "nan"]
    MDMn = MDMn.dropna(subset = ["Mn"])
    
    
    # Select data from one core
    Depth_ox = MDOx[MDOx["Core"] == Cruise_Core]['Depth'].astype(float)
    Ox = MDOx[MDOx["Core"] == Cruise_Core]['Oxygen'].astype(float)

    Depth_NO = MDNO[MDNO["Core"] == Cruise_Core]['Depth'].astype(float)
    NO = MDNO[MDNO["Core"] == Cruise_Core]['NO3'].astype(float)
    
    Depth_NH = MDNH[MDNH["Core"] == Cruise_Core]['Depth'].astype(float)
    NH = MDNH[MDNH["Core"] == Cruise_Core]['NH4'].astype(float)
    
    Depth_Fe = MDFe[MDFe["Core"] == Cruise_Core]['Depth'].astype(float)
    Fe = MDFe[MDFe["Core"] == Cruise_Core]['Fe'].astype(float)
    
    Depth_Mn = MDMn[MDMn["Core"] == Cruise_Core]['Depth'].astype(float)
    Mn = MDMn[MDMn["Core"] == Cruise_Core]['Mn'].astype(float)

    
    
    # =======================================================================
    # Add the MetaData of GEOCHEMISTRY (MDGeo) to Sequencing MetaData (MDSeq)
    # =======================================================================
    # -------- Find concentrations for each sample that is sequenced --------
    # -----------------------------------------------------------------------   
    
    for d in np.array(MDSeq[MDSeq['Core'] == str(i)].Depth):
        
        # Print data set that will be matched with GEOdata
        # print(CruiseID+"_"+CoreID+"_"+d)
        
        # -----------------------------------------------------------------------
        # ------ Append the most close concentrations to the depth +- 4 cm ------
        # -----------------------------------------------------------------------

        row_GEOdata = MDGeo[(MDGeo['CoreID'] == CoreID) & 
                        (MDGeo['CruiseID'] == CruiseID) &
                        (MDGeo.Depth.astype(float) == float(d))] 
        
        # Check if GEOdata available for this sequence
        len_row = len(row_GEOdata)
        
        # If available (len_row >= 1), CONTINUE
        if len_row >= 1:
                        
            Oxygen_value = np.nanmedian(row_GEOdata.Oxygen.astype(float))
            NO3_value = np.nanmedian(row_GEOdata.NO3.astype(float))
            Ammonia_value = np.nanmedian(row_GEOdata.NH4.astype(float))
        
            
            MDSeq.loc[(MDSeq['Core'] == str(i)) & 
                  (MDSeq['Depth'] == d),'Oxygen'] = Oxygen_value
            MDSeq.loc[(MDSeq['Core'] == str(i)) & 
                  (MDSeq['Depth'] == d),'NO3'] = NO3_value
            MDSeq.loc[(MDSeq['Core'] == str(i)) & 
                  (MDSeq['Depth'] == d),'Ammonia'] = Ammonia_value
            
        # -----------------------------------------------------------------------
        # Interpolating if data available
        # -----------------------------------------------------------------------

        if len(Ox) > 1:
            
            # Interpolate the oxygen data
            interpolate_ox = interp1d(Depth_ox, Ox)

            unique_depth =  Depth_ox.unique()
            # Interpolate from the start to end of existing data
            xdata_ox = np.arange(unique_depth[0],
                             unique_depth[len(unique_depth)-1])
            ydata_ox = interpolate_ox(xdata_ox)


            # Only predict for values inside the known values
            # If the dataset ends in zeros they have already been manually extended to the end for this dataset
            if float(d) >= Depth_ox.unique()[0] and float(d) <= Depth_ox.unique()[len(Depth_ox.unique())-1]:
                # Predict interpolated for depth (d) and insert
                Oxygen_interpolated_value = interpolate_ox(d)
                
                MDSeq.loc[(MDSeq['Core'] == str(i)) & 
                  (MDSeq['Depth'] == str(d)),'Oxygen_interpolated'] = Oxygen_interpolated_value

        if len(NO) > 1:
            # Interpolate the data
            interpolate_NO = interp1d(Depth_NO, NO)

            # Interpolate from the start to end of existing data
            unique_depth =  Depth_NO.unique()
            # Interpolate from the start to end of existing data
            xdata_NO = np.arange(unique_depth[0],
                             unique_depth[len(unique_depth)-1])
            ydata_NO = interpolate_NO(xdata_NO)

            # Only predict for values inside the known values
            # If the dataset ends in zeros they have already been manually extended to the end for this dataset
            if float(d) >= Depth_NO.unique()[0] and float(d) <= Depth_NO.unique()[len(Depth_NO.unique())-1]:
                # Predict interpolated for depth (d) and insert
                NO3_interpolated_value = interpolate_NO(d)

                MDSeq.loc[(MDSeq['Core'] == str(i)) & 
                  (MDSeq['Depth'] == str(d)),'NO3_interpolated'] = NO3_interpolated_value
        
                
        if len(NH) > 1:
            # Interpolate the data
            interpolate_NH = interp1d(Depth_NH, NH)

            # Interpolate from the start to end of existing data
            xdata_NH = np.arange(Depth_NH.unique()[0],
                              Depth_NH.unique()[len(Depth_NH)-1])
            ydata_NH = interpolate_NH(xdata_NH)

            # Only predict for values inside the known values
            # If the dataset ends in zeros they have already been manually extended to the end for this dataset
            if float(d) >= Depth_NH.unique()[0] and float(d) <= Depth_NH.unique()[len(Depth_NH.unique())-1]:
                # Predict interpolated for depth (d) and insert
                Ammonia_interpolated_value = interpolate_NH(d)

                MDSeq.loc[(MDSeq['Core'] == str(i)) & 
                  (MDSeq['Depth'] == str(d)),'Ammonia_interpolated'] = Ammonia_interpolated_value

        if len(Fe) > 1:
            # Interpolate the data
            interpolate_Fe = interp1d(Depth_Fe, Fe)

            # Interpolate from the start to end of existing data
            xdata_Fe = np.arange(Depth_Fe.unique()[0],
                              Depth_Fe.unique()[len(Depth_Fe)-1])
            ydata_Fe = interpolate_Fe(xdata_Fe)


            # Only predict for values inside the known values
            # If the dataset ends in zeros they have already been manually extended to the end for this dataset
            if float(d) >= Depth_Fe.unique()[0] and float(d) <= Depth_Fe.unique()[len(Depth_Fe.unique())-1]:
                # Predict interpolated for depth (d) and insert
                Fe_interpolated_value = interpolate_Fe(d)

                MDSeq.loc[(MDSeq['Core'] == str(i)) & 
                  (MDSeq['Depth'] == str(d)),'Fe_interpolated'] = Fe_interpolated_value

        if len(Mn) > 1:
            # Interpolate the data
            interpolate_Mn = interp1d(Depth_Mn, Mn)

            # Interpolate from the start to end of existing data
            xdata_Mn = np.arange(Depth_Mn.unique()[0],
                              Depth_Mn.unique()[len(Depth_Mn)-1])
            ydata_Mn = interpolate_Mn(xdata_Mn)


            # Only predict for values inside the known values
            # If the dataset ends in zeros they have already been manually extended to the end for this dataset
            if float(d) >= Depth_Mn.unique()[0] and float(d) <= Depth_Mn.unique()[len(Depth_Mn.unique())-1]:
                # Predict interpolated for depth (d) and insert
                Mn_interpolated_value = interpolate_Mn(d)

                MDSeq.loc[(MDSeq['Core'] == str(i)) & 
                  (MDSeq['Depth'] == str(d)),'Mn_interpolated'] = Mn_interpolated_value

        
    if CruiseID == 'SCS':
        if len(Ox) > 1:
            # Show the figures for each cores processed in the code above
            plt.figure()
            plt.title("Oxygen: "+CruiseID+"_"+CoreID)
            plt.plot(Depth_ox, Ox, 'o', xdata_ox, ydata_ox, '-')
            plt.show()

        if len(NO) > 1:  
            # Show the figures for each cores processed in the code above
            plt.figure()
            plt.title("NO3: "+CruiseID+"_"+CoreID)
            plt.plot(Depth_NO, NO, 'o', xdata_NO, ydata_NO, '-')
            plt.show()
            

        if len(NH) > 1:  
            # Show the figures for each cores processed in the code above
            plt.figure()
            plt.title("Ammonia: "+CruiseID+"_"+CoreID)
            plt.plot(Depth_NH, NH, 'o', xdata_NH, ydata_NH, '-')
            plt.show()

        if len(Fe) > 1:  
            # Show the figures for each cores processed in the code above
            plt.figure()
            plt.title("Fe: "+CruiseID+"_"+CoreID)
            plt.plot(Depth_Fe, Fe, 'o', xdata_Fe, ydata_Fe, '-')
            plt.show()

        if len(Mn) > 1:  
            # Show the figures for each cores processed in the code above
            plt.figure()
            plt.title("Mn: "+CruiseID+"_"+CoreID)
            plt.plot(Depth_Mn, Mn, 'o', xdata_Mn, ydata_Mn, '-')
            plt.show()


        
        
# %%
MDSeq.reset_index(drop=True, inplace=True)
for i in range(MDSeq.shape[0]):
    
    if (MDSeq['Oxygen'][i] != 'nan') & (MDSeq['Oxygen_interpolated'][i] == 'nan'):
        
        MDSeq['Oxygen_interpolated'][i] = MDSeq['Oxygen'][i]
        
    if (MDSeq['NO3'][i] != 'nan') & (MDSeq['NO3_interpolated'][i] == 'nan'):
        
        MDSeq['NO3_interpolated'][i] = MDSeq['NO3'][i]
        
    if (MDSeq['Ammonia'][i] != 'nan') & (MDSeq['Ammonia_interpolated'][i] == 'nan'):
        
        MDSeq['Ammonia_interpolated'][i] = MDSeq['Ammonia'][i]
        
# %%
def Header_names(df):
    names = df.columns[0]
    for i in range(len(df.columns)-1):
        names = str(names+","+df.columns[i+1])
    return(names)

# %%
# Save interpolated data
np.savetxt(os.getcwd()+"/../4)Significant_gene/MetaGeoSeqData_with_interpolation_MergedData.csv", MDSeq, delimiter=",", fmt="%s", header = Header_names(MDSeq))
