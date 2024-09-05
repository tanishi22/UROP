#in TCGA-GBM GDC data, 10% of 9p21.3 deleted patients were not CDKN2A deleted. Similarly, in TCGA-GBM CPTAC data, 14.6% of 9p21.3 deleted patients were not CDKN2A deleted. 
#what other deletions do these patients have to be classified as 9p21.3 deleted? 

#Checking for deletions in CDKN2B, MTAP, IFN, and other key chr9p21.3 genes in these exception patients.  

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np 

# TCGA-GBM GDC data 

seg_data = pd.read_csv("/Users/tanishi/Desktop/data_cna_hg38.seg", sep="\t")

# Define genomic regions to identify and extract exception patients 
chr9p21_start = 19900002
chr9p21_end = 25600000
CDKN2A_start = 21967752
CDKN2A_end = 21995324
deletion_threshold = -0.1

chr9p_seg = seg_data[(seg_data['chrom'] == "9") & (seg_data['loc.start'] >= chr9p21_start) & (seg_data['loc.end'] <= chr9p21_end)]
chr9p_seg['is_deletion'] = chr9p_seg['seg.mean'] < deletion_threshold
patients_9p_del = chr9p_seg[chr9p_seg['is_deletion']]['ID'].unique() 

CDKN2A_seg = chr9p_seg[(chr9p_seg['loc.start'] < CDKN2A_end) & (chr9p_seg['loc.end'] > CDKN2A_start)]
CDKN2A_seg['is_deletion'] = CDKN2A_seg['seg.mean'] < deletion_threshold
patients_CDKN2A_del = CDKN2A_seg[CDKN2A_seg['is_deletion']]['ID'].unique() 

print(chr9p_seg[chr9p_seg['ID'] == "TCGA-RR-A6KC-01A"])

chr9p_exception = []
for patient in patients_9p_del:
    if patient not in patients_CDKN2A_del:
        chr9p_exception.append(patient)
len(chr9p_exception) #36 patients are 9p21.3 deleted but not CDKN2A deleted
exception_df = chr9p_seg[chr9p_seg['ID'].isin(chr9p_exception)] 

#Define key chr9p21.3 genes in dictionary. defining start as the first 1kb of 9p21.3 with micellaneous genes. 
# Broadly defined IFN region as region with lots of IFN region. Specify later. starts with IFNB1 start position and ends with IFNE end position. 
ROI = {'start': (19900002, 20900002),
       'IFN' : (21077104, 21482313),
       'MTAP': (21802636, 21941115),
       'CDKN2B' : (22002903, 22009313),
       'end': (22200000, 25600000)}

#initialise results dataframe 
result_df = pd.DataFrame(columns=['ID'] + list(ROI.keys()))

# Iterate over each patient to check if any segment in the patient overlaps with the ROI and is deleted. 
for patient_id, patient_data in exception_df.groupby('ID'):
    # Initialise a dictionary to store the results for this patient
    patient_results = {'ID': patient_id}
    
    # Check each region of interest
    for region_name, (roi_start, roi_end) in ROI.items():
        deletion_found = any(
            (row['loc.start'] <= roi_end) and 
            (row['loc.end'] >= roi_start) and 
            (row['is_deletion'] == True)
            for _, row in patient_data.iterrows()
        )
        # Store the result
        patient_results[region_name] = deletion_found
    
    # Append the results for this patient to the DataFrame
    result_df = result_df.append(patient_results, ignore_index=True)
result_df = result_df.set_index('ID')
result_df
print(result_df[result_df["CDKN2B"] == True])

# Only two patients have deletions in MTAP. Both also have deletions in the IFN region and have a segment deleted from the middle of the broad IFN region to the first half of the MTAP gene. 
print(result_df[result_df["MTAP"] == True])
print(exception_df[exception_df['ID'] == "TCGA-02-0038-01A"])

# 7 patients have a deletion in the broadly defined IFN region. 
# 2/7 have deletions spanning the start region and IFN. 2/7 have deletions in IFN and MTAP (broadly). 3/7 have deletions in IFN only. 
IFN_patients = result_df[result_df["IFN"] == True]
print(exception_df[exception_df['ID'] == "TCGA-16-1048-01B"])

#Investigate segments deleted in IFN-del only patients. All of these patients have deletions in IFN genes only. 
# Should not classify them as 9p21.3 deleted because they only have one deletion in chr9p21. so probably ignore these patients. 
IFN_only_patients = result_df[(result_df["IFN"] == True) & (result_df["MTAP"] == False) & (result_df["start"] == False)]

#End region deletions
print(len(result_df[result_df["end"] == True])) #23 patients out of 36 have a del in the end region 

#Visualisations
region_counts = result_df.sum()
# Prepare a DataFrame for plotting overlapping deletions
overlap_df = result_df.apply(lambda row: ', '.join([k for k, v in row.items() if v]), axis=1).value_counts()

# Plotting the stacked bar chart for each region
plt.figure(figsize=(10, 6))
region_counts.plot(kind='bar', color='skyblue')
plt.title('Number of Patients with Deletions in Each Region')
plt.xlabel('Regions of Interest')
plt.ylabel('Number of Patients')
plt.xticks(rotation=45, ha='right')
plt.show()

# Plotting the overlapping deletions
plt.figure(figsize=(10, 6))
overlap_df.plot(kind='bar', color='salmon')
plt.title('Patients with Multiple Region Deletions')
plt.xlabel('Overlapping Regions')
plt.ylabel('Number of Patients')
plt.xticks(rotation=45, ha='right')
plt.show()

# Validate with TCGA-GBM CPTAC data. Using the same defined threshold and ROI dictionary as before
cptac_data = pd.read_csv("/Users/tanishi/Desktop/data_cna_hg19.seg", sep="\t")

chr9p_cptac = cptac_data[(cptac_data['chrom'] == "9") & (cptac_data['loc.start'] >= chr9p21_start) & (cptac_data['loc.end'] <= chr9p21_end)]
chr9p_cptac = chr9p_cptac.copy()
chr9p_cptac['is_deletion'] = chr9p_cptac['seg.mean'] < deletion_threshold
patients_9p_del_cp = chr9p_cptac[chr9p_cptac['is_deletion']]['ID'].unique()  #59 patients out of 98 are 9p21.3 deleted

CDKN2A_cptac = chr9p_cptac[(chr9p_cptac['loc.start'] <= CDKN2A_end) & (chr9p_cptac['loc.end'] >= CDKN2A_start)]
CDKN2A_cptac['is_deletion'] = CDKN2A_cptac['seg.mean'] < deletion_threshold 
patients_CDKN2A_del_cp = CDKN2A_cptac[CDKN2A_cptac['is_deletion']]['ID'].unique() #51 patients out of 98 are CDKN2A deleted

chr9p_exception_cp = []
for patient in patients_9p_del_cp:
    if patient not in patients_CDKN2A_del_cp:
        chr9p_exception_cp.append(patient)
chr9p_exception_cp
exception_df_cp = chr9p_cptac[chr9p_cptac['ID'].isin(chr9p_exception_cp)]

result_df_cp = []

for patient_id, patient_data in exception_df_cp.groupby('ID'):
    # Initialise a dictionary to store the results for this patient
    patient_results_cp = {'ID': patient_id}
    
    # Check each region of interest
    for region_name, (roi_start, roi_end) in ROI.items():
        deletion_found = any(
            (row['loc.start'] <= roi_end) and 
            (row['loc.end'] >= roi_start) and 
            (row['is_deletion'] == True)
            for _, row in patient_data.iterrows()
        )
        # Store the result
        patient_results[region_name] = deletion_found
    
    # Append the results for this patient to the DataFrame
    result_df_cp = result_df_cp.append(patient_results, ignore_index=True)
result_df_cp = result_df_cp.set_index('ID')
result_df_cp