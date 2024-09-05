# Attempting to determine proportion of patients with 9p deletions and cdkn2a/b deletions directly from cnv data. Prior analysis revealed 88 overlapping segments
# in the 9p21.3 region across patients. Cannot be a GISTIC2.0 input because that requires non-overlapping segmentation data. 

##Data from TCGA-GBM GDC accessed from cBioPortal

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np 

file_path = "/Users/tanishi/Desktop/data_cna_hg38.seg"
seg_data = pd.read_csv(file_path, sep="\t")
print(seg_data.head)

#Define and extract 9p21 and gene-specific data. Genomic coordinates obtained from the NIH Genome Data Viewer.
CDKN2A_start = 21967752
CDKN2A_end = 21995324
chr9p21_start = 19900002
chr9p21_end = 25600000

#Set deletion threshold. Obtained from prior analyses published on TCGA Copy Number Portal. 
deletion_threshold = -0.1

#Filter seg_df for chr9p data 
chr9p_seg = seg_data[(seg_data['chrom'] == "9") & (seg_data['loc.start'] >= chr9p21_start) & (seg_data['loc.end'] <= chr9p21_end)]
chr9p_seg = chr9p_seg.copy() #to avoid error when slicing from a Pandas dataframe
chr9p_seg['is_deletion'] = chr9p_seg['seg.mean'] < deletion_threshold

# Filter for CDKN2A region. 
CDKN2A_seg = chr9p_seg[(chr9p_seg['loc.start'] < CDKN2A_end) & (chr9p_seg['loc.end'] > CDKN2A_start)]
CDKN2A_seg['is_deletion'] = CDKN2A_seg['seg.mean'] < deletion_threshold

# This returns 536 segments that overlap with the CDKN2A region, which could suggest that these CNVs aren't isolated to CDKN2A but instead are a part of broader deletions affecting adjacent regions on 9p. 
# This could point to a mechanism where a larger portion of 9p, including CDKN2A, is frequently deleted in these patients. 
# Can try to find a way to visualise these segments as well. 

# Identify unique patients with 9p and CDKN2A dels 
patients_9p_del = chr9p_seg[chr9p_seg['is_deletion']]['ID'].unique() #361 patients with 9p21.3 del
patients_CDKN2A_del = CDKN2A_seg[CDKN2A_seg['is_deletion']]['ID'].unique() #328 patients with segments that overlap the cdkn2a gene
total_patients = seg_data['ID'].nunique()

# Calculate proportions
proportion_9p_del = len(patients_9p_del)/(total_patients) #60.2% of GBM patients are 9p21.3 deleted
proportion_CDKN2A_del = len(patients_CDKN2A_del)/(total_patients) #54.2% of GBM patients are CDKN2A deleted
overlap_patients = np.intersect1d(patients_CDKN2A_del, patients_9p_del)
proportion_CDKN2A_with_9p = len(overlap_patients) / len(patients_CDKN2A_del) #All CDKN2A deleted patients have a 9p21 deletion (of course)
proportion_9p_with_CDKN2A = len(overlap_patients) / len(patients_9p_del) #90% of 9p21.3 deleted patients are CDKN2A deleted

#Save results
output_df = pd.DataFrame({ 
    'Total_Patients': [len(seg_data['ID'].unique())],
    'Patients_with_CDKN2A_deletions': [len(patients_CDKN2A_del)],
    'Patients_with_9p21_deletions': [len(patients_9p_del)],
    'Proportion_CDKN2A_with_9p21': [proportion_CDKN2A_with_9p],
    'Proportion_9p21_with_CDKN2A': [proportion_9p_with_CDKN2A]
})

print(output_df)
output_df.to_csv('deletion_analysis_results.csv', index=False)

# GDC results: 90.1% of 9p21 deleted GBM patients have a CDKN2A deletion. Thus, the CDKN2A region is included when 9p deletions occur, but not exclusively. 
# What about the remaining 10% of 9p deleted GBM patients that do not have a CDKN2A deletion? 
# Note: the proportion of 9p-deleted GBM patients with CDKN2A deletions changes to 0.86 when the deletion threshold is changed to -0.3 from the original -0.1 

# What proportion of individuals with CDKN2A/B deletion are truly 9p deleted? 

# To understand this, check if CDKN2A deleted patients have substantial deletions in other segments on the chr9p21.3 arm. 
# So, define other regions in the chr9p21.3 arm and see if the patients have deletions here as well. 

# Broadly define key chr9p21.3 genes/regions in dictionary. Defining 'start' as the first 1kb of 9p21.3 with micellaneous genes. 
# Broadly defined IFN region as region with lots of IFN region. Specify later. Starts with IFNB1 start position and ends with IFNE end position. 
ROI = {'start': (19900002, 20900002),
       'IFN' : (21077104, 21482313),
       'MTAP': (21802636, 21941115),
       'CDKN2B' : (22002903, 22009313),
       'end': (22200000, 25600000)}
patients_CDKN2A_del = CDKN2A_seg[CDKN2A_seg['is_deletion']]['ID'].unique()
deletion_status = pd.DataFrame(index=patients_CDKN2A_del, columns=ROI.keys())

# Iterate over each CDKN2A-del patient to check if any other segment in the patient overlaps with the ROIs and is deleted. 
for patient_id in patients_CDKN2A_del:
    patient_data = chr9p_seg[chr9p_seg['ID'] == patient_id] #store results for 9p21 patients

    #check each ROI for deleted segments and append to deletion_status dataframe
    for region, (start, end) in ROI.items():
        region_deleted = patient_data[
            (patient_data['is_deletion'] == True) &
            (patient_data['loc.start'] <= end) &
            (patient_data['loc.end'] >= start)
        ]
        deletion_status.loc[patient_id, region] = len(region_deleted) > 0

print(deletion_status)

# Determine the proportion of patients with combinations of deleted 9p21.3 genes/regions. 
deletion_combos = []

for index, row in deletion_status.iterrows():
    combo = []
    for region in ROI.keys():
        if row[region]:
            combo.append(region)
    if combo:
        deletion_combos.append(', '.join(combo))
    else:
        deletion_combos.append('None')

deletion_status['deletion_combo'] = deletion_combos
summary_table = deletion_status['deletion_combo'].value_counts().reset_index()
summary_table.columns = ['Deletion Combination', 'Number of Patients']
summary_table['Proportion'] = summary_table['Number of Patients'] / len(deletion_status)

print(summary_table)
summary_table.to_csv('/Users/tanishi/Desktop/summary_table.csv', index=False)

# Visualise summary table 
plt.figure(figsize=(20, 6))
plt.barh(summary_table['Deletion Combination'], summary_table['Proportion'], color='salmon')
plt.title('Proportion of CDKN2A-del Patients with Additional 9p21 Deletions', fontsize=20, weight='bold', pad=20)
plt.xlabel('Proportion of CDKN2A-del Patients', fontsize=16, weight = 'bold', pad=10)
plt.ylabel('Deletion Combination', fontsize=16, weight = 'bold', pad = 10)
plt.yticks(rotation=20, ha='right', fontsize=12)
plt.gca().invert_yaxis()
plt.xticks(fontsize=14)
plt.grid(True, axis='x', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.show()


# Validating with another TCGA dataset. Using TCGA-GBM data from CPTAC accessed from cBioportal. Smaller cohort of 98 patients. 

# Focusing on CDKN2A-del patients vs 9p21 del patients 
# Using the same genomic coordinates and deletion threshold as defined before. 
 
cptac_data = pd.read_csv("/Users/tanishi/Desktop/data_cna_hg19.seg", sep="\t")
print(cptac_data.head())

#Filter cptac_data for chr9p21 data and deletion
chr9p_cptac = cptac_data[(cptac_data['chrom'] == "9") & (cptac_data['loc.start'] >= chr9p21_start) & (cptac_data['loc.end'] <= chr9p21_end)]
print(chr9p_cptac.tail)

chr9p_cptac = chr9p_cptac.copy()
chr9p_cptac['is_deletion'] = chr9p_cptac['seg.mean'] < deletion_threshold
print(chr9p_cptac.head())

#Filter for CDKN2A del only. 
# 1 patient has a deletion within the CDKN2A gene region only (code not shown). 51 patients have a segment overlapping the CDKN2A region that is deleted. 
CDKN2A_cptac = chr9p_cptac[(chr9p_cptac['loc.start'] <= CDKN2A_end) & (chr9p_cptac['loc.end'] >= CDKN2A_start)]
CDKN2A_cptac['is_deletion'] = CDKN2A_cptac['seg.mean'] < deletion_threshold
print(CDKN2A_cptac) 
print(len(CDKN2A_cptac[CDKN2A_cptac['is_deletion']]['ID'].unique()))

# Identify unique patients with 9p and CDKN2A dels 
patients_9p_del_cp = chr9p_cptac[chr9p_cptac['is_deletion']]['ID'].unique() 
patients_CDKN2A_del_cp = CDKN2A_cptac[CDKN2A_cptac['is_deletion']]['ID'].unique() 
total_patients_cp = cptac_data['ID'].nunique() 

proportion_9p_del_cp = len(patients_9p_del_cp)/(total_patients_cp) #60.2% of GBM patients are 9p21.3 deleted
proportion_CDKN2A_del_cp = len(patients_CDKN2A_del_cp)/(total_patients_cp) #52% of GBM patients have a CDKN2A/B codel
overlap_patients_cp = np.intersect1d(patients_CDKN2A_del_cp, patients_9p_del_cp)
proportion_CDKN2A_with_9p_cp = len(overlap_patients_cp) / len(patients_CDKN2A_del_cp)
proportion_9p_with_CDKN2A_cp = len(overlap_patients_cp) / len(patients_9p_del_cp)

# Results
print(f"Total number of patients: {len(cptac_data['ID'].unique())}") #98 total patients
print(f"Number of patients with 9p deletions: {len(patients_9p_del_cp)}") #59 patients out of 98 are 9p21.3 deleted
print(f"Number of patients with CDKN2A deletions: {len(patients_CDKN2A_del_cp)}") #51 patients out of 98 are CDKN2A deleted

print(f"Proportion of CDKN2A-deleted patients with 9p deletions: {proportion_CDKN2A_with_9p_cp:.2f}") #All CDKN2A deleted patients are 9p21.3 deleted 
print(f"Proportion of 9p-deleted patients with CDKN2A deletions: {proportion_9p_with_CDKN2A_cp:.2f}") #86% of 9p21.3 deleted patients are CDKN2A deleted. 

#Similarly, 86% of 9p21.3 deleted patients in this cohort are CDKN2A deleted. 
#What about the remaining 14% of 9p21.3 deleted patients? 
#Need to repeat CDKN2A-del patient analysis. 