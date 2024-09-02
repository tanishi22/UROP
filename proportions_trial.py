##Attempting to determine proportion of patients with 9p deletions and cdkn2a/b deletions directly from cnv data. Prior analysis revealed 88 overlapping segments 
##in the 9p21.3 region across patients. Cannot be a GISTIC2.0 input because that requires non-overlapping segmentation data. 

##Again, using data from TCGA-GBM GDC accessed from cBioPortal

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np 

file_path = "/Users/tanishi/Desktop/data_cna_hg38.seg"
seg_data = pd.read_csv(file_path, sep="\t")
print(seg_data.head)

#Define and extract 9p and gene-specific data. Note, chr9p refers to all of chr9p. Genomic coordinates obtained from the NIH Genome Data Viewer.
chr9p_start = 1
chr9p_end = 43000000
CDKN2A_start = 21967752
CDKN2A_end = 21995324

#Filter seg_df for chr9p data 
chr9p_seg = seg_data[(seg_data['chrom'] == "9") & (seg_data['loc.start'] >= chr9p_start) & (seg_data['loc.end'] <= chr9p_end)]
print(chr9p_seg.tail)

#Set threshold for deletion and apply to chr9p data
deletion_threshold = -0.3
chr9p_seg = chr9p_seg.copy()
chr9p_seg['is_deletion'] = chr9p_seg['seg.mean'] < deletion_threshold
print(chr9p_seg.head())

# Filter for CDKN2A region. This returns 536 segments that overlap with the CDKN2A region, which could suggest that these CNVs aren't isolated to CDKN2A
#but are a part of broader deletions affecting adjacent regions on 9p. 
#This could point to a mechanism where a larger portion of 9p, including CDKN2A, is frequently deleted in these patients. 
#Can try to find a way to visualise these segments as well. 

CDKN2A_seg = chr9p_seg[(chr9p_seg['loc.start'] < CDKN2A_end) & (chr9p_seg['loc.end'] > CDKN2A_start)]
CDKN2A_seg['is_deletion'] = CDKN2A_seg['seg.mean'] < deletion_threshold
CDKN2A_seg

#Identify unique patients with 9p and CDKN2A dels 
patients_9p_del = chr9p_seg[chr9p_seg['is_deletion']]['ID'].unique()
patients_CDKN2A_del = CDKN2A_seg[CDKN2A_seg['is_deletion']]['ID'].unique()
total_patients = seg_data['ID'].nunique()

proportion_9p_del = len(patients_9p_del)/(total_patients)
proportion_CDKN2A_del = len(patients_CDKN2A_del)/(total_patients)

##What proportion of individuals with CDKN2A/B deletion are truly 9p deleted? 
overlap_patients = np.intersect1d(patients_CDKN2A_del, patients_9p_del)
proportion_CDKN2A_with_9p = len(overlap_patients) / len(patients_CDKN2A_del)
proportion_9p_with_CDKN2A = len(overlap_patients) / len(patients_9p_del)

##Results
print(f"Total number of patients: {len(seg_data['ID'].unique())}")
print(f"Number of patients with CDKN2A deletions: {len(patients_CDKN2A_del)}")
print(f"Number of patients with 9p deletions: {len(patients_9p_del)}")
print(f"Number of patients with both CDKN2A and 9p deletions: {len(overlap_patients)}")
print(f"Proportion of CDKN2A-deleted patients with 9p deletions: {proportion_CDKN2A_with_9p:.2f}")
print(f"Proportion of 9p-deleted patients with CDKN2A deletions: {proportion_9p_with_CDKN2A:.2f}")

#Save results
output_df = pd.DataFrame({
    'Total_Patients': [len(seg_data['ID'].unique())],
    'Patients_with_CDKN2A_deletions': [len(patients_CDKN2A_del)],
    'Patients_with_9p_deletions': [len(patients_9p_del)],
    'Patients_with_both_deletions': [len(overlap_patients)],
    'Proportion_CDKN2A_with_9p': [proportion_CDKN2A_with_9p],
    'Proportion_9p_with_CDKN2A': [proportion_9p_with_CDKN2A]
})

print(output_df)
output_df.to_csv('deletion_analysis_results.csv', index=False)

#Overall results: 89% of 9p deleted GBM patients have a CDKN2A deletion. Thus, the CDKN2A region is included when 9p deletions occur, but not exclusively. 
#What about the 11% of 9p deleted GBM patients that do not have a CDKN2A deletion? 

#Note: the propotion of 9p-deleted GBM patients with CDKN2A deletions changes to 0.86 when the deletion threshold is changed to -0.3 from the original -0.1 
#The proportion of CDKN2A-del patients with 9p-del remains 1.00, as before, but the number of patients with CDKN2A and 9p deletions decreases to 429 and 498 patients, respectively
#from the original 451 and 507 patients, respectively. 