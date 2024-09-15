# Co-deletion Co-occurrence Matrix 
# Determining the co-occurence of frequent GBM CNVs (CDKN2A/B codel, MTAP del, PTPRD del, and EGFR amp)

#Install packages
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Extract patient IDs from mutation data. Only process and analyse data for which segmentation and mutation information is available to make downstream analysis simpler.  
mutation_file_path = "/Users/tanishi/Desktop/data_mutations.csv"
mutations_df = pd.read_csv(mutation_file_path, delimiter=',', skiprows=2)
mutations_df = mutations_df[(mutations_df['Annotation_Status'] == "SUCCESS") & (mutations_df["Mutation_Status"] == "Somatic")]
mutation_ids = mutations_df["Tumor_Sample_Barcode"].unique()

# Load CNV and filter segmentation data
cnv_file_path = "/Users/tanishi/Desktop/data_cna_hg38.seg"
seg_data = pd.read_csv(cnv_file_path, sep="\t")
filtered_seg_df = seg_data[seg_data['ID'].isin(mutation_ids)]
seg_filt_ids = filtered_seg_df["ID"].nunique() #Mutation and segmentation data available for 374 patients

# Define thresholds based on most recent analyses published on the TCGA Copy Number Portal.
deletion_threshold = -0.1
amp_threshold = 0.1

# Define function to extract segmentation data for each genomic region. EGFR is the only typically amplified gene from this GOI list. 
# Process data as sets to find the intersection between groups of patients later.  
def get_patient_set(start, end, chrom, threshold, is_amplification=False):
    if is_amplification:
        segment = filtered_seg_df[(filtered_seg_df['chrom'] == chrom) & 
                                  (filtered_seg_df["loc.start"] <= end) & 
                                  (filtered_seg_df['loc.end'] >= start)]
        segment['is_amplification'] = segment['seg.mean'] >= threshold
        return set(segment[segment['is_amplification']]['ID'].unique())
    else:
        segment = filtered_seg_df[(filtered_seg_df['chrom'] == chrom) & 
                                  (filtered_seg_df["loc.start"] <= end) & 
                                  (filtered_seg_df['loc.end'] >= start)]
        segment['is_deletion'] = segment['seg.mean'] < threshold
        return set(segment[segment['is_deletion']]['ID'].unique())

# Apply function to define CNV patient sets.  
CDKN2A_del_patients = get_patient_set(21967752, 21995324, "9", deletion_threshold)
CDKN2B_del_patients = get_patient_set(22002903, 22009313, "9", deletion_threshold)
codel_patients = CDKN2A_del_patients.intersection(CDKN2B_del_patients) #288 CDKN2A/b codel patients

MTAP_del_patients = get_patient_set(21802636, 21941115, "9", deletion_threshold) #288 MTAP-del patients
PTPRD_del_patients = get_patient_set(8314246, 10613002, "9", deletion_threshold) #184 PTPRD-del patients
EGFR_amp_patients = get_patient_set(55086710, 55279321, "7", amp_threshold, is_amplification=True) #353 EGFR-amp patients

total_patients = filtered_seg_df['ID'].unique() #374 total patients
total_patients_set = set(total_patients) 

#Create dictionary of patient sets for computing into co-occurence matrix 
sets_dict = {
    'CDKN2A/B codel': codel_patients,
    'MTAP del': MTAP_del_patients,
    'PTPRD del': PTPRD_del_patients,
    'EGFR amp': EGFR_amp_patients,
    'Total patients with CNV': total_patients_set
}

# Compute co-occurrence matrix
def compute_co_occurrence_matrix(sets_dict):
    results = pd.DataFrame(index=sets_dict.keys(), columns=sets_dict.keys())
    for row in sets_dict:
        for col in sets_dict:
            if row == col:
                results.at[row, col] = "NA"
            else:
                intersection_count = len(sets_dict[row].intersection(sets_dict[col]))
                results.at[row, col] = intersection_count
    return results

co_occurrence_df = compute_co_occurrence_matrix(sets_dict)
print(co_occurrence_df) #Co-occurrence count data. 

#Compute co-occurrence proportions matrix
def compute_co_occurrence_proportions(sets_dict):
    results = pd.DataFrame(index=sets_dict.keys(), columns=sets_dict.keys())
    for row in sets_dict:
        for col in sets_dict:
            if row == col:
                results.at[row, col] = 100.0
            else:
                intersection_count = len(sets_dict[row].intersection(sets_dict[col]))
                region_count = len(sets_dict[row])
                results.at[row, col] = (intersection_count/region_count) * 100
    return results

co_occurrence_props = compute_co_occurrence_proportions(sets_dict)
co_occurrence_props = co_occurrence_props.drop(index = "Total patients with CNV").drop(columns = "Total patients with CNV") #Not needed in the final matrix. 

# Visualise the co-occurrence proportions as a heatmap with labelled proportions. 
sns.heatmap(
    co_occurrence_props.astype(float), 
    annot=True,
    cmap='YlGnBu',
    fmt='.2f',  # Format as a percentage with two decimal places
    cbar=True,
    annot_kws={"size": 12}
)
plt.title('Co-occurrence Matrix of CNV Combinations', weight='bold', fontsize=14)
plt.show()

#Save output 
co_occurrence_df.to_csv('codeletion_cooccurrence_counts.csv', header = True)
co_occurrence_props.to_csv('codeletion_cooccurrence_props.csv', header = True)

# -------------------------------------------------------------------- Mutation Data Analysis -------------------------------------------------------------------- #

#Define commonly mutated genes in GBM and extract patient data for each gene
mutations_df = mutations_df[(mutations_df['Annotation_Status'] == "SUCCESS") & (mutations_df["Mutation_Status"] == "Somatic")] 
GOI = ['TP53', 'PTPRD', 'PTEN', 'EGFR', 'ATRX', 'TERT', 'IDH1', 'IDH2']
mutations_filtered = mutations_df[mutations_df['Hugo_Symbol'].isin(GOI)]

# Create sets of patient IDs for each mutation
mutation_sets_dict = {}
for gene in GOI:
    mutation_sets_dict[gene + '_mut'] = set(mutations_filtered[mutations_filtered['Hugo_Symbol'] == gene]['Tumor_Sample_Barcode'].unique())

# Create combined set for patients with CDKN2A/B and MTAP deletions. Seems redundant to keep them separate in the co-occurrence matrices if they are almost always co-deleted.
MTAP_codel_patients = codel_patients.intersection(MTAP_del_patients)

# Combine the co-deletion/amplification sets with the mutation sets dictionary
combined_sets_dict = {
    'CDKN2A/B-MTAP codel': MTAP_codel_patients,
    'PTPRD del': PTPRD_del_patients,
    'EGFR amp': EGFR_amp_patients
}
combined_sets_dict.update(mutation_sets_dict) 

# Compute the counts co-occurrence matrix 
def compute_co_occurrence_matrix(sets_dict):
    results = pd.DataFrame(index=sets_dict.keys(), columns=sets_dict.keys())
    for row in sets_dict:
        for col in sets_dict:
            if row == col:
                results.at[row, col] = len(sets_dict[row])
            else:
                intersection_count = len(sets_dict[row].intersection(sets_dict[col]))
                results.at[row, col] = intersection_count
    return results

co_occurrence_mut_df = compute_co_occurrence_matrix(combined_sets_dict)

# Compute the proportions co-occurrence matrix 
def compute_co_occurrence_proportions(sets_dict):
    results = pd.DataFrame(index=sets_dict.keys(), columns=sets_dict.keys())
    for row in sets_dict:
        for col in sets_dict:
            if row == col:
                results.at[row, col] = 100.0
            else:
                intersection_count = len(sets_dict[row].intersection(sets_dict[col]))
                region_count = len(sets_dict[row])
                results.at[row, col] = (intersection_count/region_count) * 100
    return results

co_occurrence_mut_props = compute_co_occurrence_proportions(combined_sets_dict)

#Visualise proportions co-occurrence 
plt.figure(figsize=(12, 10))
sns.heatmap(
    co_occurrence_mut_props.astype(float), 
    annot=True,
    cmap='YlGnBu',
    fmt='.2f', 
    cbar=True,
    annot_kws={"size": 12}
)
plt.title('Co-occurrence Matrix of CNV and Mutations', weight='bold', fontsize=14)
plt.savefig("co_occurrence_matrix_CNV_mut.pdf", format='pdf', bbox_inches='tight')
plt.show()

# Visualise proportions for CDKN2A/B Codel Patients only 
co_occurrence_codel_mut_props = co_occurrence_mut_props.drop(['PTPRD del', 'EGFR amp', 'TP53_mut', 'PTPRD_mut', 'PTEN_mut', 'EGFR_mut', 'ATRX_mut', 'TERT_mut', 'IDH1_mut', 'IDH2_mut'], axis = 0)
co_occurrence_codel_mut_df = co_occurrence_mut_df.drop(['PTPRD del', 'EGFR amp', 'TP53_mut', 'PTPRD_mut', 'PTEN_mut', 'EGFR_mut', 'ATRX_mut', 'TERT_mut', 'IDH1_mut', 'IDH2_mut'], axis = 0)

co_occurrence_codel_mut_df.to_csv("CDKN2A_B_codel_mut_cooccurrence.csv", header = True)

plt.figure(figsize=(12, 10))
sns.heatmap(
    co_occurrence_codel_mut_props.astype(float), 
    annot=True,
    cmap='YlGnBu',
    fmt='.2f', 
    cbar=True,
    annot_kws={"size": 12}
)
plt.title('Co-occurrence Matrix of Mutations in CDKN2A/B Codel Patients', weight='bold', fontsize=14)
plt.savefig("co_occurrence_matrix_codel_mut.pdf", format='pdf', bbox_inches='tight')
plt.show()