#Removing overlaps in segemented data. GISTIC2.0 assumes that each genomic position is covered by, at most, one segment per sample. Consequently, GISTIC2.0 can only run on 
#non-overlapping segmentation files to detect amplifications, deletions, and CNVs. Overlapping segmentation may reflect over-segmentation in the initial data collection process, where the
#process used to obtain the segmentation file detected too many boundaries between segments. 

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

file_path = "/Users/tanishi/Desktop/data_cna_hg38.seg"
seg_df = pd.read_csv(file_path, sep= "\t")
print(seg_df.head)
print(seg_df.tail())

#Two segments overlap if they belong to the same sample and chromosome, and if the Start position of the second segment is less than the End position of the first segment


#Detecting overlaps. Loop through the data to compare each segment with the previous. If the current and previous segment belong to the same chromosome
#and sample, and if the start position of the current segement [i] is less than the end position of the previous segment [i-1], then they overlap. 

overlaps = []
for i in range(1, len(seg_df)):
    if (seg_df.iloc[i]['ID'] == seg_df.iloc[i-1]['ID'] and
        seg_df.iloc[i]['chrom'] == seg_df.iloc[i-1]['chrom'] and
        seg_df.iloc[i]['loc.start'] < seg_df.iloc[i-1]['loc.end']):
        overlaps.append((seg_df.iloc[i-1], seg_df.iloc[i]))

# Print overlaps. Total overlaps detected: 2692
for overlap in overlaps[:5]:
    print(f"Overlap between:\n{overlap[0]}\nand\n{overlap[1]}\n")

print(f"Total overlaps detected: {len(overlaps)}")

##Summarising segmentation data. Chromosome 12 has the greatest no. of overlaps (428). Chr1, 3, 4, 6, 7, 9, and 11 have 100+ no. of overlaps. 
def count_overlaps_by_chromosome(overlaps):
    # Create a dictionary to count overlaps for each chromosome
    chromosome_counts = {}
    
    for overlap in overlaps:
        chromosome = overlap[0]['chrom']
        if chromosome in chromosome_counts:
            chromosome_counts[chromosome] += 1
        else:
            chromosome_counts[chromosome] = 1
    
    # Convert the dictionary to a dataframe
    counts_df = pd.DataFrame(list(chromosome_counts.items()), columns=['Chromosome', 'Number of Overlaps'])
    counts_df = counts_df.sort_values(by='Number of Overlaps', ascending=False).reset_index(drop=True)
    
    return counts_df

overlaps_df = count_overlaps_by_chromosome(overlaps)
print(overlaps_df)

# Focus on chromosome 9
#Filter overlaps for chr9. Total overlaps detected on chromosome 9: 185
def is_chr9(overlap):
    return overlap[0]['chrom'] == '9'

chr9_overlaps = list(filter(is_chr9, overlaps))

for overlap in chr9_overlaps[:5]:
    print(f"Overlap between:\n{overlap[0]}\nand\n{overlap[1]}\n")

print(f"Total overlaps detected on chromosome 9: {len(chr9_overlaps)}")  

#Visualisation for chr9
chr9_df = seg_df[seg_df['chrom'] == '9']

plt.figure(figsize=(14, 8))
sns.scatterplot(data=chr9_df, x='loc.start', y='seg.mean', hue='ID', palette='tab10', alpha=0.6)
plt.title('CNV Segments for Chromosome 9')
plt.xlabel('Genomic Start Position')
plt.ylabel('Segment Mean (Log2)')
plt.legend(title='Sample', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()

#Visualisation for chr9p21
# If you are interested in specific regions, such as 9p
region_start = 18000000
region_end = 35000000
chr9p_df = chr9_df[(chr9_df['loc.start'] < region_end) & (chr9_df['loc.end'] > region_start)]

plt.figure(figsize=(14, 8))
sns.scatterplot(data=chr9p_df, x='loc.start', y='seg.mean', hue='ID', palette='tab10', alpha=0.6)
plt.title('CNV Segments for Chromosome 9p Region')
plt.xlabel('Genomic Start Position')
plt.ylabel('Segment Mean (Log2)')
plt.legend(title='Sample', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()

###Visualising overlapping segments in 9p21.3
region_start = 19900000
region_end = 25600000

# Focus on chromosome 9
#Filter overlaps for chr9. Total overlaps detected on chromosome 9: 185
def is_chr9(overlap):
    return overlap[0]['chrom'] == '9'

chr9_overlaps = list(filter(is_chr9, overlaps))

# Filter for overlaps within chromosome 9 and the 9p21 region
chr9p21_overlaps = [
    overlap for overlap in chr9_overlaps 
    if overlap[0]['loc.start'] < region_end and overlap[0]['loc.end'] > region_start
]

# Extract the overlapping segments
overlapping_segments_9p21 = []
for overlap in chr9p21_overlaps:
    overlapping_segments_9p21.extend(overlap)

# Convert to DataFrame
overlapping_df_9p21 = pd.DataFrame(overlapping_segments_9p21)

# Filter for segments within the 9p21 region
overlapping_df_9p21 = overlapping_df_9p21[
    (overlapping_df_9p21['loc.start'] < region_end) & 
    (overlapping_df_9p21['loc.end'] > region_start)
]

# Plot CNV segments for overlapping regions in chromosome 9p21
plt.figure(figsize=(14, 8))
sns.scatterplot(data=overlapping_df_9p21, x='loc.start', y='seg.mean', hue='ID', palette='tab10', alpha=0.7)
plt.title('Overlapping CNV Segments for Chromosome 9p21.3 Region')
plt.xlabel('Genomic Start Position')
plt.ylabel('Segment Mean (Log2)')
plt.legend(title='Sample', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()