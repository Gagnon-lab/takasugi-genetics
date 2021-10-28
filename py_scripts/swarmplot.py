"""
@author: Kimberly Truong

Script to compute the avg number of edits per barcode 
in samples across different Cas systems and deployed altogether 

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re
import seaborn as sns 
import statistics 
import barcode_utils as utils

# fonts
plt.rcParams["font.family"] = "Helvetica"
plt.rc('font', size=14)

# sample names 
controls = ["Control" + str(ii) for ii in range(1,15)]
Lba_edited = ["Lba" + str(ii) for ii in range(1,9)]
Sau_edited = ["Sau" + str(ii) for ii in range(1,9)]
Spy_edited = ["Spy" + str(ii) for ii in range(1,8)]
All_edited = ["ALL" + str(ii) for ii in range(1,36)]
samples = controls + Lba_edited + Sau_edited + Spy_edited + All_edited

unedited_bar = "NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE"

def AvgNumEdits(readcounts_df):
	# create column for number of edits represented by the barcode 
	readcounts_df["numEdits"] = readcounts_df["event"].apply(utils.countEdits)

	# add new proportion col after removing unedited barcode (or not, if sample is a control)
	total_reads = readcounts_df.loc[:,'count'].sum(axis=0)
	readcounts_df["new_proportion"] = readcounts_df['count']/total_reads 

	# calculate a weighted avg number of edits per barcode 
	avgNumEdits_wt = np.sum(readcounts_df['numEdits'].values * readcounts_df['new_proportion'].values)

	return avgNumEdits_wt

list_samplerows = []
for samp_id in samples:
	readcounts_df = pd.read_csv("../data/readcounts_files/" + samp_id + ".allReadCounts", sep="\t", usecols=[0,2])

	# clean out substitutions and collapse duplicates 
	new_readcounts = utils.CleanSubstitutions(readcounts_df)

	condition = re.search(r'([A-Za-z]+)\d+', samp_id).group(1)

	# if sample is truly edited, remove the unedited barcode from computation
	# else, sample is a control and include unedited barcode as part of the weighted average 
	if condition != 'Control': 
		# remove unedited barcode 
		unedited_ind = new_readcounts[new_readcounts['event']==unedited_bar].index 
		new_readcounts.drop(unedited_ind, inplace=True)
		new_readcounts.reset_index(drop=True, inplace=True)

	editsperBar_wt = AvgNumEdits(new_readcounts)

	sample_row = [samp_id, condition, editsperBar_wt]
	list_samplerows.append(sample_row)

## repeat for previously published GESTALT datasets 

# GESTALT/scGESTALT barcode have 10 sites in total
unedited_gestalt = "NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE"
scgestalt_samples = ["hsp5-" + str(ii) for ii in range(1,9)]

for samp_id in scgestalt_samples:
	readcounts_df = pd.read_csv("../data/scGESTALT_readcounts/" + samp_id + ".allReadCounts", sep="\t", usecols=[0,2])

	# clean out substitutions and collapse duplicates 
	new_readcounts = utils.CleanSubstitutions(readcounts_df)

	# remove the unedited barcode (there may be additional reads after removing substitutions)
	unedited_ind = new_readcounts[new_readcounts['event']==unedited_gestalt].index 
	new_readcounts.drop(unedited_ind, inplace=True)
	new_readcounts.reset_index(drop=True, inplace=True)

	editsperBar_wt = AvgNumEdits(new_readcounts)

	sample_row = [samp_id, "scGESTALT", editsperBar_wt]
	list_samplerows.append(sample_row)

# load original GESTALT data 
gestalt_df = pd.read_csv("../data/embryo_all_reads_Mar_23_2016.txt", sep='\t', usecols=[0,4,6])

for name,sampdf in gestalt_df.groupby('sample'):

	# clean out substitutions and collapse duplicates 
	new_readcounts = utils.CleanSubstitutions(sampdf)

	# remove the unedited barcode (there may be additional reads after removing substitutions)
	unedited_ind = new_readcounts[new_readcounts['event']==unedited_gestalt].index 
	new_readcounts.drop(unedited_ind, inplace=True)
	new_readcounts.reset_index(drop=True, inplace=True)

	editsperBar_wt = AvgNumEdits(new_readcounts)

	sample_row = [str(name), "GESTALT", editsperBar_wt]
	list_samplerows.append(sample_row)

Cas_edits = pd.DataFrame(list_samplerows, columns=['sample_id', 'condition', 'numEditsperBarcode_wt'])

# find the median for each condition
print(Cas_edits.groupby("condition")['numEditsperBarcode_wt'].median())

# swarmplot for paper 
# default fig size in matplotlib is (6.4, 4.8)
# default marker size is 5
colors = ['#db5f57','#cce666', '#99cc66', '#66b266', '#339966', '#6a3d9a','#cab2d6']

fig = plt.figure(figsize=(6,8))
sns.swarmplot(x="condition", y="numEditsperBarcode_wt", data=Cas_edits, size=7, palette=colors)
ax = plt.gca()
ax.set(xlabel=" ", ylabel="Average number of edits per barcode")
ymax = Cas_edits['numEditsperBarcode_wt'].max()
ax.set_ylim([0, ymax+0.75])

ax_labels = [e.get_text() for e in ax.get_xticklabels()]
ax_ticks = ax.get_xticks()
w = 0.42
for idx,cond in enumerate(ax_labels):
    ax.hlines(Cas_edits[Cas_edits['condition'] == cond]['numEditsperBarcode_wt'].median(), ax_ticks[idx]-w, ax_ticks[idx]+w, color='black', linestyles='solid')

ax.set_xticklabels(['Control', 'LbaCas12a', 'SauCas9', 'SpyCas9', 'All three systems', 'scGESTALT', 'GESTALT'], rotation=45, ha='center')
plt.tight_layout()
plt.show()










