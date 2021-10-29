"""

Script to compute % of single-site edits in a typical barcode
weighted by barcode abundance within a sample 

@author: Kimberly Truong


"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns
import re 
import barcode_utils as utils 

# fonts
plt.rcParams["font.family"] = "Helvetica"
plt.rc('font', size=14)

unedited_bar = "NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE"

# samples for single-edited Cas systems 
Lba_edited = ["Lba" + str(ii) for ii in range(1,9)]
Sau_edited = ["Sau" + str(ii) for ii in range(1,9)]
Spy_edited = ["Spy" + str(ii) for ii in range(1,8)]
samples = Lba_edited + Sau_edited + Spy_edited

list_samplerows = []
for ss in samples: 
	readcounts_df = pd.read_csv("../data/readcounts_files/" + str(ss) + ".allReadCounts", sep="\t", usecols=[0,2])

	# remove substitutions and collapse duplicates 
	new_readcounts = utils.CleanSubstitutions(readcounts_df)

	condition = re.search(r'([A-Za-z]+)\d+', ss).group(1)

	# remove the unedited barcode 
	unedited_ind = new_readcounts[new_readcounts['event']==unedited_bar].index 
	new_readcounts.drop(unedited_ind, inplace=True)
	new_readcounts.reset_index(drop=True, inplace=True)

	new_readcounts['numIndels'] = new_readcounts['event'].apply(utils.countEditsbyType, args=('indel',))
	new_readcounts['numEdits'] = new_readcounts['event'].apply(utils.countEdits)

	# add new proportion col after removing unedited barcode
	total_reads = new_readcounts.loc[:,'count'].sum(axis=0)
	new_readcounts["new_proportion"] = new_readcounts['count']/total_reads 
	
	# compute % of indels in a sample
	percentIndels = np.sum((new_readcounts['numIndels'].values/new_readcounts['numEdits'].values)*new_readcounts['new_proportion'].values)
	
	sample_row = [str(ss), condition, percentIndels]

	list_samplerows.append(sample_row)

df = pd.DataFrame(list_samplerows, columns=['sample_id', 'condition', 'percentIndels'])

# convert to percentage
df['percentIndels'] = 100*df['percentIndels']

# find the median for each condition
print(df.groupby("condition")['percentIndels'].median())

# swarmplot for paper 
# default fig size in matplotlib is (6.4, 4.8)
# default marker size is 5
colors = ['#cce666', '#99cc66', '#66b266']

fig = plt.figure(figsize=(6,8))
sns.swarmplot(x="condition", y="percentIndels", data=df, size=7, palette=colors)
ax = plt.gca() 
ax.set(xlabel="CRISPR System", ylabel='Percentage of barcode edits that are single-site edits')
ax.set_ylim([0,100])

ax_labels = [e.get_text() for e in ax.get_xticklabels()]
ax_ticks = ax.get_xticks()
w = 0.25
for idx,cond in enumerate(ax_labels):
    ax.hlines(df[df['condition'] == cond]['percentIndels'].median(), ax_ticks[idx]-w, ax_ticks[idx]+w, color='black', linestyles='solid')

ax.set_xticklabels(['LbaCas12a', 'SauCas9', 'SpyCas9'])

plt.tight_layout()
plt.show()


