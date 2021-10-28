"""

Script to make CDF plot of early to late editing systems 
enabled by orthogonal Cas systems 

@author: Kimberly Truong

"""
 
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt 
import seaborn as sns 
import barcode_utils as utils 

# nice color palette 
summer_colors = sns.color_palette("summer_r", 4)

# fonts
plt.rcParams["font.family"] = "Helvetica"
plt.rc('font', size=14)

# representative samples after removing substitution noise 
samples = ["Lba2", "Sau1", "Spy3","All9"]

unedited_bar = "NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE"

# list of total barcodes in every sample
x_numBarcodes = [] 

# list of lists where every sublist is a cumulative sum of frequencies 
y_cumfreqs = []

for ss in samples: 
	readcounts_df = pd.read_csv("../data/readcounts_files/" + ss + ".allReadCounts", sep="\t", usecols=[0,2])

	# clean out substitutions and collapse duplicates 
	new_readcounts = utils.CleanSubstitutions(readcounts_df)
	
	# remove unedited barcode (there may be more reads after removing substitutions)
	unedited_ind = new_readcounts[new_readcounts['event']==unedited_bar].index 
	new_readcounts.drop(unedited_ind, inplace=True)
	new_readcounts.reset_index(drop=True, inplace=True)

	# add new proportion col after removing unedited barcode
	total_reads = new_readcounts.loc[:,'count'].sum(axis=0)
	new_readcounts["new_proportion"] = new_readcounts['count']/total_reads 

	# number of barcodes in sample
	num_bar = new_readcounts.shape[0]

	# compute the cumulative sum 
	cumfreq_arr = np.cumsum(new_readcounts['new_proportion'].values)

	# append to data structure 
	x_numBarcodes.append(num_bar)
	y_cumfreqs.append(cumfreq_arr)

labels = ['LbaCas12a', 'SauCas9', 'SpyCas9', 'All three systems']
plt.figure(figsize=(4.8,6.4))
for nbar,y_arr,label,col in zip(x_numBarcodes, y_cumfreqs, labels, summer_colors):
	plt.plot(np.arange(1,nbar+1), y_arr, label=label, linewidth=3, color=col)
plt.legend(frameon=False, loc='lower right')
plt.ylabel('Cumulative frequency')
plt.xlim([0,1000])
plt.xlabel('Unique barcodes')

plt.tight_layout()
# plt.savefig("CDF.linewidth3.legendlowerRight-v3.svg", format='svg', dpi=300)
plt.show() 




