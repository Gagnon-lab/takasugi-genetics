"""
Created on Sept 28, 2021

@author: Kimberly Truong

Script to compute cosine similarity matrix among single-edited embryos
(targeting Lba, Sau, and Spy individually) and triple-edited embryos

Key features:
-Pipeline removes all substitutions prior to calculation 
-computation takes into account of relative abundances 
-user needs to explicitly set "samples" variable to embryos of interest

"""

import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns
import barcode_utils as utils 

unedited_bar = 'NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE'

## make numpy matrix symmetrical across diagonal
def symmetrize(A):
    return A + A.T - np.diag(A.diagonal())

# sample names 
controls = ["Control" + str(ii) for ii in range(1,15)]
Lba_edited = ["Lba" + str(ii) for ii in range(1,9)]
Sau_edited = ["Sau" + str(ii) for ii in range(1,9)]
Spy_edited = ["Spy" + str(ii) for ii in range(1,8)]
All_edited = ["ALL" + str(ii) for ii in range(1,36)]

# set samples to any set of single-edited or triple-edited embryos 
samples = Lba_edited

num_samples = len(samples)
cosim_arr = np.zeros((num_samples,num_samples))

# loop through every unique sample pairwise comparison 
for ii in range(0,num_samples): 
    for jj in range(ii,num_samples): 

        # load ReadCounts files of pairs to be compared 
        df1 = pd.read_csv("../data/readcounts_files/" + str(samples[ii]) + ".allReadCounts", sep="\t", usecols=[0,2])
        df2 = pd.read_csv("../data/readcounts_files/" + str(samples[jj]) + ".allReadCounts", sep="\t", usecols=[0,2])

        # clean out substitutions and collapse duplicates
        newdf1 = utils.CleanSubstitutions(df1)
        newdf2 = utils.CleanSubstitutions(df2)

        # load event column as sets; this is more efficient for rapid lookup
        eventsetA = set(newdf1['event'])
        eventsetB = set(newdf2['event'])
        
        # remove unedited barcode
        if unedited_bar in eventsetA: 
            newdf1 = newdf1[newdf1['event']!=unedited_bar]
            eventsetA.remove(unedited_bar)

        if unedited_bar in eventsetB: 
            newdf2 = newdf2[newdf2['event']!=unedited_bar]
            eventsetB.remove(unedited_bar)
        
        total_reads1 = newdf1.loc[:,'count'].sum(axis=0)
        newdf1["new_proportion"] = newdf1['count']/total_reads1

        total_reads2 = newdf2.loc[:,'count'].sum(axis=0)
        newdf2["new_proportion"] = newdf2['count']/total_reads2

        # form dictionaries mapping barcodes to proportions post-filtering for rapid lookup
        dictA = dict(newdf1[['event', 'new_proportion']].values)
        dictB = dict(newdf2[['event', 'new_proportion']].values)

        # number of unique barcodes across 2 samples 
        uniqBarcodes = eventsetA.union(eventsetB)
        uniqBarcodes = list(uniqBarcodes) # form a list so that every position corresponds to a unique barcode 
        
        # create two column vectors to represent barcodes present in samples A and B
        vA = np.zeros(len(uniqBarcodes), dtype=float)
        vB = np.zeros(len(uniqBarcodes), dtype=float)
        
        # this may take a while due to O(n) for list searching 
        for bb in eventsetA: 
            k = uniqBarcodes.index(bb)
            vA[k] = dictA[bb]

        for bb in eventsetB:
            k = uniqBarcodes.index(bb)
            vB[k] = dictB[bb]
            
        cosim = np.dot(vA,vB) /  (np.sqrt(np.dot(vA,vA)) * np.sqrt(np.dot(vB,vB)))
        cosim_arr[ii][jj] = cosim 

# we are computing half of the matrix so we need to reflect entries across y=x
cosim_arr = symmetrize(cosim_arr)
cosim_uppertri = np.triu(cosim_arr, k=1)

# make heatmap
cosim_heatmap = sns.heatmap(cosim_arr, mask=cosim_uppertri, cmap='viridis_r', square=True, annot=True, fmt='.2f', vmin=0, center=0.5, vmax=1, linewidths=0.5)
plt.xticks(np.arange(0,len(samples))+0.5, samples)
plt.yticks(np.arange(0,len(samples))+0.5, samples)
cosim_heatmap.set_xticklabels(cosim_heatmap.get_xticklabels(), rotation=45)
cosim_heatmap.set_yticklabels(cosim_heatmap.get_yticklabels(), rotation=0)

plt.show()
