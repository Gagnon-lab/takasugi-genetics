"""

Useful functions for barcode preprocessing
@author: Kimberly Truong

Sept 24, 2021 

"""

from collections import Counter 

def MismatchExists(siteEdit):
    if 'S' in siteEdit: 
        return True
    else: 
        return False

def is_multiallelic(siteEdit): 
	if '&' in siteEdit: 
		return True 
	else:
		return False
    
def RemoveSubstitutions(barcode):
    edits_arr = barcode.split("_")

    valid_edits_arr = []
    for x in edits_arr: 
        
        # edit at site is a valid edit i.e. deletion or insertion 
        if not MismatchExists(x):
            valid_edits_arr.append(x)
        
        # edit at site contains at least one substitution 
        else:
            # substitution is found at site alone 
            if not is_multiallelic(x): 
                valid_edits_arr.append('NONE')
            # substitution is found alongside other potentially valid edits 
            else: 
                alleles = x.split('&')
                # make a copy of the list to ensure going through each element of the list 
                alleles_copy = alleles.copy()
                for aa in alleles_copy:
                    if MismatchExists(aa):
                        alleles.remove(aa)
                valid_edits_arr.append('&'.join(alleles)) 
    
    return '_'.join(valid_edits_arr)

def CleanSubstitutions(dataframe): 
    
    dfcopy = dataframe.copy()
    
    # find subset of dataframe containing substitutions and clean those rows
    df_toclean = dfcopy[dfcopy['event'].apply(MismatchExists)]
    
    for index,row in df_toclean.iterrows(): 
        newEvent = RemoveSubstitutions(row['event'])

        # modify original dataframe (the deep copy) at given row  
        dfcopy.at[index,'event'] = newEvent 

    # collapse duplicates; this sets the 'event' column as the index
    result = dfcopy.groupby('event').agg('sum')
    result.reset_index(drop=False, inplace=True) # change index back to numerical
    result.sort_values(by=['count'], axis=0, ascending=False, inplace=True, na_position='last')
        
    return result

def countEdits(barcodeString): 
	edits_by_site = barcodeString.split("_")
	edits_set = set(edits_by_site)

	edits_set.discard('NONE')
	uniqEdits_list = list(edits_set)

	allEdits = []
	for x in uniqEdits_list:
		if is_multiallelic(x):
			allEdits.extend(x.split("&"))
		else: 
			allEdits.append(x)
  	
	return len(allEdits)

def countEditsbyType(barcode, editType):
    siteEdits = barcode.split("_")
    
    # remove non-edits
    while 'NONE' in siteEdits:
        siteEdits.remove('NONE')

    # returns dictionary where keys are unique edits and values are essentially number of sites spanned by the edit 
    counter_dict = Counter(siteEdits)
    
    # identify edits by whether they are multisite or indels
    indels = []
    multisite_edits = []
    for key,value in counter_dict.items():
        if value > 1:
            multisite_edits.append(key)
        elif value == 1: 
            indels.append(key)
        else:
            pass 
    
    numIndels, numMultiSites = len(indels), len(multisite_edits)
    
    if editType == 'indel': 
        return numIndels
    elif editType == 'multisite':
        return numMultiSites
    else: 
        print('Not a valid edit type! Specify either \'indel\' or \'multisite\'' )