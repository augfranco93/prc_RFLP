#--------------------
#Project: PCR-RFLP simulation
#Name: Augusto Franco
#--------------------

from Bio import *
from Bio import SeqIO
from Bio.Seq import Seq
import math
from Bio.Restriction import *
from pandas import DataFrame
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np
#Data
data = pd.read_csv("/home/augfranco/Documents/Methylobacterium/PCR-RFLP/16S_sequences/arb-silva.de_testprime_taxlist_966407.csv")
# data as a data frame
sequences = pd.DataFrame(data)

#Filtering for only Methylobacterium species. 
df =sequences.filter(["primaryAccession", "path", "organismName","length"])
df2 = df[df["path"].str.contains("Methylobacterium")]

# Dataframe to a csv file
df2.to_csv("/home/augfranco/Documents/Methylobacterium/PCR-RFLP/16S_sequences/filter_data.csv")

# Filter data
Obtain Methylobacterium species that were present 3 or more times at the PCR results

#---------Filter by primaryAccession >=3 ---------------------

df_counts = df2.groupby(['primaryAccession']).size().reset_index(name='counts')
df_filter= df_counts[df_counts.counts>2]

# Dataframe to a csv file
df_filter.to_csv("/home/augfranco/Documents/Methylobacterium/PCR-RFLP/16S_sequences/Data_filter_Counts.csv")


# Restriction digest fuction from Bactome library
----Function-----------------------------------
def restriction_digest(seq, enzyme, max_band=23130, min_band=2000, 
                       linear=False, p='yes'):
    """
    Performs restriction endonuclease digestion on a sequence and group the
    resulting fragments into 3 groups so simulate different agarose gel
    electrophoresis:
    1. fragment length more than the maximum size
    2. fragment length between the maximum and minimum size
    3. fragment length less than the minimum size
    
    Parameters:
        seq = DNA sequence for restriction endonuclease digestion
        enzyme = Restriction endonuclease object from Bio.Restriction
            package
        max_band = size of maximum band in basepairs. Default = 23130
        min_band = size of minimum band in basepairs. Default = 2000
        linear = flag to define if DNA sequence is linear.
            Default = False (DNA is circular)
        p = flag to determine if the data is to be printed. Default = yes
        
    Result:
        (Number of fragments after digestion, 
        List of fragments with molecular size above max_band, 
        List of fragments with molecular size between max_band and min_band, 
        List of fragments with molecular size below min_band)
    """
    digest = enzyme.search(seq, linear=linear)
    digest.sort()
    fragment = [digest[x+1] - digest[x]
                for x in range(len(digest) - 1)]
    fragment.sort()
    ogel = [x for x in fragment if x > max_band]
    gel = [x for x in fragment if x <= max_band and x >= min_band]
    ugel = [x for x in fragment if x < min_band]
    ogel.sort()
    gel.sort()
    ugel.sort()
    if p == 'yes':
        print('Enzyme: ' + str(enzyme))
        print('Restriction site: ' + enzyme.site)
        print('Number of fragments: ' + str(len(fragment)))
        print('Number of fragments (x > ' + str(max_band) + '): ' + \
            str(len(ogel)))
        print('Number of fragments (' + str(max_band) + ' < x < ' + \
            str(min_band) + '): ' + str(len(gel)))
        print('Number of fragments (x < ' + str(min_band) + '): ' + \
            str(len(ugel)))
        
    #return (len(fragment), ogel, gel, ugel)
    return (str(enzyme),len(fragment), fragment)

# PCR-RFLP analysis

#----------------------------------PCR-RFLP analysis---------------------------------------
results = [] # empty list for restriction digest results
ids = [] # empty list for Ids
test = [BstUI,CfoI,RsaI,AluI] # list of endonucleases to test
accessions = []    
for i in (df_filter['primaryAccession']):
    
    for sequences in SeqIO.parse("C:/Users/XPS 9590/OneDrive/Documentos/Personal Augusto Franco/Documentos personales/InPP/Methylobacterium/16S_sequences/Filter_sequences.fasta","fasta"):
        
        if i in (sequences.id):
            
              
            for enzymes in test:
                #(restriction_digest(sequences.seq,enzymes))
                x =  restriction_digest(sequences.seq,enzymes,p=False, max_band=3000, min_band=0)
                results.append(x)
                y = sequences.description
                ids.append(y)
                accessions.append(i)

#---------------------------Results as DataFrames------------------------------------------
df = DataFrame(results, columns=["Enzyme","N_fragments","fragments_Size"])
df2 = DataFrame(ids, columns=["Organism"])
df3 = DataFrame(accessions, columns=['Accession'])
df_row = pd.concat([df3,df2,df],axis=1)
df_final = df_row.drop_duplicates(subset=['Accession','Organism','Enzyme'])
df_final.to_csv("/home/augfranco/Documents/Methylobacterium/PCR-RFLP/16S_sequences/PCR_results.csv") # DataFrame to a CSV
df_separated = df_final.explode('fragments_Size').reset_index(drop=True)
df_separated.to_csv("/home/augfranco/Documents/Methylobacterium/PCR-RFLP/16S_sequences/Final_PCR-Results.csv") # DataFrame to a CSV
