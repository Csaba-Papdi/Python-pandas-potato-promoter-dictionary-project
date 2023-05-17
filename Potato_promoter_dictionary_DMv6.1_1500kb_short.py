#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np


# ## I. Find the location of CDS start sites (ATGs) on the chromosomes
# Genome annotation file spuddb.uga.edu/dataDM_1-3_516_R44_potato.v6.1.hc_gene_models.gff3.gz was converted to csv.

# ### Task I/1. Create a dataframe that contains the location of all CDS start site.

# In[2]:


""" import csv file and convert it to a dataframe """

data = pd.read_csv("DM_1-3_516_R44_potato.v6.1.hc_gene_models.csv")
df1 = pd.DataFrame(data)
df1.head(2)


# In[3]:


""" ATGs are on the first CDSs, drop rows if the Feature column contains the "exon" string  """

df1 = df1[~df1['Feature'].isin(['exon'])]
df1.head(2)


# In[4]:


""" split the "Attributes" column """

df1[["transcript_ID","gene_ID"]] = df1["Attributes"].str.split("; ",expand=True)  
df1.head(2)


# In[5]:


""" drop the unnecessary columns """

df2= df1.drop(['Attributes', 'transcript_ID','Source', 'Feature', 'Score', 'Frame'], axis=1)        #default is axis = 0 which is rows columns is axis = 1.
df2.head(2)


# In[6]:


""" clean the gene_ID column """

df2['gene_ID'] = df2['gene_ID'].map(lambda x: x.lstrip('gene_id "').rstrip('";'))
df2.head(2)


# ### Task I/2. Split the dataframe in two (based on the genes are located on the positive or negative DNA strand). 

# In[7]:


""" define positive strand dataframe (drop rows with -), drop Strand column """

df_pos_strand = df2.drop(df2[df2['Strand'] == '-'].index)
df_pos_strand = df_pos_strand.drop(['End','Strand'], axis=1)  
df_pos_strand.head(2)


# In[8]:


""" define negative strand dataframe (drop rows with +), drop Strand column """

df_neg_strand = df2.drop(df2[df2['Strand'] == '+'].index)
df_neg_strand = df_neg_strand.drop(['Start', 'Strand'], axis=1) 
df_neg_strand.head(2)


# ### Task I/3. Format the two new dataframes: for positive DNA strand keep the Start column, for negative DNA strand dataframe keep the End column.

# In[9]:


""" drop duplicated gene_IDs, keep the first row of each - those should be the ones starting with the ATGs """

df_pos_strand_ = df_pos_strand.drop_duplicates(subset='gene_ID', keep="first")
df_pos_strand_ = df_pos_strand_.reset_index(drop=True)
df_pos_strand.head(2)


# In[10]:


""" substract a number from each value to fix indexing problem later """

df_pos_strand_["Start"] = df_pos_strand_["Start"] - 1  
df_pos_strand_.head(2)


# In[11]:


""" drop duplicated gene_IDs, keep the last row of each those should end with CAT (ATG), drop Strand column """

df_neg_strand_ = df_neg_strand.drop_duplicates(subset='gene_ID', keep="last")
df_neg_strand_= df_neg_strand_.reset_index(drop=True)
df_neg_strand_.head(2)


# In[12]:


""" add a number to each value to fix indexing problem that would come later """

df_neg_strand_["End"] = df_neg_strand_["End"] + 1  
df_neg_strand_.head(2)


# ### Task I/4. Both DNA strand dataframes are to be separated to chromosomes.

# In[13]:


""" separate chrosomes in both strands -> 24 new DFs with the mapped ATGs """

# chromosome separation for positive DNA strand datframes

chr1_pos_mapped = df_pos_strand_[df_pos_strand_['Seqname'].str.contains("chr01")]
chr2_pos_mapped = df_pos_strand_[df_pos_strand_['Seqname'].str.contains("chr02")]
chr3_pos_mapped = df_pos_strand_[df_pos_strand_['Seqname'].str.contains("chr03")]
chr4_pos_mapped = df_pos_strand_[df_pos_strand_['Seqname'].str.contains("chr04")]
chr5_pos_mapped = df_pos_strand_[df_pos_strand_['Seqname'].str.contains("chr05")]
chr6_pos_mapped = df_pos_strand_[df_pos_strand_['Seqname'].str.contains("chr06")]
chr7_pos_mapped = df_pos_strand_[df_pos_strand_['Seqname'].str.contains("chr07")]
chr8_pos_mapped = df_pos_strand_[df_pos_strand_['Seqname'].str.contains("chr08")]
chr9_pos_mapped = df_pos_strand_[df_pos_strand_['Seqname'].str.contains("chr09")]
chr10_pos_mapped = df_pos_strand_[df_pos_strand_['Seqname'].str.contains("chr10")]
chr11_pos_mapped = df_pos_strand_[df_pos_strand_['Seqname'].str.contains("chr11")]
chr12_pos_mapped = df_pos_strand_[df_pos_strand_['Seqname'].str.contains("chr12")]

# chromosome separation for negative DNA strand datframes

chr1_neg_mapped = df_neg_strand_[df_neg_strand_['Seqname'].str.contains("chr01")]
chr2_neg_mapped = df_neg_strand_[df_neg_strand_['Seqname'].str.contains("chr02")]
chr3_neg_mapped = df_neg_strand_[df_neg_strand_['Seqname'].str.contains("chr03")]
chr4_neg_mapped = df_neg_strand_[df_neg_strand_['Seqname'].str.contains("chr04")]
chr5_neg_mapped = df_neg_strand_[df_neg_strand_['Seqname'].str.contains("chr05")]
chr6_neg_mapped = df_neg_strand_[df_neg_strand_['Seqname'].str.contains("chr06")]
chr7_neg_mapped = df_neg_strand_[df_neg_strand_['Seqname'].str.contains("chr07")]
chr8_neg_mapped = df_neg_strand_[df_neg_strand_['Seqname'].str.contains("chr08")]
chr9_neg_mapped = df_neg_strand_[df_neg_strand_['Seqname'].str.contains("chr09")]
chr10_neg_mapped = df_neg_strand_[df_neg_strand_['Seqname'].str.contains("chr10")]
chr11_neg_mapped = df_neg_strand_[df_neg_strand_['Seqname'].str.contains("chr11")]
chr12_neg_mapped = df_neg_strand_[df_neg_strand_['Seqname'].str.contains("chr12")]


# In[14]:


""" drop Seqname columns in all the 24 datafarmes """

for i in range(1, 13):
    chr_pos_mapped_ = globals().get(f"chr{i}_pos_mapped")
    chr_neg_mapped_ = globals().get(f"chr{i}_neg_mapped")
    chr_pos_mapped_ = chr_pos_mapped_.drop(['Seqname'], axis=1)
    chr_neg_mapped_ = chr_neg_mapped_.drop(['Seqname'], axis=1)
    globals()[f"chr{i}_pos_mapped_"] = chr_pos_mapped_
    globals()[f"chr{i}_neg_mapped_"] = chr_neg_mapped_
# test one of them
chr9_neg_mapped_.head(2)


# ## II. Assign chromosome sequences to variables
# Genome sequence file DM_1-3_516_R44_potato_genome_assembly.v6.1.fa is downloaded from SpuDB and uncompressed on usegalaxy.org site.

# In[15]:


""" open fasta file as a string and assign to a variable """    

filename = "DM_1-3_516_R44_potato_genome_assembly.v6.1.fa_uncompressed.fasta"
genome = open(filename, "r")
whole_seq_withLineBreaks = genome.read()
genome.close()
whole_seq = whole_seq_withLineBreaks.replace("\n", "") 


# In[16]:


""" Find chromosomes on the genomic DNA sequence by indexing the first occurrence of chr01-chr12 Substrings """

whole_seq_index1 = whole_seq.find(">chr01")
whole_seq_index2 = whole_seq.find(">chr02")
whole_seq_index3 = whole_seq.find(">chr03")
whole_seq_index4 = whole_seq.find(">chr04")
whole_seq_index5 = whole_seq.find(">chr05")
whole_seq_index6 = whole_seq.find(">chr06")
whole_seq_index7 = whole_seq.find(">chr07")
whole_seq_index8 = whole_seq.find(">chr08")
whole_seq_index9 = whole_seq.find(">chr09")
whole_seq_index10 = whole_seq.find(">chr10")
whole_seq_index11 = whole_seq.find(">chr11")
whole_seq_index12 = whole_seq.find(">chr12")
whole_seq_index13 = whole_seq.find(">sc")


# In[17]:


""" separate the chromosomes, based on the first occurance of Substrings  """

chr_1 = whole_seq[0:88591692]
chr_2 = whole_seq[88591692:134694613]
chr_3 = whole_seq[134694613:195402189]
chr_4 = whole_seq[195402189:264638526]
chr_5 = whole_seq[264638526:320238229]
chr_6 = whole_seq[320238229:379329813]
chr_7 = whole_seq[379329813:436969136]
chr_8 = whole_seq[436969136:496195142]
chr_9 = whole_seq[496195142:563795448]
chr_10 = whole_seq[563795448:624839605]
chr_11 = whole_seq[624839605:671616998]
chr_12 = whole_seq[671616998:731287759]


# In[18]:


""" delete chr labels of each string """

chr1 = chr_1[7:]
chr2 = chr_2[7:]
chr3 = chr_3[7:]
chr4 = chr_4[7:]
chr5 = chr_5[7:]
chr6 = chr_6[7:]
chr7 = chr_7[7:]
chr8 = chr_8[7:]
chr9 = chr_9[7:]
chr10 = chr_10[7:]
chr11 = chr_11[7:]
chr12 = chr_12[7:]


# ## III. Mapping CDS starting positions on the chromosomes

# ### Task III/1. Create lists of gene positions on the chromosomes.

# In[19]:


""" the Start/End columns are to be converted to a list (of integers) """

# on positive DNA strand chromosomes

chr1_pos_genes_pos = chr1_pos_mapped_['Start'].tolist()    
chr2_pos_genes_pos = chr2_pos_mapped_['Start'].tolist()
chr3_pos_genes_pos = chr3_pos_mapped_['Start'].tolist()
chr4_pos_genes_pos = chr4_pos_mapped_['Start'].tolist()
chr5_pos_genes_pos = chr5_pos_mapped_['Start'].tolist()
chr6_pos_genes_pos = chr6_pos_mapped_['Start'].tolist()
chr7_pos_genes_pos = chr7_pos_mapped_['Start'].tolist()
chr8_pos_genes_pos = chr8_pos_mapped_['Start'].tolist()
chr9_pos_genes_pos = chr9_pos_mapped_['Start'].tolist()
chr10_pos_genes_pos = chr10_pos_mapped_['Start'].tolist()
chr11_pos_genes_pos = chr11_pos_mapped_['Start'].tolist()
chr12_pos_genes_pos = chr12_pos_mapped_['Start'].tolist()

# on negative DNA strand chromosomes

chr1_neg_genes_pos = chr1_neg_mapped_['End'].tolist()    
chr2_neg_genes_pos = chr2_neg_mapped_['End'].tolist()
chr3_neg_genes_pos = chr3_neg_mapped_['End'].tolist()
chr4_neg_genes_pos = chr4_neg_mapped_['End'].tolist()
chr5_neg_genes_pos = chr5_neg_mapped_['End'].tolist()
chr6_neg_genes_pos = chr6_neg_mapped_['End'].tolist()
chr7_neg_genes_pos = chr7_neg_mapped_['End'].tolist()
chr8_neg_genes_pos = chr8_neg_mapped_['End'].tolist()
chr9_neg_genes_pos = chr9_neg_mapped_['End'].tolist()
chr10_neg_genes_pos = chr10_neg_mapped_['End'].tolist()
chr11_neg_genes_pos = chr11_neg_mapped_['End'].tolist()
chr12_neg_genes_pos = chr12_neg_mapped_['End'].tolist()


# ### Task III/2. Making the 2 x 12 lists of promoters (1500nt in length).

# Substracting 1500 from ATG positions on the positive strand

# In[20]:


""" chr 1 + """

gene_index = 0 
chr1_pos_promoters = [] # it creates an empty list
for element in chr1_pos_genes_pos:      # for the every elements of the list of integers (n = 2142)
    xld_1 = (chr1[(chr1_pos_genes_pos[gene_index]-1500):chr1_pos_genes_pos[gene_index]])
    chr1_pos_promoters.append(xld_1)
    gene_index += 1


# In[21]:


""" chr 2 + """

gene_index = 0 
chr2_pos_promoters = [] 
for element in chr2_pos_genes_pos:     
    xld_2 = (chr2[(chr2_pos_genes_pos[gene_index]-1500):chr2_pos_genes_pos[gene_index]])
    chr2_pos_promoters.append(xld_2)
    gene_index += 1


# In[22]:


""" chr 3 + """

gene_index = 0 
chr3_pos_promoters = [] # it creates an empty list
for element in chr3_pos_genes_pos:      # for the every elements of the list of integers (n = 2142)
    xld_3 = (chr3[(chr3_pos_genes_pos[gene_index]-1500):chr3_pos_genes_pos[gene_index]])
    chr3_pos_promoters.append(xld_3)
    gene_index += 1


# In[23]:


""" chr 4 + """

gene_index = 0 
chr4_pos_promoters = [] # it creates an empty list
for element in chr4_pos_genes_pos:      # for the every elements of the list of integers (n = 2142)
    xld_4 = (chr4[(chr4_pos_genes_pos[gene_index]-1500):chr4_pos_genes_pos[gene_index]])
    chr4_pos_promoters.append(xld_4)
    gene_index += 1


# In[24]:


""" chr 5 + """

gene_index = 0 
chr5_pos_promoters = [] # it creates an empty list
for element in chr5_pos_genes_pos:      # for the every elements of the list of integers (n = 2142)
    xld_5 = (chr5[(chr5_pos_genes_pos[gene_index]-1500):chr5_pos_genes_pos[gene_index]])
    chr5_pos_promoters.append(xld_5)
    gene_index += 1


# In[25]:


""" chr 6 + """

gene_index = 0 
chr6_pos_promoters = [] # it creates an empty list
for element in chr6_pos_genes_pos:      # for the every elements of the list of integers (n = 2142)
    xld_6 = (chr6[(chr6_pos_genes_pos[gene_index]-1500):chr6_pos_genes_pos[gene_index]])
    chr6_pos_promoters.append(xld_6)
    gene_index += 1


# In[26]:


""" chr 7 + """

gene_index = 0 
chr7_pos_promoters = [] # it creates an empty list
for element in chr7_pos_genes_pos:      # for the every elements of the list of integers (n = 2142)
    xld_7 = (chr7[(chr7_pos_genes_pos[gene_index]-1500):chr7_pos_genes_pos[gene_index]])
    chr7_pos_promoters.append(xld_7)
    gene_index += 1


# In[27]:


""" chr 8 + """

gene_index = 0 
chr8_pos_promoters = [] # it creates an empty list
for element in chr8_pos_genes_pos:      # for the every elements of the list of integers (n = 2142)
    xld_8 = (chr8[(chr8_pos_genes_pos[gene_index]-1500):chr8_pos_genes_pos[gene_index]])
    chr8_pos_promoters.append(xld_8)
    gene_index += 1


# In[28]:


""" chr 9 + """

gene_index = 0 
chr9_pos_promoters = [] # it creates an empty list
for element in chr9_pos_genes_pos:      # for the every elements of the list of integers (n = 2142)
    xld_9 = (chr9[(chr9_pos_genes_pos[gene_index]-1500):chr9_pos_genes_pos[gene_index]])
    chr9_pos_promoters.append(xld_9)
    gene_index += 1


# In[29]:


""" chr 10 + """

gene_index = 0 
chr10_pos_promoters = [] # it creates an empty list
for element in chr10_pos_genes_pos:      # for the every elements of the list of integers (n = 2142)
    xld_10 = (chr10[(chr10_pos_genes_pos[gene_index]-1500):chr10_pos_genes_pos[gene_index]])
    chr10_pos_promoters.append(xld_10)
    gene_index += 1


# In[30]:


""" chr 11 + """

gene_index = 0 
chr11_pos_promoters = [] # it creates an empty list
for element in chr11_pos_genes_pos:      # for the every elements of the list of integers (n = 2142)
    xld_11 = (chr11[(chr11_pos_genes_pos[gene_index]-1500):chr11_pos_genes_pos[gene_index]])
    chr11_pos_promoters.append(xld_11)
    gene_index += 1


# In[31]:


""" chr 12 + """

gene_index = 0 
chr12_pos_promoters = [] # it creates an empty list
for element in chr12_pos_genes_pos:      # for the every elements of the list of integers (n = 2142)
    xld_12 = (chr12[(chr12_pos_genes_pos[gene_index]-1500):chr12_pos_genes_pos[gene_index]])
    chr12_pos_promoters.append(xld_12)
    gene_index += 1


# Adding 1500 to ATG positions on the negative (reverse complement) strand

# In[32]:


""" chr 1 - """

gene_index = 0 
chr1_neg_promoters = [] 
for element in chr1_neg_genes_pos:     
    xlf_1 = (chr1[(chr1_neg_genes_pos[gene_index]):chr1_neg_genes_pos[gene_index]+1500])
    chr1_neg_promoters.append(xlf_1)
    gene_index += 1


# In[33]:


""" chr 2 - """

gene_index = 0 
chr2_neg_promoters = [] 
for element in chr2_neg_genes_pos:     
    xlf_2 = (chr2[(chr2_neg_genes_pos[gene_index]):chr2_neg_genes_pos[gene_index]+1500])
    chr2_neg_promoters.append(xlf_2)
    gene_index += 1


# In[34]:


""" chr 3 - """

gene_index = 0 
chr3_neg_promoters = [] 
for element in chr3_neg_genes_pos:      
    xlf_3 = (chr3[(chr3_neg_genes_pos[gene_index]):chr3_neg_genes_pos[gene_index]+1500])
    chr3_neg_promoters.append(xlf_3)
    gene_index += 1


# In[35]:


""" chr 4 - """

gene_index = 0 
chr4_neg_promoters = [] 
for element in chr4_neg_genes_pos:      
    xlf_4 = (chr4[(chr4_neg_genes_pos[gene_index]):chr4_neg_genes_pos[gene_index]+1500])
    chr4_neg_promoters.append(xlf_4)
    gene_index += 1


# In[36]:


""" chr 5 - """

gene_index = 0 
chr5_neg_promoters = [] 
for element in chr5_neg_genes_pos:      
    xlf_5 = (chr5[(chr5_neg_genes_pos[gene_index]):chr5_neg_genes_pos[gene_index]+1500])
    chr5_neg_promoters.append(xlf_5)
    gene_index += 1


# In[37]:


""" chr 6 - """

gene_index = 0 
chr6_neg_promoters = [] 
for element in chr6_neg_genes_pos:      
    xlf_6 = (chr6[(chr6_neg_genes_pos[gene_index]):chr6_neg_genes_pos[gene_index]+1500])
    chr6_neg_promoters.append(xlf_6)
    gene_index += 1


# In[38]:


""" chr 7 - """

gene_index = 0 
chr7_neg_promoters = [] 
for element in chr7_neg_genes_pos:      
    xlf_7 = (chr7[(chr7_neg_genes_pos[gene_index]):chr7_neg_genes_pos[gene_index]+1500])
    chr7_neg_promoters.append(xlf_7)
    gene_index += 1


# In[39]:


""" chr 8 - """

gene_index = 0 
chr8_neg_promoters = [] 
for element in chr8_neg_genes_pos:      
    xlf_8 = (chr8[(chr8_neg_genes_pos[gene_index]):chr8_neg_genes_pos[gene_index]+1500])
    chr8_neg_promoters.append(xlf_8)
    gene_index += 1


# In[40]:


""" chr 9 - """

gene_index = 0 
chr9_neg_promoters = [] 
for element in chr9_neg_genes_pos:      
    xlf_9 = (chr9[(chr9_neg_genes_pos[gene_index]):chr9_neg_genes_pos[gene_index]+1500])
    chr9_neg_promoters.append(xlf_9)
    gene_index += 1


# In[41]:


""" chr 10 - """

gene_index = 0 
chr10_neg_promoters = [] 
for element in chr10_neg_genes_pos:      
    xlf_10 = (chr10[(chr10_neg_genes_pos[gene_index]):chr10_neg_genes_pos[gene_index]+1500])
    chr10_neg_promoters.append(xlf_10)
    gene_index += 1


# In[42]:


""" chr 11 - """

gene_index = 0 
chr11_neg_promoters = [] 
for element in chr11_neg_genes_pos:      
    xlf_11 = (chr11[(chr11_neg_genes_pos[gene_index]):chr11_neg_genes_pos[gene_index]+1500])
    chr11_neg_promoters.append(xlf_11)
    gene_index += 1


# In[43]:


""" chr 12 - """

gene_index = 0 
chr12_neg_promoters = [] 
for element in chr12_neg_genes_pos:      
    xlf_12 = (chr12[(chr12_neg_genes_pos[gene_index]):chr12_neg_genes_pos[gene_index]+1500])
    chr12_neg_promoters.append(xlf_12)
    gene_index += 1


# ### Task III/3. Create a dictionary of gene IDs.

# In[44]:


""" the gene ID columns are to be converted to a list of strings """

# chromosome_positive-strand_genes_positions:

chr1_pos_gene_ID = chr1_pos_mapped_['gene_ID'].tolist()  
chr2_pos_gene_ID = chr2_pos_mapped_['gene_ID'].tolist()
chr3_pos_gene_ID = chr3_pos_mapped_['gene_ID'].tolist()
chr4_pos_gene_ID = chr4_pos_mapped_['gene_ID'].tolist()
chr5_pos_gene_ID = chr5_pos_mapped_['gene_ID'].tolist()
chr6_pos_gene_ID = chr6_pos_mapped_['gene_ID'].tolist()
chr7_pos_gene_ID = chr7_pos_mapped_['gene_ID'].tolist()
chr8_pos_gene_ID = chr8_pos_mapped_['gene_ID'].tolist()
chr9_pos_gene_ID = chr9_pos_mapped_['gene_ID'].tolist()
chr10_pos_gene_ID = chr10_pos_mapped_['gene_ID'].tolist()
chr11_pos_gene_ID = chr11_pos_mapped_['gene_ID'].tolist()
chr12_pos_gene_ID = chr12_pos_mapped_['gene_ID'].tolist()

# chromosome_negative-strand_genes_positions:

chr1_neg_gene_ID = chr1_neg_mapped_['gene_ID'].tolist()    
chr2_neg_gene_ID = chr2_neg_mapped_['gene_ID'].tolist()
chr3_neg_gene_ID = chr3_neg_mapped_['gene_ID'].tolist()
chr4_neg_gene_ID = chr4_neg_mapped_['gene_ID'].tolist()
chr5_neg_gene_ID = chr5_neg_mapped_['gene_ID'].tolist()
chr6_neg_gene_ID = chr6_neg_mapped_['gene_ID'].tolist()
chr7_neg_gene_ID = chr7_neg_mapped_['gene_ID'].tolist()
chr8_neg_gene_ID = chr8_neg_mapped_['gene_ID'].tolist()
chr9_neg_gene_ID = chr9_neg_mapped_['gene_ID'].tolist()
chr10_neg_gene_ID = chr10_neg_mapped_['gene_ID'].tolist()
chr11_neg_gene_ID = chr11_neg_mapped_['gene_ID'].tolist()
chr12_neg_gene_ID = chr12_neg_mapped_['gene_ID'].tolist()


# In[45]:


""" Adding the list of promoters (from Task III/2) to the dictionaries """

# positive strand

chr1_pos_gene_ID_promoters_dict = dict(zip(chr1_pos_gene_ID, chr1_pos_promoters))
chr2_pos_gene_ID_promoters_dict = dict(zip(chr2_pos_gene_ID, chr2_pos_promoters))
chr3_pos_gene_ID_promoters_dict = dict(zip(chr3_pos_gene_ID, chr3_pos_promoters))
chr4_pos_gene_ID_promoters_dict = dict(zip(chr4_pos_gene_ID, chr4_pos_promoters))
chr5_pos_gene_ID_promoters_dict = dict(zip(chr5_pos_gene_ID, chr5_pos_promoters))
chr6_pos_gene_ID_promoters_dict = dict(zip(chr6_pos_gene_ID, chr6_pos_promoters))
chr7_pos_gene_ID_promoters_dict = dict(zip(chr7_pos_gene_ID, chr7_pos_promoters))
chr8_pos_gene_ID_promoters_dict = dict(zip(chr8_pos_gene_ID, chr8_pos_promoters))
chr9_pos_gene_ID_promoters_dict = dict(zip(chr9_pos_gene_ID, chr9_pos_promoters))
chr10_pos_gene_ID_promoters_dict = dict(zip(chr10_pos_gene_ID, chr10_pos_promoters))
chr11_pos_gene_ID_promoters_dict = dict(zip(chr11_pos_gene_ID, chr11_pos_promoters))
chr12_pos_gene_ID_promoters_dict = dict(zip(chr12_pos_gene_ID, chr12_pos_promoters))


# negative strand

chr1_neg_gene_ID_promoters_dict = dict(zip(chr1_neg_gene_ID, chr1_neg_promoters))
chr2_neg_gene_ID_promoters_dict = dict(zip(chr2_neg_gene_ID, chr2_neg_promoters))
chr3_neg_gene_ID_promoters_dict = dict(zip(chr3_neg_gene_ID, chr3_neg_promoters))
chr4_neg_gene_ID_promoters_dict = dict(zip(chr4_neg_gene_ID, chr4_neg_promoters))
chr5_neg_gene_ID_promoters_dict = dict(zip(chr5_neg_gene_ID, chr5_neg_promoters))
chr6_neg_gene_ID_promoters_dict = dict(zip(chr6_neg_gene_ID, chr6_neg_promoters))
chr7_neg_gene_ID_promoters_dict = dict(zip(chr7_neg_gene_ID, chr7_neg_promoters))
chr8_neg_gene_ID_promoters_dict = dict(zip(chr8_neg_gene_ID, chr8_neg_promoters))
chr9_neg_gene_ID_promoters_dict = dict(zip(chr9_neg_gene_ID, chr9_neg_promoters))
chr10_neg_gene_ID_promoters_dict = dict(zip(chr10_neg_gene_ID, chr10_neg_promoters))
chr11_neg_gene_ID_promoters_dict = dict(zip(chr11_neg_gene_ID, chr11_neg_promoters))
chr12_neg_gene_ID_promoters_dict = dict(zip(chr12_neg_gene_ID, chr12_neg_promoters))


# In[46]:


""" Merging positive strand dictionaries """

POSITIVE_STRAND_GENE_ID_PROMOTERS_DICTIONARY = {**chr1_pos_gene_ID_promoters_dict, **chr2_pos_gene_ID_promoters_dict, 
**chr3_pos_gene_ID_promoters_dict, **chr4_pos_gene_ID_promoters_dict, **chr5_pos_gene_ID_promoters_dict,
**chr6_pos_gene_ID_promoters_dict, **chr7_pos_gene_ID_promoters_dict, **chr8_pos_gene_ID_promoters_dict,
**chr9_pos_gene_ID_promoters_dict, **chr10_pos_gene_ID_promoters_dict, **chr11_pos_gene_ID_promoters_dict,
**chr12_pos_gene_ID_promoters_dict}


# In[47]:


""" Merging negative strand dictionaries """

Negative_strand_gene_ID_promoters_dictionary = {**chr1_neg_gene_ID_promoters_dict, **chr2_neg_gene_ID_promoters_dict, 
**chr3_neg_gene_ID_promoters_dict, **chr4_neg_gene_ID_promoters_dict, **chr5_neg_gene_ID_promoters_dict,
**chr6_neg_gene_ID_promoters_dict, **chr7_neg_gene_ID_promoters_dict, **chr8_neg_gene_ID_promoters_dict,
**chr9_neg_gene_ID_promoters_dict, **chr10_neg_gene_ID_promoters_dict, **chr11_neg_gene_ID_promoters_dict,
**chr12_neg_gene_ID_promoters_dict}


# ### Task III/4. Make the reverse complement of negative DNA strand sequences.

# In[48]:


""" A list that contains only the values of the dictionary """

neg_dict_vals = Negative_strand_gene_ID_promoters_dictionary.values() # it returns a dict-value object which is a set-like object 
neg_dict_vals_list = list(neg_dict_vals)                                  #  it converted to a list


# In[49]:


""" first concatenate items in a list to a single string """

neg_dict_vals_string = ''.join(neg_dict_vals)


# In[50]:


""" make the reverse complement of the string with Biopython"""

from Bio.Seq import Seq  
seq = Seq(neg_dict_vals_string)
seq_rev_compl = seq.reverse_complement()


# In[51]:


""" convert Bio.Seq.Seq output to a string """

seq_rev_compl_str = str(seq_rev_compl)
print(type(seq_rev_compl_str))


# In[52]:


""" split string to list based on character number n=1500 """

import re
resplit = [seq_rev_compl_str for seq_rev_compl_str in re.split(r'(\w{1500})', seq_rev_compl_str) if seq_rev_compl_str]
resplit [:3]


# In[53]:


""" first, make a list that contains only the keys of the dictionary - negative strand """

neg_dict_keys = Negative_strand_gene_ID_promoters_dictionary.keys() # it returns a dict-value object which is a set-like object 
neg_dict_keys_list = list(neg_dict_keys)
ned_dict_keys_list_rev = neg_dict_keys_list.reverse ()


# In[54]:


""" Convert Two Lists into a Dictionary """

Negative_strand_gene_ID_promoters_dictionary_v2 = dict(zip(neg_dict_keys_list, resplit))


# In[55]:


""" Reverse the Order of Dictionary Keys """

NEGATIVE_STRAND_GENE_ID_PROMOTERS_DICTIONARY = dict(reversed(list(Negative_strand_gene_ID_promoters_dictionary_v2.items())))


# In[56]:


### Task III/6. Create the final dictionary.


# In[57]:


""" Merging two dictionaries """

GENE_ID_PROMOTERS_DICTIONARY = {**POSITIVE_STRAND_GENE_ID_PROMOTERS_DICTIONARY,
**NEGATIVE_STRAND_GENE_ID_PROMOTERS_DICTIONARY}


# In[58]:


# test: print 1st + 2nd element of the dictionary

print(list(GENE_ID_PROMOTERS_DICTIONARY.keys())[:2])


# In[59]:


# test: print last + second to last elements of the dictionary

print(list(GENE_ID_PROMOTERS_DICTIONARY.keys())[-2:])


# In[60]:


""" Export dictionary to csv """

with open('GENE_ID_1500BP_PROMOTER_DICTIONARY.csv', 'w') as f:
    for key in GENE_ID_PROMOTERS_DICTIONARY.keys():
        f.write("%s,%s\n"%(key,GENE_ID_PROMOTERS_DICTIONARY[key]))

