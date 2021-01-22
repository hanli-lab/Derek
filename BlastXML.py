from Bio.Blast import NCBIXML as xml
import pandas as pd
import collections as c
import matplotlib.pyplot as plt
import numpy as np
from Bio import Entrez 

def blast_dataframe(blast_records): 
    diction = {}
    Source_Database = []
    Accessions = []
    Description = []
    length = []
    hsp_qseq = []
    hsp_hseq = []
    hsp_mseq = []
    hsp_score = []
    hsp_query_start = []
    hsp_align_length = []
    hsp_hit_start = []
    hsp_expect = []
    hsp_bits = []
    hsp_identities =[]
    hsp_positives = []
    hsp_query_end = []
    hsp_sbjct_end = []
    hsp_gaps = []
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            rawtitle = alignment.title
            parts = rawtitle.strip().split('|')
            Source_Database.append( parts[0] )
            Accessions.append( parts[1] )
            Description.append( parts[2])
            length.append( alignment.length)
            for h,hsp in enumerate(alignment.hsps):
                #print(hsp)
                #print(h)
                if h >= 1:  #Very important this removes worse alignments from the same gene, multiple alignments per gene will cause there to be too many alignments to pair with genes (array error)
                    break
                hsp_qseq.append(hsp.query)
                hsp_hseq.append(hsp.sbjct)
                hsp_mseq.append(hsp.match)
                hsp_score.append(hsp.score)
                hsp_query_start.append(hsp.query_start)
                hsp_align_length.append(hsp.align_length)
                hsp_hit_start.append(hsp.sbjct_start)
                hsp_expect.append(hsp.expect)
                hsp_bits.append(hsp.bits)
                hsp_identities.append(hsp.identities)
                hsp_positives.append(hsp.positives)
                hsp_query_end.append(hsp.query_end)
                hsp_sbjct_end.append(hsp.sbjct_end)
                hsp_gaps.append(hsp.gaps)
    diction =    {'Source_Database' : Source_Database,
                  'Accessions' : Accessions,
                  'Description' : Description,
                  'length': length,
                  'hsp_qseq' : hsp_qseq, 
                  'hsp_hseq' : hsp_hseq,
                  'hsp_mseq' : hsp_mseq,
                  'hsp_score' : hsp_score,
                  'hsp_query_start' : hsp_query_start,
                  'hsp_align_length' :hsp_align_length,
                  'hsp_hit_start' : hsp_hit_start,
                  'hsp_expect': hsp_expect,
                  'hsp_bits' : hsp_bits,
                  'hsp_identities' : hsp_identities,
                  'hsp_positives' : hsp_positives,
                  'hsp_query_end' : hsp_query_end,
                  'hsp_sbjct_end' : hsp_sbjct_end,
                  'hsp_gaps' : hsp_gaps
                 } 
    keys = [key for key in diction]
    df = pd.DataFrame(diction,columns = keys) 
    return df.set_index('Accessions')

def add_to_df_and_plot(orig_dataframe, dictionary_of_matches_,collabel,dom_q):
    lis = list(dictionary_of_matches_.keys())
    df_ = orig_dataframe.loc[lis,:]
    df_[collabel] = pd.Series(dictionary_of_matches_)
    df_[collabel+'_query'] = pd.Series(dom_q)
    Hit_res_ = list(df_[collabel])
    Hit_res_count = c.Counter(Hit_res_)
    Sum = sum(Hit_res_count.values())
    percentfreqs = {}
    for key in Hit_res_count:
        percentfreqs[key] = (Hit_res_count[key] / Sum)
    Hit_res_count.values
    print(Hit_res_count)
    y_pos = np.arange(len(percentfreqs))
    fig = plt.figure()
    plt.bar(y_pos,percentfreqs.values(), align='center', alpha=1, color = ['#A5A5A5','#FFC000','#5B9BD5','#9751CB','#ED7D31','#70AD47'])
    plt.xticks(y_pos, percentfreqs.keys())
    plt.ylabel('Freq')
    plt.title('Frequency of sequences with {}'.format(collabel))
    plt.show()
    return df_

def intersector(list_of_dictionary):
    intersect = set(list_of_dictionary[0].keys())
    for _list in list_of_dictionary:
        intersect = intersect & set(_list)
    return intersect

def get_position(df, posi,Expected_WT_Residue):
    print('the length of the blast dataframe is: ' + str(len(df)))
    output= {}
    output2 = {}
    error = {}
    gaploc= {}
    hseq = df['hsp_hseq']
    qseq = df['hsp_qseq']
    qi = df['hsp_query_start']
    for number,gene in enumerate(hseq.index):
        gaploc[gene] = [n for n in range(len(qseq[gene])) if qseq[gene].find('-', n) == n]
        shift = 0
        pos = int(posi)
        h = hseq[gene]
        q = qseq[gene]
        qindex = int(qi[gene])
        for gap in gaploc[gene]:
            if ((pos - qindex) >= (gap)):
                shift +=1   
        q_pos = pos - qindex + shift
        if  len(q) > q_pos and q_pos > 0:
            if Expected_WT_Residue == q[q_pos]:
                output[gene] = h[q_pos]
                output2[gene] = q[q_pos]
            else:
                print('Gap Error: '+gene[0:100]+'\nwas not included in dataframe because the position requested was not at expected index\n Please check this sequence manually\n')
                error[gene] = q[q_pos]
        elif (q_pos < 0):
            print('Residue {} outside alignment for gene: \n{}'.format(pos,gene))
            
    return output, output2, error

def append_percent_id(dataframe):
    dataframe.loc[:,'Percent Identity'] = 100*dataframe['hsp_identities'] / (dataframe['hsp_query_end']-dataframe['hsp_query_start'])
    return dataframe

def df_append_accesions(dataframe):
    genes = list(dataframe.index)
    accessions = {}
    for gene in genes:
        genesplit = gene.split('|')
        accessions[gene] = genesplit[1] #the 0 is the database code (ref,gb,emd), 1 is the accession
    dataframe.loc[:,'Accessions'] = pd.Series(accessions) 
    return dataframe   

def writeAccessionsFasta(dataframe,filename): #requires import of Entrez from Bio module
    Entrez.email = 'daspacio@uci.edu' #You will get an email if you excessively use the Entrez webserver for NCBI database lookups
    handle = Entrez.efetch(db="protein", id= dataframe.index, rettype="fasta", retmode="text")
    fasta = open(filename, "w")
    print('writing accessions to '+ filename + '...')
    for line in handle:
        fasta.write(line)
    fasta.close()
    print('file write completed')
    return