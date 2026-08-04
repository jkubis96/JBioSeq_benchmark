import json
import pandas as pd

with open("random_CDS.json", "r", encoding="utf-8") as f:
    CDS = json.load(f)


condb = pd.DataFrame()

for k in CDS.keys():
    df = pd.json_normalize(CDS[k])
    df['type'] = k
    
    condb = pd.concat([condb,df])
    
condb.to_excel('random_CDS.xlsx', index=False)




import pandas as pd

condb = pd.read_excel('random_CDS.xlsx')
from jbst import seq_tools as st

metadata = st.load_metadata(linkers = False, 
                                    loops = False, 
                                    regulators = False, 
                                    fluorescent_tag = False, 
                                    promoters = False, 
                                    polya = False, 
                                    marker = False, 
                                    utr5 = False, 
                                    utr3 = False) 



from tqdm import tqdm
import re
condb = condb.reset_index(drop = True)

condb['aa'] = None
condb['frequence'] = None
condb['mfe'] = None
condb['gc'] = None
condb['jbst_sequence'] = None
condb['jbst_mfe'] = None
condb['jbst_gc'] = None
condb['jbst_frequence'] = None
condb['jbst_G_max[n]'] = None
condb['jbst_A_max[n]'] = None
condb['jbst_C_max[n]'] = None
condb['jbst_T_max[n]'] = None
condb['jbst_codon_change'] = None
condb['jbst_nucleotide_change'] = None



for i in tqdm(condb.index):
    
    tmp = st.codon_optimization(condb.loc[i,'seq'], metadata, species = 'human')
    
    condb.loc[i,'frequence'] = tmp['frequence'][0]
    condb.loc[i,'aa'] = tmp['sequence_aa'][0]
    condb.loc[i,'jbst_sequence'] = tmp['sequence_na'][1]
    condb.loc[i,'jbst_mfe'] = tmp['MFE'][1]
    condb.loc[i,'jbst_gc'] = tmp['GC%'][1]
    condb.loc[i,'mfe'] = tmp['MFE'][0]
    condb.loc[i,'gc'] = tmp['GC%'][0]
    condb.loc[i,'jbst_frequence'] = tmp['frequence'][1]
    for r in list(range(1,20, 1)):
        for n in ['A', 'C', 'T', 'G']:
            if n*r in tmp['sequence_na'][1]:
                condb.loc[i,f'jbst_{n}_max[n]'] = r
    
    results_1 = st.compare_sequences(tmp['sequence_na'][0], 
                        'native',
                        tmp['sequence_na'][1], 
                        'jbst',
                        sep = 1)
    pct = float(re.sub(r".*Changed positions percent \[%\]: ([0-9.]+).*", r"\1", results_1, flags=re.S))
    pct = round(pct, 2)
    condb.loc[i,'jbst_nucleotide_change'] = pct

    
    
    results_2 = st.compare_sequences(tmp['sequence_na'][0], 
                        'native',
                        tmp['sequence_na'][1], 
                        'jbst',
                        sep = 3)
    
    pct2 = float(re.sub(r".*Changed positions percent \[%\]: ([0-9.]+).*", r"\1", results_2, flags=re.S))
    pct2 = round(pct2, 2)
    
    condb.loc[i,'jbst_codon_change'] = pct2

                






condb.to_excel('random_CDS_jbst_extended.xlsx', index=False)





# RNAI

import pandas as pd

condb = pd.read_excel('random_CDS_jbst_RNAi_top1.xlsx')
from jbst import seq_tools as st

metadata = st.load_metadata(linkers = False, 
                                    loops = False, 
                                    regulators = False, 
                                    fluorescent_tag = False, 
                                    promoters = False, 
                                    polya = False, 
                                    marker = False, 
                                    utr5 = False, 
                                    utr3 = False) 


from tqdm import tqdm

cols = [c for c in condb.columns if c not in ['target_seq', 'source']]

df_distinct = condb.drop_duplicates(subset=cols).reset_index(drop = True)

for s in tqdm(df_distinct.index):
    print(s)
    tmp_df = pd.DataFrame(df_distinct.iloc[s:s+1,:])
    sequence = st.clear_sequence(tmp_df['full_seq'][s])
    RNAi_data_21 =  st.FindRNAi(sequence, metadata, length = 21, n = 1000, max_repeat_len = 3, max_off = 1, species = 'human', output = None, database_name = "refseq_select_rna",  evalue = 1e-3, outfmt =  5, word_size = 7, max_hsps = 20, reward = 1, penalty = -3, gapopen = 5, gapextend = 2, dust = "no", extension = 'xml')    
    RNAi_data_21 = RNAi_data_21[(RNAi_data_21["GC%"] > 30) & (RNAi_data_21["GC%"] < 60)].reset_index(drop = True)
    RNAi_data_21 = RNAi_data_21.sort_values(
    by=["specificity", "repeated_motif_pct", "complemenatry_pct", "score"],
    ascending=[True, True, True, False],
    ).head(10).reset_index(drop = True)
    tmp_df['source'] = 'JBST_21'
    tmp_df['target_seq'] = RNAi_data_21['RNAi_sense'][0]
    condb = pd.concat([condb, tmp_df])
    RNAi_data_19 =  st.FindRNAi(sequence, metadata, length = 19, n = 1000, max_repeat_len = 3, max_off = 1, species = 'human', output = None, database_name = "refseq_select_rna",  evalue = 1e-3, outfmt =  5, word_size = 7, max_hsps = 20, reward = 1, penalty = -3, gapopen = 5, gapextend = 2, dust = "no", extension = 'xml')    
    RNAi_data_19 = RNAi_data_19[(RNAi_data_19["GC%"] > 30) & (RNAi_data_19["GC%"] < 60)].reset_index(drop = True)
    RNAi_data_19 = RNAi_data_19.sort_values(
    by=["specificity", "repeated_motif_pct", "complemenatry_pct", "score"],
    ascending=[True, True, True, False],
    ).head(10).reset_index(drop = True)
    tmp_df['source'] = 'JBST_19'
    tmp_df['target_seq'] = RNAi_data_19['RNAi_sense'][0]
    condb = pd.concat([condb, tmp_df])


condb = condb.reset_index(drop = True)

condb.to_excel('random_CDS_jbst_RNAi_top1_with_jbst.xlsx', index=False)





