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






seq = st.get_sequences_gene('KIT', species = 'human', max_results = 20)
    

sequence = load_sequence()
sequence = clear_sequence(sequence)




seq = st.get_sequences_gene('PAX3', species = 'human', max_results = 20)
    


sequence = st.load_sequence()
sequence = st.clear_sequence(sequence)



# JBioSeq

## KIT

sequence = st.load_sequence()

sequence = st.clear_sequence(sequence)

st.check_coding(sequence)

optimized = st.codon_otymization(sequence, metadata, species = 'human')



## PAX3

sequence = st.load_sequence()

sequence = st.clear_sequence(sequence)

st.check_coding(sequence)

optimized = st.codon_otymization(sequence, metadata, species = 'human')








# CodonTransformer

## KIT

sequence = st.load_sequence()

sequence = st.clear_sequence(sequence)

st.check_coding(sequence)

optimized = st.codon_otymization(sequence, metadata, species = 'human')



## PAX3

sequence = st.load_sequence()

sequence = st.clear_sequence(sequence)

st.check_coding(sequence)

optimized = st.codon_otymization(sequence, metadata, species = 'human')




# GenScript

## KIT

sequence = st.load_sequence()

sequence = st.clear_sequence(sequence)

st.check_coding(sequence)

optimized = st.codon_otymization(sequence, metadata, species = 'human')



## PAX3

sequence = st.load_sequence()

sequence = st.clear_sequence(sequence)

st.check_coding(sequence)

optimized = st.codon_otymization(sequence, metadata, species = 'human')




# VectorBuilder

## KIT

sequence = st.load_sequence()

sequence = st.clear_sequence(sequence)

st.check_coding(sequence)

optimized = st.codon_otymization(sequence, metadata, species = 'human')



## PAX3

sequence = st.load_sequence()

sequence = st.clear_sequence(sequence)

st.check_coding(sequence)

optimized = st.codon_otymization(sequence, metadata, species = 'human')







################################################################################


# Native

## KIT

sequence = st.load_sequence()

sequence = st.clear_sequence(sequence)

_, dot = st.predict_structure(sequence, 
                  anty_sequence = '',
                  height=None, 
                  width=None, 
                  dis_alpha = 0.15, 
                  seq_force = 27, 
                  pair_force = 3, 
                  show_plot = True)




## PAX3

sequence = st.load_sequence()

sequence = st.clear_sequence(sequence)


_, dot = st.predict_structure(sequence, 
                  anty_sequence = '',
                  height=None, 
                  width=None, 
                  dis_alpha = 0.15, 
                  seq_force = 27, 
                  pair_force = 3, 
                  show_plot = True)




# JBioSeq

## KIT

sequence = st.load_sequence()

sequence = st.clear_sequence(sequence)

_, dot = st.predict_structure(sequence, 
                  anty_sequence = '',
                  height=None, 
                  width=None, 
                  dis_alpha = 0.15, 
                  seq_force = 27, 
                  pair_force = 3, 
                  show_plot = True)



## PAX3

sequence = st.load_sequence()

sequence = st.clear_sequence(sequence)

_, dot = st.predict_structure(sequence, 
                  anty_sequence = '',
                  height=None, 
                  width=None, 
                  dis_alpha = 0.15, 
                  seq_force = 27, 
                  pair_force = 3, 
                  show_plot = True)



# CodonTransformer

## KIT

sequence = st.load_sequence()

sequence = st.clear_sequence(sequence)

_, dot = st.predict_structure(sequence, 
                  anty_sequence = '',
                  height=None, 
                  width=None, 
                  dis_alpha = 0.15, 
                  seq_force = 27, 
                  pair_force = 3, 
                  show_plot = True)



## PAX3

sequence = st.load_sequence()

sequence = st.clear_sequence(sequence)

_, dot = st.predict_structure(sequence, 
                  anty_sequence = '',
                  height=None, 
                  width=None, 
                  dis_alpha = 0.15, 
                  seq_force = 27, 
                  pair_force = 3, 
                  show_plot = True)




# GenScript

## KIT

sequence = st.load_sequence()

sequence = st.clear_sequence(sequence)

_, dot = st.predict_structure(sequence, 
                  anty_sequence = '',
                  height=None, 
                  width=None, 
                  dis_alpha = 0.15, 
                  seq_force = 27, 
                  pair_force = 3, 
                  show_plot = True)




## PAX3

sequence = st.load_sequence()

sequence = st.clear_sequence(sequence)

_, dot = st.predict_structure(sequence, 
                  anty_sequence = '',
                  height=None, 
                  width=None, 
                  dis_alpha = 0.15, 
                  seq_force = 27, 
                  pair_force = 3, 
                  show_plot = True)





# VectorBuilder

## KIT

sequence = st.load_sequence()

sequence = st.clear_sequence(sequence)

_, dot = st.predict_structure(sequence, 
                  anty_sequence = '',
                  height=None, 
                  width=None, 
                  dis_alpha = 0.15, 
                  seq_force = 27, 
                  pair_force = 3, 
                  show_plot = True)



## PAX3

sequence = st.load_sequence()

sequence = st.clear_sequence(sequence)

_, dot = st.predict_structure(sequence, 
                  anty_sequence = '',
                  height=None, 
                  width=None, 
                  dis_alpha = 0.15, 
                  seq_force = 27, 
                  pair_force = 3, 
                  show_plot = True)



################################################################################





# Native

## KIT

native = st.load_sequence()

native = st.clear_sequence(native)

sequence_name_1 = 'KIT_native'

sequence_name_2 = 'KIT_optimized'





# JBioSeq

## KIT

sequence = st.load_sequence()

sequence = st.clear_sequence(sequence)


results_1_2 = st.compare_sequences(native, 
                    sequence_name_1,
                    sequence, 
                    sequence_name_2,
                    sep = 1)

results_1_2 = st.compare_sequences(native, 
                    sequence_name_1,
                    sequence, 
                    sequence_name_2,
                    sep = 3)





# CodonTransformer

## KIT

sequence = st.load_sequence()

sequence = st.clear_sequence(sequence)



results_1_2 = st.compare_sequences(native, 
                    sequence_name_1,
                    sequence, 
                    sequence_name_2,
                    sep = 1)

results_1_2 = st.compare_sequences(native, 
                    sequence_name_1,
                    sequence, 
                    sequence_name_2,
                    sep = 3)







# GenScript

## KIT

sequence = st.load_sequence()

sequence = st.clear_sequence(sequence)



results_1_2 = st.compare_sequences(native, 
                    sequence_name_1,
                    sequence, 
                    sequence_name_2,
                    sep = 1)

results_1_2 = st.compare_sequences(native, 
                    sequence_name_1,
                    sequence, 
                    sequence_name_2,
                    sep = 3)





# VectorBuilder

## KIT

sequence = st.load_sequence()

sequence = st.clear_sequence(sequence)



results_1_2 = st.compare_sequences(native, 
                    sequence_name_1,
                    sequence, 
                    sequence_name_2,
                    sep = 1)

results_1_2 = st.compare_sequences(native, 
                    sequence_name_1,
                    sequence, 
                    sequence_name_2,
                    sep = 3)





# Native




## PAX3

native = st.load_sequence()

native = st.clear_sequence(native)





# JBioSeq



## PAX3

sequence = st.load_sequence()

sequence = st.clear_sequence(sequence)


results_1_2 = st.compare_sequences(native, 
                    sequence_name_1,
                    sequence, 
                    sequence_name_2,
                    sep = 1)

results_1_2 = st.compare_sequences(native, 
                    sequence_name_1,
                    sequence, 
                    sequence_name_2,
                    sep = 3)




# CodonTransformer


## PAX3

sequence = st.load_sequence()

sequence = st.clear_sequence(sequence)



results_1_2 = st.compare_sequences(native, 
                    sequence_name_1,
                    sequence, 
                    sequence_name_2,
                    sep = 1)

results_1_2 = st.compare_sequences(native, 
                    sequence_name_1,
                    sequence, 
                    sequence_name_2,
                    sep = 3)



# GenScript




## PAX3

sequence = st.load_sequence()

sequence = st.clear_sequence(sequence)



results_1_2 = st.compare_sequences(native, 
                    sequence_name_1,
                    sequence, 
                    sequence_name_2,
                    sep = 1)

results_1_2 = st.compare_sequences(native, 
                    sequence_name_1,
                    sequence, 
                    sequence_name_2,
                    sep = 3)


# VectorBuilder




## PAX3

sequence = st.load_sequence()

sequence = st.clear_sequence(sequence)


results_1_2 = st.compare_sequences(native, 
                    sequence_name_1,
                    sequence, 
                    sequence_name_2,
                    sep = 1)

results_1_2 = st.compare_sequences(native, 
                    sequence_name_1,
                    sequence, 
                    sequence_name_2,
                    sep = 3)




##################################################################################



import RNA
import pandas as pd


data = pd.read_excel('optimization.xlsx',    engine='openpyxl')
data['MEF'] = None

for i in data.index:
    print(i)

    seq = data.loc[i,'sequence']
    structure = data.loc[i,'dot']
    
    energy = RNA.eval_structure_simple(seq, structure)
    
    data.loc[i,'MEF'] = energy
    print(energy)
    



# stats



import pandas as pd
sequences = pd.read_excel('optimization.xlsx')

sequences['G_max[n]'] = None
sequences['A_max[n]'] = None
sequences['C_max[n]'] = None
sequences['T_max[n]'] = None

for i in sequences.index:
    print(i)
    
    for r in list(range(3,30, 1)):
        for n in ['A', 'C', 'T', 'G']:
            if n*r in sequences.loc[i, 'sequence']:
                sequences.loc[i,f'{n}_max[n]'] = r
                
            
sequences.to_excel('opt_results.xlsx', index=False, engine='openpyxl')

