from jbst import seq_tools as st



# Test1

metadata = st.load_metadata() 

    


seq = st.get_sequences_gene('KIT', species = 'human', max_results = 20)
    

sequence = st.load_sequence()
sequence = st.clear_sequence(sequence)




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

