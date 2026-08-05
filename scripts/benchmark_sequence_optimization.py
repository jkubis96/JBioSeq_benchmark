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
    )
    tmp_df['source'] = 'JBST_21'
    tmp_df['target_seq'] = RNAi_data_21['RNAi_sense'][0]
    condb = pd.concat([condb, tmp_df])
    RNAi_data_19 =  st.FindRNAi(sequence, metadata, length = 19, n = 1000, max_repeat_len = 3, max_off = 1, species = 'human', output = None, database_name = "refseq_select_rna",  evalue = 1e-3, outfmt =  5, word_size = 7, max_hsps = 20, reward = 1, penalty = -3, gapopen = 5, gapextend = 2, dust = "no", extension = 'xml')    
    RNAi_data_19 = RNAi_data_19[(RNAi_data_19["GC%"] > 30) & (RNAi_data_19["GC%"] < 60)].reset_index(drop = True)
    RNAi_data_19 = RNAi_data_19.sort_values(
    by=["specificity", "repeated_motif_pct", "complemenatry_pct", "score"],
    ascending=[True, True, True, False],
    )
    tmp_df['source'] = 'JBST_19'
    tmp_df['target_seq'] = RNAi_data_19['RNAi_sense'][0]
    condb = pd.concat([condb, tmp_df])


condb = condb.reset_index(drop = True)

condb.to_excel('random_CDS_jbst_RNAi_top1_with_jbst.xlsx', index=False)


###############################################################################

condb = pd.read_excel('random_CDS_jbst_RNAi_top1_with_jbst.xlsx')


from jbst.seq_tools import *

metadata = load_metadata() 


def rnai_scroing_base(sequence):
    sequence = reverse(sequence=complement(sequence=sequence))

    scoring = metadata["rnai"]

    score = 0

    for i in sequence[:3]:
        if i in ["G", "C"]:
            score = score + 1
        elif i in ["A", "T"]:
            score = score - 1

    for i in sequence[-3:]:
        if i in ["G", "C"]:
            score = score - 1
        elif i in ["A", "T"]:
            score = score + 1

    for j in scoring.index:
        if scoring["position"][j] != "last":
            if sequence[scoring["position"][j]] == scoring["element"][j]:
                if "+" in scoring["operation"][j]:
                    score = score + float(scoring["score"][j])
                elif "-" in scoring["operation"][j]:
                    score = score - float(scoring["score"][j])
        elif scoring["position"][j] == "last":
            if sequence[-1] == scoring["element"][j]:
                if "+" in scoring["operation"][j]:
                    score = score + float(scoring["score"][j])
                elif "-" in scoring["operation"][j]:
                    score = score - float(scoring["score"][j])

    return score, sequence



def find_self_complementarity(sequence, min_length=5):
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    self_complementary_regions = []
    
    max_range =  len(sequence)

    while True:
        
        min_range = min_length
    
        intervals = list(range(min_range,max_range + 1,  1))
        
        if min_range + min_length  <= max(intervals):     
            for i in intervals:
                if i + min_length <= max_range:
                    pre_seq = sequence[0:i]
                    post_seq = sequence[i:i + min_length][::-1]
                    post_seq_complement = ''.join([complement[x] for x in post_seq])

                    
                    if post_seq_complement in pre_seq:
                        self_complementary_regions.append((post_seq_complement,post_seq))
                    
                else:
                    break
                
            min_length += 1
        
        else:
            break
        
    # unification
    self_complementary_regions = [
        x for x in self_complementary_regions
        if all(
            x[0] not in y[0] and x[1] not in y[1]
            for y in self_complementary_regions
            if y != x
        )
    ]
    
    return self_complementary_regions



def repeat_scoring(seq, max_repeat_len):
    repeat_list = []
    i = 0
    while i < len(seq):
        repeat_char = seq[i]
        count = 1
        while i + 1 < len(seq) and seq[i + 1] == repeat_char:
            count += 1
            i += 1
        if count > max_repeat_len:
            repeat_list.append(repeat_char * count)
        i += 1

    full_len = sum(len(rep) for rep in repeat_list)
    pct = round(full_len / len(seq), 2) if len(seq) > 0 else 0.0

    return repeat_list, pct


    
# args

max_repeat_len = 3 
max_off: int = 1
species: str = "human"
database_name: str = "refseq_select_rna"
evalue=1e-2
outfmt=5
word_size: int = 7
max_hsps: int = 20
reward=1
penalty=-2
gapopen=5
gapextend=2
dust="no"



import jbst
source = str(jbst.__file__)
source = re.sub(r'\\__init__.py', '', source)


system = platform.system()

if system == "Windows":
    print("\nWindows operating system")
    blast_executable = os.path.join(
        source, "blast/windows/ncbi-blast-2.14.1+/bin/"
    )
    command = "blastn.exe"

elif system == "Linux":
    print("\nLinux operating system")
    blast_executable = os.path.join(
        source, "blast/linux/ncbi-blast-2.14.1+/bin/"
    )
    command = "./blastn"


output_file = os.getcwd()
output_file = output_file + "\\tmp_blast_out.xml"


    
    
    
# put the RNAi list


# Eurofinsgenomics_siRNA_KIT
 
predicted_rnai = condb['target_seq'] 
predicted_rnai = [rna_to_dna(x) for x in predicted_rnai]

predicted_rnai = [
    reverse(sequence=complement(sequence=x)) for x in predicted_rnai
]

names = ["RNAi"] * len(predicted_rnai)
fasta_string = ""
unique_names = [f"{name}_{i}" for i, name in enumerate(names, start=1)]

condb['id'] = unique_names

for name, seq in zip(unique_names, predicted_rnai):
    fasta_string += f">{re.sub(' ', '_', name)}\n{seq}\n"

random_prefix = random_name(length=30)

input_file = os.path.join(source, "tmp", random_prefix + "_tmp_rnai.fasta")

with open(input_file, "w") as fasta_file:
    fasta_file.write(fasta_string)

command_list = [
    command,
    "-query",
    input_file,
    "-db",
    database_name,
    "-out",
    output_file,
    "-evalue",
    str(evalue),
    "-outfmt",
    str(outfmt),
    "-word_size",
    str(word_size),
    "-max_hsps",
    str(max_hsps),
    "-reward",
    str(reward),
    "-penalty",
    str(penalty),
    "-gapopen",
    str(gapopen),
    "-gapextend",
    str(gapextend),
    "-dust",
    str(dust),
]

if system == "Windows":
    subprocess.run(command_list, cwd=blast_executable, shell=True)

elif system == "Linux":
    subprocess.run(command_list, cwd=blast_executable, shell=False)

try:
    os.remove(input_file)
    print(f"{input_file} successfully deleted.")
except OSError as e:
    print(f"Error: {input_file} - {e.strerror}")

tree = ET.parse(output_file)

try:
    os.remove(output_file)
    print(f"{output_file} successfully deleted.")
except OSError as e:
    print(f"Error: {output_file} - {e.strerror}")

root = tree.getroot()

# Create lists to store data
query_ids = []
subject_ids = []
e_values = []
bit_scores = []
alignment_lengths = []
query_sequences = []
subject_sequences = []

# Iterate through the XML tree and extract relevant data
for iteration in root.findall(".//Iteration"):
    query_id = iteration.find(".//Iteration_query-def").text

    if len(iteration.findall(".//Hit")) > 0:
        for hit in iteration.findall(".//Hit"):
            subject_id = hit.find(".//Hit_def").text
            e_value = hit.find(".//Hsp_evalue").text
            bit_score = hit.find(".//Hsp_bit-score").text
            alignment_length = hit.find(".//Hsp_align-len").text
            query_sequence = hit.find(".//Hsp_qseq").text
            subject_sequence = hit.find(".//Hsp_hseq").text

            query_ids.append(query_id)
            subject_ids.append(subject_id)
            e_values.append(float(e_value))
            bit_scores.append(float(bit_score))
            alignment_lengths.append(int(alignment_length))
            query_sequences.append(query_sequence)
            subject_sequences.append(subject_sequence)
    else:
        query_ids.append(query_id)
        subject_ids.append(None)
        e_values.append(None)
        bit_scores.append(None)
        alignment_lengths.append(0)
        query_sequences.append(None)
        subject_sequences.append(None)

# Create a DataFrame
data = {
    "target": subject_ids,
    "e-value": e_values,
    "bit_score": bit_scores,
    "alignment_length": alignment_lengths,
    "target_seq": subject_sequences,
    "RNAi_name": query_ids,
}

df = pd.DataFrame(data)

name_mapping = dict(zip(unique_names, predicted_rnai))
df["RNAi_seq"] = df["RNAi_name"].map(name_mapping)

df["target_gene_name"] = [
    (
        re.sub(r"\).*", "", re.sub(r".*\(", "", x)).upper()
        if x is not None
        else None
    )
    for x in df["target"]
]

df["species"] = [
    (
        " ".join(re.sub(r"^PREDICTED: ", "", x).split()[:2])
        if x is not None
        else None
    )
    for x in df["target"]
]

try:
    species_map = {
        "human": ["Homo sapiens"],
        "mouse": ["Mus musculus"],
        "rat": ["Rattus norvegicus"],
        "both": ["Mus musculus", "Homo sapiens"],
        "both2": ["Rattus norvegicus", "Homo sapiens"],
        "multi": ["Mus musculus", "Rattus norvegicus", "Homo sapiens"],
    }

    species_lower = species.lower()
    if species_lower in species_map:
        allowed_species = species_map[species_lower]
        df = df[df["species"].isin(allowed_species) | df["species"].isna()]

except:
    None

df = (
    df.groupby(["RNAi_name", "RNAi_seq"])[
        [
            "target",
            "e-value",
            "bit_score",
            "alignment_length",
            "target_seq",
            "target_gene_name",
            "species",
        ]
    ]
    .agg(list)
    .reset_index()
)

df["specificity"] = None
df["complemenatry_regions"] = None
df["complemenatry_pct"] = None
df["RNAi_sense"] = None
df["repeated_motif"] = None
df["repeated_motif_pct"] = None
df["score"] = None
df["GC%"] = None

for i in df.index:

    if None in df["target"][i]:
        df["target"][i] = [y for y in df["target"][i] if y is not None]
        df["target_seq"][i] = [
            y for y in df["target_seq"][i] if y is not None
        ]
        df["target_gene_name"][i] = [
            y for y in df["target_gene_name"][i] if y is not None
        ]
        df["species"][i] = [y for y in df["species"][i] if y is not None]
        df["e-value"][i] = [y for y in df["e-value"][i] if y == y]
        df["bit_score"][i] = [y for y in df["bit_score"][i] if y == y]

    df["specificity"][i] = len(
        set([x.upper() for x in df["target_gene_name"][i]])
    )
    df["complemenatry_regions"][i] = list(
        set(find_self_complementarity(df["RNAi_seq"][i], min_length=5))
    )
    amount = 0
    for l in df["complemenatry_regions"][i]:
        amount = amount + len(l)

    try:
        df["complemenatry_pct"][i] = amount / len(df["RNAi_seq"][i])
    except:
        df["complemenatry_pct"][i] = 0

    df["RNAi_sense"][i] = rnai_scroing_base(df["RNAi_seq"][i])[1]
    df["score"][i] = rnai_scroing_base(df["RNAi_seq"][i])[0]

    df["repeated_motif"][i] = repeat_scoring(
        df["RNAi_seq"][i], max_repeat_len
    )[0]
    df["repeated_motif_pct"][i] = repeat_scoring(
        df["RNAi_seq"][i], max_repeat_len
    )[1]
    df["GC%"][i] = round(
        (df["RNAi_seq"][i].count("C")
        + df["RNAi_seq"][i].count("G")) / len(df["RNAi_seq"][i]) * 100,
        2,
    )


df = df.sort_values(
    by=["specificity", "repeated_motif_pct", "complemenatry_pct", "score"],
    ascending=[True, True, True, False],
)

df = df.reset_index(drop=True)

condb = condb.merge(
    df,
    left_on="id",
    right_on = "RNAi_name",
    how="left"
)

condb.to_excel('random_CDS_jbst_RNAi_top1_with_jbst_with_statistics.xlsx', index=False, engine='openpyxl')




