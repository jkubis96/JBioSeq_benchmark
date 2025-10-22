from jbst import seq_tools as st
import pandas as pd



metadata = st.load_metadata() 

    
seq = st.get_sequences_gene('KIT', species = 'human', max_results = 100)

data_dict = st.get_sequences_accesion(accesion_list= ['NM_001385286'])

fasta_string = st.generate_fasta_string(data_dict)

alignment_file = st.MuscleMultipleSequenceAlignment(fasta_string, output = None, gapopen = 10, gapextend = 0.5)

decoded_alignment_file = st.decode_alignments(alignment_file)

consensuse = st.ExtractConsensuse(alignment_file, refseq_sequences = data_dict)

res_mut = pd.DataFrame()

for i, _ in enumerate(consensuse['side']):
    print(i)
    
    res_rnai =  st.FindRNAi(consensuse['sequence'][i], 
                            metadata, 
                            length = 21, 
                            n = 100, 
                            max_repeat_len = 3, 
                            max_off = 1, 
                            species = 'human', 
                            output = None, 
                            database_name = "refseq_select_rna",  
                            evalue = 1e-2, 
                            outfmt =  5, 
                            word_size = 7, 
                            max_hsps = 20, 
                            reward = 1, 
                            penalty = -2, 
                            gapopen = 5, 
                            gapextend = 2, 
                            dust = "no", 
                            extension = 'xml')    
    
    res_rnai['side'] = consensuse['side'][i]
    res_rnai = res_rnai[(res_rnai["GC%"] > 30) & (res_rnai["GC%"] < 60)]

 
    res_mut = pd.concat([res_mut, res_rnai])
    
      
res_mut = res_mut.sort_values(
    by=["specificity", "repeated_motif_pct", "complemenatry_pct", "score"],
    ascending=[True, True, True, False],
).head(10)
    



res_mut.to_excel('data/siRNA_bench/JBioSeq_siRNA_KIT.xlsx', index=False, engine='openpyxl')





seq = st.get_sequences_gene('PAX3', species = 'human', max_results = 100)

data_dict = st.get_sequences_accesion(accesion_list= ['NM_181461.4'])

fasta_string = st.generate_fasta_string(data_dict)

alignment_file = st.MuscleMultipleSequenceAlignment(fasta_string, output = None, gapopen = 10, gapextend = 0.5)

decoded_alignment_file = st.decode_alignments(alignment_file)

consensuse = st.ExtractConsensuse(alignment_file, refseq_sequences = data_dict)

res_mut = pd.DataFrame()

for i, _ in enumerate(consensuse['side']):
    print(i)
    
    res_rnai =  st.FindRNAi(consensuse['sequence'][i], 
                            metadata, 
                            length = 21, 
                            n = 100, 
                            max_repeat_len = 3, 
                            max_off = 1, 
                            species = 'human', 
                            output = None, 
                            database_name = "refseq_select_rna",  
                            evalue = 1e-2, 
                            outfmt =  5, 
                            word_size = 7, 
                            max_hsps = 20, 
                            reward = 1, 
                            penalty = -2, 
                            gapopen = 5, 
                            gapextend = 2, 
                            dust = "no", 
                            extension = 'xml')    
    
    res_rnai['side'] = consensuse['side'][i]
    res_rnai = res_rnai[(res_rnai["GC%"] > 30) & (res_rnai["GC%"] < 60)]

 
    res_mut = pd.concat([res_mut, res_rnai])
    
      
res_mut = res_mut.sort_values(
    by=["specificity", "repeated_motif_pct", "complemenatry_pct", "score"],
    ascending=[True, True, True, False],
).head(10)
    



res_mut.to_excel('data/siRNA_bench/JBioSeq_siRNA_PAX3.xlsx', index=False, engine='openpyxl')



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

def find_self_complementarity(sequence, min_length=3):
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    self_complementary_regions = []

    for i in range(len(sequence) - min_length + 1):
        for j in range(i + min_length, len(sequence) + 1):
            subsequence = sequence[i:j]
            reverse_complement = "".join(
                complement[base] for base in subsequence[::-1]
            )

            if subsequence == reverse_complement:
                self_complementary_regions.append(subsequence)

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

max_repeat_len = 10 # to pass all results
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

file_1 = pd.read_excel('data/siRNA_bench/Eurofinsgenomics_siRNA_KIT.xls')   
 
predicted_rnai = file_1['siRNA Sequence (Sense)'] 
predicted_rnai = [rna_to_dna(x) for x in predicted_rnai]

predicted_rnai = [
    reverse(sequence=complement(sequence=x)) for x in predicted_rnai
]
fasta_string = ""
names = ["RNAi"] * len(predicted_rnai)
unique_names = [f"{name}_{i}" for i, name in enumerate(names, start=1)]
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
        set(find_self_complementarity(df["RNAi_seq"][i], min_length=3))
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
        df["RNAi_seq"][i].count("C")
        + df["RNAi_seq"][i].count("G") / len(df["RNAi_seq"][i]) * 100,
        2,
    )

df = df.sort_values(
    by=["specificity", "repeated_motif_pct", "complemenatry_pct", "score"],
    ascending=[True, True, True, False],
)

df = df.reset_index(drop=True)


df.to_excel('data/siRNA_bench/Eurofinsgenomics_siRNA_KIT_jbst.xlsx', index=False, engine='openpyxl')



#################################################################################


# Eurofinsgenomics_siRNA_PAX3

file_1 = pd.read_excel('data/siRNA_bench/Eurofinsgenomics_siRNA_PAX3.xls')   
 
predicted_rnai = file_1['siRNA Sequence (Sense)'] 
predicted_rnai = [rna_to_dna(x) for x in predicted_rnai]

predicted_rnai = [
    reverse(sequence=complement(sequence=x)) for x in predicted_rnai
]
fasta_string = ""
names = ["RNAi"] * len(predicted_rnai)
unique_names = [f"{name}_{i}" for i, name in enumerate(names, start=1)]
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
        set(find_self_complementarity(df["RNAi_seq"][i], min_length=3))
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
        df["RNAi_seq"][i].count("C")
        + df["RNAi_seq"][i].count("G") / len(df["RNAi_seq"][i]) * 100,
        2,
    )

df = df.sort_values(
    by=["specificity", "repeated_motif_pct", "complemenatry_pct", "score"],
    ascending=[True, True, True, False],
)

df = df.reset_index(drop=True)


df.to_excel('data/siRNA_bench/Eurofinsgenomics_siRNA_PAX3_jbst.xlsx', index=False, engine='openpyxl')




#################################################################################


# GenScript_siRNA_KIT

file_1 = pd.read_excel('data/siRNA_bench/GenScript_siRNA_KIT.xls')   
 
predicted_rnai = file_1['Sequence'] 
# predicted_rnai = [rna_to_dna(x) for x in predicted_rnai]

predicted_rnai = [
    reverse(sequence=complement(sequence=x)) for x in predicted_rnai
]
fasta_string = ""
names = ["RNAi"] * len(predicted_rnai)
unique_names = [f"{name}_{i}" for i, name in enumerate(names, start=1)]
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
        set(find_self_complementarity(df["RNAi_seq"][i], min_length=3))
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
        df["RNAi_seq"][i].count("C")
        + df["RNAi_seq"][i].count("G") / len(df["RNAi_seq"][i]) * 100,
        2,
    )

df = df.sort_values(
    by=["specificity", "repeated_motif_pct", "complemenatry_pct", "score"],
    ascending=[True, True, True, False],
)

df = df.reset_index(drop=True)


df.to_excel('data/siRNA_bench/GenScript_siRNA_KIT_jbst.xlsx', index=False, engine='openpyxl')





#################################################################################


# GenScript_siRNA_PAX3

file_1 = pd.read_excel('data/siRNA_bench/GenScript_siRNA_PAX3.xls')   
 
predicted_rnai = file_1['Sequence'] 
# predicted_rnai = [rna_to_dna(x) for x in predicted_rnai]

predicted_rnai = [
    reverse(sequence=complement(sequence=x)) for x in predicted_rnai
]
fasta_string = ""
names = ["RNAi"] * len(predicted_rnai)
unique_names = [f"{name}_{i}" for i, name in enumerate(names, start=1)]
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
        set(find_self_complementarity(df["RNAi_seq"][i], min_length=3))
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
        df["RNAi_seq"][i].count("C")
        + df["RNAi_seq"][i].count("G") / len(df["RNAi_seq"][i]) * 100,
        2,
    )

df = df.sort_values(
    by=["specificity", "repeated_motif_pct", "complemenatry_pct", "score"],
    ascending=[True, True, True, False],
)

df = df.reset_index(drop=True)


df.to_excel('data/siRNA_bench/GenScript_siRNA_PAX3_jbst.xlsx', index=False, engine='openpyxl')





#################################################################################


# Invivogen_siRNA_KIT

file_1 = pd.read_excel('data/siRNA_bench/Invivogen_siRNA_KIT.xls')   
 
predicted_rnai = file_1['Sequence'] 
# predicted_rnai = [rna_to_dna(x) for x in predicted_rnai]

predicted_rnai = [
    reverse(sequence=complement(sequence=x)) for x in predicted_rnai
]
fasta_string = ""
names = ["RNAi"] * len(predicted_rnai)
unique_names = [f"{name}_{i}" for i, name in enumerate(names, start=1)]
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
        set(find_self_complementarity(df["RNAi_seq"][i], min_length=3))
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
        df["RNAi_seq"][i].count("C")
        + df["RNAi_seq"][i].count("G") / len(df["RNAi_seq"][i]) * 100,
        2,
    )

df = df.sort_values(
    by=["specificity", "repeated_motif_pct", "complemenatry_pct", "score"],
    ascending=[True, True, True, False],
)

df = df.reset_index(drop=True)


df.to_excel('data/siRNA_bench/Invivogen_siRNA_KIT_jbst.xlsx', index=False, engine='openpyxl')




#################################################################################


# Invivogen_siRNA_PAX3

file_1 = pd.read_excel('data/siRNA_bench/Invivogen_siRNA_PAX3.xls')   
 
predicted_rnai = file_1['Sequence'] 
# predicted_rnai = [rna_to_dna(x) for x in predicted_rnai]

predicted_rnai = [
    reverse(sequence=complement(sequence=x)) for x in predicted_rnai
]
fasta_string = ""
names = ["RNAi"] * len(predicted_rnai)
unique_names = [f"{name}_{i}" for i, name in enumerate(names, start=1)]
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
        set(find_self_complementarity(df["RNAi_seq"][i], min_length=3))
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
        df["RNAi_seq"][i].count("C")
        + df["RNAi_seq"][i].count("G") / len(df["RNAi_seq"][i]) * 100,
        2,
    )

df = df.sort_values(
    by=["specificity", "repeated_motif_pct", "complemenatry_pct", "score"],
    ascending=[True, True, True, False],
)

df = df.reset_index(drop=True)


df.to_excel('data/siRNA_bench/Invivogen_siRNA_PAX3_jbst.xlsx', index=False, engine='openpyxl')





#################################################################################


# VectorBuilder_siRNA_KIT

file_1 = pd.read_excel('data/siRNA_bench/VectorBuilder_siRNA_KIT.xls')   
 
predicted_rnai = file_1['Target Sequence'] 
# predicted_rnai = [rna_to_dna(x) for x in predicted_rnai]

predicted_rnai = [
    reverse(sequence=complement(sequence=x)) for x in predicted_rnai
]
fasta_string = ""
names = ["RNAi"] * len(predicted_rnai)
unique_names = [f"{name}_{i}" for i, name in enumerate(names, start=1)]
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
        set(find_self_complementarity(df["RNAi_seq"][i], min_length=3))
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
        df["RNAi_seq"][i].count("C")
        + df["RNAi_seq"][i].count("G") / len(df["RNAi_seq"][i]) * 100,
        2,
    )

df = df.sort_values(
    by=["specificity", "repeated_motif_pct", "complemenatry_pct", "score"],
    ascending=[True, True, True, False],
)

df = df.reset_index(drop=True)


df.to_excel('data/siRNA_bench/VectorBuilder_siRNA_KIT_jbst.xlsx', index=False, engine='openpyxl')




#################################################################################


# VectorBuilder_siRNA_PAX3

file_1 = pd.read_excel('data/siRNA_bench/VectorBuilder_siRNA_PAX3.xls')   
 
predicted_rnai = file_1['Target Sequence'] 
# predicted_rnai = [rna_to_dna(x) for x in predicted_rnai]

predicted_rnai = [
    reverse(sequence=complement(sequence=x)) for x in predicted_rnai
]
fasta_string = ""
names = ["RNAi"] * len(predicted_rnai)
unique_names = [f"{name}_{i}" for i, name in enumerate(names, start=1)]
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
        set(find_self_complementarity(df["RNAi_seq"][i], min_length=3))
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
        df["RNAi_seq"][i].count("C")
        + df["RNAi_seq"][i].count("G") / len(df["RNAi_seq"][i]) * 100,
        2,
    )

df = df.sort_values(
    by=["specificity", "repeated_motif_pct", "complemenatry_pct", "score"],
    ascending=[True, True, True, False],
)

df = df.reset_index(drop=True)


df.to_excel('data/siRNA_bench/VectorBuilder_siRNA_PAX3_jbst.xlsx', index=False, engine='openpyxl')



