import random
from Bio import Entrez, SeqIO
from Bio.SeqUtils import gc_fraction

# Pamiętaj podać swój e-mail (NCBI tego wymaga)
Entrez.email = "jseq.info@gmail.com"

def get_grouped_cds():
    print("1. Searching NCBI database for mRNA transcripts...")
    
    # Query for Human RefSeq mRNA transcripts
    query = "Homo sapiens[Organism] AND RefSeq[Filter] AND biomol_mrna[PROP]"
    search_handle = Entrez.esearch(db="nucleotide", term=query, retmax=1000)
    search_results = Entrez.read(search_handle)
    search_handle.close()
    
    id_list = search_results["IdList"]
    print(f"Retrieved {len(id_list)} record IDs. Scanning CDS regions...\n")
    
    # Initialize 6 specific categories (5 genes per group = 30 total)
    groups = {
        "Long | High GC": [],
        "Long | Low GC": [],
        "Medium | High GC": [],
        "Medium | Low GC": [],
        "Short | High GC": [],
        "Short | Low GC": []
    }
    
    # Zbiór śledzący użyte nazwy genów, aby uniknąć duplikatów
    seen_gene_names = set()
    
    batch_size = 50
    for i in range(0, len(id_list), batch_size):
        if all(len(genes) >= 5 for genes in groups.values()):
            print("Successfully collected enough genes for all 6 categories!")
            break

        batch_ids = id_list[i:i+batch_size]
        fetch_handle = Entrez.efetch(db="nucleotide", id=batch_ids, rettype="gb", retmode="text")
        records = SeqIO.parse(fetch_handle, "genbank")
        
        for record in records:
            # Próba pobrania nazwy genu z adnotacji rekordów
            gene_name = None
            if "gene" in record.annotations:
                gene_name = record.annotations["gene"]
            
            for feature in record.features:
                if feature.type == "CDS":
                    # Jeśli nie znaleziono w annotations, szukamy w kwalifikatorach cechy CDS
                    if not gene_name and "gene" in feature.qualifiers:
                        gene_name = feature.qualifiers["gene"][0]
                    
                    # Jeśli nadal brak nazwy genu lub nazwa już wystąpiła, pomijamy ten rekord
                    if not gene_name or gene_name in seen_gene_names:
                        break
                        
                    try:
                        cds_seq = feature.extract(record.seq).upper()
                        
                        # Filter out invalid or incomplete sequences
                        if set(cds_seq) - set("ACGT") or len(cds_seq) < 100:
                            break
                            
                        length = len(cds_seq)
                        gc = round(gc_fraction(cds_seq) * 100, 2)
                        
                        # Classification logic based on defined thresholds
                        key = None
                        if length > 2000:  # Long
                            if gc > 55.0: key = "Long | High GC"
                            elif gc < 42.0: key = "Long | Low GC"
                        elif 800 <= length <= 2000:  # Medium
                            if gc > 55.0: key = "Medium | High GC"
                            elif gc < 42.0: key = "Medium | Low GC"
                        else:  # Short (< 800 bp)
                            if gc > 55.0: key = "Short | High GC"
                            elif gc < 42.0: key = "Short | Low GC"
                        
                        if key and len(groups[key]) < 10:  
                            # Dodajemy nazwę do zbioru unikanych genów
                            seen_gene_names.add(gene_name)
                            groups[key].append({
                                "id": record.id,
                                "gene_name": gene_name,
                                "length": length,
                                "gc": gc,
                                "seq": str(cds_seq)
                            })
                            break
                    except Exception:
                        continue
        fetch_handle.close()

    # Randomly select exactly 5 genes from each category
    final_30 = {}
    print("\n2. Sampling Results (5 genes per group):\n")
    
    for group_name, genes in groups.items():
        count_to_sample = min(len(genes), 5)
        sampled = random.sample(genes, count_to_sample)
        final_30[group_name] = sampled
        
        print(f"=== GROUP: {group_name} (Found: {len(genes)}) ===")
        print(f"{'Gene Name':<12} | {'Accession ID':<15} | {'CDS Length (bp)':<18} | {'GC Content (%)':<15}")
        print("-" * 70)
        for g in sampled:
            print(f"{g['gene_name']:<12} | {g['id']:<15} | {g['length']:<18} | {g['gc']:<15}%")
        print()

    return final_30

# Run script
results = get_grouped_cds()

import json
with open('random_CDS.json', "w", encoding="utf-8") as f:
        json.dump(results, f, indent=4)
