import pandas as pd
import os
import re 

list_jbst = os.listdir('data/siRNA_bench')


list_jbst = [x for x in list_jbst if 'jbst' in x]


list_jbst_pax = [x for x in list_jbst if 'PAX' in x]


df = pd.DataFrame()
for p in list_jbst_pax:
    tmp = pd.read_excel(os.path.join('data', 'siRNA_bench', p))
    tmp['source'] = re.sub('_siRNA.*', '', p)
    df = pd.concat([df, tmp])
    
    
df.to_excel('data/siRNA_bench/mutual_PAX3.xlsx', index=False, engine='openpyxl')


list_jbst_kit = [x for x in list_jbst if 'KIT' in x]


df = pd.DataFrame()
for p in list_jbst_kit:
    tmp = pd.read_excel(os.path.join('data', 'siRNA_bench', p))
    tmp['source'] = re.sub('_siRNA.*', '', p)
    df = pd.concat([df, tmp])
    
    
df.to_excel('data/siRNA_bench/mutual_KIT.xlsx', index=False, engine='openpyxl')



################################KIT####################################################


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import ast

df = pd.read_excel('data/siRNA_bench/mutual_KIT.xlsx')



optimal_df = df[(df['GC%'] >= 30) & (df['GC%'] <= 60)]
outside_df = df[(df['GC%'] < 30) | (df['GC%'] > 60)]

fig = plt.figure(figsize=(10, 6))

sns.boxplot(x='source', y='GC%', data=df, palette='Set2', fliersize=0)

sns.stripplot(x='source', y='GC%', data=optimal_df, color='green', 
              size=7, marker='o')

sns.stripplot(x='source', y='GC%', data=outside_df, color='red', 
              size=7, marker='X')

plt.axhline(30, color='red', linestyle='--', label='30%')
plt.axhline(60, color='blue', linestyle='--', label='60%')

plt.title('GC% - KIT', fontsize=13)
plt.xlabel('Tool')
plt.ylabel('GC (%)')
plt.tight_layout()
plt.show()

fig.savefig('KIT_GC.svg', dpi = 300)



fig = plt.figure(figsize=(10, 6))
sns.boxplot(x='source', y='score', data=df, palette='Set2')
plt.title('Score - KIT')
plt.xlabel('Tool')
plt.ylabel('Score')
plt.show()

fig.savefig('KIT_score.svg', dpi = 300)



##############################################################################

max_score_df = df.loc[df.groupby('source')['score'].idxmax()]

fig = plt.figure(figsize=(8, 6))
sns.scatterplot(data=max_score_df, x='GC%', y='score', hue='source', s=150, edgecolor='black')

plt.axvline(30, color='red', linestyle='--', label='Lower GC threshold')
plt.axvline(60, color='blue', linestyle='--', label='Upper GC threshold')

plt.title('Top detected siRNA [score ~ GC%] - KIT', fontsize=13)
plt.xlabel('GC (%)')
plt.ylabel('Score')

plt.legend(title='', loc='upper right', bbox_to_anchor=(0.95, 1))

plt.tight_layout()
plt.show()

fig.savefig('KIT_score~gc.svg', dpi = 300)


######################################################################################


import ast
import seaborn as sns
import matplotlib.pyplot as plt

def count_list_elements(x):
    if isinstance(x, str):
        try:
            parsed = ast.literal_eval(x)
            if isinstance(parsed, list):
                return len(parsed)
        except:
            return 0
    elif isinstance(x, list):
        return len(x)
    return 0

df['complemenatry_count'] = df['complemenatry_regions'].apply(count_list_elements)

df['complemenatry_pct_percent'] = df['complemenatry_pct'] * 100

fig = plt.figure(figsize=(10, 6))
sns.boxplot(x='source', y='complemenatry_count', data=df, palette='Set2')
plt.title('Self-complementary regions in siRNA - KIT')
plt.xlabel('Tool')
plt.ylabel('Count')
plt.tight_layout()
plt.show()

fig.savefig('KIT_complementary.svg', dpi = 300)




################################PAX3####################################################


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import ast

df = pd.read_excel('data/siRNA_bench/mutual_PAX3.xlsx')



optimal_df = df[(df['GC%'] >= 30) & (df['GC%'] <= 60)]
outside_df = df[(df['GC%'] < 30) | (df['GC%'] > 60)]

fig = plt.figure(figsize=(10, 6))

sns.boxplot(x='source', y='GC%', data=df, palette='Set2', fliersize=0)

sns.stripplot(x='source', y='GC%', data=optimal_df, color='green', 
              size=7, marker='o')

sns.stripplot(x='source', y='GC%', data=outside_df, color='red', 
              size=7, marker='X')

plt.axhline(30, color='red', linestyle='--', label='30%')
plt.axhline(60, color='blue', linestyle='--', label='60%')

plt.title('GC% - PAX3', fontsize=13)
plt.xlabel('Tool')
plt.ylabel('GC (%)')
plt.tight_layout()
plt.show()

fig.savefig('PAX3_GC.svg', dpi = 300)



fig = plt.figure(figsize=(10, 6))
sns.boxplot(x='source', y='score', data=df, palette='Set2')
plt.title('Score - PAX3')
plt.xlabel('Tool')
plt.ylabel('Score')
plt.show()

fig.savefig('PAX3_score.svg', dpi = 300)



##############################################################################

max_score_df = df.loc[df.groupby('source')['score'].idxmax()]

fig = plt.figure(figsize=(8, 6))
sns.scatterplot(data=max_score_df, x='GC%', y='score', hue='source', s=150, edgecolor='black')

plt.axvline(30, color='red', linestyle='--', label='Lower GC threshold')
plt.axvline(60, color='blue', linestyle='--', label='Upper GC threshold')

plt.title('Top detected siRNA [score ~ GC%] - PAX3', fontsize=13)
plt.xlabel('GC (%)')
plt.ylabel('Score')

plt.legend(title='', loc='upper right', bbox_to_anchor=(0.95, 1))

plt.tight_layout()
plt.show()

fig.savefig('PAX3_score~gc.svg', dpi = 300)


######################################################################################


import ast
import seaborn as sns
import matplotlib.pyplot as plt

def count_list_elements(x):
    if isinstance(x, str):
        try:
            parsed = ast.literal_eval(x)
            if isinstance(parsed, list):
                return len(parsed)
        except:
            return 0
    elif isinstance(x, list):
        return len(x)
    return 0

df['complemenatry_count'] = df['complemenatry_regions'].apply(count_list_elements)

df['complemenatry_pct_percent'] = df['complemenatry_pct'] * 100

fig = plt.figure(figsize=(10, 6))
sns.boxplot(x='source', y='complemenatry_count', data=df, palette='Set2')
plt.title('Self-complementary regions in siRNA - PAX3')
plt.xlabel('Tool')
plt.ylabel('Count')
plt.tight_layout()
plt.show()

fig.savefig('PAX3_complementary.svg', dpi = 300)


################################################################################


## codon optimization


import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

df = pd.read_excel("data/opt_results.xlsx", sheet_name="Sheet1")

genes = ["KIT", "PAX3"]

sns.set(style="whitegrid", context="talk")

for gene in genes:
    sub = df[df["gene"] == gene]


    fig = plt.figure(figsize=(10,6))
    sns.barplot(x="tool", y="GC", data=sub, palette="Set2", ci="sd", edgecolor="black")
    sns.stripplot(x="tool", y="GC", data=sub, color="black", size=6, jitter=True)
    plt.title(f"{gene} — GC%")
    plt.ylabel("GC [%]")
    plt.xlabel("Sequence/Tool")
    plt.tight_layout()
    plt.show()
    
    fig.savefig(f'{gene}_optimization_GC.svg', dpi = 300)
    
    
    fig = plt.figure(figsize=(10,6))
    sns.barplot(x="tool", y="freq", data=sub, palette="Set2", ci="sd", edgecolor="black")
    sns.stripplot(x="tool", y="freq", data=sub, color="black", size=6, jitter=True)
    plt.title(f"{gene} — Codon Frequency")
    plt.ylabel("Codon Frequency")
    plt.xlabel("Sequence/Tool")
    plt.tight_layout()
    plt.show()
    
    fig.savefig(f'{gene}_optimization_freq.svg', dpi = 300)
    
     
    fig = plt.figure(figsize=(10,6))
    sns.barplot(x="tool", y="changed_bases %", data=sub, palette="Set2", ci="sd", edgecolor="black")
    sns.stripplot(x="tool", y="changed_bases %", data=sub, color="black", size=6, jitter=True)
    plt.title(f"{gene} - Changed bases")
    plt.ylabel("Changed bases [%]")
    plt.xlabel("Sequence/Tool")
    plt.tight_layout()
    plt.show()
    
    fig.savefig(f'{gene}_optimization_bases.svg', dpi = 300)
    
    
    fig = plt.figure(figsize=(10,6))
    sns.barplot(x="tool", y="changed_codons %", data=sub, palette="Set2", ci="sd", edgecolor="black")
    sns.stripplot(x="tool", y="changed_codons %", data=sub, color="black", size=6, jitter=True)
    plt.title(f"{gene} — Changed codons")
    plt.ylabel("Changed codons [%]")
    plt.xlabel("Sequence/Tool")
    plt.tight_layout()
    plt.show()
    
    fig.savefig(f'{gene}_optimization_codons.svg', dpi = 300)


    fig = plt.figure(figsize=(10,6))
    sns.barplot(x="tool", y="MFE [kcal/mol]", data=sub, palette="Set2", ci="sd", edgecolor="black")
    sns.stripplot(x="tool", y="MFE [kcal/mol]", data=sub, color="black", size=6, jitter=True)
    plt.title(f"{gene} - Minimum Free Energy (MFE)")
    plt.ylabel("MFE [kcal/mol]")
    plt.xlabel("Sequence/Tool")
    plt.tight_layout()
    plt.show()
    
    fig.savefig(f'{gene}_optimization_MEF.svg', dpi = 300)


    # --- GC vs MFE (scatter) ---
    fig = plt.figure(figsize=(10,6))
    sns.scatterplot(
        data=sub, x="GC", y="MFE [kcal/mol]",
        hue="tool", s=150, palette="coolwarm", edgecolor="black"
    )
    
    plt.title(f"{gene} — GC% ~ MFE", fontsize=16)
    plt.xlabel("GC [%]", fontsize=14)
    plt.ylabel("MFE [kcal/mol]", fontsize=14)
    
    leg = plt.legend(title=None, fontsize=12)
    plt.setp(leg.get_texts(), fontweight="bold")
    
    plt.grid(False)
    
    plt.tight_layout()
    plt.show()
    
    fig.savefig(f'{gene}_optimization_MEF_GC.svg', dpi = 300)


native = df[df["tool"] == "Native"].set_index("gene")
opt = df[df["tool"] != "Native"].set_index("gene")

delta = opt.copy()
delta["ΔMFE"] = opt["MFE [kcal/mol]"] - native.loc[opt.index, "MFE [kcal/mol]"].values

for gene in genes:


    fig = plt.figure(figsize=(10,6))
    sns.barplot(x="tool", y="ΔMFE", data=delta[delta.index == gene], palette="rocket", edgecolor="black")
    plt.axhline(0, color="gray", linestyle="--")
    plt.title(f"{gene} - Change in MFE relative to the native sequence")
    plt.ylabel("ΔMFE [kcal/mol]")
    plt.tight_layout()
    plt.show()
    
    fig.savefig(f'{gene}_optimization_deltaMEF.svg', dpi = 300)






import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import re


cols = ["gene", "tool", "G_max[n]", "A_max[n]", "C_max[n]", "T_max[n]"]

genes = df["gene"].unique()

sns.set(style="white", context="talk")

for gene in genes:
    sub = df[df["gene"] == gene][cols].set_index("tool")

    data_heat = sub[["G_max[n]", "A_max[n]", "C_max[n]", "T_max[n]"]]
    data_heat.columns = [re.sub(r'_max\[n\]', '', x) for x in data_heat.columns]

    fig = plt.figure(figsize=(10,8))
    sns.heatmap(
        data_heat,
        annot=True, fmt=".0f",
        cmap="viridis",
        cbar_kws={'label': 'max_count'},
        linewidths=1, linecolor="white"
    )
    plt.title(f"{gene} - Max repeat tract of G/A/C/T", fontsize=16)
    plt.xlabel("")
    plt.ylabel("Sequence/Tool")
    plt.tight_layout()
    plt.show()
    
    fig.savefig(f'{gene}_optimization_heatmap.svg', dpi = 300)





