import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from scipy.stats import ttest_ind, pearsonr
import re

# --------------------------
# DATASET 1
# --------------------------
data1 = {
    'Catalyst': ['WT', 'E7A', 'A11L', 'L18A', 'N19A', 'K22A', 'N88A', 
                 'M89A', 'A92E', 'F93A', 'S95A', 'S97A', 'D100A'],
    'Rep1': [38, 38, 38, 28, 35, 49, 44, 34, 36, 33, 45, 50, 15],
    'Rep2': [37, 40, 40, 33, 36, 51, 39, 32, 37, 34, 38, 49, 12],
    'Rep3': [47, 40, 43, 40, 48, 42, 26, 40, 35, 26, 34, 43, 21],
    'Rep4': [46, 43, 43, 42, 47, 43, 37, 42, 34, 26, 36, 44, 41]
}
df1 = pd.DataFrame(data1)
df1['Mean'] = df1[['Rep1','Rep2','Rep3','Rep4']].mean(axis=1)
df1['SEM'] = df1[['Rep1','Rep2','Rep3','Rep4']].sem(axis=1)

# Sort: WT first, then by mean
wt1 = df1[df1['Catalyst']=='WT']
mutants1 = df1[df1['Catalyst']!='WT'].sort_values('Mean').reset_index(drop=True)
df1 = pd.concat([wt1, mutants1]).reset_index(drop=True)
x1 = np.arange(len(df1))

# Compute p-values against WT
wt_values1 = wt1[['Rep1','Rep2','Rep3','Rep4']].values.flatten().astype(float)
p_values1 = []
for i, row in df1.iterrows():
    if row['Catalyst'] == 'WT':
        p_values1.append(np.nan)
    else:
        mutant_values = row[['Rep1','Rep2','Rep3','Rep4']].values.astype(float)
        _, p_val = ttest_ind(mutant_values, wt_values1, equal_var=False)
        p_values1.append(p_val)
df1['p_value'] = p_values1

def format_pvalue(p):
    if pd.isna(p):
        return ''
    else:
        return f'{p:.3f}'

df1['p_label'] = df1['p_value'].apply(format_pvalue)

# --------------------------
# DATASET 2
# --------------------------
data2 = {
    'Catalyst': ['WT', 'W96A', 'E7A', 'A11L', 'L18A', 'N19A', 'K22A', 'N88A', 
                 'M89A', 'A92E', 'F93A', 'S95A', 'S97A', 'D100A', 'V99A', 
                 'I16A', 'N14A', 'L9A', 'R10A', 'K101A', 'E104A'],
    'Replicate_a': [17, 5, 17, 17, 16, 17, 14, 14, 13, 14, 7, 14, 17, 4, 19, 13, 16, 18, 16, 12, 11],
    'Replicate_b': [17, 5, 18, 19, 19, 17, 14, 13, 12, 14, 7, 16, 15, 3, 19, 9, 14, 17, 16, 11, 16],
    'Replicate_c': [17, 5, 16, 16, 17, 14, 13, 15, 13, 14, 7, 16, 18, 3, 13, 9, 10, 12, 11, 18, 15]
}
df2 = pd.DataFrame(data2)
df2['Average'] = df2[['Replicate_a','Replicate_b','Replicate_c']].mean(axis=1)
df2['Std_Error'] = df2[['Replicate_a','Replicate_b','Replicate_c']].sem(axis=1)

# Sort: WT first, then by average
wt2 = df2[df2['Catalyst']=='WT']
mutants2 = df2[df2['Catalyst']!='WT'].sort_values('Average').reset_index(drop=True)
df2 = pd.concat([wt2, mutants2]).reset_index(drop=True)
x2 = np.arange(len(df2))

# Compute p-values against WT
wt_values2 = wt2[['Replicate_a','Replicate_b','Replicate_c']].values.flatten().astype(float)
p_values2 = []
for i, row in df2.iterrows():
    if row['Catalyst'] == 'WT':
        p_values2.append(np.nan)
    else:
        variant_values = [row['Replicate_a'], row['Replicate_b'], row['Replicate_c']]
        _, p_val = ttest_ind(wt_values2, variant_values, equal_var=False)
        p_values2.append(p_val)
df2['p_value'] = p_values2
df2['p_label'] = df2['p_value'].apply(format_pvalue)

# Color coding for Dataset 2
def classify_color(value):
    if value < 10: return '#1f77b4'
    elif value < 15: return '#ff7f0e'
    else: return '#2ca02c'
bar_colors2 = df2['Average'].apply(classify_color).tolist()
bar_colors2[0] = '#800080'  # WT purple

# --------------------------
# DELTA YIELD DATA
# --------------------------
previous_data = {
    'Variant': ['L18A','K22A','F93A','S95A','S97A','E7A','A11L','N19A','N88A','M89A','A92E','D100A'],
    'DeltaYield': [-6.4,4.3,-11,-3.7,4.7,-1.7,-0.9,-3.5,-7.8,-4.8,-6.3,-19]
}
current_data = {
    'Variant': ['W96A','E7A','A11L','L18A','N19A','K22A','N88A','M89A','A92E','F93A','S95A','S97A','D100A','V99A','I16A','N14A','L9A','R10A','K101A','E104A'],
    'DeltaYield': [-12,0,0.5,-0.5,0,-3.5,-3.5,-4.5,-3,-10,-1.5,-0.5,-13.5,0,-6.5,-4,-1.5,-1,-3.5,-3]
}
df_prev = pd.DataFrame(previous_data)
df_curr = pd.DataFrame(current_data)

def extract_pos(variant):
    match = re.search(r'(\d+)', variant)
    return int(match.group(1)) if match else np.inf

df_prev['Position'] = df_prev['Variant'].apply(extract_pos)
df_curr['Position'] = df_curr['Variant'].apply(extract_pos)

overlap_variants = set(df_prev['Variant']) & set(df_curr['Variant'])
unique_prev = set(df_prev['Variant']) - set(df_curr['Variant'])
unique_curr = set(df_curr['Variant']) - set(df_prev['Variant'])

df_overlap = pd.merge(
    df_prev[['Variant','DeltaYield','Position']],
    df_curr[['Variant','DeltaYield']],
    on='Variant', suffixes=('_Old','_New')
).sort_values('Position')
df_unique_prev = df_prev[df_prev['Variant'].isin(unique_prev)].sort_values('Position')
df_unique_prev['DeltaYield_New'] = 0
df_unique_curr = df_curr[df_curr['Variant'].isin(unique_curr)].sort_values('Position')
df_unique_curr['DeltaYield_Old'] = 0
df_unique_curr = df_unique_curr.rename(columns={'DeltaYield':'DeltaYield_New'})

# --------------------------
# FIGURE SETUP
# --------------------------
fig, axs = plt.subplots(3,2,figsize=(24,20))
axs = axs.flatten()
fig.suptitle('Comprehensive Yield Analysis', fontsize=16, fontweight='bold', y=0.995)

# --- Plot 1: Dataset 1 all replicates ---
w1 = 0.18
axs[0].bar(x1 - 1.5*w1, df1['Rep1'], w1, label='Rep a')
axs[0].bar(x1 - 0.5*w1, df1['Rep2'], w1, label='Rep b')
axs[0].bar(x1 + 0.5*w1, df1['Rep3'], w1, label='Rep c')
axs[0].bar(x1 + 1.5*w1, df1['Rep4'], w1, label='Rep d')
axs[0].set_xticks(x1)
axs[0].set_xticklabels(df1['Catalyst'], rotation=45, ha='right', fontsize=6)
axs[0].set_ylabel('Yield (%)')
axs[0].set_title('Dataset 1: All Replicates', fontweight='bold')
axs[0].grid(axis='y', linestyle='--', alpha=0.3)
axs[0].legend(fontsize=6)

# --- Plot 2: Dataset 1 mean ± SEM with p-values ---
bar_colors1 = ['#8b5cf6' if c=='WT' else '#10b981' if m>=40 else '#f59e0b' if m>=30 else '#ef4444'
               for c,m in zip(df1['Catalyst'], df1['Mean'])]
axs[1].bar(x1, df1['Mean'], color=bar_colors1, edgecolor='black', linewidth=0.5)
axs[1].errorbar(x1, df1['Mean'], yerr=df1['SEM'], fmt='none', ecolor='black', capsize=2)

# Add p-values
for i, (idx, row) in enumerate(df1.iterrows()):
    if row['p_label']:
        y_pos = row['Mean'] + row['SEM'] + 2
        axs[1].text(i, y_pos, row['p_label'], ha='center', va='bottom', fontsize=6)

axs[1].set_xticks(x1)
axs[1].set_xticklabels(df1['Catalyst'], rotation=45, ha='right', fontsize=6)
axs[1].set_ylabel('Average Yield (%)')
axs[1].set_title('Dataset 1: Mean ± SEM', fontweight='bold')
axs[1].grid(axis='y', linestyle='--', alpha=0.3)
legend_patches1 = [Patch(facecolor='#8b5cf6', label='WT'),
                   Patch(facecolor='#10b981', label='High ≥40'),
                   Patch(facecolor='#f59e0b', label='Medium 30–39'),
                   Patch(facecolor='#ef4444', label='Low <30')]
axs[1].legend(handles=legend_patches1, fontsize=5, loc='upper left')

# --- Plot 3: Dataset 2 all replicates ---
width2 = 0.25
axs[2].bar(x2 - width2, df2['Replicate_a'], width2, label='Rep a')
axs[2].bar(x2, df2['Replicate_b'], width2, label='Rep b')
axs[2].bar(x2 + width2, df2['Replicate_c'], width2, label='Rep c')
axs[2].set_xticks(x2)
axs[2].set_xticklabels(df2['Catalyst'], rotation=45, ha='right', fontsize=5)
axs[2].set_ylabel('Yield (%)')
axs[2].set_title('Dataset 2: All Replicates', fontweight='bold')
axs[2].grid(axis='y', linestyle='--', alpha=0.2)
axs[2].legend(fontsize=5)

# --- Plot 4: Dataset 2 mean ± SEM with p-values ---
axs[3].bar(x2, df2['Average'], color=bar_colors2, edgecolor='black', linewidth=0.5)
axs[3].errorbar(x2, df2['Average'], yerr=df2['Std_Error'], fmt='none', ecolor='black', capsize=4)

# Add p-values
for i, (idx, row) in enumerate(df2.iterrows()):
    if row['p_label']:
        y_pos = row['Average'] + row['Std_Error'] + 0.5
        axs[3].text(i, y_pos, row['p_label'], ha='center', va='bottom', fontsize=5)

axs[3].set_xticks(x2)
axs[3].set_xticklabels(df2['Catalyst'], rotation=45, ha='right', fontsize=5)
axs[3].set_ylabel('Average Yield (%)')
axs[3].set_title('Dataset 2: Mean ± SEM', fontweight='bold')
axs[3].grid(axis='y', linestyle='--', alpha=0.2)
legend_patches2 = [Patch(facecolor='#800080', label='WT'),
                   Patch(facecolor='#1f77b4', label='Low <10'),
                   Patch(facecolor='#ff7f0e', label='Medium 10–15'),
                   Patch(facecolor='#2ca02c', label='High ≥15')]
axs[3].legend(handles=legend_patches2, fontsize=5, loc='upper left')

# --- Plot 5: Overlapping ΔYield ---
variants_overlap = df_overlap['Variant'].tolist()
x_overlap = np.arange(len(variants_overlap))
width = 0.35
axs[4].bar(x_overlap - width/2, df_overlap['DeltaYield_Old'], width, color='#2ca02c', label='Old')
axs[4].bar(x_overlap + width/2, df_overlap['DeltaYield_New'], width, color='#1f77b4', label='New')
axs[4].set_xticks(x_overlap)
axs[4].set_xticklabels(variants_overlap, rotation=45, ha='right', fontsize=9)
axs[4].set_ylabel('ΔYield (%)')
axs[4].set_title('Alanine Scan: Overlapping Variants', fontweight='bold')
axs[4].grid(axis='y', alpha=0.3)
axs[4].axhline(0, color='black', linewidth=0.5)
axs[4].legend(fontsize=8)
corr, _ = pearsonr(df_overlap['DeltaYield_Old'], df_overlap['DeltaYield_New'])
axs[4].text(0.98, 0.98, f'R² = {corr**2:.3f}', transform=axs[4].transAxes,
            ha='right', va='top', fontsize=8,
            bbox=dict(boxstyle='round', facecolor='white', edgecolor='black', alpha=0.8))

# --- Plot 6: Unique variants ---
unique_variants = list(df_unique_prev['Variant']) + list(df_unique_curr['Variant'])
x_unique = np.arange(len(unique_variants))
y_new = list(df_unique_prev['DeltaYield_New']) + list(df_unique_curr['DeltaYield_New'])
axs[5].bar(x_unique, y_new, color='#FF6B35', label='Unique')
axs[5].set_xticks(x_unique)
axs[5].set_xticklabels(unique_variants, rotation=45, ha='right', fontsize=9)
axs[5].set_ylabel('ΔYield (%)')
axs[5].set_title('Alanine Scan: Unique Variants', fontweight='bold')
axs[5].axhline(0, color='black', linewidth=0.5)
axs[5].grid(axis='y', alpha=0.3)
axs[5].legend(fontsize=8)

plt.tight_layout()
plt.savefig("Comprehensive_Yield_Analysis_with_R2.png", dpi=400)
plt.show()