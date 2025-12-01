import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from matplotlib.gridspec import GridSpec
from scipy.stats import ttest_ind, ttest_rel, pearsonr
import re
import os

# --------------------------
# FONT SETTINGS
# --------------------------
plt.rcParams.update({
    "font.size": 6,
    "xtick.labelsize": 5,
    "ytick.labelsize": 5,
    "legend.fontsize": 6
})

# =============================================================================
# DATASET 1: First experiment data
# =============================================================================
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

# Separate WT and sort mutants
wt1 = df1[df1['Catalyst']=='WT']
mutants1 = df1[df1['Catalyst']!='WT'].sort_values('Mean').reset_index(drop=True)
df1 = pd.concat([wt1, mutants1]).reset_index(drop=True)
x1 = np.arange(len(df1))

# p-values vs WT
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
df1['p_label'] = df1['p_value'].apply(lambda p: f"{p:.3f}" if pd.notna(p) else "")

# =============================================================================
# DATASET 2: Second experiment data
# =============================================================================
data2 = {
    'Catalyst': ['WT', 'W96A', 'E7A', 'A11L', 'L18A', 'N19A', 'K22A', 'N88A', 
                 'M89A', 'A92E', 'F93A', 'S95A', 'S97A', 'D100A', 'V99A', 
                 'I16A', 'N14A', 'L9A', 'R10A', 'K101A', 'E104A'],
    'Replicate_a': [17, 5, 17, 17, 16, 17, 14, 14, 13, 14, 7, 14, 17, 4, 19, 13, 16, 18, 16, 12, 11],
    'Replicate_b': [17, 5, 18, 19, 19, 17, 14, 13, 12, 14, 7, 16, 15, 3, 19, 9, 14, 17, 16, 11, 16],
    'Replicate_c': [17, 5, 16, 16, 17, 14, 13, 15, 13, 14, 7, 16, 18, 3, 13, 9, 10, 12, 11, 18, 15]
}

df2 = pd.DataFrame(data2)
df2['Average'] = df2[['Replicate_a', 'Replicate_b', 'Replicate_c']].mean(axis=1)
df2['Std_Error'] = df2[['Replicate_a', 'Replicate_b', 'Replicate_c']].sem(axis=1)

# Separate WT and sort mutants
wt2 = df2[df2['Catalyst']=='WT']
mutants2 = df2[df2['Catalyst']!='WT'].sort_values('Average').reset_index(drop=True)
df2 = pd.concat([wt2, mutants2]).reset_index(drop=True)
x2 = np.arange(len(df2))

# p-values vs WT
wt_values2 = wt2[['Replicate_a','Replicate_b','Replicate_c']].values.flatten().astype(float)
p_values2 = []
for i, row in df2.iterrows():
    if row['Catalyst'] == 'WT':
        p_values2.append(1.0)
    else:
        variant_values = [row['Replicate_a'], row['Replicate_b'], row['Replicate_c']]
        _, p_val = ttest_ind(wt_values2, variant_values, equal_var=False)
        p_values2.append(p_val)
df2['p_value'] = p_values2

# Color classification for Dataset 2
def classify_color(value):
    if value < 10:
        return '#1f77b4'  # BLUE = low
    elif value < 15:
        return '#ff7f0e'  # ORANGE = medium
    else:
        return '#2ca02c'  # GREEN = high
bar_colors2 = df2['Average'].apply(classify_color).tolist()
bar_colors2[0] = '#800080'  # WT = purple

# =============================================================================
# DELTA YIELD COMPARISON DATA
# =============================================================================
previous_data = {
    'Variant': ['L18A', 'K22A', 'F93A', 'S95A', 'S97A',
                'E7A', 'A11L', 'N19A', 'N88A', 'M89A', 'A92E', 'D100A'],
    'DeltaYield': [-6.4, 4.3, -11, -3.7, 4.7, -1.7, -0.9, -3.5, -7.8, -4.8, -6.3, -19]
}
df_previous = pd.DataFrame(previous_data)

current_data = {
    'Variant': ['W96A', 'E7A', 'A11L', 'L18A', 'N19A', 'K22A',
                'N88A', 'M89A', 'A92E', 'F93A', 'S95A', 'S97A',
                'D100A', 'V99A', 'I16A', 'N14A', 'L9A', 'R10A',
                'K101A', 'E104A'],
    'DeltaYield': [-12.0, 0.0, 0.5, -0.5, 0.0, -3.5,
                   -3.5, -4.5, -3.0, -10.0, -1.5, -0.5,
                   -13.5, 0.0, -6.5, -4.0, -1.5, -1.0,
                   -3.5, -3.0]
}
df_current = pd.DataFrame(current_data)

def extract_position(variant):
    match = re.search(r'(\d+)', variant)
    return int(match.group(1)) if match else np.inf

df_previous['Position'] = df_previous['Variant'].apply(extract_position)
df_current['Position'] = df_current['Variant'].apply(extract_position)

previous_variants = set(df_previous['Variant'])
current_variants = set(df_current['Variant'])
overlapping_variants = previous_variants & current_variants
unique_to_previous = previous_variants - current_variants
unique_to_current = current_variants - previous_variants

df_overlap = pd.merge(
    df_previous[['Variant', 'DeltaYield', 'Position']],
    df_current[['Variant', 'DeltaYield']],
    on='Variant',
    suffixes=('_Old', '_New')
).sort_values('Position')

df_unique_old = df_previous[df_previous['Variant'].isin(unique_to_previous)].sort_values('Position')
df_unique_old['DeltaYield_New'] = 0

df_unique_new = df_current[df_current['Variant'].isin(unique_to_current)].sort_values('Position')
df_unique_new['DeltaYield_Old'] = 0
df_unique_new = df_unique_new.rename(columns={'DeltaYield': 'DeltaYield_New'})

# =============================================================================
# FIGURE SETUP (6-panel 3x2)
# =============================================================================
fig = plt.figure(figsize=(24, 20))
fig.suptitle('Comprehensive Experimental Yield Analysis', fontsize=16, fontweight='bold', y=0.995)

gs = GridSpec(3, 2, figure=fig, hspace=0.40, wspace=0.25, top=0.96)
ax1 = fig.add_subplot(gs[0,0])
ax2 = fig.add_subplot(gs[0,1])
ax3 = fig.add_subplot(gs[1,0])
ax4 = fig.add_subplot(gs[1,1])
ax5 = fig.add_subplot(gs[2,0])
ax6 = fig.add_subplot(gs[2,1])

# --------------------------
# PLOT 1: Dataset 1 - All Replicates
# --------------------------
w1 = 0.18
ax1.bar(x1 - 1.5*w1, df1['Rep1'], w1, label='Replicate a', alpha=0.85)
ax1.bar(x1 - 0.5*w1, df1['Rep2'], w1, label='Replicate b', alpha=0.85)
ax1.bar(x1 + 0.5*w1, df1['Rep3'], w1, label='Replicate c', alpha=0.85)
ax1.bar(x1 + 1.5*w1, df1['Rep4'], w1, label='Replicate d', alpha=0.85)

for i, catalyst in enumerate(df1['Catalyst']):
    top = max(df1.loc[i, ['Rep1','Rep2','Rep3','Rep4']])
    ax1.text(i, top + 1.5, catalyst, ha='center', fontsize=8, fontweight='bold')

ax1.set_title('Dataset 1: All Replicates', fontsize=10, fontweight='bold')
ax1.set_ylabel('Yield (%)')
ax1.set_xticks(x1)
ax1.set_xticklabels(df1['Catalyst'], rotation=45, ha='right', fontsize=6)
ax1.legend(fontsize=6)
ax1.grid(axis='y', alpha=0.3, linestyle='--')
ax1.set_ylim(0, 60)

# --------------------------
# PLOT 2: Dataset 1 - Mean ± SEM
# --------------------------
bar_colors1 = ['#8b5cf6' if c=='WT' else 
               '#10b981' if m >= 40 else
               '#f59e0b' if m >= 30 else '#ef4444'
               for c,m in zip(df1['Catalyst'], df1['Mean'])]

ax2.bar(x1, df1['Mean'], color=bar_colors1, alpha=0.75, edgecolor='black', linewidth=0.5)
ax2.errorbar(x1, df1['Mean'], yerr=df1['SEM'], fmt='none', ecolor='black', elinewidth=1, capsize=2)
for i, (mean, sem, p_label) in enumerate(zip(df1['Mean'], df1['SEM'], df1['p_label'])):
    if p_label: ax2.text(i, mean + sem + 2, p_label, ha='center', fontsize=6, color='black')

ax2.set_title('Dataset 1: Mean ± SEM', fontsize=10, fontweight='bold')
ax2.set_xticks(x1)
ax2.set_xticklabels(df1['Catalyst'], rotation=45, ha='right', fontsize=6)
ax2.set_ylabel('Average Yield (%)')
ax2.grid(axis='y', alpha=0.3, linestyle='--')
ax2.set_ylim(0, 60)
legend_patches1 = [
    Patch(facecolor='#8b5cf6', label='Wild Type'),
    Patch(facecolor='#10b981', label='High (≥40%)'),
    Patch(facecolor='#f59e0b', label='Medium (30–39%)'),
    Patch(facecolor='#ef4444', label='Low (<30%)')
]
ax2.legend(handles=legend_patches1, fontsize=5, loc='upper left')

# --------------------------
# PLOT 3: Dataset 2 - All Replicates
# --------------------------
width2 = 0.25
ax3.bar(x2 - width2, df2['Replicate_a'], width2, label='Replicate a', alpha=0.8)
ax3.bar(x2, df2['Replicate_b'], width2, label='Replicate b', alpha=0.8)
ax3.bar(x2 + width2, df2['Replicate_c'], width2, label='Replicate c', alpha=0.8)
for i, variant in enumerate(df2['Catalyst']):
    max_h = max(df2.loc[i, ['Replicate_a','Replicate_b','Replicate_c']])
    ax3.text(i, max_h + 0.5, variant, ha='center', va='bottom', fontsize=6, fontweight='bold')
ax3.set_title('Dataset 2: All Replicates', fontsize=10, fontweight='bold')
ax3.set_xticks(x2)
ax3.set_xticklabels(df2['Catalyst'], rotation=45, ha='right', fontsize=5)
ax3.set_ylabel('Yield (%)')
ax3.set_ylim(0, 27)
ax3.set_xlabel('Variants')
ax3.grid(axis='y', alpha=0.2, linestyle='--')
ax3.legend(fontsize=5, loc='upper left')

# --------------------------
# PLOT 4: Dataset 2 - Average ± SEM
# --------------------------
ax4.bar(x2, df2['Average'], color=bar_colors2, alpha=0.7, edgecolor='black', linewidth=0.5)
ax4.errorbar(x2, df2['Average'], yerr=df2['Std_Error'], fmt='none', ecolor='black', capsize=4, elinewidth=1.3)
for i, row in df2.iterrows():
    ax4.text(i, row['Average'] + row['Std_Error'] + 0.5, f"{row['p_value']:.3f}", ha='center', fontsize=5)

ax4.set_title('Dataset 2: Average ± SEM', fontsize=10, fontweight='bold')
ax4.set_xticks(x2)
ax4.set_xticklabels(df2['Catalyst'], rotation=45, ha='right', fontsize=5)
ax4.set_ylim(0, 27)
ax4.set_xlabel('Variants')
ax4.set_ylabel('Average Yield (%)')
ax4.grid(axis='y', alpha=0.2, linestyle='--')

legend_patches2 = [
    Patch(facecolor='#800080', label='Wild Type'),
    Patch(facecolor='#1f77b4', label='Low yield (<10%)'),
    Patch(facecolor='#ff7f0e', label='Medium (10–15%)'),
    Patch(facecolor='#2ca02c', label='High yield (≥15%)')
]
ax4.legend(handles=legend_patches2, fontsize=5, loc='upper left')

# --------------------------
# PLOT 5: Alanine Scanning - Overlapping
# --------------------------
variants_overlap = df_overlap['Variant'].tolist()
x_overlap = np.arange(len(variants_overlap))
width = 0.35
ax5.bar(x_overlap - width/2, df_overlap['DeltaYield_Old'], width, color='#2ca02c', edgecolor='black', linewidth=0.5, label='Old Dataset (WT: 40%)')
ax5.bar(x_overlap + width/2, df_overlap['DeltaYield_New'], width, color='#1f77b4', edgecolor='black', linewidth=0.5, label='New Dataset (WT: 17%)')
ax5.set_xticks(x_overlap)
ax5.set_xticklabels(variants_overlap, rotation=45, ha='right', fontsize=9)
ax5.set_ylabel('ΔYield (%)', fontsize=10)
ax5.legend(fontsize=8, loc='lower left')
ax5.grid(axis='y', alpha=0.3)
ax5.axhline(y=0, color='black', linestyle='-', linewidth=0.5)

# R² removed from this plot

# --------------------------
# PLOT 6: Alanine Scanning - Unique
# --------------------------
variants_unique_old = df_unique_old['Variant'].tolist()
variants_unique_new = df_unique_new['Variant'].tolist()
all_unique = variants_unique_old + variants_unique_new
x_unique = np.arange(len(all_unique))
y_new = [0]*len(variants_unique_old) + list(df_unique_new['DeltaYield_New'])
ax6.bar(x_unique + width/2, y_new, width, color='#FF6B35', edgecolor='black', linewidth=0.5, label='New Mutants Tested (WT: 17%)')
ax6.set_xticks(x_unique)
ax6.set_xticklabels(all_unique, rotation=45, ha='right', fontsize=9)
ax6.set_ylabel('ΔYield (%)', fontsize=10)
ax6.set_title('Unique Mutants', fontsize=10, fontweight='bold')
ax6.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
ax6.grid(axis='y', alpha=0.3)
ax6.legend(fontsize=8, loc='lower left')

# R² removed from this plot

# --------------------------
# SAVE FIGURE
# --------------------------
plt.tight_layout()
plt.savefig("Comprehensive_Yield_Analysis.png", dpi=400)
plt.show()