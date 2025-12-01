import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
import os

# --------------------------
# DATA
# --------------------------
data = {
    'Catalyst': ['WT', 'W96A', 'E7A', 'A11L', 'L18A', 'N19A', 'K22A', 'N88A', 
                 'M89A', 'A92E', 'F93A', 'S95A', 'S97A', 'D100A', 'V99A', 
                 'I16A', 'N14A', 'L9A', 'R10A', 'K101A', 'E104A'],
    'Replicate_a': [17, 5, 17, 17, 16, 17, 14, 14, 13, 14, 7, 14, 17, 4, 19, 13, 16, 18, 16, 12, 11],
    'Replicate_b': [17, 5, 18, 19, 19, 17, 14, 13, 12, 14, 7, 16, 15, 3, 19, 9, 14, 17, 16, 11, 16],
    'Replicate_c': [17, 5, 16, 16, 17, 14, 13, 15, 13, 14, 7, 16, 18, 3, 13, 9, 10, 12, 11, 18, 15]
}

df = pd.DataFrame(data)
df['Average'] = df[['Replicate_a', 'Replicate_b', 'Replicate_c']].mean(axis=1)
df['Std_Error'] = df[['Replicate_a', 'Replicate_b', 'Replicate_c']].sem(axis=1)

# --------------------------
# CALCULATE P-VALUES VS WT
# --------------------------
wt_values = df.loc[df['Catalyst'] == 'WT', ['Replicate_a', 'Replicate_b', 'Replicate_c']].values.flatten()

df['p_value'] = df.apply(
    lambda row: 1.0 if row['Catalyst'] == 'WT' 
    else stats.ttest_ind(wt_values, [row['Replicate_a'], row['Replicate_b'], row['Replicate_c']], equal_var=False)[1], 
    axis=1
)

# --------------------------
# COLOR RULES FOR RIGHT PLOT
# --------------------------
def classify_color(value):
    if value < 10:
        return '#1f77b4'  # BLUE
    elif value < 15:
        return '#ff7f0e'  # ORANGE
    else:
        return '#2ca02c'  # GREEN

bar_colors = df['Average'].apply(classify_color).tolist()
bar_colors[0] = '#800080'  # WT = purple

# --------------------------
# PLOTTING
# --------------------------
fig, axes = plt.subplots(1, 2, figsize=(16, 6))
fig.suptitle('Experimental Yield Analysis (vs WT)', fontsize=14, fontweight='bold')

x = np.arange(len(df))
width = 0.25
p_fontsize = 5

# --------------------------
# LEFT PLOT: All replicates (CLEANER)
# --------------------------
ax1 = axes[0]
replicates = ['Replicate_a', 'Replicate_b', 'Replicate_c']
for i, rep in enumerate(replicates):
    ax1.bar(x + (i-1)*width, df[rep], width, label=rep, alpha=0.8)

# Variant labels above the highest bar for each variant
for i, variant in enumerate(df['Catalyst']):
    max_height = max(df.loc[i, 'Replicate_a'], df.loc[i, 'Replicate_b'], df.loc[i, 'Replicate_c'])
    ax1.text(i, max_height + 0.5, variant, ha='center', va='bottom', fontsize=6, fontweight='bold')

ax1.set_xticks(x)
ax1.set_xticklabels(df['Catalyst'], rotation=45, ha='right', fontsize=7)
ax1.set_ylim(0, 27)
ax1.set_xlabel('Variants', fontsize=10)
ax1.set_ylabel('Yield (%)', fontsize=10)
ax1.grid(axis='y', alpha=0.2, linestyle='--')
ax1.legend(fontsize=8, loc='upper right')
ax1.set_title('All Replicates', fontsize=11, fontweight='bold', pad=10)

# --------------------------
# RIGHT PLOT: Average yields
# --------------------------
ax2 = axes[1]
bars = ax2.bar(x, df['Average'], color=bar_colors, alpha=0.7, edgecolor='black')
ax2.errorbar(x, df['Average'], yerr=df['Std_Error'], fmt='none', ecolor='black', capsize=4, elinewidth=1.3)

# P-values above error bars
for i, row in df.iterrows():
    h = row['Average'] + row['Std_Error']
    ax2.text(i, h + 0.5, f"{row['p_value']:.3f}", ha='center', va='bottom', fontsize=p_fontsize, fontweight='normal')

ax2.set_xticks(x)
ax2.set_xticklabels(df['Catalyst'], rotation=45, ha='right', fontsize=7)
ax2.set_ylim(0, 27)
ax2.set_xlabel('Variants', fontsize=10)
ax2.set_ylabel('Average Yield (%)', fontsize=10)
ax2.grid(axis='y', alpha=0.2, linestyle='--')
ax2.set_title('Average Yield ± SEM', fontsize=11, fontweight='bold', pad=10)

# Legend with Wild Type included
wt_patch = plt.Rectangle((0,0),1,1,color='#800080', label='Wild Type')
low_patch = plt.Rectangle((0,0),1,1,color='#1f77b4', label='Low yield (<10%)')
med_patch = plt.Rectangle((0,0),1,1,color='#ff7f0e', label='Medium (10-15%)')
high_patch = plt.Rectangle((0,0),1,1,color='#2ca02c', label='High yield (≥15%)')
ax2.legend(handles=[wt_patch, low_patch, med_patch, high_patch], fontsize=8, loc='upper right')

plt.tight_layout()

# --------------------------
# SAVE FIGURES
# --------------------------
save_paths = [
    r"D:\PhD_Thesis\LmrR_EVB\concerted\analysis\yield_analysis_variants.png",
    r"D:\PhD_Thesis\LmrR_Paper\figures\yield_analysis_variants.png"
]

for path in save_paths:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches='tight')
    print(f"Saved figure to: {path}")

plt.show()