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
wt_values = df[df['Catalyst'] == 'WT'][['Replicate_a', 'Replicate_b', 'Replicate_c']].values.flatten()

p_values = []
for idx, row in df.iterrows():
    if row['Catalyst'] == 'WT':
        p_values.append(1.0)
    else:
        variant_values = [row['Replicate_a'], row['Replicate_b'], row['Replicate_c']]
        _, p_val = stats.ttest_ind(wt_values, variant_values, equal_var=False)
        p_values.append(p_val)

df['p_value'] = p_values

# --------------------------
# PLOTTING
# --------------------------
fig, axes = plt.subplots(1, 2, figsize=(16, 6))
fig.suptitle('Experimental Yield Analysis with p-values (vs WT)', fontsize=14, fontweight='bold')

x = np.arange(len(df))
width = 0.25
p_fontsize = 4

# Plot 1: All replicates
ax1 = axes[0]
ax1.bar(x - width, df['Replicate_a'], width, label='Replicate a', alpha=0.8)
ax1.bar(x, df['Replicate_b'], width, label='Replicate b', alpha=0.8)
ax1.bar(x + width, df['Replicate_c'], width, label='Replicate c', alpha=0.8)
for i, row in df.iterrows():
    max_h = max(row['Replicate_a'], row['Replicate_b'], row['Replicate_c'])
    ax1.text(i, max_h + 0.3, f"p={row['p_value']:.3f}", ha='center', fontsize=p_fontsize)
ax1.set_xticks(x)
ax1.set_xticklabels(df['Catalyst'], rotation=45, ha='right', fontsize=7)
ax1.set_ylim(0, 27)
ax1.set_xlabel('Variants')
ax1.set_ylabel('Yield (%)')
ax1.grid(axis='y', alpha=0.2)
ax1.legend(fontsize=7)

# Plot 2: Average ± SEM
ax2 = axes[1]
bars = ax2.bar(x, df['Average'], alpha=0.7, edgecolor='black')
ax2.errorbar(x, df['Average'], yerr=df['Std_Error'], fmt='none', ecolor='black', capsize=4, elinewidth=1.3)
for i, row in df.iterrows():
    h = row['Average'] + row['Std_Error']
    ax2.text(i, h + 0.3, f"p={row['p_value']:.3f}", ha='center', fontsize=p_fontsize)
ax2.set_xticks(x)
ax2.set_xticklabels(df['Catalyst'], rotation=45, ha='right', fontsize=7)
ax2.set_ylim(0, 27)
ax2.set_xlabel('Variants')
ax2.set_ylabel('Average Yield (%)')
ax2.grid(axis='y', alpha=0.2)

plt.tight_layout()

# --------------------------
# SAVE FIGURES (ENSURE PATH EXISTS)
# --------------------------
save_paths = [
    r"D:\PhD_Thesis\LmrR_EVB\concerted\analysis\yield_analysis_variants.png",
    r"D:\PhD_Thesis\LmrR_Paper\figures\yield_analysis_variants.png"
]

for path in save_paths:
    os.makedirs(os.path.dirname(path), exist_ok=True)  # create directory if it doesn't exist
    fig.savefig(path, dpi=300, bbox_inches='tight')
    print(f"Saved figure to: {path}")

plt.show()
