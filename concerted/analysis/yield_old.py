import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from scipy.stats import ttest_ind
import os

# --------------------------
# SMALL FONT SETTINGS
# --------------------------
plt.rcParams.update({
    "font.size": 6,
    "xtick.labelsize": 5,
    "ytick.labelsize": 5,
    "legend.fontsize": 6
})

# --------------------------
# DATA
# --------------------------
data = {
    'Catalyst': ['WT', 'E7A', 'A11L', 'L18A', 'N19A', 'K22A', 'N88A', 
                 'M89A', 'A92E', 'F93A', 'S95A', 'S97A', 'D100A'],
    'Rep1': [38, 38, 38, 28, 35, 49, 44, 34, 36, 33, 45, 50, 15],
    'Rep2': [37, 40, 40, 33, 36, 51, 39, 32, 37, 34, 38, 49, 12],
    'Rep3': [47, 40, 43, 40, 48, 42, 26, 40, 35, 26, 34, 43, 21],
    'Rep4': [46, 43, 43, 42, 47, 43, 37, 42, 34, 26, 36, 44, 41]
}

df = pd.DataFrame(data)
df['Mean'] = df[['Rep1', 'Rep2', 'Rep3', 'Rep4']].mean(axis=1)
df['SEM'] = df[['Rep1', 'Rep2', 'Rep3', 'Rep4']].sem(axis=1)

# --------------------------
# CALCULATE P-VALUES VS WT
# --------------------------
wt_values = df.loc[df['Catalyst']=='WT', ['Rep1','Rep2','Rep3','Rep4']].values.flatten().astype(float)

p_values = []
for i, row in df.iterrows():
    if row['Catalyst'] == 'WT':
        p_values.append(np.nan)
    else:
        mutant_values = row[['Rep1','Rep2','Rep3','Rep4']].values.astype(float)
        t_stat, p_val = ttest_ind(mutant_values, wt_values, equal_var=False)
        p_values.append(p_val)

df['p_value'] = p_values

# --------------------------
# FORMAT P-VALUES AS NUMERIC
# --------------------------
def format_pval_numeric(p):
    if pd.isna(p):
        return ""
    else:
        return f"{p:.3f}"  # shows 3 decimal places

df['p_label'] = df['p_value'].apply(format_pval_numeric)

# --------------------------
# FIGURE
# --------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
fig.suptitle('Experimental Yield Analysis', fontsize=9)

x = np.arange(len(df))
w = 0.18

# --------------------------
# Plot 1: All Replicates
# --------------------------
ax1.bar(x - 1.5*w, df['Rep1'], w, label='Replicate a', alpha=0.85)
ax1.bar(x - 0.5*w, df['Rep2'], w, label='Replicate b', alpha=0.85)
ax1.bar(x + 0.5*w, df['Rep3'], w, label='Replicate c', alpha=0.85)
ax1.bar(x + 1.5*w, df['Rep4'], w, label='Replicate d', alpha=0.85)

# Add catalyst labels
for i, catalyst in enumerate(df['Catalyst']):
    top = max(df.loc[i, ['Rep1', 'Rep2', 'Rep3', 'Rep4']])
    ax1.text(i, top + 1.2, catalyst, ha='center', fontsize=5)
    ax1.text(i, -3.0, catalyst, ha='center', fontsize=5)

ax1.set_title('All Replicates', fontsize=8)
ax1.set_ylabel('Yield (%)')
ax1.set_xticks(x)
ax1.set_xticklabels([])
ax1.legend(loc='upper right', fontsize=5)
ax1.grid(axis='y', alpha=0.3, linestyle='--')
ax1.set_ylim(-5, 60)

# --------------------------
# Plot 2: Mean + SEM
# --------------------------
# Color logic
bar_colors = []
for _, row in df.iterrows():
    if row['Catalyst'] == 'WT':
        bar_colors.append('#8b5cf6')
    elif row['Mean'] >= 40:
        bar_colors.append('#10b981')
    elif row['Mean'] >= 30:
        bar_colors.append('#f59e0b')
    else:
        bar_colors.append('#ef4444')

bars = ax2.bar(x, df['Mean'], color=bar_colors, alpha=0.75, edgecolor='black', linewidth=0.5)

ax2.errorbar(
    x, df['Mean'], yerr=df['SEM'], fmt='none',
    ecolor='black', elinewidth=1, capsize=2
)

# Add catalyst labels and numeric p-values
for i, (mean, sem, catalyst, p_label) in enumerate(zip(df['Mean'], df['SEM'], df['Catalyst'], df['p_label'])):
    ax2.text(i, mean + sem + 1.5, catalyst, ha='center', fontsize=5)
    if p_label:
        ax2.text(i, mean + sem + 5, p_label, ha='center', fontsize=6, color='red')
    ax2.text(i, -3.0, catalyst, ha='center', fontsize=5)

ax2.set_title('Mean Yield ± SEM', fontsize=8)
ax2.set_ylabel('Average Yield (%)')
ax2.set_xticks(x)
ax2.set_xticklabels([])
ax2.grid(axis='y', alpha=0.3, linestyle='--')
ax2.set_ylim(-5, 60)

# Color legend
legend_patches = [
    Patch(facecolor='#8b5cf6', label='Wild Type'),
    Patch(facecolor='#10b981', label='High (≥40%)'),
    Patch(facecolor='#f59e0b', label='Medium (30–39%)'),
    Patch(facecolor='#ef4444', label='Low (<30%)')
]
ax2.legend(handles=legend_patches, fontsize=6, loc='upper right')

plt.tight_layout()

# --------------------------
# SAVE FIGURES
# --------------------------
save_paths = [
    r"D:\PhD_Thesis\LmrR_EVB\concerted\analysis\yield_analysis.png",
    r"D:\PhD_Thesis\LmrR_Paper\figures\yield_analysis.png"
]

for path in save_paths:
    # Ensure directory exists
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches='tight')

plt.show()
