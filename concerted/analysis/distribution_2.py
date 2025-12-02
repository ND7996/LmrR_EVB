import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch
from scipy.stats import ttest_ind

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
df1['SEM']   = df1[['Rep1','Rep2','Rep3','Rep4']].sem(axis=1)

# --------------------------
# DATASET 2
# --------------------------
data2 = {
    'Catalyst': ['WT','W96A','E7A','A11L','L18A','N19A','K22A','N88A',
                 'M89A','A92E','F93A','S95A','S97A','D100A','V99A',
                 'I16A','N14A','L9A','R10A','K101A','E104A'],
    'Replicate_a': [17,5,17,17,16,17,14,14,13,14,7,14,17,4,19,13,16,18,16,12,11],
    'Replicate_b': [17,5,18,19,19,17,14,13,12,14,7,16,15,3,19,9,14,17,16,11,16],
    'Replicate_c': [17,5,16,16,17,14,13,15,13,14,7,16,18,3,13,9,10,12,11,18,15]
}

df2 = pd.DataFrame(data2)
df2['Mean'] = df2[['Replicate_a','Replicate_b','Replicate_c']].mean(axis=1)
df2['SEM']  = df2[['Replicate_a','Replicate_b','Replicate_c']].sem(axis=1)

# --------------------------
# Sort by mean for violin plots (ascending) with WT first
# --------------------------
def sort_wt_first_then_ascending(df):
    wt = df[df['Catalyst'] == 'WT']
    others = df[df['Catalyst'] != 'WT'].sort_values('Mean')
    return pd.concat([wt, others]).reset_index(drop=True)

df1_sorted = sort_wt_first_then_ascending(df1)
df2_sorted = sort_wt_first_then_ascending(df2)

# --------------------------
# VIOLIN PLOT FUNCTION
# --------------------------
def plot_violin(ax, df, rep_cols, wt_color, mut_color,
                special_variants=None, special_color='green', title=""):

    if special_variants is None:
        special_variants = []

    violin_data = [df.loc[i, rep_cols].values.astype(float) for i in range(len(df))]
    positions = np.arange(len(df))

    parts = ax.violinplot(violin_data, positions=positions, showmeans=True,
                          showextrema=True, widths=0.7)

    # Color violins
    for i, body in enumerate(parts['bodies']):
        cat = df.loc[i, 'Catalyst']
        if cat == 'WT':
            color = wt_color
        elif cat in special_variants:
            color = special_color
        else:
            color = mut_color
        body.set_facecolor(color)
        body.set_alpha(0.7)
        body.set_edgecolor('black')

    # Scatter on top
    for i, vals in enumerate(violin_data):
        x = np.random.normal(i, 0.04, size=len(vals))
        cat = df.loc[i, 'Catalyst']
        if cat == 'WT':
            color = wt_color
        elif cat in special_variants:
            color = special_color
        else:
            color = mut_color
        ax.scatter(x, vals, s=40, color=color, edgecolor='black', zorder=3)

    ax.set_xticks(positions)
    ax.set_xticklabels(df['Catalyst'], rotation=45, ha='right')
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.set_ylabel("Yield (%)")
    ax.grid(axis='y', linestyle='--', alpha=0.3)


# --------------------------
# P-VALUE PANEL (3rd & 4th)
# --------------------------
def pvalue_plot(ax, df, rep_cols, wt_values, title):

    pvals = []

    # Force WT numeric
    wt_values = pd.to_numeric(wt_values, errors='coerce').astype(float)
    wt_values = wt_values[~np.isnan(wt_values)]

    for _, row in df.iterrows():
        vals = row[rep_cols].values

        # Force numeric
        vals = pd.to_numeric(vals, errors='coerce').astype(float)
        vals = vals[~np.isnan(vals)]

        if len(vals) < 2 or len(wt_values) < 2:
            pvals.append(np.nan)
            continue

        p = ttest_ind(vals, wt_values, equal_var=False)[1]
        pvals.append(p)

    ax.bar(df['Catalyst'], pvals, color="teal")
    ax.axhline(0.05, color='red', linestyle='--')
    ax.set_yscale('log')
    ax.set_ylabel("p-value (log scale)")
    ax.set_title(title)
    ax.set_xticklabels(df['Catalyst'], rotation=45)


# --------------------------
# CREATE 4-PANEL FIGURE (2x2)
# --------------------------
fig, axes = plt.subplots(2, 2, figsize=(18, 10))
ax1, ax2, ax3, ax4 = axes.flatten()

fig.suptitle("Yield Analysis - Violin Plots and Distributions", fontsize=18, fontweight='bold')

# ---- Panels 1 & 2: Violin ----
plot_violin(ax1, df1_sorted, ['Rep1','Rep2','Rep3','Rep4'],
            wt_color='purple', mut_color='red',
            title="Dataset 1 – Violin (WT first, then Low to High Yield)")

special_list = ['W96A','V99A','I16A','N14A','L9A','R10A','K101A','E104A']

plot_violin(ax2, df2_sorted,
            ['Replicate_a','Replicate_b','Replicate_c'],
            wt_color='purple', mut_color='blue',
            special_variants=special_list, special_color='yellow',
            title="Dataset 2 – Violin (WT first, then Low to High Yield)")

# ---- Panels 3 & 4: Individual Histograms ----
data1_flat = df1[['Rep1','Rep2','Rep3','Rep4']].values.flatten()
data2_flat = df2[['Replicate_a','Replicate_b','Replicate_c']].values.flatten()

sns.histplot(data1_flat, bins=10, kde=True, ax=ax3, color='red', alpha=0.6)
ax3.set_title("Dataset 1 – Yield Distribution", fontsize=12, fontweight='bold')
ax3.set_xlabel("Yield (%)", fontsize=11)
ax3.set_ylabel("Frequency (Number of Observations)", fontsize=11)

sns.histplot(data2_flat, bins=10, kde=True, ax=ax4, color='blue', alpha=0.6)
ax4.set_title("Dataset 2 – Yield Distribution", fontsize=12, fontweight='bold')
ax4.set_xlabel("Yield (%)", fontsize=11)
ax4.set_ylabel("Frequency (Number of Observations)", fontsize=11)

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig('yield_analysis_4panel.png', dpi=300, bbox_inches='tight')
plt.savefig('yield_analysis_4panel.pdf', bbox_inches='tight')
print("Figures saved as 'yield_analysis_4panel.png' and 'yield_analysis_4panel.pdf'")
plt.show()