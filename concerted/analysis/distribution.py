import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from matplotlib.gridspec import GridSpec
from scipy.stats import ttest_ind, ttest_rel, pearsonr, norm
import seaborn as sns
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
df1['Mean'] = df1[['Rep1', 'Rep2', 'Rep3', 'Rep4']].mean(axis=1)
df1['SEM'] = df1[['Rep1', 'Rep2', 'Rep3', 'Rep4']].sem(axis=1)

# Calculate p-values and delta yield vs WT for dataset 1
wt_values1 = df1.loc[df1['Catalyst']=='WT', ['Rep1','Rep2','Rep3','Rep4']].values.flatten().astype(float)
wt_mean1 = wt_values1.mean()

p_values1 = []
delta_yield1 = []
for i, row in df1.iterrows():
    if row['Catalyst'] == 'WT':
        p_values1.append(np.nan)
        delta_yield1.append(0)
    else:
        mutant_values = row[['Rep1','Rep2','Rep3','Rep4']].values.astype(float)
        t_stat, p_val = ttest_ind(mutant_values, wt_values1, equal_var=False)
        p_values1.append(p_val)
        delta_yield1.append(row['Mean'] - wt_mean1)

df1['p_value'] = p_values1
df1['delta_yield'] = delta_yield1
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

# Calculate p-values and delta yield vs WT for dataset 2
wt_values2 = df2[df2['Catalyst'] == 'WT'][['Replicate_a', 'Replicate_b', 'Replicate_c']].values.flatten()
wt_mean2 = wt_values2.mean()

p_values2 = []
delta_yield2 = []
for idx, row in df2.iterrows():
    if row['Catalyst'] == 'WT':
        p_values2.append(1.0)
        delta_yield2.append(0)
    else:
        variant_values = [row['Replicate_a'], row['Replicate_b'], row['Replicate_c']]
        _, p_val = ttest_ind(wt_values2, variant_values, equal_var=False)
        p_values2.append(p_val)
        delta_yield2.append(row['Average'] - wt_mean2)

df2['p_value'] = p_values2
df2['delta_yield'] = delta_yield2

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
# CREATE FIGURE WITH DISTRIBUTION PLOTS
# =============================================================================
fig = plt.figure(figsize=(26, 24))
fig.suptitle('Comprehensive Experimental Yield Analysis with Distributions', 
             fontsize=14, fontweight='bold', y=0.995)

gs = GridSpec(4, 2, figure=fig, hspace=0.35, wspace=0.25, top=0.96, bottom=0.05)

ax1 = fig.add_subplot(gs[0, 0])  # Dataset 1: Mean ± SEM
ax2 = fig.add_subplot(gs[0, 1])  # Dataset 2: Average ± SEM
ax3 = fig.add_subplot(gs[1, 0])  # Dataset 1: Delta Yield vs -log10(p-value)
ax4 = fig.add_subplot(gs[1, 1])  # Dataset 2: Delta Yield vs -log10(p-value)
ax5 = fig.add_subplot(gs[2, 0])  # Dataset 1: Distribution histograms
ax6 = fig.add_subplot(gs[2, 1])  # Dataset 2: Distribution histograms
ax7 = fig.add_subplot(gs[3, 0])  # Dataset 1: Violin plot
ax8 = fig.add_subplot(gs[3, 1])  # Dataset 2: Violin plot

# --------------------------
# PLOT 1: Dataset 1 - Mean + SEM (Improved)
# --------------------------
x1 = np.arange(len(df1))
bar_colors1 = []
for _, row in df1.iterrows():
    if row['Catalyst'] == 'WT':
        bar_colors1.append('#8b5cf6')
    elif row['Mean'] >= 40:
        bar_colors1.append('#10b981')
    elif row['Mean'] >= 30:
        bar_colors1.append('#f59e0b')
    else:
        bar_colors1.append('#ef4444')

bars1 = ax1.bar(x1, df1['Mean'], color=bar_colors1, alpha=0.75, edgecolor='black', linewidth=0.5)
ax1.errorbar(x1, df1['Mean'], yerr=df1['SEM'], fmt='none', ecolor='black', elinewidth=1, capsize=2)

for i, (mean, sem, p_label) in enumerate(zip(df1['Mean'], df1['SEM'], df1['p_label'])):
    if p_label:
        ax1.text(i, mean + sem + 2, p_label, ha='center', fontsize=6, color='black')

ax1.set_title('Dataset 1: Mean Yield ± SEM', fontsize=10, fontweight='bold')
ax1.set_ylabel('Average Yield (%)', fontsize=8)
ax1.set_xticks(x1)
ax1.set_xticklabels(df1['Catalyst'], rotation=45, ha='right', fontsize=7)
ax1.grid(axis='y', alpha=0.3, linestyle='--')
ax1.set_ylim(0, 60)

legend_patches1 = [
    Patch(facecolor='#8b5cf6', label='Wild Type'),
    Patch(facecolor='#10b981', label='High (≥40%)'),
    Patch(facecolor='#f59e0b', label='Medium (30–39%)'),
    Patch(facecolor='#ef4444', label='Low (<30%)')
]
ax1.legend(handles=legend_patches1, fontsize=6, loc='upper right')

# --------------------------
# PLOT 2: Dataset 2 - Average ± SEM
# --------------------------
x2 = np.arange(len(df2))
bars2 = ax2.bar(x2, df2['Average'], color=bar_colors2, alpha=0.7, edgecolor='black', linewidth=0.5)
ax2.errorbar(x2, df2['Average'], yerr=df2['Std_Error'], fmt='none', ecolor='black', capsize=4, elinewidth=1.3)

for i, row in df2.iterrows():
    h = row['Average'] + row['Std_Error']
    ax2.text(i, h + 0.5, f"{row['p_value']:.3f}", ha='center', va='bottom', fontsize=5, color='black')

ax2.set_title('Dataset 2: Average Yield ± SEM', fontsize=10, fontweight='bold')
ax2.set_xticks(x2)
ax2.set_xticklabels(df2['Catalyst'], rotation=45, ha='right', fontsize=7)
ax2.set_ylim(0, 27)
ax2.set_xlabel('Variants', fontsize=8)
ax2.set_ylabel('Average Yield (%)', fontsize=8)
ax2.grid(axis='y', alpha=0.2, linestyle='--')

wt_patch = Patch(facecolor='#800080', label='Wild Type')
low_patch = Patch(facecolor='#1f77b4', label='Low yield (<10%)')
med_patch = Patch(facecolor='#ff7f0e', label='Medium (10-15%)')
high_patch = Patch(facecolor='#2ca02c', label='High yield (≥15%)')
ax2.legend(handles=[wt_patch, low_patch, med_patch, high_patch], fontsize=6, loc='upper right')

# --------------------------
# PLOT 3: Dataset 1 - Delta Yield vs -log10(p-value) VOLCANO PLOT
# --------------------------
df1_filtered = df1[df1['Catalyst'] != 'WT'].copy()
# Handle log10 safely - replace inf/nan values
with np.errstate(divide='ignore', invalid='ignore'):
    df1_filtered['-log10(p)'] = -np.log10(df1_filtered['p_value'])
    df1_filtered['-log10(p)'] = df1_filtered['-log10(p)'].replace([np.inf, -np.inf], np.nan)

# Color by significance and effect
colors_volcano1 = []
for _, row in df1_filtered.iterrows():
    if row['p_value'] < 0.05 and abs(row['delta_yield']) > 5:
        colors_volcano1.append('#ef4444')  # Significant & large effect
    elif row['p_value'] < 0.05:
        colors_volcano1.append('#f59e0b')  # Significant
    else:
        colors_volcano1.append('#94a3b8')  # Not significant

ax3.scatter(df1_filtered['delta_yield'], df1_filtered['-log10(p)'], 
           c=colors_volcano1, s=100, alpha=0.7, edgecolors='black', linewidth=0.5)

# Add labels for significant variants
for _, row in df1_filtered.iterrows():
    if row['p_value'] < 0.05:
        ax3.annotate(row['Catalyst'], 
                    (row['delta_yield'], row['-log10(p)']),
                    fontsize=6, ha='center', va='bottom')

# Add reference lines
ax3.axhline(y=-np.log10(0.05), color='red', linestyle='--', linewidth=1, alpha=0.5, label='p=0.05')
ax3.axvline(x=5, color='blue', linestyle='--', linewidth=1, alpha=0.5, label='ΔYield=±5%')
ax3.axvline(x=-5, color='blue', linestyle='--', linewidth=1, alpha=0.5)

ax3.set_xlabel('ΔYield vs WT (%)', fontsize=8)
ax3.set_ylabel('-log₁₀(p-value)', fontsize=8)
ax3.set_title('Dataset 1: Volcano Plot (Effect Size vs Significance)', fontsize=10, fontweight='bold')
ax3.grid(alpha=0.3)
ax3.legend(fontsize=6)

# --------------------------
# PLOT 4: Dataset 2 - Delta Yield vs -log10(p-value) VOLCANO PLOT
# --------------------------
df2_filtered = df2[df2['Catalyst'] != 'WT'].copy()
# Handle log10 safely
with np.errstate(divide='ignore', invalid='ignore'):
    df2_filtered['-log10(p)'] = -np.log10(df2_filtered['p_value'])
    df2_filtered['-log10(p)'] = df2_filtered['-log10(p)'].replace([np.inf, -np.inf], np.nan)

colors_volcano2 = []
for _, row in df2_filtered.iterrows():
    if row['p_value'] < 0.05 and abs(row['delta_yield']) > 3:
        colors_volcano2.append('#ef4444')
    elif row['p_value'] < 0.05:
        colors_volcano2.append('#f59e0b')
    else:
        colors_volcano2.append('#94a3b8')

ax4.scatter(df2_filtered['delta_yield'], df2_filtered['-log10(p)'], 
           c=colors_volcano2, s=100, alpha=0.7, edgecolors='black', linewidth=0.5)

for _, row in df2_filtered.iterrows():
    if row['p_value'] < 0.05:
        ax4.annotate(row['Catalyst'], 
                    (row['delta_yield'], row['-log10(p)']),
                    fontsize=6, ha='center', va='bottom')

ax4.axhline(y=-np.log10(0.05), color='red', linestyle='--', linewidth=1, alpha=0.5, label='p=0.05')
ax4.axvline(x=3, color='blue', linestyle='--', linewidth=1, alpha=0.5, label='ΔYield=±3%')
ax4.axvline(x=-3, color='blue', linestyle='--', linewidth=1, alpha=0.5)

ax4.set_xlabel('ΔYield vs WT (%)', fontsize=8)
ax4.set_ylabel('-log₁₀(p-value)', fontsize=8)
ax4.set_title('Dataset 2: Volcano Plot (Effect Size vs Significance)', fontsize=10, fontweight='bold')
ax4.grid(alpha=0.3)
ax4.legend(fontsize=6)

# --------------------------
# PLOT 5: Dataset 1 - Distribution Overlay (WT vs all mutants)
# --------------------------
# Create histogram data
all_mutant_values1 = []
for _, row in df1[df1['Catalyst'] != 'WT'].iterrows():
    all_mutant_values1.extend(row[['Rep1', 'Rep2', 'Rep3', 'Rep4']].values)

ax5.hist(wt_values1, bins=8, alpha=0.6, color='#8b5cf6', edgecolor='black', 
         linewidth=1, label=f'WT (μ={wt_mean1:.1f}%)', density=True)
ax5.hist(all_mutant_values1, bins=15, alpha=0.6, color='#3b82f6', edgecolor='black', 
         linewidth=1, label=f'All Mutants (μ={np.mean(all_mutant_values1):.1f}%)', density=True)

# Add density curves safely
x_range = np.linspace(10, 55, 100)
if len(wt_values1) > 1 and wt_values1.std() > 0:
    wt_kde = norm.pdf(x_range, wt_values1.mean(), wt_values1.std())
    ax5.plot(x_range, wt_kde, color='#8b5cf6', linewidth=2, linestyle='--')

if len(all_mutant_values1) > 1 and np.std(all_mutant_values1) > 0:
    mutant_kde = norm.pdf(x_range, np.mean(all_mutant_values1), np.std(all_mutant_values1))
    ax5.plot(x_range, mutant_kde, color='#3b82f6', linewidth=2, linestyle='--')

ax5.set_xlabel('Yield (%)', fontsize=8)
ax5.set_ylabel('Density', fontsize=8)
ax5.set_title('Dataset 1: Yield Distribution (WT vs Mutants)', fontsize=10, fontweight='bold')
ax5.legend(fontsize=7)
ax5.grid(alpha=0.3)

# --------------------------
# PLOT 6: Dataset 2 - Distribution Overlay
# --------------------------
all_mutant_values2 = []
for _, row in df2[df2['Catalyst'] != 'WT'].iterrows():
    all_mutant_values2.extend([row['Replicate_a'], row['Replicate_b'], row['Replicate_c']])

ax6.hist(wt_values2, bins=5, alpha=0.6, color='#800080', edgecolor='black', 
         linewidth=1, label=f'WT (μ={wt_mean2:.1f}%)', density=True)
ax6.hist(all_mutant_values2, bins=15, alpha=0.6, color='#3b82f6', edgecolor='black', 
         linewidth=1, label=f'All Mutants (μ={np.mean(all_mutant_values2):.1f}%)', density=True)

x_range2 = np.linspace(2, 22, 100)
if len(wt_values2) > 1 and wt_values2.std() > 0:
    wt_kde2 = norm.pdf(x_range2, wt_values2.mean(), wt_values2.std())
    ax6.plot(x_range2, wt_kde2, color='#800080', linewidth=2, linestyle='--')

if len(all_mutant_values2) > 1 and np.std(all_mutant_values2) > 0:
    mutant_kde2 = norm.pdf(x_range2, np.mean(all_mutant_values2), np.std(all_mutant_values2))
    ax6.plot(x_range2, mutant_kde2, color='#3b82f6', linewidth=2, linestyle='--')

ax6.set_xlabel('Yield (%)', fontsize=8)
ax6.set_ylabel('Density', fontsize=8)
ax6.set_title('Dataset 2: Yield Distribution (WT vs Mutants)', fontsize=10, fontweight='bold')
ax6.legend(fontsize=7)
ax6.grid(alpha=0.3)

# --------------------------
# PLOT 7: Dataset 1 - Violin Plot with Individual Points
# --------------------------
# Prepare data for violin plot - must be list of arrays
violin_data1 = []
violin_labels1 = []

for _, row in df1.iterrows():
    values = row[['Rep1', 'Rep2', 'Rep3', 'Rep4']].values.astype(float)
    violin_data1.append(values)
    violin_labels1.append(row['Catalyst'])

# Create violin plot with proper error handling
try:
    parts1 = ax7.violinplot(violin_data1, positions=range(len(df1)), 
                            widths=0.7, showmeans=True, showextrema=True)
    
    # Color violins
    for i, pc in enumerate(parts1['bodies']):
        pc.set_facecolor(bar_colors1[i])
        pc.set_alpha(0.6)
except Exception as e:
    print(f"Warning: Violin plot failed for Dataset 1: {e}")
    # Fallback to box plot
    ax7.boxplot(violin_data1, positions=range(len(df1)), widths=0.7)

# Overlay individual points
for i, values in enumerate(violin_data1):
    y = values
    x = np.random.normal(i, 0.04, size=len(y))
    ax7.scatter(x, y, alpha=0.6, s=30, color='black', zorder=3)

ax7.set_xticks(range(len(df1)))
ax7.set_xticklabels(df1['Catalyst'], rotation=45, ha='right', fontsize=7)
ax7.set_ylabel('Yield (%)', fontsize=8)
ax7.set_title('Dataset 1: Yield Distribution by Variant (Violin Plot)', fontsize=10, fontweight='bold')
ax7.grid(axis='y', alpha=0.3)

# --------------------------
# PLOT 8: Dataset 2 - Violin Plot with Individual Points
# --------------------------
violin_data2 = []
violin_labels2 = []

for _, row in df2.iterrows():
    values = np.array([row['Replicate_a'], row['Replicate_b'], row['Replicate_c']], dtype=float)
    violin_data2.append(values)
    violin_labels2.append(row['Catalyst'])

try:
    parts2 = ax8.violinplot(violin_data2, positions=range(len(df2)), 
                            widths=0.7, showmeans=True, showextrema=True)
    
    for i, pc in enumerate(parts2['bodies']):
        pc.set_facecolor(bar_colors2[i])
        pc.set_alpha(0.6)
except Exception as e:
    print(f"Warning: Violin plot failed for Dataset 2: {e}")
    # Fallback to box plot
    ax8.boxplot(violin_data2, positions=range(len(df2)), widths=0.7)

for i, values in enumerate(violin_data2):
    y = values
    x = np.random.normal(i, 0.04, size=len(y))
    ax8.scatter(x, y, alpha=0.6, s=30, color='black', zorder=3)

ax8.set_xticks(range(len(df2)))
ax8.set_xticklabels(df2['Catalyst'], rotation=45, ha='right', fontsize=7)
ax8.set_ylabel('Yield (%)', fontsize=8)
ax8.set_title('Dataset 2: Yield Distribution by Variant (Violin Plot)', fontsize=10, fontweight='bold')
ax8.grid(axis='y', alpha=0.3)

plt.tight_layout()

# --------------------------
# SAVE FIGURE
# --------------------------
save_path = "comprehensive_yield_analysis_with_distributions.png"
fig.savefig(save_path, dpi=300, bbox_inches='tight')
print(f"✅ Saved comprehensive figure with distributions to: {save_path}")

# =============================================================================
# ENHANCED STATISTICAL SUMMARY
# =============================================================================
print("\n" + "="*80)
print("ENHANCED STATISTICAL ANALYSIS")
print("="*80)

print("\n--- DATASET 1 ---")
print(f"WT Mean: {wt_mean1:.2f} ± {wt_values1.std():.2f}%")
print(f"All Mutants Mean: {np.mean(all_mutant_values1):.2f} ± {np.std(all_mutant_values1):.2f}%")
print(f"\nSignificant variants (p < 0.05):")
sig_variants1 = df1_filtered[df1_filtered['p_value'] < 0.05]
for _, row in sig_variants1.iterrows():
    print(f"  {row['Catalyst']}: ΔYield = {row['delta_yield']:+.2f}%, p = {row['p_value']:.4f}")

print("\n--- DATASET 2 ---")
print(f"WT Mean: {wt_mean2:.2f} ± {wt_values2.std():.2f}%")
print(f"All Mutants Mean: {np.mean(all_mutant_values2):.2f} ± {np.std(all_mutant_values2):.2f}%")
print(f"\nSignificant variants (p < 0.05):")
sig_variants2 = df2_filtered[df2_filtered['p_value'] < 0.05]
for _, row in sig_variants2.iterrows():
    print(f"  {row['Catalyst']}: ΔYield = {row['delta_yield']:+.2f}%, p = {row['p_value']:.4f}")

print("\n" + "="*80)

plt.show()