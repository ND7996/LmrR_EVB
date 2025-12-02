import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# =============================================================================
# DATA PREPARATION
# =============================================================================
data1 = {
    'Catalyst': ['WT', 'E7A', 'A11L', 'L18A', 'N19A', 'K22A', 'N88A', 
                 'M89A', 'A92E', 'F93A', 'S95A', 'S97A', 'D100A'],
    'Rep1': [38, 38, 38, 28, 35, 49, 44, 34, 36, 33, 45, 50, 15],
    'Rep2': [37, 40, 40, 33, 36, 51, 39, 32, 37, 34, 38, 49, 12],
    'Rep3': [47, 40, 43, 40, 48, 42, 26, 40, 35, 26, 34, 43, 21],
    'Rep4': [46, 43, 43, 42, 47, 43, 37, 42, 34, 26, 36, 44, 41]
}

data2 = {
    'Catalyst': ['WT', 'W96A', 'E7A', 'A11L', 'L18A', 'N19A', 'K22A', 'N88A', 
                 'M89A', 'A92E', 'F93A', 'S95A', 'S97A', 'D100A', 'V99A', 
                 'I16A', 'N14A', 'L9A', 'R10A', 'K101A', 'E104A'],
    'Replicate_a': [17, 5, 17, 17, 16, 17, 14, 14, 13, 14, 7, 14, 17, 4, 19, 13, 16, 18, 16, 12, 11],
    'Replicate_b': [17, 5, 18, 19, 19, 17, 14, 13, 12, 14, 7, 16, 15, 3, 19, 9, 14, 17, 16, 11, 16],
    'Replicate_c': [17, 5, 16, 16, 17, 14, 13, 15, 13, 14, 7, 16, 18, 3, 13, 9, 10, 12, 11, 18, 15]
}

df1 = pd.DataFrame(data1)
df1['Mean'] = df1[['Rep1','Rep2','Rep3','Rep4']].mean(axis=1)

df2 = pd.DataFrame(data2)
df2['Mean'] = df2[['Replicate_a', 'Replicate_b', 'Replicate_c']].mean(axis=1)

# Get WT values
WT1 = df1[df1['Catalyst']=='WT']['Mean'].values[0]
WT2 = df2[df2['Catalyst']=='WT']['Mean'].values[0]
DELTA_WT = WT1 - WT2

# Calculate metrics for all mutants
df1['Percent_Change'] = ((df1['Mean'] - WT1) / WT1) * 100
df2['Percent_Change'] = ((df2['Mean'] - WT2) / WT2) * 100

# Find common mutants (excluding WT)
common_mutants = sorted(list(set(df1['Catalyst']) & set(df2['Catalyst']) - {'WT'}))

# =============================================================================
# CREATE HEATMAP DATA
# =============================================================================
heatmap_data = []

for mutant in common_mutants:
    raw1 = df1[df1['Catalyst']==mutant]['Mean'].values[0]
    raw2 = df2[df2['Catalyst']==mutant]['Mean'].values[0]
    pc1 = df1[df1['Catalyst']==mutant]['Percent_Change'].values[0]
    pc2 = df2[df2['Catalyst']==mutant]['Percent_Change'].values[0]
    norm_diff = (raw1 - raw2) / DELTA_WT
    
    heatmap_data.append({
        'Mutant': mutant,
        'Dataset 1\nYield (%)': raw1,
        'Dataset 2\nYield (%)': raw2,
        'Change from WT\n(Dataset 1, %)': pc1,
        'Change from WT\n(Dataset 2, %)': pc2,
        'Normalized\nDifference': norm_diff
    })

df_heatmap = pd.DataFrame(heatmap_data)
df_heatmap = df_heatmap.set_index('Mutant')

# Sort by absolute difference from WT behavior (normalized difference from 1.0)
# Sort DESCENDING - most deviation first
df_heatmap['sort_key'] = abs(df_heatmap['Normalized\nDifference'] - 1.0)
df_heatmap = df_heatmap.sort_values('sort_key', ascending=False)
df_heatmap = df_heatmap.drop('sort_key', axis=1)

# =============================================================================
# CREATE FIGURE WITH TWO HEATMAPS
# =============================================================================
fig = plt.figure(figsize=(14, 10))
gs = fig.add_gridspec(2, 1, height_ratios=[1, 1], hspace=0.4)

# Title
fig.suptitle('Mutant Yield Comparison: Dataset 1 vs Dataset 2', 
             fontsize=16, fontweight='bold', y=0.98)

# Add info box
info_text = f'WT Old Dataset (D1) = 42% yield  |  WT New Dataset (D2) = 17% yield  |  Δ = {DELTA_WT:.1f}%'
fig.text(0.5, 0.94, info_text, ha='center', fontsize=11, 
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

# =============================================================================
# TOP HEATMAP: Raw yields and percent changes
# =============================================================================
ax1 = fig.add_subplot(gs[0])

data_top = df_heatmap[['Dataset 1\nYield (%)', 'Dataset 2\nYield (%)', 
                        'Change from WT\n(Dataset 1, %)', 'Change from WT\n(Dataset 2, %)']].T

sns.heatmap(data_top, annot=True, fmt='.1f', cmap='RdYlGn', center=0,
            cbar_kws={'label': 'Value'}, ax=ax1, linewidths=0.5, 
            linecolor='gray', vmin=-100, vmax=100)

ax1.set_title('Raw Yields and Percent Changes from WT', 
              fontweight='bold', fontsize=13, pad=15)
ax1.set_ylabel('Metric', fontweight='bold', fontsize=11)
ax1.set_xlabel('')

# =============================================================================
# BOTTOM HEATMAP: Normalized difference (key metric)
# =============================================================================
ax2 = fig.add_subplot(gs[1])

data_bottom = df_heatmap[['Normalized\nDifference']].T

# Custom colormap: green around 1.0, red for deviations
sns.heatmap(data_bottom, annot=True, fmt='.2f', cmap='RdYlGn', center=1.0,
            cbar_kws={'label': 'Normalized Difference'}, ax=ax2, 
            linewidths=0.5, linecolor='gray', vmin=0, vmax=2)

ax2.set_title('Normalized Difference: (Mutant_D1 - Mutant_D2) / ΔWT\n' + 
              'Green ≈1.0 = Mutant behaves like WT | Red = Different behavior',
              fontweight='bold', fontsize=13, pad=15)
ax2.set_ylabel('Key Metric', fontweight='bold', fontsize=11)
ax2.set_xlabel('Mutants (sorted by deviation from WT behavior: MOST → LEAST)', 
               fontweight='bold', fontsize=11)

plt.tight_layout()
plt.savefig("Mutant_Comparison_Heatmap.png", dpi=300, bbox_inches='tight')
plt.show()

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================
print("\n" + "="*80)
print("MUTANT BEHAVIOR SUMMARY")
print("="*80)
print("\nMutants sorted by deviation from WT behavior:")
print("(MOST deviation → LEAST deviation)\n")

summary = df_heatmap.copy()
summary['Deviation from WT'] = abs(summary['Normalized\nDifference'] - 1.0)
summary = summary.sort_values('Deviation from WT', ascending=False)

for idx, row in summary.iterrows():
    norm_diff = row['Normalized\nDifference']
    deviation = row['Deviation from WT']
    
    if deviation < 0.2:
        behavior = "✓ WT-like"
    elif norm_diff < 0.8:
        behavior = "↓ More affected in Dataset 2"
    else:
        behavior = "↑ More affected in Dataset 1"
    
    print(f"{idx:8s} | Norm.Diff: {norm_diff:5.2f} | Deviation: {deviation:5.2f} | {behavior}")

print("\n" + "="*80)