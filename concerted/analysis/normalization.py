import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

# =============================================================================
# YOUR ACTUAL DATA
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

print("="*80)
print("NORMALIZATION EXPLANATION")
print("="*80)
print(f"\nDataset 1 WT: {WT1:.2f}%")
print(f"Dataset 2 WT: {WT2:.2f}%")
print(f"ΔWT (WT1 - WT2): {DELTA_WT:.2f}%")
print("\nThe baseline is ~25% higher in Dataset 1")
print("We need to account for this when comparing mutants!")
print("="*80)

# =============================================================================
# METHOD 1: SIMPLE NORMALIZATION - Percent Change from WT
# =============================================================================
df1['Percent_Change'] = ((df1['Mean'] - WT1) / WT1) * 100
df2['Percent_Change'] = ((df2['Mean'] - WT2) / WT2) * 100

# Find common mutants
common_mutants = sorted(list(set(df1['Catalyst']) & set(df2['Catalyst']) - {'WT'}))

# Create comparison
comparison = []
for mutant in common_mutants:
    pc1 = df1[df1['Catalyst']==mutant]['Percent_Change'].values[0]
    pc2 = df2[df2['Catalyst']==mutant]['Percent_Change'].values[0]
    raw1 = df1[df1['Catalyst']==mutant]['Mean'].values[0]
    raw2 = df2[df2['Catalyst']==mutant]['Mean'].values[0]
    
    comparison.append({
        'Mutant': mutant,
        'Raw_D1': raw1,
        'Raw_D2': raw2,
        'PercentChange_D1': pc1,
        'PercentChange_D2': pc2,
        'Agreement': abs(pc1 - pc2)
    })

df_comp = pd.DataFrame(comparison)

# =============================================================================
# METHOD 2: NORMALIZED DIFFERENCE
# =============================================================================
for mutant in common_mutants:
    raw1 = df1[df1['Catalyst']==mutant]['Mean'].values[0]
    raw2 = df2[df2['Catalyst']==mutant]['Mean'].values[0]
    normalized_diff = (raw1 - raw2) / DELTA_WT
    df_comp.loc[df_comp['Mutant']==mutant, 'Normalized_Diff'] = normalized_diff

# =============================================================================
# VISUALIZATION — ONLY RIGHT PLOT
# =============================================================================
fig, ax2 = plt.subplots(figsize=(10, 6))
fig.suptitle('Yield Normalization: Normalized Difference', fontsize=14, fontweight='bold')

x_pos = np.arange(len(common_mutants))

colors = ['#10b981' if abs(n-1.0) < 0.2 else '#f59e0b' if n < 1.0 else '#ef4444' 
          for n in df_comp['Normalized_Diff']]

ax2.bar(x_pos, df_comp['Normalized_Diff'], color=colors, alpha=0.85, 
        edgecolor='black', linewidth=1)

ax2.axhline(y=1.0, color='black', linestyle='-', linewidth=2,
            label='Expected if mutant = WT behavior')
ax2.axhline(y=0.8, color='gray', linestyle='--', linewidth=1, alpha=0.5)
ax2.axhline(y=1.2, color='gray', linestyle='--', linewidth=1, alpha=0.5)

ax2.set_xlabel('Mutants', fontweight='bold', fontsize=12)
ax2.set_ylabel('(Mutant_D1 - Mutant_D2) / ΔWT', fontweight='bold', fontsize=12)
ax2.set_title('METHOD 2: Normalized Difference\n(Does mutant match WT difference?)',
              fontweight='bold', fontsize=11)

ax2.set_xticks(x_pos)
ax2.set_xticklabels(common_mutants, rotation=45, ha='right')
ax2.grid(axis='y', alpha=0.3)

from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='#10b981', label='≈WT behavior (0.8–1.2)'),
    Patch(facecolor='#f59e0b', label='More affected in D2 (<0.8)'),
    Patch(facecolor='#ef4444', label='More affected in D1 (>1.2)')
]
ax2.legend(handles=legend_elements, fontsize=9, loc='upper right')

plt.tight_layout()
plt.savefig("Normalization_Method2_Only.png", dpi=300, bbox_inches='tight')
plt.show()
