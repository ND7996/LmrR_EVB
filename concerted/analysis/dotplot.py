import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

# Data from yield measurements
yield_data = pd.DataFrame({
    'Variant': ['W96A','E7A','A11L','L18A','N19A','K22A','N88A','M89A','A92E','F93A','S95A','S97A','D100A','V99A','I16A','N14A','L9A','R10A','K101A','E104A'],
    'DeltaYield': [-12,0,0.5,-0.5,0,-3.5,-3.5,-4.5,-3,-10,-1.5,-0.5,-13.5,0,-6.5,-4,-1.5,-1,-3.5,-3]
})

# Data from kcal measurements (ΔG*)
kcal_data = pd.DataFrame({
    "Variant": ['WT', 'ARG10A', 'ASN18A', 'ASP99A', 'GLU7A', 'LEU17A', 
                'LYS100A', 'MET88A', 'SER94A', 'TRP95A', 'ASN14A', 'ASN87A',
                'ILE15A', 'GLU103A', 'LEU9A', 'PHE92A', 'LYS21A', 'SER96A',
                'VAL98A', 'ALA91E', 'ALA11L'],
    "ΔG*": [20.800, 19.329, 19.385, 20.328, 19.329, 18.328, 
            19.328, 20.529, 18.330, 20.331, 21.329, 18.328,
            21.385, 18.880, 21.127, 22.329, 19.828, 21.330,
            20.385, 19.349, 18.125]
})

# Calculate ΔΔG = ΔG*(variant) - ΔG*(WT)
wt_dg = kcal_data[kcal_data['Variant'] == 'WT']['ΔG*'].values[0]
kcal_data['ΔΔG'] = kcal_data['ΔG*'] - wt_dg

# Create a mapping function to standardize variant names
def standardize_variant(name):
    """Convert variant names to standard format"""
    mapping = {
        'W96A': 'TRP95A', 'E7A': 'GLU7A', 'A11L': 'ALA11L', 'L18A': 'ASN18A',
        'N19A': 'ASN18A', 'K22A': 'LYS21A', 'N88A': 'ASN87A', 'M89A': 'MET88A',
        'A92E': 'ALA91E', 'F93A': 'PHE92A', 'S95A': 'SER94A', 'S97A': 'SER96A',
        'D100A': 'ASP99A', 'V99A': 'VAL98A', 'I16A': 'ILE15A', 'N14A': 'ASN14A',
        'L9A': 'LEU9A', 'R10A': 'ARG10A', 'K101A': 'LYS100A', 'E104A': 'GLU103A'
    }
    return mapping.get(name, name)

# Standardize variant names in yield data
yield_data['Variant_std'] = yield_data['Variant'].apply(standardize_variant)

# Merge datasets
merged_data = pd.merge(yield_data, kcal_data, 
                       left_on='Variant_std', right_on='Variant', 
                       suffixes=('_yield', '_kcal'))

# Calculate statistics
pearson_r, p_value = stats.pearsonr(merged_data['DeltaYield'], merged_data['ΔG*'])
spearman_r, spearman_p = stats.spearmanr(merged_data['DeltaYield'], merged_data['ΔG*'])

# Create figure with single scatter plot
fig, ax = plt.subplots(figsize=(10, 8))

# Define discrete color categories based on ΔG* values
def categorize_dg(dg_value):
    if dg_value < 19.0:
        return 'Low ΔG* (<19)', '#d73027'  # Dark red
    elif dg_value < 20.0:
        return 'Medium ΔG* (19-20)', '#fc8d59'  # Orange
    elif dg_value < 21.0:
        return 'High ΔG* (20-21)', '#91bfdb'  # Light blue
    else:
        return 'Very High ΔG* (>21)', '#4575b4'  # Dark blue

# Add categories to dataframe
merged_data = merged_data.copy()
merged_data['category'], merged_data['color'] = zip(*merged_data['ΔG*'].apply(categorize_dg))

# Plot all points as circles with discrete colors
for category in merged_data['category'].unique():
    subset = merged_data[merged_data['category'] == category]
    ax.scatter(subset['DeltaYield'], subset['ΔG*'], 
               c=subset['color'], s=100, alpha=0.7, 
               edgecolors='black', linewidth=1.5, 
               label=category)

# Add WT point in purple (WT has DeltaYield = 0 by definition)
wt_row = kcal_data[kcal_data['Variant'] == 'WT']
ax.scatter([0], [wt_row['ΔG*'].values[0]], 
           s=200, alpha=0.9, color='purple', edgecolors='black', 
           linewidth=2, marker='o', label='WT', zorder=5)

# Add labels for ALL points
for idx, row in merged_data.iterrows():
    ax.annotate(row['Variant_yield'], 
                (row['DeltaYield'], row['ΔG*']),
                xytext=(5, 5), textcoords='offset points',
                fontsize=9, alpha=0.75, color='black', weight='normal')

# Add WT label
ax.annotate('WT', (0, wt_row['ΔG*'].values[0]),
            xytext=(5, 5), textcoords='offset points',
            fontsize=9, alpha=0.75, color='purple', weight='bold')

# Regression line
z = np.polyfit(merged_data['DeltaYield'], merged_data['ΔG*'], 1)
p_fit = np.poly1d(z)
x_line = np.linspace(merged_data['DeltaYield'].min(), merged_data['DeltaYield'].max(), 100)
ax.plot(x_line, p_fit(x_line), "r--", alpha=0.8, linewidth=2, label='Linear fit')

ax.set_xlabel('ΔYield (%)', fontsize=14, fontweight='bold')
ax.set_ylabel('ΔG* (kcal/mol)', fontsize=14, fontweight='bold')
ax.set_title(f'ΔYield vs ΔG*\nPearson r = {pearson_r:.3f} (R² = {pearson_r**2:.3f}), p = {p_value:.4f}', 
             fontsize=15, fontweight='bold', pad=20)
ax.grid(True, alpha=0.3, linestyle='--')
ax.legend(fontsize=11, loc='best')

plt.tight_layout()

# Save the figure
plt.savefig('protein_variant_scatter.png', dpi=300, bbox_inches='tight')
plt.savefig('protein_variant_scatter.pdf', bbox_inches='tight')
print("\n✓ Figures saved as:")
print("  - protein_variant_scatter.png (high resolution)")
print("  - protein_variant_scatter.pdf (vector format)")

plt.show()

# Print detailed analysis
print("="*80)
print("CORRELATION ANALYSIS: ΔYield vs ΔG*")
print("="*80)
print(f"\n1. CORRELATION STATISTICS:")
print(f"   Pearson correlation (r): {pearson_r:.4f}")
print(f"   R-squared (R²): {pearson_r**2:.4f} ({pearson_r**2*100:.1f}%)")
print(f"   P-value: {p_value:.6f}")
print(f"   Spearman correlation (ρ): {spearman_r:.4f}")
print(f"   Spearman p-value: {spearman_p:.6f}")

print(f"\n2. INTERPRETATION:")
print(f"   R² = {pearson_r**2:.3f} means {pearson_r**2*100:.1f}% of variance in ΔG* is explained by ΔYield")
print(f"   This is typical for biological systems where multiple factors contribute to outcomes")

print(f"\n3. DATA SUMMARY:")
print(f"   Number of variants analyzed: {len(merged_data)}")
print(f"   ΔYield range: [{merged_data['DeltaYield'].min():.1f}, {merged_data['DeltaYield'].max():.1f}] %")
print(f"   ΔG* range: [{merged_data['ΔG*'].min():.3f}, {merged_data['ΔG*'].max():.3f}] kcal/mol")
print(f"   WT ΔG*: {wt_dg:.3f} kcal/mol")

print("="*80)