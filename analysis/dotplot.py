import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import seaborn as sns

# Data from yield measurements
yield_data = pd.DataFrame({
    "Variant": ['W96A', 'E7A', 'A11L', 'L18A', 'N19A', 'K22A',
                'N88A', 'M89A', 'A92E', 'F93A', 'S95A', 'S97A',
                'D100A', 'V99A', 'I16A', 'N14A', 'L9A', 'R10A',
                'K101A', 'E104A'],
    "ΔYield":  [-12.0, 0.0, 0.5, -0.5, 0.0, -3.5,
                -3.5, -4.5, -10.0, -3.0, -1.5, -0.5,
                -13.5, 0.0, -6.5, -4.0, -1.5, -1.0,
                -3.5, -3.0]
})

# Data from kcal measurements (ΔΔG)
kcal_data = pd.DataFrame({
    "Variant": ['WT', 'ARG10A', 'ASN18A', 'ASP99A', 'GLU7A', 'LEU17A', 
                'LYS100A', 'MET88A', 'SER94A', 'TRP95A', 'ASN14A', 'ASN87A',
                'ILE15A', 'GLU103A', 'LEU9A', 'PHE92A', 'LYS21A', 'SER96A',
                'VAL98A', 'ALA91E', 'ALA11L'],
    "ΔΔG": [0.000, -1.471, -1.415, -0.472, -1.471, -2.472, 
            -1.472, -0.271, -2.470, -0.469, 0.529, -2.472,
            0.585, -1.920, 0.327, 1.529, -0.972, 0.530,
            -0.415, -1.451, -2.675]
})

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
pearson_r, p_value = stats.pearsonr(merged_data['ΔYield'], merged_data['ΔΔG'])
spearman_r, spearman_p = stats.spearmanr(merged_data['ΔYield'], merged_data['ΔΔG'])

# Identify outliers using residuals
z = np.polyfit(merged_data['ΔYield'], merged_data['ΔΔG'], 1)
p_fit = np.poly1d(z)
merged_data['predicted_ddG'] = p_fit(merged_data['ΔYield'])
merged_data['residuals'] = merged_data['ΔΔG'] - merged_data['predicted_ddG']
merged_data['abs_residuals'] = np.abs(merged_data['residuals'])

# Z-score for outlier detection
merged_data['z_score'] = np.abs(stats.zscore(merged_data['residuals']))
outliers = merged_data[merged_data['z_score'] > 1.5]

# Create comprehensive figure
fig = plt.figure(figsize=(16, 12))
gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)

# 1. Main scatter plot with outliers highlighted
ax1 = fig.add_subplot(gs[0:2, 0:2])
normal_points = merged_data[merged_data['z_score'] <= 1.5]
outlier_points = merged_data[merged_data['z_score'] > 1.5]

ax1.scatter(normal_points['ΔYield'], normal_points['ΔΔG'], 
           s=100, alpha=0.6, color='steelblue', edgecolors='black', 
           linewidth=1, label='Normal points')
if len(outlier_points) > 0:
    ax1.scatter(outlier_points['ΔYield'], outlier_points['ΔΔG'], 
               s=150, alpha=0.8, color='red', edgecolors='darkred', 
               linewidth=2, marker='^', label='Potential outliers')

# Add labels for ALL points
for idx, row in merged_data.iterrows():
    color = 'red' if row['z_score'] > 1.5 else 'black'
    weight = 'bold' if row['z_score'] > 1.5 else 'normal'
    ax1.annotate(row['Variant_yield'], 
                (row['ΔYield'], row['ΔΔG']),
                xytext=(5, 5), textcoords='offset points',
                fontsize=7, alpha=0.75, color=color, weight=weight)

# Regression line
x_line = np.linspace(merged_data['ΔYield'].min(), merged_data['ΔYield'].max(), 100)
ax1.plot(x_line, p_fit(x_line), "r--", alpha=0.8, linewidth=2, label='Linear fit')

ax1.set_xlabel('ΔYield (%)', fontsize=12, fontweight='bold')
ax1.set_ylabel('ΔΔG (kcal/mol)', fontsize=12, fontweight='bold')
ax1.set_title(f'ΔYield vs ΔΔG\nPearson r = {pearson_r:.3f} (R² = {pearson_r**2:.3f}), p = {p_value:.4f}', 
             fontsize=13, fontweight='bold')
ax1.grid(True, alpha=0.3, linestyle='--')
ax1.legend(fontsize=9)

# 2. Residuals plot
ax2 = fig.add_subplot(gs[0, 2])
ax2.scatter(merged_data['ΔYield'], merged_data['residuals'], 
           c=merged_data['z_score'], cmap='RdYlBu_r', s=80, edgecolors='black', linewidth=0.5)
ax2.axhline(y=0, color='black', linestyle='--', linewidth=2)
ax2.axhline(y=merged_data['residuals'].std(), color='red', linestyle=':', alpha=0.5, label='±1 SD')
ax2.axhline(y=-merged_data['residuals'].std(), color='red', linestyle=':', alpha=0.5)
ax2.set_xlabel('ΔYield (%)', fontsize=10, fontweight='bold')
ax2.set_ylabel('Residuals', fontsize=10, fontweight='bold')
ax2.set_title('Residuals vs ΔYield', fontsize=11, fontweight='bold')
ax2.grid(True, alpha=0.3)
ax2.legend(fontsize=8)

# 3. Histogram of residuals
ax3 = fig.add_subplot(gs[1, 2])
ax3.hist(merged_data['residuals'], bins=10, color='steelblue', 
         edgecolor='black', alpha=0.7)
ax3.axvline(x=0, color='red', linestyle='--', linewidth=2)
ax3.set_xlabel('Residuals', fontsize=10, fontweight='bold')
ax3.set_ylabel('Frequency', fontsize=10, fontweight='bold')
ax3.set_title('Distribution of Residuals', fontsize=11, fontweight='bold')
ax3.grid(True, alpha=0.3, axis='y')

# 4. ΔYield distribution
ax4 = fig.add_subplot(gs[2, 0])
ax4.hist(merged_data['ΔYield'], bins=12, color='coral', 
         edgecolor='black', alpha=0.7)
ax4.axvline(x=merged_data['ΔYield'].mean(), color='red', 
           linestyle='--', linewidth=2, label='Mean')
ax4.set_xlabel('ΔYield (%)', fontsize=10, fontweight='bold')
ax4.set_ylabel('Frequency', fontsize=10, fontweight='bold')
ax4.set_title('Distribution of ΔYield', fontsize=11, fontweight='bold')
ax4.legend(fontsize=8)
ax4.grid(True, alpha=0.3, axis='y')

# 5. ΔΔG distribution
ax5 = fig.add_subplot(gs[2, 1])
ax5.hist(merged_data['ΔΔG'], bins=12, color='lightgreen', 
         edgecolor='black', alpha=0.7)
ax5.axvline(x=merged_data['ΔΔG'].mean(), color='red', 
           linestyle='--', linewidth=2, label='Mean')
ax5.set_xlabel('ΔΔG (kcal/mol)', fontsize=10, fontweight='bold')
ax5.set_ylabel('Frequency', fontsize=10, fontweight='bold')
ax5.set_title('Distribution of ΔΔG', fontsize=11, fontweight='bold')
ax5.legend(fontsize=8)
ax5.grid(True, alpha=0.3, axis='y')

# 6. Statistics table
ax6 = fig.add_subplot(gs[2, 2])
ax6.axis('tight')
ax6.axis('off')

stats_data = [
    ['Metric', 'Value'],
    ['─'*20, '─'*15],
    ['N variants', f'{len(merged_data)}'],
    ['Pearson r', f'{pearson_r:.4f}'],
    ['R²', f'{pearson_r**2:.4f}'],
    ['P-value', f'{p_value:.6f}'],
    ['Spearman ρ', f'{spearman_r:.4f}'],
    ['Spearman p', f'{spearman_p:.6f}'],
    ['RMSE', f'{np.sqrt(np.mean(merged_data["residuals"]**2)):.3f}'],
    ['Outliers', f'{len(outliers)}'],
]

table = ax6.table(cellText=stats_data, cellLoc='left', loc='center',
                  colWidths=[0.6, 0.4])
table.auto_set_font_size(False)
table.set_fontsize(9)
table.scale(1, 2)

# Style header row
for i in range(2):
    table[(0, i)].set_facecolor('#4472C4')
    table[(0, i)].set_text_props(weight='bold', color='white')

plt.suptitle('Comprehensive Analysis: Protein Variant ΔYield vs ΔΔG', 
             fontsize=16, fontweight='bold', y=0.995)

plt.tight_layout()

# Save the figure
plt.savefig('protein_variant_analysis.png', dpi=300, bbox_inches='tight')
plt.savefig('protein_variant_analysis.pdf', bbox_inches='tight')
print("\n✓ Figures saved as:")
print("  - protein_variant_analysis.png (high resolution)")
print("  - protein_variant_analysis.pdf (vector format)")

plt.show()

# Print detailed analysis
print("="*80)
print("COMPREHENSIVE CORRELATION ANALYSIS")
print("="*80)
print(f"\n1. CORRELATION STATISTICS:")
print(f"   Pearson correlation (r): {pearson_r:.4f}")
print(f"   R-squared (R²): {pearson_r**2:.4f} ({pearson_r**2*100:.1f}%)")
print(f"   P-value: {p_value:.6f}")
print(f"   Spearman correlation (ρ): {spearman_r:.4f} (rank-based, robust to outliers)")
print(f"   Spearman p-value: {spearman_p:.6f}")

print(f"\n2. MODEL FIT QUALITY:")
print(f"   RMSE (Root Mean Square Error): {np.sqrt(np.mean(merged_data['residuals']**2)):.3f} kcal/mol")
print(f"   Mean Absolute Error: {merged_data['abs_residuals'].mean():.3f} kcal/mol")
print(f"   Regression equation: ΔΔG = {z[0]:.4f} × ΔYield + {z[1]:.4f}")

print(f"\n3. DATA DISTRIBUTION:")
print(f"   ΔYield range: [{merged_data['ΔYield'].min():.1f}, {merged_data['ΔYield'].max():.1f}] %")
print(f"   ΔYield mean ± SD: {merged_data['ΔYield'].mean():.2f} ± {merged_data['ΔYield'].std():.2f} %")
print(f"   ΔΔG range: [{merged_data['ΔΔG'].min():.3f}, {merged_data['ΔΔG'].max():.3f}] kcal/mol")
print(f"   ΔΔG mean ± SD: {merged_data['ΔΔG'].mean():.3f} ± {merged_data['ΔΔG'].std():.3f} kcal/mol")

print(f"\n4. OUTLIER ANALYSIS (|z-score| > 1.5):")
if len(outliers) > 0:
    print(f"   Number of outliers detected: {len(outliers)}")
    print(f"\n   {'Variant':<12} {'ΔYield':<10} {'ΔΔG':<10} {'Residual':<12} {'Z-score':<10}")
    print(f"   {'-'*60}")
    for _, row in outliers.sort_values('abs_residuals', ascending=False).iterrows():
        print(f"   {row['Variant_yield']:<12} {row['ΔYield']:>8.1f} {row['ΔΔG']:>9.3f} {row['residuals']:>11.3f} {row['z_score']:>9.2f}")
else:
    print("   No major outliers detected")

print(f"\n5. INTERPRETATION OF R² = {pearson_r**2:.3f} ({pearson_r**2*100:.1f}%):")
print(f"   - {pearson_r**2*100:.1f}% of variance in ΔΔG is explained by ΔYield")
print(f"   - {(1-pearson_r**2)*100:.1f}% of variance is due to other factors")
print("\n   Possible reasons for low R²:")
print("   • Protein stability (ΔΔG) and expression yield are largely independent")
print("   • Yield depends on multiple factors: folding kinetics, aggregation, degradation")
print("   • ΔΔG measures thermodynamic stability, not folding efficiency")
print("   • Other factors: expression level, chaperone interactions, post-translational modifications")
print("   • Measurement uncertainty in both ΔYield and ΔΔG")

print(f"\n6. TOP 5 MOST DELETERIOUS VARIANTS (by ΔYield):")
top_deleterious = merged_data.nsmallest(5, 'ΔYield')[['Variant_yield', 'ΔYield', 'ΔΔG']]
print(top_deleterious.to_string(index=False))

print(f"\n7. TOP 5 MOST DESTABILIZING VARIANTS (by ΔΔG):")
top_stabilizing = merged_data.nsmallest(5, 'ΔΔG')[['Variant_yield', 'ΔYield', 'ΔΔG']]
print(top_stabilizing.to_string(index=False))

print("\n" + "="*80)
print("\nFull dataset with predictions and residuals:")
print("="*80)
# FIXED: Include abs_residuals in display_cols OR sort by a column that IS in display_cols
display_cols = ['Variant_yield', 'ΔYield', 'ΔΔG', 'predicted_ddG', 'residuals', 'abs_residuals']
print(merged_data[display_cols].sort_values('abs_residuals', ascending=False).to_string(index=False))
print("="*80)