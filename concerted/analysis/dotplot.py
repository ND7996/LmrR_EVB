import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

# -----------------------------
# 1. INPUT DATA
# -----------------------------
# Raw replicate data - FIXED: Added L17A as the 21st entry
data2 = {
    'Catalyst': ['WT', 'W96A', 'E7A', 'A11L', 'L18A', 'N19A', 'K22A', 'N88A', 
                 'M89A', 'A92E', 'F93A', 'S95A', 'S97A', 'D100A', 'V99A', 
                 'I16A', 'L17A', 'L9A', 'R10A', 'K101A', 'E104A'],
    'Replicate_a': [17, 5, 17, 17, 16, 17, 14, 14, 13, 14, 7, 14, 17, 4, 19, 13, 16, 18, 16, 12, 11],
    'Replicate_b': [17, 5, 18, 19, 19, 17, 14, 13, 12, 14, 7, 16, 15, 3, 19, 9, 14, 17, 16, 11, 16],
    'Replicate_c': [17, 5, 16, 16, 17, 14, 13, 15, 13, 14, 7, 16, 18, 3, 13, 9, 10, 12, 11, 18, 15]
}

yield_data = pd.DataFrame(data2)
# Calculate average yield across replicates
yield_data['AvgYield'] = yield_data[['Replicate_a', 'Replicate_b', 'Replicate_c']].mean(axis=1)

kcal_data = pd.DataFrame({
    "Variant": ['WT','ARG10A','ASN18A','ASP99A','GLU7A','LEU17A',
                'LYS100A','MET88A','SER94A','TRP95A','ASN19A','ASN87A',
                'ILE15A','GLU103A','LEU9A','PHE92A','LYS21A','SER96A',
                'VAL98A','ALA91E','ALA11L'],
    "ΔG*": [20.800,19.329,19.385,20.328,19.329,18.328,
            19.328,20.529,18.330,20.331,21.329,18.328,
            21.385,18.880,21.127,22.329,19.828,21.330,
            20.385,19.349,18.125]
})

# ΔΔG calculation (not used for R² plot, but kept)
wt_dg = kcal_data[kcal_data['Variant'] == 'WT']['ΔG*'].values[0]
kcal_data['ΔΔG'] = kcal_data['ΔG*'] - wt_dg

# -----------------------------
# 2. STANDARDIZE VARIANT NAMES
# -----------------------------
def standardize_variant(name):
    mapping = {
        'W96A': 'TRP95A', 'E7A': 'GLU7A', 'A11L': 'ALA11L', 'L18A': 'ASN18A',
        'N19A': 'ASN19A', 'K22A': 'LYS21A', 'N88A': 'ASN87A', 'M89A': 'MET88A',
        'A92E': 'ALA91E', 'F93A': 'PHE92A', 'S95A': 'SER94A', 'S97A': 'SER96A',
        'D100A': 'ASP99A', 'V99A': 'VAL98A', 'I16A': 'ILE15A', 'L17A': 'LEU17A',
        'L9A': 'LEU9A', 'R10A': 'ARG10A', 'K101A': 'LYS100A', 'E104A': 'GLU103A'
    }
    return mapping.get(name, name)

yield_data['Variant_std'] = yield_data['Catalyst'].apply(standardize_variant)

merged_data = pd.merge(
        yield_data, kcal_data,
        left_on='Variant_std', right_on='Variant',
        suffixes=('_yield','_kcal')
)

# -----------------------------
# 3. CORRELATION (R² only)
# -----------------------------
pearson_r, _ = stats.pearsonr(merged_data['AvgYield'], merged_data['ΔG*'])
r2 = pearson_r**2

# -----------------------------
# 4. COLOR CATEGORIES
# -----------------------------
def categorize_dg(dg_value):
    if dg_value < 19.0:
        return 'Low ΔG*', '#d73027'
    elif dg_value < 20.0:
        return 'Medium ΔG*', '#fc8d59'
    elif dg_value < 21.0:
        return 'High ΔG*', '#91bfdb'
    else:
        return 'Very High ΔG*', '#4575b4'

merged_data['category'], merged_data['color'] = zip(*merged_data['ΔG*'].apply(categorize_dg))

# -----------------------------
# 5. PLOT (R² ONLY)
# -----------------------------
fig, ax = plt.subplots(figsize=(10, 8))

# scatter groups
for category in merged_data['category'].unique():
    subset = merged_data[merged_data['category'] == category]
    ax.scatter(subset['AvgYield'], subset['ΔG*'],
               c=subset['color'], s=100, alpha=0.7,
               edgecolors='black', linewidth=1.5, label=category)

# WT point
wt_row = merged_data[merged_data['Catalyst'] == 'WT']
ax.scatter(wt_row['AvgYield'], wt_row['ΔG*'],
           s=200, alpha=0.9, color='purple',
           edgecolors='black', linewidth=2,
           marker='o', label='WT')

# variant labels
for idx, row in merged_data.iterrows():
    label = row['Catalyst']
    ax.annotate(label,
                (row['AvgYield'], row['ΔG*']),
                xytext=(5,5), textcoords='offset points',
                fontsize=9,
                color='purple' if label == 'WT' else 'black',
                weight='bold' if label == 'WT' else 'normal')

# regression line
z = np.polyfit(merged_data['AvgYield'], merged_data['ΔG*'], 1)
p_fit = np.poly1d(z)
x_line = np.linspace(merged_data['AvgYield'].min(), merged_data['AvgYield'].max(), 100)
ax.plot(x_line, p_fit(x_line), "r--", linewidth=2)

# -----------------------------
# Updated labels, title, and annotation
# -----------------------------
ax.set_xlabel("Experimental Average Yield (%)", fontsize=14, fontweight="bold")
ax.set_ylabel("Computational ΔG* (kcal/mol)", fontsize=14, fontweight="bold")
ax.set_title(f"Experimental Average Yield vs Computational ΔG*  (R² = {r2:.3f})",
             fontsize=16, fontweight="bold", pad=20)

ax.grid(True, alpha=0.3)
ax.legend(fontsize=11)

plt.tight_layout()
plt.savefig("avg_yield_vs_dg_plot.png", dpi=300)
plt.show()

print(f"\nTotal variants analyzed: {len(merged_data)}")
print(f"R² = {r2:.3f}")