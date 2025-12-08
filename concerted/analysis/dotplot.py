import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

# -----------------------------
# 1. INPUT DATA
# -----------------------------
yield_data = pd.DataFrame({
    'Variant': ['W96A','E7A','A11L','L18A','N19A','K22A','N88A','M89A',
                'A92E','F93A','S95A','S97A','D100A','V99A','I16A',
                'N14A','L9A','R10A','K101A','E104A'],
    'DeltaYield': [-12,0,0.5,-0.5,0,-3.5,-3.5,-4.5,-3,-10,-1.5,-0.5,
                   -13.5,0,-6.5,-4,-1.5,-1,-3.5,-3]
})

kcal_data = pd.DataFrame({
    "Variant": ['WT','ARG10A','ASN18A','ASP99A','GLU7A','LEU17A',
                'LYS100A','MET88A','SER94A','TRP95A','ASN14A','ASN87A',
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
        'N19A': 'ASN18A', 'K22A': 'LYS21A', 'N88A': 'ASN87A', 'M89A': 'MET88A',
        'A92E': 'ALA91E', 'F93A': 'PHE92A', 'S95A': 'SER94A', 'S97A': 'SER96A',
        'D100A': 'ASP99A', 'V99A': 'VAL98A', 'I16A': 'ILE15A', 'N14A': 'ASN14A',
        'L9A': 'LEU9A', 'R10A': 'ARG10A', 'K101A': 'LYS100A', 'E104A': 'GLU103A'
    }
    return mapping.get(name, name)

yield_data['Variant_std'] = yield_data['Variant'].apply(standardize_variant)

merged_data = pd.merge(
        yield_data, kcal_data,
        left_on='Variant_std', right_on='Variant',
        suffixes=('_yield','_kcal')
)

# -----------------------------
# 3. CORRELATION (R² only)
# -----------------------------
pearson_r, _ = stats.pearsonr(merged_data['DeltaYield'], merged_data['ΔG*'])
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
    ax.scatter(subset['DeltaYield'], subset['ΔG*'],
               c=subset['color'], s=100, alpha=0.7,
               edgecolors='black', linewidth=1.5, label=category)

# WT point
wt_row = kcal_data[kcal_data['Variant'] == 'WT']
ax.scatter([0], [wt_row['ΔG*'].values[0]],
           s=200, alpha=0.9, color='purple',
           edgecolors='black', linewidth=2,
           marker='o', label='WT')

# variant labels
for idx, row in merged_data.iterrows():
    ax.annotate(row['Variant_yield'],
                (row['DeltaYield'], row['ΔG*']),
                xytext=(5,5), textcoords='offset points',
                fontsize=9)

ax.annotate("WT", (0, wt_row['ΔG*'].values[0]),
            xytext=(5,5), textcoords="offset points",
            fontsize=9, color='purple', weight='bold')

# regression line
z = np.polyfit(merged_data['DeltaYield'], merged_data['ΔG*'], 1)
p_fit = np.poly1d(z)
x_line = np.linspace(merged_data['DeltaYield'].min(), merged_data['DeltaYield'].max(), 100)
ax.plot(x_line, p_fit(x_line), "r--", linewidth=2)

# -----------------------------
# Updated labels, title, and annotation
# -----------------------------
ax.set_xlabel("Experimental ΔYield (%)", fontsize=14, fontweight="bold")
ax.set_ylabel("Computational ΔG* (kcal/mol)", fontsize=14, fontweight="bold")
ax.set_title(f"Experimental ΔYield vs Computational ΔG*  (R² = {r2:.3f})",
             fontsize=16, fontweight="bold", pad=20)

ax.grid(True, alpha=0.3)
ax.legend(fontsize=11)

plt.tight_layout()
plt.savefig("R2_only_plot.png", dpi=300)
plt.show()
