import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import pearsonr, spearmanr, shapiro
import matplotlib.pyplot as plt
import re

# =============================================================================
# DATASET 1 — Previous
# =============================================================================
previous_data = {
    'Variant': ['L18A', 'K22A', 'F93A', 'S95A', 'S97A',
                'E7A', 'A11L', 'N19A', 'N88A', 'N89A', 'A92E', 'D100A'],
    'DeltaYield': [-22, -25, -30, -18, -20, 0, 5, 2, 0, 0, -5, 0]
}
df_previous = pd.DataFrame(previous_data)

# =============================================================================
# DATASET 2 — Current
# =============================================================================
current_data = {
    'Variant': ['W96A', 'E7A', 'A11L', 'L18A', 'N19A', 'K22A',
                'N88A', 'M89A', 'A92E', 'F93A', 'S95A', 'S97A',
                'D100A', 'V99A', 'I16A', 'N14A', 'L9A', 'R10A',
                'K101A', 'E104A'],
    'DeltaYield': [-12.0, 0.0, 0.5, -0.5, 0.0, -3.5,
                   -3.5, -4.5, -10.0, -3.0, -1.5, -0.5,
                   -13.5, 0.0, -6.5, -4.0, -1.5, -1.0,
                   -3.5, -3.0]
}
df_current = pd.DataFrame(current_data)

# =============================================================================
# HELPER — Extract residue position
# =============================================================================
def extract_position(variant):
    match = re.search(r'(\d+)', variant)
    return int(match.group(1)) if match else np.inf

df_previous['Position'] = df_previous['Variant'].apply(extract_position)
df_current['Position'] = df_current['Variant'].apply(extract_position)
df_previous = df_previous.sort_values('Position')
df_current = df_current.sort_values('Position')

# =============================================================================
# MERGE
# =============================================================================
df_combined = pd.concat([df_previous, df_current], ignore_index=True)
df_merged = pd.merge(df_previous[['Variant', 'DeltaYield']],
                     df_current[['Variant', 'DeltaYield']],
                     on='Variant', how='outer', suffixes=('_Old', '_New'))
df_merged['Position'] = df_merged['Variant'].apply(extract_position)
df_merged = df_merged.sort_values('Position')

# =============================================================================
# REGRESSION FUNCTION
# =============================================================================
def perform_regression(df):
    x = np.arange(1, len(df) + 1)
    y = df['DeltaYield'].values
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    y_pred = slope * x + intercept
    return {
        'x': x, 'y': y, 'y_pred': y_pred,
        'r': r_value, 'r2': r_value**2,
        'p': p_value, 'slope': slope, 'intercept': intercept
    }

results_prev = perform_regression(df_previous)
results_curr = perform_regression(df_current)
results_comb = perform_regression(df_combined)

# =============================================================================
# CORRELATION & EXTRA ANALYSES
# =============================================================================
# 1. Correlation old vs new (shared variants)
shared = df_merged.dropna(subset=['DeltaYield_Old', 'DeltaYield_New'])
r_pearson, p_pearson = pearsonr(shared['DeltaYield_Old'], shared['DeltaYield_New'])
r_spear, p_spear = spearmanr(shared['DeltaYield_Old'], shared['DeltaYield_New'])

# 2. ΔΔYield
df_merged['ΔΔYield'] = df_merged['DeltaYield_New'] - df_merged['DeltaYield_Old']

# 3. Z-score outlier detection
df_combined['Zscore'] = stats.zscore(df_combined['DeltaYield'])
outliers = df_combined[np.abs(df_combined['Zscore']) > 2]

# 4. Position correlation
r_pos, p_pos = pearsonr(df_combined['Position'], df_combined['DeltaYield'])

# 5. Hydrophobicity-based grouping
hydrophobic = ['A','V','L','I','M','F','W','Y']
df_combined['Original'] = df_combined['Variant'].str[0]
df_combined['Hydrophobicity'] = np.where(df_combined['Original'].isin(hydrophobic), 'Hydrophobic', 'Polar')
hydro_mean = df_combined.groupby('Hydrophobicity')['DeltaYield'].mean()

# 6. Normality test
stat, p_norm = shapiro(df_combined['DeltaYield'])

# =============================================================================
# PRINT RESULTS
# =============================================================================
print("="*90)
print("ALANINE SCANNING ΔYIELD — EXTENDED ANALYSIS")
print("="*90)
print(f"Correlation (Old vs New): Pearson r={r_pearson:.3f} (p={p_pearson:.3e}), Spearman r={r_spear:.3f}")
print(f"Correlation (Position vs ΔYield): r={r_pos:.3f} (p={p_pos:.3e})")
print(f"Hydrophobic mean ΔYield:\n{hydro_mean}")
print(f"Normality test (Shapiro): stat={stat:.3f}, p={p_norm:.3f}")
print(f"Outliers (|Z|>2):\n{outliers[['Variant','DeltaYield','Zscore']]}")
print("="*90)

# =============================================================================
# PLOTS
# =============================================================================
fig, axes = plt.subplots(3, 2, figsize=(12, 10))
fig.suptitle('Alanine Scanning ΔYield — Comparative Analysis', fontsize=16, fontweight='bold', y=1.03)

# (1) Previous Regression
axes[0,0].scatter(results_prev['x'], results_prev['y'], c='green', label='Data')
axes[0,0].plot(results_prev['x'], results_prev['y_pred'], 'r--', label='Regression')
axes[0,0].set_xticks(results_prev['x'])
axes[0,0].set_xticklabels(df_previous['Variant'], rotation=70, fontsize=8)
axes[0,0].set_ylabel('ΔYield (%)')
axes[0,0].set_title('Previous Dataset')
axes[0,0].legend(); axes[0,0].grid(alpha=0.3)

# (2) Current Regression
axes[0,1].scatter(results_curr['x'], results_curr['y'], c='blue', label='Data')
axes[0,1].plot(results_curr['x'], results_curr['y_pred'], 'r--', label='Regression')
axes[0,1].set_xticks(results_curr['x'])
axes[0,1].set_xticklabels(df_current['Variant'], rotation=70, fontsize=8)
axes[0,1].set_ylabel('ΔYield (%)')
axes[0,1].set_title('Current Dataset')
axes[0,1].legend(); axes[0,1].grid(alpha=0.3)

# (3) Combined Regression
axes[1,0].scatter(results_comb['x'], results_comb['y'], c='purple', alpha=0.7)
axes[1,0].plot(results_comb['x'], results_comb['y_pred'], 'r--')
axes[1,0].set_xticks(results_comb['x'])
axes[1,0].set_xticklabels(df_combined['Variant'], rotation=70, fontsize=8)
axes[1,0].set_ylabel('ΔYield (%)')
axes[1,0].set_title('Combined Dataset')
axes[1,0].grid(alpha=0.3)

# (4) ΔΔYield comparison
axes[1,1].bar(df_merged['Variant'], df_merged['ΔΔYield'], color='gray')
axes[1,1].axhline(0, color='red', linestyle='--')
axes[1,1].set_title('ΔΔYield (New - Old)')
axes[1,1].set_ylabel('ΔΔYield (%)')
axes[1,1].tick_params(axis='x', rotation=70)
axes[1,1].grid(alpha=0.3)

# (5) Hydrophobic vs Polar
bars = axes[2,0].bar(hydro_mean.index, hydro_mean.values, color=['#ffcc66','#6699ff'])
axes[2,0].set_title('Hydrophobicity vs Mean ΔYield')
axes[2,0].set_ylabel('Mean ΔYield (%)')
for i, v in enumerate(hydro_mean.values):
    axes[2,0].text(i, v + 0.5, f"{v:.1f}", ha='center', fontsize=8)
axes[2,0].grid(alpha=0.3)

# (6) ΔYield Distribution + Outliers
axes[2,1].hist(df_combined['DeltaYield'], bins=10, color='teal', edgecolor='black', alpha=0.7)
axes[2,1].set_title('Distribution of ΔYield')
axes[2,1].set_xlabel('ΔYield (%)')
axes[2,1].set_ylabel('Frequency')
for _, row in outliers.iterrows():
    axes[2,1].text(row['DeltaYield'], 0.5, row['Variant'], rotation=45, color='red')

plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.savefig('alanine_scanning_extended_analysis.png', dpi=300, bbox_inches='tight')
plt.show()

print("\n✅ Extended figure saved as 'alanine_scanning_extended_analysis.png'")
