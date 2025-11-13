<<<<<<< HEAD
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import re

# =============================================================================
# DATASET 1 — Previous
# =============================================================================
previous_data = {
    'Variant': ['L18A', 'K22A', 'F93A', 'S95A', 'S97A',
                'E7A', 'A11L', 'N19A', 'N88A', 'M89A', 'A92E', 'D100A'],
    'DeltaYield': [-22, -25, -30, -18, -20, 0, 5, 2, 0, 0, -5, 0]
}
df_previous = pd.DataFrame(previous_data)

# =============================================================================
# DATASET 2 — New Dataset
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
# SORT VARIANTS NUMERICALLY
# =============================================================================
def extract_position(variant):
    match = re.search(r'(\d+)', variant)
    return int(match.group(1)) if match else np.inf

df_previous['Position'] = df_previous['Variant'].apply(extract_position)
df_current['Position'] = df_current['Variant'].apply(extract_position)
df_previous = df_previous.sort_values('Position')
df_current = df_current.sort_values('Position')

# =============================================================================
# MERGE DATASETS
# =============================================================================
df_combined = pd.concat([df_previous, df_current], ignore_index=True)

# =============================================================================
# REGRESSION FUNCTION
# =============================================================================
def perform_regression(df, label):
    x = np.arange(1, len(df) + 1)
    y = df['DeltaYield'].values
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    y_pred = slope * x + intercept
    residuals = y - y_pred
    r_squared = r_value ** 2
    return {
        'label': label,
        'n': len(x),
        'r': r_value,
        'r_squared': r_squared,
        'p_value': p_value,
        'slope': slope,
        'intercept': intercept,
        'x': x,
        'y': y,
        'y_pred': y_pred,
        'residuals': residuals
    }

previous_results = perform_regression(df_previous, 'Previous')
current_results = perform_regression(df_current, 'Current')
combined_results = perform_regression(df_combined, 'Combined')

# =============================================================================
# PRINT RESULTS
# =============================================================================
def print_results(name, res):
    print(f"\n{name} Dataset Results")
    print("-" * 80)
    print(f"Sample size (n):        {res['n']}")
    print(f"Pearson r:              {res['r']:.4f}")
    print(f"R-squared:              {res['r_squared']:.4f}")
    print(f"P-value:                {res['p_value']:.4e}")
    print(f"Slope:                  {res['slope']:.4f}")
    print(f"Intercept:              {res['intercept']:.4f}")
    print(f"Equation: ΔYield = {res['slope']:.3f} × index + {res['intercept']:.3f}")

print("="*80)
print("ALANINE SCANNING ΔYIELD REGRESSION ANALYSIS")
print("="*80)
print_results("Previous", previous_results)
print_results("Current", current_results)
print_results("Combined", combined_results)

# =============================================================================
# MERGE ALL VARIANTS (INCLUDE UNIQUE ONES)
# =============================================================================
df_merged = pd.merge(df_previous[['Variant', 'DeltaYield']],
                     df_current[['Variant', 'DeltaYield']],
                     on='Variant', how='outer', suffixes=('_Old', '_New'))
df_merged['Position'] = df_merged['Variant'].apply(extract_position)
df_merged = df_merged.sort_values('Position')

# =============================================================================
# PLOTTING — COMPACT, NON-OVERLAPPING
# =============================================================================
plt.figure(figsize=(12, 7))  # smaller figure
plt.suptitle('Alanine Scanning ΔYield Comparison', fontsize=14, fontweight='bold', y=1.02)

# --- (1) Previous dataset
plt.subplot(2, 2, 1)
plt.scatter(previous_results['x'], previous_results['y'], color='#2ca02c', s=40, label='Old Data')
plt.plot(previous_results['x'], previous_results['y_pred'], 'r--', linewidth=1)
plt.xticks(previous_results['x'], df_previous['Variant'], rotation=70, ha='right', fontsize=7)
plt.xlabel('Variant', fontsize=8); plt.ylabel('ΔYield (%)', fontsize=8)
plt.title(f'Previous Dataset (n={previous_results["n"]})', fontsize=10)
plt.legend(fontsize=7); plt.grid(alpha=0.3)

# --- (2) Current dataset
plt.subplot(2, 2, 2)
plt.scatter(current_results['x'], current_results['y'], color='#1f77b4', s=40, label='New Data')
plt.plot(current_results['x'], current_results['y_pred'], 'r--', linewidth=1)
plt.xticks(current_results['x'], df_current['Variant'], rotation=70, ha='right', fontsize=7)
plt.xlabel('Variant', fontsize=8); plt.ylabel('ΔYield (%)', fontsize=8)
plt.title(f'Current Dataset (n={current_results["n"]})', fontsize=10)
plt.legend(fontsize=7); plt.grid(alpha=0.3)

# --- (3) Combined dataset
plt.subplot(2, 2, 3)
plt.scatter(combined_results['x'], combined_results['y'], color='purple', s=35, alpha=0.7, label='Combined')
plt.plot(combined_results['x'], combined_results['y_pred'], 'r--', linewidth=1)
plt.xticks(combined_results['x'], df_combined['Variant'], rotation=70, ha='right', fontsize=7)
plt.xlabel('Variant', fontsize=8); plt.ylabel('ΔYield (%)', fontsize=8)
plt.title(f'Combined Dataset (n={combined_results["n"]})', fontsize=10)
plt.legend(fontsize=7); plt.grid(alpha=0.3)

# --- (4) ΔYield comparison barplot (NO SEABORN)
plt.subplot(2, 2, 4)
variants = df_merged['Variant'].tolist()
x = np.arange(len(variants))
width = 0.35
y_old = df_merged['DeltaYield_Old'].fillna(0).values
y_new = df_merged['DeltaYield_New'].fillna(0).values

plt.bar(x - width/2, y_old, width, label='Old', color='#2ca02c', edgecolor='black', linewidth=0.5)
plt.bar(x + width/2, y_new, width, label='New', color='#1f77b4', edgecolor='black', linewidth=0.5)

plt.xticks(x, variants, rotation=75, ha='right', fontsize=7)
plt.ylabel('ΔYield (%)', fontsize=8)
plt.title('ΔYield Comparison Across All Variants', fontsize=10)
plt.legend(fontsize=7, title='Dataset')
plt.grid(alpha=0.3, axis='y')

# Adjust layout tightly for smaller figure
plt.tight_layout(rect=[0, 0, 1, 0.95], h_pad=1.5, w_pad=1.5)
plt.savefig('alanine_scanning_deltayield_analysis_compact.png', dpi=300, bbox_inches='tight')
plt.show()

# =============================================================================
# SUMMARY
# =============================================================================
print("\n✅ Compact figure saved as 'alanine_scanning_deltayield_analysis_compact.png'")
print(f"Total unique variants: {df_combined['Variant'].nunique()}")
print(f"Variants analyzed (sorted): {', '.join(df_merged['Variant'])}")
print("="*80)
print("Analysis complete and compact figure generated!")
print("="*80)
=======
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import re

# =============================================================================
# DATASET 1 — Previous
# =============================================================================
previous_data = {
    'Variant': ['L18A', 'K22A', 'F93A', 'S95A', 'S97A',
                'E7A', 'A11L', 'N19A', 'N88A', 'M89A', 'A92E', 'D100A'],
    'DeltaYield': [-22, -25, -30, -18, -20, 0, 5, 2, 0, 0, -5, 0]
}
df_previous = pd.DataFrame(previous_data)

# =============================================================================
# DATASET 2 — New Dataset
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
# SORT VARIANTS NUMERICALLY
# =============================================================================
def extract_position(variant):
    match = re.search(r'(\d+)', variant)
    return int(match.group(1)) if match else np.inf

df_previous['Position'] = df_previous['Variant'].apply(extract_position)
df_current['Position'] = df_current['Variant'].apply(extract_position)
df_previous = df_previous.sort_values('Position')
df_current = df_current.sort_values('Position')

# =============================================================================
# MERGE DATASETS
# =============================================================================
df_combined = pd.concat([df_previous, df_current], ignore_index=True)

# =============================================================================
# REGRESSION FUNCTION
# =============================================================================
def perform_regression(df, label):
    x = np.arange(1, len(df) + 1)
    y = df['DeltaYield'].values
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    y_pred = slope * x + intercept
    residuals = y - y_pred
    r_squared = r_value ** 2
    return {
        'label': label,
        'n': len(x),
        'r': r_value,
        'r_squared': r_squared,
        'p_value': p_value,
        'slope': slope,
        'intercept': intercept,
        'x': x,
        'y': y,
        'y_pred': y_pred,
        'residuals': residuals
    }

previous_results = perform_regression(df_previous, 'Previous')
current_results = perform_regression(df_current, 'Current')
combined_results = perform_regression(df_combined, 'Combined')

# =============================================================================
# PRINT RESULTS
# =============================================================================
def print_results(name, res):
    print(f"\n{name} Dataset Results")
    print("-" * 80)
    print(f"Sample size (n):        {res['n']}")
    print(f"Pearson r:              {res['r']:.4f}")
    print(f"R-squared:              {res['r_squared']:.4f}")
    print(f"P-value:                {res['p_value']:.4e}")
    print(f"Slope:                  {res['slope']:.4f}")
    print(f"Intercept:              {res['intercept']:.4f}")
    print(f"Equation: ΔYield = {res['slope']:.3f} × index + {res['intercept']:.3f}")

print("="*80)
print("ALANINE SCANNING ΔYIELD REGRESSION ANALYSIS")
print("="*80)
print_results("Previous", previous_results)
print_results("Current", current_results)
print_results("Combined", combined_results)

# =============================================================================
# MERGE ALL VARIANTS (INCLUDE UNIQUE ONES)
# =============================================================================
df_merged = pd.merge(df_previous[['Variant', 'DeltaYield']],
                     df_current[['Variant', 'DeltaYield']],
                     on='Variant', how='outer', suffixes=('_Old', '_New'))
df_merged['Position'] = df_merged['Variant'].apply(extract_position)
df_merged = df_merged.sort_values('Position')

# =============================================================================
# PLOTTING — COMPACT, NON-OVERLAPPING
# =============================================================================
plt.figure(figsize=(12, 7))  # smaller figure
plt.suptitle('Alanine Scanning ΔYield Comparison', fontsize=14, fontweight='bold', y=1.02)

# --- (1) Previous dataset
plt.subplot(2, 2, 1)
plt.scatter(previous_results['x'], previous_results['y'], color='#2ca02c', s=40, label='Old Data')
plt.plot(previous_results['x'], previous_results['y_pred'], 'r--', linewidth=1)
plt.xticks(previous_results['x'], df_previous['Variant'], rotation=70, ha='right', fontsize=7)
plt.xlabel('Variant', fontsize=8); plt.ylabel('ΔYield (%)', fontsize=8)
plt.title(f'Previous Dataset (n={previous_results["n"]})', fontsize=10)
plt.legend(fontsize=7); plt.grid(alpha=0.3)

# --- (2) Current dataset
plt.subplot(2, 2, 2)
plt.scatter(current_results['x'], current_results['y'], color='#1f77b4', s=40, label='New Data')
plt.plot(current_results['x'], current_results['y_pred'], 'r--', linewidth=1)
plt.xticks(current_results['x'], df_current['Variant'], rotation=70, ha='right', fontsize=7)
plt.xlabel('Variant', fontsize=8); plt.ylabel('ΔYield (%)', fontsize=8)
plt.title(f'Current Dataset (n={current_results["n"]})', fontsize=10)
plt.legend(fontsize=7); plt.grid(alpha=0.3)

# --- (3) Combined dataset
plt.subplot(2, 2, 3)
plt.scatter(combined_results['x'], combined_results['y'], color='purple', s=35, alpha=0.7, label='Combined')
plt.plot(combined_results['x'], combined_results['y_pred'], 'r--', linewidth=1)
plt.xticks(combined_results['x'], df_combined['Variant'], rotation=70, ha='right', fontsize=7)
plt.xlabel('Variant', fontsize=8); plt.ylabel('ΔYield (%)', fontsize=8)
plt.title(f'Combined Dataset (n={combined_results["n"]})', fontsize=10)
plt.legend(fontsize=7); plt.grid(alpha=0.3)

# --- (4) ΔYield comparison barplot (NO SEABORN)
plt.subplot(2, 2, 4)
variants = df_merged['Variant'].tolist()
x = np.arange(len(variants))
width = 0.35
y_old = df_merged['DeltaYield_Old'].fillna(0).values
y_new = df_merged['DeltaYield_New'].fillna(0).values

plt.bar(x - width/2, y_old, width, label='Old', color='#2ca02c', edgecolor='black', linewidth=0.5)
plt.bar(x + width/2, y_new, width, label='New', color='#1f77b4', edgecolor='black', linewidth=0.5)

plt.xticks(x, variants, rotation=75, ha='right', fontsize=7)
plt.ylabel('ΔYield (%)', fontsize=8)
plt.title('ΔYield Comparison Across All Variants', fontsize=10)
plt.legend(fontsize=7, title='Dataset')
plt.grid(alpha=0.3, axis='y')

# Adjust layout tightly for smaller figure
plt.tight_layout(rect=[0, 0, 1, 0.95], h_pad=1.5, w_pad=1.5)
plt.savefig('alanine_scanning_deltayield_analysis_compact.png', dpi=300, bbox_inches='tight')
plt.show()

# =============================================================================
# SUMMARY
# =============================================================================
print("\n✅ Compact figure saved as 'alanine_scanning_deltayield_analysis_compact.png'")
print(f"Total unique variants: {df_combined['Variant'].nunique()}")
print(f"Variants analyzed (sorted): {', '.join(df_merged['Variant'])}")
print("="*80)
print("Analysis complete and compact figure generated!")
print("="*80)
>>>>>>> 5ab00af8f2b3a8fa0c8224ac7705e9b8d1653770
