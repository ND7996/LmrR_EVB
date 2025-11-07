import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import re

# =============================================================================
# DATASET 1 — PREVIOUS
# =============================================================================
previous_data = {
    'Variant': ['L18A', 'K22A', 'F93A', 'S95A', 'S97A',
                'E7A', 'A11L', 'N19A', 'N88A', 'N89A', 'A92E', 'D100A'],
    'DeltaYield': [-22, -25, -30, -18, -20, 0, 5, 2, 0, 0, -5, 0]
}
df_prev = pd.DataFrame(previous_data)

# =============================================================================
# DATASET 2 — NEW
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
df_curr = pd.DataFrame(current_data)

# =============================================================================
# DATASET 3 — LATEST (from image, variants 99–104)
# =============================================================================
latest_data = {
    'Variant': ['V99A', 'D100A', 'K101A', 'G102A', 'G103A', 'E104A'],
    'DeltaYield': [-0.5, -1.0, -3.0, -2.0, -1.5, -3.0]
}
df_latest = pd.DataFrame(latest_data)

# =============================================================================
# PROCESSING — add numeric position, sort, merge
# =============================================================================
def pos(v):
    m = re.search(r'(\d+)', v)
    return int(m.group(1)) if m else np.nan

for df in [df_prev, df_curr, df_latest]:
    df["Position"] = df["Variant"].apply(pos)
    df.sort_values("Position", inplace=True)

df_combined = pd.concat([df_prev, df_curr, df_latest], ignore_index=True)
df_combined = df_combined.dropna(subset=['DeltaYield']).sort_values("Position").reset_index(drop=True)

# =============================================================================
# REGRESSION FUNCTION
# =============================================================================
def regression(df, label):
    x = np.arange(1, len(df) + 1)
    y = df["DeltaYield"].values
    slope, intercept, r, p, std_err = stats.linregress(x, y)
    return {
        "label": label,
        "x": x,
        "y": y,
        "slope": slope,
        "intercept": intercept,
        "r": r,
        "r2": r**2,
        "p": p,
        "y_pred": slope * x + intercept
    }

results = regression(df_combined, "Combined")

# =============================================================================
# PRINT STATS
# =============================================================================
print("="*70)
print("ΔYIELD ANALYSIS SUMMARY")
print("="*70)
print(f"Total Variants: {df_combined['Variant'].nunique()}")
print(f"Mean ΔYield: {df_combined['DeltaYield'].mean():.3f}")
print(f"Std ΔYield: {df_combined['DeltaYield'].std():.3f}")
print(f"Min ΔYield: {df_combined['DeltaYield'].min():.3f}")
print(f"Max ΔYield: {df_combined['DeltaYield'].max():.3f}")
print("\n--- Regression ---")
print(f"R: {results['r']:.3f} | R²: {results['r2']:.3f} | p={results['p']:.3e}")
print(f"Slope: {results['slope']:.3f} | Intercept: {results['intercept']:.3f}")
print("="*70)

# =============================================================================
# PLOTTING (NO SEABORN)
# =============================================================================
plt.figure(figsize=(15, 9))
plt.suptitle("Alanine Scanning ΔYield Analysis", fontsize=16, fontweight='bold')

# (1) Regression scatter
plt.subplot(2, 3, 1)
plt.scatter(results['x'], results['y'], color='teal', s=50, label='Data')
plt.plot(results['x'], results['y_pred'], 'r--', label='Regression')
plt.xlabel("Variant Index")
plt.ylabel("ΔYield (%)")
plt.title("ΔYield Regression Fit")
plt.legend()
plt.grid(alpha=0.3)

# (2) Histogram
plt.subplot(2, 3, 2)
plt.hist(df_combined['DeltaYield'], bins=12, color='cornflowerblue', edgecolor='black')
plt.xlabel("ΔYield (%)")
plt.ylabel("Frequency")
plt.title("ΔYield Distribution")
plt.grid(alpha=0.3)

# (3) Residuals
residuals = results['y'] - results['y_pred']
plt.subplot(2, 3, 3)
plt.scatter(results['x'], residuals, color='indianred', s=40)
plt.axhline(0, color='black', linestyle='--')
plt.xlabel("Variant Index")
plt.ylabel("Residuals")
plt.title("Residuals of Regression")
plt.grid(alpha=0.3)

# (4) ΔYield per Variant
plt.subplot(2, 3, 4)
x = np.arange(len(df_combined))
plt.bar(x, df_combined['DeltaYield'], color='mediumpurple', alpha=0.8)
plt.xticks(x, df_combined['Variant'], rotation=90, fontsize=7)
plt.xlabel("Variant")
plt.ylabel("ΔYield (%)")
plt.title("ΔYield per Variant (Sorted)")
plt.grid(alpha=0.3, axis='y')

# (5) Compare old vs new vs latest (average per dataset)
plt.subplot(2, 3, 5)
means = [df_prev['DeltaYield'].mean(), df_curr['DeltaYield'].mean(), df_latest['DeltaYield'].mean()]
labels = ['Previous', 'New', 'Latest']
plt.bar(labels, means, color=['#2ca02c', '#1f77b4', '#9467bd'], alpha=0.8)
plt.ylabel("Average ΔYield (%)")
plt.title("Mean ΔYield per Dataset")
plt.grid(alpha=0.3, axis='y')

# (6) Top vs bottom 5
plt.subplot(2, 3, 6)
top5 = df_combined.nlargest(5, 'DeltaYield')
bottom5 = df_combined.nsmallest(5, 'DeltaYield')
plt.bar(top5['Variant'], top5['DeltaYield'], color='darkgreen', label='Top 5')
plt.bar(bottom5['Variant'], bottom5['DeltaYield'], color='darkred', label='Bottom 5')
plt.xticks(rotation=60, fontsize=8)
plt.ylabel("ΔYield (%)")
plt.title("Top & Bottom 5 Variants")
plt.legend()
plt.grid(alpha=0.3, axis='y')

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig("combined_deltayield_analysis.png", dpi=300, bbox_inches='tight')
plt.show()

print("\n✅ Figure saved as 'combined_deltayield_analysis.png'")
print("✅ Complete statistical & graphical analysis done.")
