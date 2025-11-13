<<<<<<< HEAD
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import re

# =============================================================================
# 🧬 1️⃣ Define datasets manually
# =============================================================================
df_prev = pd.DataFrame({
    "Variant": ['L18A', 'K22A', 'F93A', 'S95A', 'S97A',
                'E7A', 'A11L', 'N19A', 'N88A', 'N89A', 'A92E', 'D100A'],
    "ΔYield":  [-22, -25, -30, -18, -20, 0, 5, 2, 0, 0, -5, 0]
})

df_curr = pd.DataFrame({
    "Variant": ['W96A', 'E7A', 'A11L', 'L18A', 'N19A', 'K22A',
                'N88A', 'M89A', 'A92E', 'F93A', 'S95A', 'S97A',
                'D100A', 'V99A', 'I16A', 'N14A', 'L9A', 'R10A',
                'K101A', 'E104A'],
    "ΔYield":  [-12.0, 0.0, 0.5, -0.5, 0.0, -3.5,
                -3.5, -4.5, -10.0, -3.0, -1.5, -0.5,
                -13.5, 0.0, -6.5, -4.0, -1.5, -1.0,
                -3.5, -3.0]
})

df_latest = pd.DataFrame({
    "Variant": ['V99A', 'I16A', 'N14A', 'L9A', 'R10A',
                'K101A', 'E104A', 'Q49A', 'A50L', 'Y56A'],
    "ΔYield":  [0.5, -0.2, 1.0, -2.5, -1.8, 0.8, -0.9, -3.5, 1.5, -2.0]
})

# =============================================================================
# 🧮 2️⃣ Helper functions
# =============================================================================
def get_pos(v):
    m = re.search(r'\d+', v)
    return int(m.group()) if m else np.nan

for df in [df_prev, df_curr, df_latest]:
    df["Pos"] = df["Variant"].apply(get_pos)

df_all = pd.concat([df_prev, df_curr, df_latest], ignore_index=True)

# =============================================================================
# 📊 3️⃣ Descriptive Statistics
# =============================================================================
print("="*80)
print("📊 Descriptive Statistics for ΔYield")
print("="*80)
for name, df in zip(["Previous", "Current", "Latest"], [df_prev, df_curr, df_latest]):
    print(f"\n{name} dataset:")
    print(df["ΔYield"].describe())

# =============================================================================
# 📈 4️⃣ Distribution Plots
# =============================================================================
plt.figure(figsize=(10, 5))
for df, label, color in zip([df_prev, df_curr, df_latest],
                            ["Previous", "Current", "Latest"],
                            ["#2ca02c", "#1f77b4", "purple"]):
    plt.hist(df["ΔYield"], bins=8, alpha=0.5, label=label, color=color)
plt.xlabel("ΔYield (%)")
plt.ylabel("Frequency")
plt.title("Distribution of ΔYield Across Datasets")
plt.legend()
plt.tight_layout()
plt.savefig("distribution_deltayield.png", dpi=300)
plt.show()

# =============================================================================
# 🔗 5️⃣ Correlation Analysis
# =============================================================================
df_merge = df_prev.merge(df_curr, on="Variant", how="outer", suffixes=("_prev", "_curr"))
df_merge = df_merge.merge(df_latest, on="Variant", how="outer")
corr = df_merge[["ΔYield_prev", "ΔYield_curr", "ΔYield"]].corr()
print("\n🔗 Correlation Matrix:")
print(corr)

plt.figure(figsize=(5,4))
plt.matshow(corr, cmap="coolwarm", fignum=1)
plt.colorbar(label="Correlation")
plt.xticks(range(len(corr.columns)), corr.columns, rotation=45)
plt.yticks(range(len(corr.columns)), corr.columns)
plt.title("Correlation between Datasets", pad=20)
plt.tight_layout()
plt.savefig("correlation_matrix.png", dpi=300)
plt.show()

# =============================================================================
# 🧬 6️⃣ Identify Stabilizing/Destabilizing Mutants
# =============================================================================
print("\n🧬 Most Stabilizing and Destabilizing Mutants")
print("="*80)
for name, df in zip(["Previous", "Current", "Latest"], [df_prev, df_curr, df_latest]):
    best = df.loc[df["ΔYield"].idxmax()]
    worst = df.loc[df["ΔYield"].idxmin()]
    print(f"\n{name} dataset:")
    print(f"  ↑ Most stabilizing: {best['Variant']} ({best['ΔYield']}%)")
    print(f"  ↓ Most destabilizing: {worst['Variant']} ({worst['ΔYield']}%)")

# =============================================================================
# 🧩 7️⃣ Variants common across datasets
# =============================================================================
common = set(df_prev["Variant"]).intersection(df_curr["Variant"]).intersection(df_latest["Variant"])
print("\n🧩 Variants appearing in all datasets:", common if common else "None")

# =============================================================================
# 📊 8️⃣ Average ΔYield across all datasets
# =============================================================================
avg_yield = df_all.groupby("Variant")["ΔYield"].mean().sort_values()
plt.figure(figsize=(10, 5))
plt.bar(avg_yield.index, avg_yield.values, color="teal")
plt.xticks(rotation=70, ha="right")
plt.ylabel("Average ΔYield (%)")
plt.title("Average ΔYield per Variant (All Datasets Combined)")
plt.tight_layout()
plt.savefig("average_deltayield_barplot.png", dpi=300)
plt.show()

# =============================================================================
# 🧭 9️⃣ Position-Based Analysis
# =============================================================================
plt.figure(figsize=(8, 5))
plt.scatter(df_all["Pos"], df_all["ΔYield"], color="darkred", alpha=0.7)
plt.xlabel("Residue Position")
plt.ylabel("ΔYield (%)")
plt.title("ΔYield by Residue Position")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("deltayield_by_position.png", dpi=300)
plt.show()

# Trendline
x = df_all["Pos"].dropna().values
y = df_all["ΔYield"].dropna().values
coef, intercept = np.polyfit(x, y, 1)
r = np.corrcoef(x, y)[0, 1]
plt.figure(figsize=(8,5))
plt.scatter(x, y, alpha=0.7, label="Variants")
plt.plot(x, coef*x + intercept, "r--", label=f"y={coef:.2f}x+{intercept:.2f}, r={r:.2f}")
plt.xlabel("Residue Position")
plt.ylabel("ΔYield (%)")
plt.title("Trend of ΔYield Across Sequence Positions")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("deltayield_trendline.png", dpi=300)
plt.show()

# =============================================================================
# ⚠️ 10️⃣ Outlier Detection
# =============================================================================
z = np.abs(stats.zscore(df_all["ΔYield"]))
outliers = df_all[z > 2]
print("\n⚠️ Potential outlier variants (|z| > 2):")
print(outliers[["Variant", "ΔYield"]])

# =============================================================================
# ✅ Done
# =============================================================================
print("\n✅ Analysis complete — all figures saved in the current folder!")
print("="*80)
=======
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import re

# =============================================================================
# 🧬 1️⃣ Define datasets manually
# =============================================================================
df_prev = pd.DataFrame({
    "Variant": ['L18A', 'K22A', 'F93A', 'S95A', 'S97A',
                'E7A', 'A11L', 'N19A', 'N88A', 'N89A', 'A92E', 'D100A'],
    "ΔYield":  [-22, -25, -30, -18, -20, 0, 5, 2, 0, 0, -5, 0]
})

df_curr = pd.DataFrame({
    "Variant": ['W96A', 'E7A', 'A11L', 'L18A', 'N19A', 'K22A',
                'N88A', 'M89A', 'A92E', 'F93A', 'S95A', 'S97A',
                'D100A', 'V99A', 'I16A', 'N14A', 'L9A', 'R10A',
                'K101A', 'E104A'],
    "ΔYield":  [-12.0, 0.0, 0.5, -0.5, 0.0, -3.5,
                -3.5, -4.5, -10.0, -3.0, -1.5, -0.5,
                -13.5, 0.0, -6.5, -4.0, -1.5, -1.0,
                -3.5, -3.0]
})

df_latest = pd.DataFrame({
    "Variant": ['V99A', 'I16A', 'N14A', 'L9A', 'R10A',
                'K101A', 'E104A', 'Q49A', 'A50L', 'Y56A'],
    "ΔYield":  [0.5, -0.2, 1.0, -2.5, -1.8, 0.8, -0.9, -3.5, 1.5, -2.0]
})

# =============================================================================
# 🧮 2️⃣ Helper functions
# =============================================================================
def get_pos(v):
    m = re.search(r'\d+', v)
    return int(m.group()) if m else np.nan

for df in [df_prev, df_curr, df_latest]:
    df["Pos"] = df["Variant"].apply(get_pos)

df_all = pd.concat([df_prev, df_curr, df_latest], ignore_index=True)

# =============================================================================
# 📊 3️⃣ Descriptive Statistics
# =============================================================================
print("="*80)
print("📊 Descriptive Statistics for ΔYield")
print("="*80)
for name, df in zip(["Previous", "Current", "Latest"], [df_prev, df_curr, df_latest]):
    print(f"\n{name} dataset:")
    print(df["ΔYield"].describe())

# =============================================================================
# 📈 4️⃣ Distribution Plots
# =============================================================================
plt.figure(figsize=(10, 5))
for df, label, color in zip([df_prev, df_curr, df_latest],
                            ["Previous", "Current", "Latest"],
                            ["#2ca02c", "#1f77b4", "purple"]):
    plt.hist(df["ΔYield"], bins=8, alpha=0.5, label=label, color=color)
plt.xlabel("ΔYield (%)")
plt.ylabel("Frequency")
plt.title("Distribution of ΔYield Across Datasets")
plt.legend()
plt.tight_layout()
plt.savefig("distribution_deltayield.png", dpi=300)
plt.show()

# =============================================================================
# 🔗 5️⃣ Correlation Analysis
# =============================================================================
df_merge = df_prev.merge(df_curr, on="Variant", how="outer", suffixes=("_prev", "_curr"))
df_merge = df_merge.merge(df_latest, on="Variant", how="outer")
corr = df_merge[["ΔYield_prev", "ΔYield_curr", "ΔYield"]].corr()
print("\n🔗 Correlation Matrix:")
print(corr)

plt.figure(figsize=(5,4))
plt.matshow(corr, cmap="coolwarm", fignum=1)
plt.colorbar(label="Correlation")
plt.xticks(range(len(corr.columns)), corr.columns, rotation=45)
plt.yticks(range(len(corr.columns)), corr.columns)
plt.title("Correlation between Datasets", pad=20)
plt.tight_layout()
plt.savefig("correlation_matrix.png", dpi=300)
plt.show()

# =============================================================================
# 🧬 6️⃣ Identify Stabilizing/Destabilizing Mutants
# =============================================================================
print("\n🧬 Most Stabilizing and Destabilizing Mutants")
print("="*80)
for name, df in zip(["Previous", "Current", "Latest"], [df_prev, df_curr, df_latest]):
    best = df.loc[df["ΔYield"].idxmax()]
    worst = df.loc[df["ΔYield"].idxmin()]
    print(f"\n{name} dataset:")
    print(f"  ↑ Most stabilizing: {best['Variant']} ({best['ΔYield']}%)")
    print(f"  ↓ Most destabilizing: {worst['Variant']} ({worst['ΔYield']}%)")

# =============================================================================
# 🧩 7️⃣ Variants common across datasets
# =============================================================================
common = set(df_prev["Variant"]).intersection(df_curr["Variant"]).intersection(df_latest["Variant"])
print("\n🧩 Variants appearing in all datasets:", common if common else "None")

# =============================================================================
# 📊 8️⃣ Average ΔYield across all datasets
# =============================================================================
avg_yield = df_all.groupby("Variant")["ΔYield"].mean().sort_values()
plt.figure(figsize=(10, 5))
plt.bar(avg_yield.index, avg_yield.values, color="teal")
plt.xticks(rotation=70, ha="right")
plt.ylabel("Average ΔYield (%)")
plt.title("Average ΔYield per Variant (All Datasets Combined)")
plt.tight_layout()
plt.savefig("average_deltayield_barplot.png", dpi=300)
plt.show()

# =============================================================================
# 🧭 9️⃣ Position-Based Analysis
# =============================================================================
plt.figure(figsize=(8, 5))
plt.scatter(df_all["Pos"], df_all["ΔYield"], color="darkred", alpha=0.7)
plt.xlabel("Residue Position")
plt.ylabel("ΔYield (%)")
plt.title("ΔYield by Residue Position")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("deltayield_by_position.png", dpi=300)
plt.show()

# Trendline
x = df_all["Pos"].dropna().values
y = df_all["ΔYield"].dropna().values
coef, intercept = np.polyfit(x, y, 1)
r = np.corrcoef(x, y)[0, 1]
plt.figure(figsize=(8,5))
plt.scatter(x, y, alpha=0.7, label="Variants")
plt.plot(x, coef*x + intercept, "r--", label=f"y={coef:.2f}x+{intercept:.2f}, r={r:.2f}")
plt.xlabel("Residue Position")
plt.ylabel("ΔYield (%)")
plt.title("Trend of ΔYield Across Sequence Positions")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("deltayield_trendline.png", dpi=300)
plt.show()

# =============================================================================
# ⚠️ 10️⃣ Outlier Detection
# =============================================================================
z = np.abs(stats.zscore(df_all["ΔYield"]))
outliers = df_all[z > 2]
print("\n⚠️ Potential outlier variants (|z| > 2):")
print(outliers[["Variant", "ΔYield"]])

# =============================================================================
# ✅ Done
# =============================================================================
print("\n✅ Analysis complete — all figures saved in the current folder!")
print("="*80)
>>>>>>> 5ab00af8f2b3a8fa0c8224ac7705e9b8d1653770
