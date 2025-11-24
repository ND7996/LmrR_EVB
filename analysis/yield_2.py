import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score, mean_squared_error
import re

# =============================================================================
# 📊 Load Data
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
                -3.5, -4.5, -3.0, -10.0, -1.5, -0.5,
                -13.5, 0.0, -6.5, -4.0, -1.5, -1.0,
                -3.5, -3.0]
})

def get_pos(v):
    m = re.search(r'\d+', v)
    return int(m.group()) if m else np.nan

for df, name in zip([df_prev, df_curr], ["Old", "New"]):
    df["Pos"] = df["Variant"].apply(get_pos)
    df["Dataset"] = name

df_all = pd.concat([df_prev, df_curr], ignore_index=True)

# =============================================================================
# 1️⃣ STATISTICAL TESTS - Compare distributions across datasets
# =============================================================================
print("\n" + "="*80)
print("📊 STATISTICAL TESTS")
print("="*80)

# Normality tests
print("\n1. Shapiro-Wilk Normality Test:")
for name, df in zip(["Old", "New"], [df_prev, df_curr]):
    stat, p = stats.shapiro(df["ΔYield"])
    print(f"   {name}: W={stat:.4f}, p={p:.4f} {'(Normal)' if p > 0.05 else '(Non-normal)'}")

# Mann-Whitney U test (2 datasets)
stat, p = stats.mannwhitneyu(df_prev["ΔYield"], df_curr["ΔYield"])
print(f"\n2. Mann-Whitney U Test (Old vs New): U={stat:.2f}, p={p:.4f}")
if p < 0.05:
    print("   ✓ Significant difference between datasets")
else:
    print("   ✗ No significant difference between datasets")

# Levene's test for homogeneity of variance
stat, p = stats.levene(df_prev["ΔYield"], df_curr["ΔYield"])
print(f"\n3. Levene's Test (variance homogeneity): W={stat:.4f}, p={p:.4f}")

# =============================================================================
# 2️⃣ MUTATION TYPE ANALYSIS
# =============================================================================
print("\n" + "="*80)
print("🧬 MUTATION TYPE ANALYSIS")
print("="*80)

def extract_aa_info(variant):
    match = re.match(r'([A-Z])(\d+)([A-Z])', variant)
    if match:
        return match.group(1), int(match.group(2)), match.group(3)
    return None, None, None

df_all["WT_AA"], df_all["Position"], df_all["Mut_AA"] = zip(*df_all["Variant"].apply(extract_aa_info))

# Amino acid properties
hydrophobic = set("AVILMFYW")
charged = set("DEKR")
polar = set("STNQ")
aromatic = set("FYW")

def classify_mutation(row):
    wt, mut = row["WT_AA"], row["Mut_AA"]
    if wt is None or mut is None:
        return "Unknown"
    if mut == "A":
        if wt in hydrophobic:
            return "Hydrophobic→Ala"
        elif wt in charged:
            return "Charged→Ala"
        elif wt in polar:
            return "Polar→Ala"
        elif wt in aromatic:
            return "Aromatic→Ala"
    elif mut == "L":
        return "→Leucine"
    return "Other"

df_all["Mutation_Type"] = df_all.apply(classify_mutation, axis=1)

print("\n5. Average ΔYield by Mutation Type:")
mut_summary = df_all.groupby("Mutation_Type")["ΔYield"].agg(['mean', 'std', 'count'])
print(mut_summary)

# =============================================================================
# 3️⃣ REPRODUCIBILITY ANALYSIS - Variants in multiple datasets
# =============================================================================
print("\n" + "="*80)
print("🔄 REPRODUCIBILITY ANALYSIS")
print("="*80)

df_wide = df_all.pivot_table(index="Variant", columns="Dataset", values="ΔYield")
df_replicated = df_wide.dropna(thresh=2)  # At least 2 datasets

print(f"\n6. Variants measured in multiple datasets: {len(df_replicated)}")
print(df_replicated)

# Calculate coefficient of variation for replicated variants
df_replicated["Mean"] = df_replicated.mean(axis=1)
df_replicated["Std"] = df_replicated.std(axis=1)
df_replicated["CV%"] = (df_replicated["Std"] / abs(df_replicated["Mean"]) * 100).round(2)
print("\n7. Coefficient of Variation for Replicated Variants:")
print(df_replicated[["Mean", "Std", "CV%"]].sort_values("CV%", ascending=False))

# Concordance analysis
def classify_effect(val):
    if pd.isna(val):
        return "NA"
    elif val > 1:
        return "Beneficial"
    elif val < -1:
        return "Deleterious"
    else:
        return "Neutral"

for col in df_wide.columns:
    df_wide[f"{col}_Class"] = df_wide[col].apply(classify_effect)

print("\n8. Effect Classification Concordance (for variants in both datasets):")
df_both = df_wide.dropna()
if len(df_both) > 0:
    for variant in df_both.index:
        classes = [df_wide.loc[variant, f"{col}_Class"] for col in ["Old", "New"]]
        concordant = len(set(classes)) == 1
        print(f"   {variant}: {classes} {'✓ Concordant' if concordant else '✗ Discordant'}")
else:
    print("   No variants present in both datasets")

# =============================================================================
# 4️⃣ POSITIONAL HOTSPOT ANALYSIS
# =============================================================================
print("\n" + "="*80)
print("📍 POSITIONAL HOTSPOT ANALYSIS")
print("="*80)

# Divide sequence into regions
df_all["Region"] = pd.cut(df_all["Pos"], bins=[0, 25, 50, 75, 100, 125], 
                           labels=["1-25", "26-50", "51-75", "76-100", "101+"])

print("\n9. Average ΔYield by Sequence Region:")
region_summary = df_all.groupby("Region")["ΔYield"].agg(['mean', 'std', 'count'])
print(region_summary)

# =============================================================================
# 5️⃣ PREDICTIVE MODELING - Predict New from Old
# =============================================================================
print("\n" + "="*80)
print("🤖 PREDICTIVE MODELING")
print("="*80)

# Predict New from Old
overlap_prev_curr = df_wide[["Old", "New"]].dropna()
if len(overlap_prev_curr) > 2:
    X = overlap_prev_curr["Old"].values.reshape(-1, 1)
    y = overlap_prev_curr["New"].values
    
    model = LinearRegression()
    model.fit(X, y)
    y_pred = model.predict(X)
    
    r2 = r2_score(y, y_pred)
    rmse = np.sqrt(mean_squared_error(y, y_pred))
    
    print(f"\n10. Linear Regression: Old → New")
    print(f"    R² = {r2:.4f}")
    print(f"    RMSE = {rmse:.4f}")
    print(f"    Equation: New = {model.coef_[0]:.3f} × Old + {model.intercept_:.3f}")

# =============================================================================
# 6️⃣ CLUSTERING ANALYSIS
# =============================================================================
print("\n" + "="*80)
print("🎯 CLUSTERING ANALYSIS")
print("="*80)

# K-means clustering on variants with data in multiple datasets
# Only use numeric columns for clustering
numeric_cols = ["Old", "New"]
df_cluster_numeric = df_wide[numeric_cols].fillna(df_wide[numeric_cols].mean())

if len(df_cluster_numeric) > 3:
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(df_cluster_numeric)
    
    df_cluster = df_cluster_numeric.copy()
    
    kmeans = KMeans(n_clusters=3, random_state=42, n_init=10)
    df_cluster["Cluster"] = kmeans.fit_predict(X_scaled)
    
    print("\n11. K-Means Clustering (k=3):")
    for cluster in range(3):
        variants = df_cluster[df_cluster["Cluster"] == cluster].index.tolist()
        cluster_data = df_cluster[df_cluster["Cluster"] == cluster][numeric_cols]
        mean_yield = cluster_data.mean().mean()
        print(f"    Cluster {cluster}: {len(variants)} variants, Avg ΔYield = {mean_yield:.2f}%")
        print(f"       {', '.join(variants[:5])}{'...' if len(variants) > 5 else ''}")

# =============================================================================
# 7️⃣ BOOTSTRAP CONFIDENCE INTERVALS
# =============================================================================
print("\n" + "="*80)
print("📊 BOOTSTRAP CONFIDENCE INTERVALS")
print("="*80)

def bootstrap_ci(data, n_bootstrap=1000, ci=95):
    means = []
    for _ in range(n_bootstrap):
        sample = np.random.choice(data, size=len(data), replace=True)
        means.append(np.mean(sample))
    lower = np.percentile(means, (100-ci)/2)
    upper = np.percentile(means, 100-(100-ci)/2)
    return lower, upper

print(f"\n12. Bootstrap 95% CI for Mean ΔYield:")
for name, df in zip(["Old", "New"], [df_prev, df_curr]):
    lower, upper = bootstrap_ci(df["ΔYield"].values)
    mean = df["ΔYield"].mean()
    print(f"    {name}: {mean:.2f}% [{lower:.2f}, {upper:.2f}]")

# =============================================================================
# 8️⃣ RANK CORRELATION (Spearman)
# =============================================================================
print("\n" + "="*80)
print("📈 RANK CORRELATION ANALYSIS")
print("="*80)

print("\n13. Spearman Rank Correlation:")
overlap = df_wide[["Old", "New"]].dropna()
if len(overlap) > 2:
    rho, p = stats.spearmanr(overlap["Old"], overlap["New"])
    print(f"    Old vs New: ρ={rho:.4f}, p={p:.4f}")

# =============================================================================
# 📊 VISUALIZATIONS
# =============================================================================

# 1. Heatmap of all variants across datasets
fig, ax = plt.subplots(figsize=(8, 12))
# Only use numeric columns for heatmap
numeric_cols_heatmap = ["Old", "New"]
df_heatmap = df_wide[numeric_cols_heatmap].fillna(0)
variants = df_heatmap.index.tolist()
datasets = df_heatmap.columns.tolist()

# Create heatmap manually
im = ax.imshow(df_heatmap.values, cmap='RdYlGn', aspect='auto', vmin=-30, vmax=5)
ax.set_xticks(range(len(datasets)))
ax.set_yticks(range(len(variants)))
ax.set_xticklabels(datasets)
ax.set_yticklabels(variants)
plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

# Add text annotations
for i in range(len(variants)):
    for j in range(len(datasets)):
        text = ax.text(j, i, f'{df_heatmap.values[i, j]:.1f}',
                      ha="center", va="center", color="black", fontsize=8)

cbar = plt.colorbar(im, ax=ax)
cbar.set_label('ΔYield (%)', rotation=270, labelpad=15)
plt.title("ΔYield Heatmap: Variants × Datasets", fontweight="bold", fontsize=14, pad=20)
plt.tight_layout()
plt.savefig("heatmap_variants_datasets.png", dpi=300)
plt.show()

# 2. Violin plots comparing distributions
fig, ax = plt.subplots(figsize=(8, 6))
data_for_violin = [df_prev["ΔYield"].values, df_curr["ΔYield"].values]
positions = [1, 2]
colors = ['#2ca02c', '#1f77b4']

for i, (data, pos, color) in enumerate(zip(data_for_violin, positions, colors)):
    parts = ax.violinplot([data], positions=[pos], showmeans=True, showmedians=True)
    for pc in parts['bodies']:
        pc.set_facecolor(color)
        pc.set_alpha(0.7)

ax.set_xticks([1, 2])
ax.set_xticklabels(["Old", "New"])
ax.set_ylabel("ΔYield (%)", fontweight="bold")
ax.set_title("Distribution Comparison (Violin Plots)", fontweight="bold", fontsize=14)
ax.grid(axis='y', alpha=0.3)
plt.tight_layout()
plt.savefig("violin_plot_datasets.png", dpi=300)
plt.show()

# 3. Box plots with statistical annotations
fig, ax = plt.subplots(figsize=(8, 6))
bp = ax.boxplot([df_prev["ΔYield"], df_curr["ΔYield"]], 
                 labels=["Old", "New"], patch_artist=True)
for patch, color in zip(bp['boxes'], ['#2ca02c', '#1f77b4']):
    patch.set_facecolor(color)
    patch.set_alpha(0.6)
ax.set_ylabel("ΔYield (%)", fontweight="bold")
ax.set_title("ΔYield Distribution with Outliers", fontweight="bold", fontsize=14)
ax.grid(axis='y', alpha=0.3)
plt.tight_layout()
plt.savefig("boxplot_with_outliers.png", dpi=300)
plt.show()

# 4. Mutation type comparison
fig, ax = plt.subplots(figsize=(10, 6))
mut_data = df_all.groupby("Mutation_Type")["ΔYield"].apply(list)
mut_labels = mut_data.index.tolist()
bp = ax.boxplot([mut_data[m] for m in mut_labels], labels=mut_labels, patch_artist=True)
for patch in bp['boxes']:
    patch.set_facecolor('#3b82f6')
    patch.set_alpha(0.6)
plt.xticks(rotation=45, ha='right')
ax.set_ylabel("ΔYield (%)", fontweight="bold")
ax.set_title("ΔYield by Mutation Type", fontweight="bold", fontsize=14)
ax.grid(axis='y', alpha=0.3)
plt.tight_layout()
plt.savefig("mutation_type_comparison.png", dpi=300)
plt.show()

# 5. Regional hotspot visualization
fig, ax = plt.subplots(figsize=(10, 6))
region_means = df_all.groupby("Region")["ΔYield"].mean()
colors = ['#ef4444' if x < -5 else '#f59e0b' if x < 0 else '#10b981' for x in region_means]
ax.bar(range(len(region_means)), region_means.values, color=colors, alpha=0.7, edgecolor='black')
ax.set_xticks(range(len(region_means)))
ax.set_xticklabels(region_means.index)
ax.set_ylabel("Average ΔYield (%)", fontweight="bold")
ax.set_xlabel("Sequence Region", fontweight="bold")
ax.set_title("Regional Hotspot Analysis", fontweight="bold", fontsize=14)
ax.axhline(y=0, color='black', linestyle='--', linewidth=0.8)
ax.grid(axis='y', alpha=0.3)
plt.tight_layout()
plt.savefig("regional_hotspots.png", dpi=300)
plt.show()

# 6. Hierarchical clustering dendrogram
if len(df_cluster) > 2:
    fig, ax = plt.subplots(figsize=(12, 6))
    linkage_matrix = linkage(X_scaled, method='ward')
    dendrogram(linkage_matrix, labels=df_cluster.index, leaf_rotation=90, ax=ax)
    ax.set_title("Hierarchical Clustering of Variants", fontweight="bold", fontsize=14)
    ax.set_xlabel("Variant", fontweight="bold")
    ax.set_ylabel("Distance", fontweight="bold")
    plt.tight_layout()
    plt.savefig("dendrogram_clustering.png", dpi=300)
    plt.show()

# 7. Trajectory plot for replicated variants
if len(df_replicated) > 0:
    fig, ax = plt.subplots(figsize=(10, 6))
    for variant in df_replicated.index:
        values = [df_wide.loc[variant, col] for col in ["Old", "New"]]
        values = [v if not pd.isna(v) else None for v in values]
        x_positions = [i for i, v in enumerate(values) if v is not None]
        y_values = [v for v in values if v is not None]
        ax.plot(x_positions, y_values, marker='o', label=variant, alpha=0.7, linewidth=2, markersize=8)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Old", "New"])
    ax.set_ylabel("ΔYield (%)", fontweight="bold")
    ax.set_title("Variant Trajectories Across Datasets", fontweight="bold", fontsize=14)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
    ax.grid(alpha=0.3)
    ax.axhline(y=0, color='red', linestyle='--', linewidth=0.8)
    plt.tight_layout()
    plt.savefig("variant_trajectories.png", dpi=300, bbox_inches='tight')
    plt.show()

# 8. Bland-Altman plot for Old vs New
overlap_pc = df_wide[["Old", "New"]].dropna()
if len(overlap_pc) > 2:
    mean_vals = overlap_pc.mean(axis=1)
    diff_vals = overlap_pc["New"] - overlap_pc["Old"]
    mean_diff = diff_vals.mean()
    std_diff = diff_vals.std()
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(mean_vals, diff_vals, alpha=0.6, s=80)
    ax.axhline(mean_diff, color='red', linestyle='--', label=f'Mean diff: {mean_diff:.2f}')
    ax.axhline(mean_diff + 1.96*std_diff, color='gray', linestyle='--', label=f'+1.96 SD: {mean_diff + 1.96*std_diff:.2f}')
    ax.axhline(mean_diff - 1.96*std_diff, color='gray', linestyle='--', label=f'-1.96 SD: {mean_diff - 1.96*std_diff:.2f}')
    for variant, x, y in zip(overlap_pc.index, mean_vals, diff_vals):
        ax.annotate(variant, (x, y), fontsize=8, alpha=0.7)
    ax.set_xlabel("Mean of Old and New ΔYield", fontweight="bold")
    ax.set_ylabel("Difference (New - Old)", fontweight="bold")
    ax.set_title("Bland-Altman Plot: Agreement Between Datasets", fontweight="bold", fontsize=14)
    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("bland_altman_plot.png", dpi=300)
    plt.show()

print("\n" + "="*80)
print("✅ EXTENDED ANALYSIS COMPLETE")
print("="*80)
print("📊 All additional figures saved!")
print("="*80)