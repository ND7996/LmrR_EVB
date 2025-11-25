import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Data for all catalysts with their four replicates
data = {
    'Catalyst': ['WT', 'E7A', 'A11L', 'L18A', 'N19A', 'K22A', 'N88A', 
                 'M89A', 'A92E', 'F93A', 'S95A', 'S97A', 'D100A'],
    'Rep1': [38, 38, 38, 28, 35, 49, 44, 34, 36, 33, 45, 50, 15],
    'Rep2': [37, 40, 40, 33, 36, 51, 39, 32, 37, 34, 38, 49, 12],
    'Rep3': [47, 40, 43, 40, 48, 42, 26, 40, 35, 26, 34, 43, 21],
    'Rep4': [46, 43, 43, 42, 47, 43, 37, 42, 34, 26, 36, 44, 41]
}

df = pd.DataFrame(data)
df['Mean'] = df[['Rep1', 'Rep2', 'Rep3', 'Rep4']].mean(axis=1)
df['SEM'] = df[['Rep1', 'Rep2', 'Rep3', 'Rep4']].sem(axis=1)

# Create figure with 2 subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
fig.suptitle('Experimental Yield Analysis', fontsize=16, fontweight='bold')

# --- Plot 1: Individual Replicates ---
x = np.arange(len(df))
w = 0.2

bars1 = ax1.bar(x - 1.5*w, df['Rep1'], w, label='Replicate a', color='#3b82f6', alpha=0.85)
bars2 = ax1.bar(x - 0.5*w, df['Rep2'], w, label='Replicate b', color='#10b981', alpha=0.85)
bars3 = ax1.bar(x + 0.5*w, df['Rep3'], w, label='Replicate c', color='#f59e0b', alpha=0.85)
bars4 = ax1.bar(x + 1.5*w, df['Rep4'], w, label='Replicate d', color='#ec4899', alpha=0.85)

# Add catalyst labels on top and bottom of bars in Plot 1
for i, catalyst in enumerate(df['Catalyst']):
    # Top label
    max_height = max(df.loc[i, ['Rep1', 'Rep2', 'Rep3', 'Rep4']])
    ax1.text(i, max_height + 2, catalyst, ha='center', va='bottom', 
             fontsize=9, fontweight='bold', rotation=0)
    # Bottom label
    ax1.text(i, -3, catalyst, ha='center', va='top', 
             fontsize=9, fontweight='bold', rotation=0)

ax1.set_xlabel('Variants', fontweight='bold', fontsize=12)
ax1.set_ylabel('Yield (%)', fontweight='bold', fontsize=12)
ax1.set_title('All Replicates', fontweight='bold', fontsize=13)
ax1.set_xticks(x)
ax1.set_xticklabels([])
ax1.legend(loc='upper right', fontsize=10)
ax1.grid(axis='y', alpha=0.3, linestyle='--')
ax1.set_ylim(-8, 65)  # Extended range for top and bottom labels

# --- Plot 2: Mean with Error Bars ---
bar_colors = []
for i, row in df.iterrows():
    if row['Catalyst'] == 'WT':
        bar_colors.append('#8b5cf6')  # Purple for WT
    elif row['Mean'] >= 40:
        bar_colors.append('#10b981')  # Green for high
    elif row['Mean'] >= 30:
        bar_colors.append('#f59e0b')  # Orange for medium
    else:
        bar_colors.append('#ef4444')  # Red for low

bars = ax2.bar(x, df['Mean'], color=bar_colors, alpha=0.75, edgecolor='black', linewidth=0.8)
ax2.errorbar(x, df['Mean'], yerr=df['SEM'], fmt='none', 
             ecolor='black', elinewidth=2, capsize=5, capthick=2)

# Add catalyst labels on top and bottom of bars in Plot 2
for i, (bar, mean, sem, catalyst) in enumerate(zip(bars, df['Mean'], df['SEM'], df['Catalyst'])):
    # Top label
    height = mean + sem
    ax2.text(i, height + 2, catalyst, ha='center', va='bottom', 
             fontsize=9, fontweight='bold', rotation=0)
    # Bottom label
    ax2.text(i, -3, catalyst, ha='center', va='top', 
             fontsize=9, fontweight='bold', rotation=0)

ax2.set_xlabel('Variants', fontweight='bold', fontsize=12)
ax2.set_ylabel('Average Yield (%)', fontweight='bold', fontsize=12)
ax2.set_title('Average Yield with Standard Error', fontweight='bold', fontsize=13)
ax2.set_xticks(x)
ax2.set_xticklabels([])
ax2.grid(axis='y', alpha=0.3, linestyle='--')
ax2.set_ylim(-8, 65)  # Extended range for top and bottom labels

# Legend for color coding
from matplotlib.patches import Patch
legend = [Patch(facecolor='#8b5cf6', label='Wild Type (WT)'),
          Patch(facecolor='#10b981', label='High (≥40%)'),
          Patch(facecolor='#f59e0b', label='Medium (30-39%)'),
          Patch(facecolor='#ef4444', label='Low (<30%)')]
ax2.legend(handles=legend, loc='upper right', fontsize=10)

plt.tight_layout()
plt.savefig('yield_analysis_labeled.png', dpi=300, bbox_inches='tight')
plt.show()

# Summary statistics
print("\n" + "="*90)
print("CATALYST YIELD ANALYSIS SUMMARY")
print("="*90)
df_sorted = df.sort_values('Mean', ascending=False)
print(df_sorted[['Catalyst', 'Rep1', 'Rep2', 'Rep3', 'Rep4', 'Mean', 'SEM']].to_string(index=False))
print("\n" + "="*90)
print(f"Highest Mean: {df_sorted.iloc[0]['Catalyst']} ({df_sorted.iloc[0]['Mean']:.2f}%)")
print(f"Lowest Mean:  {df_sorted.iloc[-1]['Catalyst']} ({df_sorted.iloc[-1]['Mean']:.2f}%)")
print("="*90)