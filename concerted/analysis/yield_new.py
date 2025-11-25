import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Data for all catalysts with their three replicates (a, b, c)
data = {
    'Catalyst': ['WT', 'W96A', 'E7A', 'A11L', 'L18A', 'N19A', 'K22A', 'N88A', 
                 'M89A', 'A92E', 'F93A', 'S95A', 'S97A', 'D100A', 'V99A', 
                 'I16A', 'N14A', 'L9A', 'R10A', 'K101A', 'E104A'],
    'Replicate_a': [17, 5, 17, 17, 16, 17, 14, 14, 13, 14, 7, 14, 17, 4, 19, 13, 16, 18, 16, 12, 11],
    'Replicate_b': [17, 5, 18, 19, 19, 17, 14, 13, 12, 14, 7, 16, 15, 3, 19, 9, 14, 17, 16, 11, 16],
    'Replicate_c': [17, 5, 16, 16, 17, 14, 13, 15, 13, 14, 7, 16, 18, 3, 13, 9, 10, 12, 11, 18, 15]
}

df = pd.DataFrame(data)
df['Average'] = df[['Replicate_a', 'Replicate_b', 'Replicate_c']].mean(axis=1)
df['Std_Error'] = df[['Replicate_a', 'Replicate_b', 'Replicate_c']].sem(axis=1)

# Create figure with 2 subplots
fig, axes = plt.subplots(1, 2, figsize=(16, 6))
fig.suptitle('Experimental Yield Analysis', fontsize=14, fontweight='bold')

# --- Plot 1: Grouped Bar Chart (All Replicates) ---
x = np.arange(len(df))
width = 0.25

ax1 = axes[0]
bars1 = ax1.bar(x - width, df['Replicate_a'], width, label='Replicate a', color='#3b82f6', alpha=0.8)
bars2 = ax1.bar(x, df['Replicate_b'], width, label='Replicate b', color='#10b981', alpha=0.8)
bars3 = ax1.bar(x + width, df['Replicate_c'], width, label='Replicate c', color='#f59e0b', alpha=0.8)

# Add catalyst labels on top of bars in Plot 1
for i, catalyst in enumerate(df['Catalyst']):
    # Top label
    max_height = max(df.loc[i, ['Replicate_a', 'Replicate_b', 'Replicate_c']])
    ax1.text(i, max_height + 0.5, catalyst, ha='center', va='bottom', 
             fontsize=6, fontweight='bold', rotation=0)

ax1.set_xlabel('Variants', fontweight='bold', fontsize=9)
ax1.set_ylabel('Yield (%)', fontweight='bold', fontsize=9)
ax1.set_title('All Replicates', fontweight='bold', fontsize=10)
ax1.set_xticks(x)
ax1.set_xticklabels(df['Catalyst'], rotation=45, ha='right', fontsize=7)
ax1.legend(fontsize=7)
ax1.grid(axis='y', alpha=0.3)
ax1.set_ylim(0, 25)
ax1.tick_params(axis='y', labelsize=7)

# --- Plot 2: Average Yield with Standard Error Bars ---
ax2 = axes[1]
# Assign colors: WT gets purple, others get performance-based colors
colors = []
for i, (catalyst, avg) in enumerate(zip(df['Catalyst'], df['Average'])):
    if catalyst == 'WT':
        colors.append('#8b5cf6')  # Purple for wild type
    elif avg < 8:
        colors.append('#ef4444')  # Red for low
    elif avg < 14:
        colors.append('#f59e0b')  # Orange for medium
    else:
        colors.append('#10b981')  # Green for high

bars = ax2.bar(x, df['Average'], color=colors, alpha=0.7, edgecolor='black', linewidth=0.5)

# Add error bars (standard error)
ax2.errorbar(x, df['Average'], yerr=df['Std_Error'], 
             fmt='none', ecolor='black', elinewidth=1.5, capsize=4, capthick=1.5)

# Add catalyst labels on top of bars in Plot 2
for i, (avg, std_err, catalyst) in enumerate(zip(df['Average'], df['Std_Error'], df['Catalyst'])):
    # Top label
    height = avg + std_err
    ax2.text(i, height + 0.5, catalyst, ha='center', va='bottom', 
             fontsize=6, fontweight='bold', rotation=0)

ax2.set_xlabel('Variants', fontweight='bold', fontsize=9)
ax2.set_ylabel('Average Yield (%)', fontweight='bold', fontsize=9)
ax2.set_title('Average Yield with Standard Error', fontweight='bold', fontsize=10)
ax2.set_xticks(x)
ax2.set_xticklabels(df['Catalyst'], rotation=45, ha='right', fontsize=7)
ax2.grid(axis='y', alpha=0.3)
ax2.set_ylim(0, 25)
ax2.tick_params(axis='y', labelsize=7)

# Add color legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor='#8b5cf6', label='Wild Type (WT)'),
                   Patch(facecolor='#10b981', label='High (≥14%)'),
                   Patch(facecolor='#f59e0b', label='Medium (8-13%)'),
                   Patch(facecolor='#ef4444', label='Low (<8%)')]
ax2.legend(handles=legend_elements, loc='upper right', fontsize=7)

plt.tight_layout()
plt.savefig('yield_analysis_new.png', dpi=300, bbox_inches='tight')
plt.show()

# Print summary statistics
print("\n" + "="*80)
print("CATALYST YIELD ANALYSIS SUMMARY")
print("="*80)
print("\nAll Catalysts (sorted by average yield):")
df_sorted = df.sort_values('Average', ascending=False)
print(df_sorted[['Catalyst', 'Replicate_a', 'Replicate_b', 'Replicate_c', 'Average', 'Std_Error']].to_string(index=False))
print("\n" + "="*80)