import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import make_interp_spline

rc = np.linspace(0, 100, 2000)
x_pts = [0, 12, 25, 35, 45, 55, 65, 78, 84, 92, 100]
y_pts = [0, -16, -18, 16.5, -10, 9.0, -6, 20.8, 8, -22, -32]
spl = make_interp_spline(x_pts, y_pts, k=3)
y = spl(rc)

plt.figure(figsize=(13, 8))
plt.plot(rc, y, color='#000080', linewidth=6)

# Baseline
plt.axhline(0, color='gray', linestyle='--', linewidth=1.2, alpha=0.8)

# Transition states
plt.plot(35, 16.5, 's', color='#D32F2F', markersize=12, markeredgecolor='black', markeredgewidth=1.2, zorder=10)
plt.plot(55,  9.0, 's', color='#388E3C', markersize=12, markeredgecolor='black', markeredgewidth=1.2, zorder=10)
plt.plot(78, 20.8, 's', color='#D32F2F', markersize=12, markeredgecolor='black', markeredgewidth=1.2, zorder=10)

# Labels (reduced sizes)
plt.text(8,   4, 'reactants',   ha='center', fontsize=12, fontweight='bold')
plt.text(25, -14, '1',          ha='center', fontsize=18, fontweight='bold')
plt.text(50,  -7, '2',          ha='center', fontsize=18, fontweight='bold')
plt.text(94, -27, 'products',   ha='center', fontsize=12, fontweight='bold')
plt.text(23, 13, 'X₁', fontsize=13)
plt.text(48,  7, 'X₂', fontsize=13)
plt.text(72, 16, 'X₃', fontsize=13)

# Ea values (smaller and slightly repositioned)
plt.text(35, 23, 'Ea₁ = 16.5', ha='center', fontsize=10, color='#D32F2F', fontweight='bold')
plt.text(55, 14, 'Ea₂ = 9.0',  ha='center', fontsize=10, color='#388E3C', fontweight='bold')
plt.text(78, 27, 'Ea₃ = 20.8', ha='center', fontsize=10, color='#D32F2F', fontweight='bold')

# 1. Experimental ΔG‡ (purple)
plt.plot([5, 78], [20.8, 20.8], color='#8E24AA', linewidth=4)
plt.plot([78, 78], [0, 20.8],   color='#8E24AA', linewidth=4)
plt.text(41, 28, 'ΔG‡ = kcat/KM = 20.8 kcal/mol\n(experimental activation barrier)', 
         fontsize=11, fontweight='bold', color='#8E24AA', ha='center', va='bottom',
         bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="#8E24AA", linewidth=1.5))

# 2. Intrinsic kcat barrier (orange)
plt.plot([25, 78], [20.8, 20.8], color='#FF8F00', linewidth=4)
plt.plot([25, 25], [-18, 20.8], color='#FF8F00', linewidth=4)
plt.text(52, 33, 'kcat barrier ≈ 38.8 kcal/mol', 
         fontsize=11, fontweight='bold', color='#FF8F00', ha='center', va='bottom',
         bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="#FF8F00", linewidth=1.5))

# 3. Overall ΔG°
plt.annotate('', xy=(98, -32), xytext=(98, 0),
             arrowprops=dict(arrowstyle='<->', lw=4, color='#00695C'))
plt.text(99.5, -16, 'ΔG° ≈ –32 kcal/mol', 
         ha='left', va='center', fontsize=11, fontweight='bold', color='#00695C',
         bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="#00695C", linewidth=1.5))

# Title (split into two lines to save vertical space)
plt.title('Organocatalytic Cycle — Experimental ΔG‡ = 20.8 kcal/mol (kcat/KM)\n'
          'Rate-limiting step is turnover from deep intermediate', 
          fontsize=14, fontweight='bold', pad=20)

plt.ylabel('Free Energy (kcal/mol)', fontsize=13, fontweight='bold')
plt.xlabel('Reaction coordinate', fontsize=13, fontweight='bold')
plt.xlim(0, 100)
plt.ylim(-40, 45)
plt.xticks([])
plt.yticks(np.arange(-40, 50, 10))
plt.grid(True, alpha=0.2, linestyle='--')

plt.tight_layout()
plt.savefig("energy_profile_smaller_fonts.png", dpi=500, bbox_inches='tight')
plt.savefig("energy_profile_smaller_fonts.pdf", bbox_inches='tight')
plt.show()