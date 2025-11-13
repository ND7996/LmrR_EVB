<<<<<<< HEAD
import matplotlib.pyplot as plt

# Kyte–Doolittle hydrophobicity scale
hydrophobicity_scale = {
    'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9, 
    'A': 1.8, 'G': -0.4, 'T': -0.7, 'S': -0.8, 'W': -0.9, 'Y': -1.3, 
    'P': -1.6, 'H': -3.2, 'E': -3.5, 'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5
}

# Separate residues and their scores
residues = list(hydrophobicity_scale.keys())
scores = [hydrophobicity_scale[r] for r in residues]

# Sort by hydrophobicity
residues, scores = zip(*sorted(zip(residues, scores), key=lambda x: x[1], reverse=True))

# Plot
plt.figure(figsize=(8, 5))
bars = plt.bar(residues, scores, color=['darkorange' if s > 0 else 'skyblue' for s in scores])

plt.axhline(0, color='black', linewidth=0.8)
plt.title("Hydrophobicity–Hydrophilicity Scale (Kyte–Doolittle)")
plt.xlabel("Amino Acid")
plt.ylabel("Hydrophobicity Score")
plt.tight_layout()
plt.show()
=======
import matplotlib.pyplot as plt

# Kyte–Doolittle hydrophobicity scale
hydrophobicity_scale = {
    'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9, 
    'A': 1.8, 'G': -0.4, 'T': -0.7, 'S': -0.8, 'W': -0.9, 'Y': -1.3, 
    'P': -1.6, 'H': -3.2, 'E': -3.5, 'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5
}

# Separate residues and their scores
residues = list(hydrophobicity_scale.keys())
scores = [hydrophobicity_scale[r] for r in residues]

# Sort by hydrophobicity
residues, scores = zip(*sorted(zip(residues, scores), key=lambda x: x[1], reverse=True))

# Plot
plt.figure(figsize=(8, 5))
bars = plt.bar(residues, scores, color=['darkorange' if s > 0 else 'skyblue' for s in scores])

plt.axhline(0, color='black', linewidth=0.8)
plt.title("Hydrophobicity–Hydrophilicity Scale (Kyte–Doolittle)")
plt.xlabel("Amino Acid")
plt.ylabel("Hydrophobicity Score")
plt.tight_layout()
plt.show()
>>>>>>> 5ab00af8f2b3a8fa0c8224ac7705e9b8d1653770
