import csv

def filter_and_sort_atoms():
    input_file = '/home/hp/nayanika/github/LmrR_EVB/structures/nearby_atoms.csv'
    output_10A = '/home/hp/nayanika/github/LmrR_EVB/structures/atoms_within_10A.csv'
    output_10_15A = '/home/hp/nayanika/github/LmrR_EVB/structures/atoms_10_to_15A.csv'
    
    # Residues to exclude
    exclude_resi = {7, 11, 17, 18, 21, 87, 88, 91, 92, 94, 96, 99}
    
    atoms_10A = []
    atoms_10_15A = []
    
    # Read the CSV file
    with open(input_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            resi = int(row['resi'])
            distance = float(row['distance'])
            
            # Skip excluded residues
            if resi in exclude_resi:
                continue
            
            # Filter by distance
            if distance <= 10.0:
                atoms_10A.append(row)
            elif distance <= 15.0:
                atoms_10_15A.append(row)
    
    # Sort by distance
    atoms_10A.sort(key=lambda x: float(x['distance']))
    atoms_10_15A.sort(key=lambda x: float(x['distance']))
    
    # Write atoms within 10A
    with open(output_10A, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['chain', 'resn', 'resi', 'name', 'distance'])
        writer.writeheader()
        writer.writerows(atoms_10A)
    
    # Write atoms between 10-15A
    with open(output_10_15A, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['chain', 'resn', 'resi', 'name', 'distance'])
        writer.writeheader()
        writer.writerows(atoms_10_15A)
    
    print(f'Saved {len(atoms_10A)} atoms within 10Å to {output_10A}')
    print(f'Saved {len(atoms_10_15A)} atoms between 10-15Å to {output_10_15A}')
    print(f'Excluded residues: {sorted(exclude_resi)}')

if __name__ == '__main__':
    filter_and_sort_atoms()