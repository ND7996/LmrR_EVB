#!/usr/bin/env python3
"""
Simple script to make two separate single mutations A11L and A91E from WT structure
Creates two output PDB files: one for each mutation
"""

import os
import sys

def make_mutation_pymol(input_pdb, output_pdb, position, from_aa, to_aa):
    """
    Make a single mutation using PyMOL
    """
    try:
        import pymol
        from pymol import cmd
        
        # Initialize PyMOL
        pymol.finish_launching(['pymol', '-c'])  # -c for command line mode
        
        # Load the structure
        cmd.load(input_pdb, "protein")
        
        # Make the mutation
        print(f"Making mutation {from_aa}{position}{to_aa}...")
        cmd.wizard("mutagenesis")
        cmd.get_wizard().set_mode(from_aa)  # Current residue
        cmd.get_wizard().do_select(f"protein and resi {position} and name CA")
        cmd.get_wizard().set_mode(to_aa)  # Target residue
        cmd.get_wizard().apply()
        
        # Save the mutated structure
        cmd.save(output_pdb, "protein")
        
        print(f"✅ Mutant {from_aa}{position}{to_aa} saved to {output_pdb}")
        
        # Clean up
        cmd.reinitialize()
        
    except ImportError:
        print("❌ PyMOL not available. Install PyMOL or use BioPython method.")
        return False
    except Exception as e:
        print(f"❌ Error with PyMOL method: {e}")
        return False
    
    return True

def make_mutation_biopython(input_pdb, output_pdb, position, from_aa, to_aa):
    """
    Make a single mutation using BioPython (simple backbone substitution)
    Note: This only changes residue names, doesn't build side chains properly
    """
    try:
        from Bio import PDB
        
        # Parse the structure
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("protein", input_pdb)
        
        mutation_made = False
        
        # Find and mutate residue
        for model in structure:
            for chain in model:
                for residue in chain:
                    res_id = residue.get_id()[1]  # Get residue number
                    res_name = residue.get_resname()
                    
                    if res_id == position and res_name == from_aa:
                        residue.resname = to_aa
                        print(f"✅ Mutated position {position}: {from_aa} → {to_aa}")
                        mutation_made = True
                        break
        
        if not mutation_made:
            print(f"⚠️  Warning: No {from_aa} residue found at position {position}")
            print("   Check your residue numbering or input structure")
            
        # Save the structure
        io = PDB.PDBIO()
        io.set_structure(structure)
        io.save(output_pdb)
        
        print(f"✅ Structure saved to {output_pdb}")
        if mutation_made:
            print("⚠️  Note: BioPython method only changes residue names.")
            print("   You'll need to rebuild side chains with a modeling program.")
        
    except ImportError:
        print("❌ BioPython not available. Install with: pip install biopython")
        return False
    except Exception as e:
        print(f"❌ Error with BioPython method: {e}")
        return False
    
    return True

def main():
    # File paths - UPDATE THESE!
    input_pdb = "LMRR_WT2.pdb"  # Your WT structure
    output_a11l = "A11L_mutant.pdb"
    output_a91e = "A91E_mutant.pdb"
    
    # Check if input file exists
    if not os.path.exists(input_pdb):
        print(f"❌ Input file {input_pdb} not found!")
        print("Please update the 'input_pdb' variable with your WT structure file.")
        sys.exit(1)
    
    print(f"🧬 Making single mutations from {input_pdb}")
    print("=" * 50)
    
    # Mutation 1: A11L
    print("\n🔄 Making first mutation: A11L")
    success1 = make_mutation_pymol(input_pdb, output_a11l, 11, "ALA", "LEU")
    if not success1:
        print("🔄 Trying BioPython method for A11L...")
        success1 = make_mutation_biopython(input_pdb, output_a11l, 11, "ALA", "LEU")
    
    # Mutation 2: A91E
    print("\n🔄 Making second mutation: A91E")
    success2 = make_mutation_pymol(input_pdb, output_a91e, 91, "ALA", "GLU")
    if not success2:
        print("🔄 Trying BioPython method for A91E...")
        success2 = make_mutation_biopython(input_pdb, output_a91e, 91, "ALA", "GLU")
    
    # Summary
    print("\n" + "=" * 50)
    print("📊 Summary:")
    if success1:
        print(f"✅ A11L mutant: {output_a11l}")
    else:
        print("❌ A11L mutation failed")
        
    if success2:
        print(f"✅ A91E mutant: {output_a91e}")
    else:
        print("❌ A91E mutation failed")
    
    if success1 or success2:
        print("\n📝 Next steps:")
        print("   1. Check the mutant structures visually")
        print("   2. Consider energy minimization")
        print("   3. Run MD equilibration if needed")

if __name__ == "__main__":
    main()