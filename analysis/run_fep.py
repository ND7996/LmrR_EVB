#!/usr/bin/env python3
"""
Fixed Q FEP .en file parser for Qdyn 5.10.27 format
"""

import glob
import struct
import numpy as np
from pathlib import Path

def read_en_file_fixed(filename, skip=10, max_frames=10000):
    """
    Read Q FEP .en files in Qdyn 5.10.27 format
    """
    gaps = []
    frame_count = 0
    
    try:
        with open(filename, 'rb') as f:
            print(f"\nReading {filename}...")
            
            # Read and parse header
            canary = struct.unpack('<i', f.read(4))[0]  # 0x6c = 108
            arrays = struct.unpack('<i', f.read(4))[0]  # 0x0539 = 1337
            totresid = struct.unpack('<i', f.read(4))[0]  # 0x01
            
            print(f"  Header: canary={canary}, arrays={arrays}, totresid={totresid}")
            
            # Skip arrays data (types, numres, gcnum)
            f.seek(arrays * 12, 1)
            # Skip resid array
            f.seek(totresid * 4, 1)
            
            # Read version string (80 bytes)
            version = f.read(80).decode('ascii', errors='ignore').strip()
            print(f"  Version: {version}")
            
            # Now read energy records
            while frame_count < max_frames:
                pos = f.tell()
                
                # Read record header (8 bytes: rec_len + unknown)
                header = f.read(8)
                if len(header) < 8:
                    break
                
                rec_len = struct.unpack('<i', header[0:4])[0]
                unknown = struct.unpack('<i', header[4:8])[0]
                
                # Validate record length - should be 124 based on hexdump
                if rec_len != 124:
                    print(f"  Warning: Unexpected record length {rec_len} at frame {frame_count}, expected 124")
                    # Try to resync - look for next 0x7c000000 pattern
                    f.seek(pos)
                    sync_bytes = f.read(4)
                    while sync_bytes:
                        if len(sync_bytes) == 4 and struct.unpack('<i', sync_bytes)[0] == 124:
                            f.seek(-4, 1)  # Back to start of valid record
                            break
                        sync_bytes = f.read(4)
                    continue
                
                # Read the 124-byte record
                record_data = f.read(rec_len)
                if len(record_data) < rec_len:
                    break
                
                # Parse record - structure based on hexdump:
                # Bytes 0-3: record length (124)
                # Bytes 4-7: unknown (0x00000001 in your example)
                # Bytes 8-11: n_states? (0x00000001)
                # Then energy data for states
                
                n_states = struct.unpack('<i', record_data[8:12])[0]
                
                if n_states == 1:
                    # Single state - read EQtot from bytes 12-20
                    e1 = struct.unpack('<d', record_data[12:20])[0]
                    e2 = 0.0  # No second state
                    gap = e1
                elif n_states == 2:
                    # Two states - read both EQtot values
                    # State 1: bytes 12-20
                    e1 = struct.unpack('<d', record_data[12:20])[0]
                    # State 2: bytes 124-132 in next record? Wait for pattern
                    
                    # Actually, based on hexdump, after the first state data,
                    # there's another record marker, then second state
                    # Let's read the trailing marker first
                    trail_marker = f.read(4)
                    if len(trail_marker) < 4:
                        break
                    
                    # Now read second state record
                    header2 = f.read(8)
                    if len(header2) < 8:
                        break
                    
                    rec_len2 = struct.unpack('<i', header2[0:4])[0]
                    if rec_len2 != 124:
                        print(f"  Warning: Second record length {rec_len2} != 124")
                        break
                    
                    record_data2 = f.read(rec_len2)
                    if len(record_data2) < rec_len2:
                        break
                    
                    # Second state EQtot is at bytes 12-20 of second record
                    e2 = struct.unpack('<d', record_data2[12:20])[0]
                    
                    gap = e1 - e2
                    gaps.append(gap)
                    frame_count += 1
                    
                    # Read trailing marker of second record
                    f.read(4)
                else:
                    print(f"  Warning: Unexpected number of states: {n_states}")
                    continue
                
                # Skip to next record (we already read the data)
            
            print(f"  Successfully read {frame_count} frames")
            return gaps[skip:] if len(gaps) > skip else gaps
            
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        import traceback
        traceback.print_exc()
        return []

def read_en_file_simple(filename, skip=10, max_frames=10000):
    """
    Simplified parser that looks for the specific pattern in your files
    """
    gaps = []
    
    try:
        with open(filename, 'rb') as f:
            print(f"\nReading {filename} (simple parser)...")
            
            # Read entire file
            data = f.read()
            pos = 0
            frame_count = 0
            
            # Look for the pattern: 0x7c000000 followed by state data
            while pos < len(data) - 132 and frame_count < max_frames:
                # Look for record marker 0x7c000000 (124 in little-endian)
                if (pos <= len(data) - 4 and 
                    struct.unpack('<i', data[pos:pos+4])[0] == 124):
                    
                    # Check if this looks like a valid energy record
                    if pos + 132 <= len(data):
                        # Try to extract energy value (8 bytes after record header)
                        try:
                            energy = struct.unpack('<d', data[pos+12:pos+20])[0]
                            
                            # Check if this is a reasonable energy value
                            if -10000 < energy < 10000:
                                # Look for the next state in sequence
                                next_pos = pos + 128  # Skip current record + marker
                                if (next_pos <= len(data) - 4 and 
                                    struct.unpack('<i', data[next_pos:next_pos+4])[0] == 124):
                                    
                                    # Extract second energy
                                    energy2 = struct.unpack('<d', data[next_pos+12:next_pos+20])[0]
                                    
                                    if -10000 < energy2 < 10000:
                                        gap = energy - energy2
                                        gaps.append(gap)
                                        frame_count += 1
                                        pos = next_pos + 128  # Move to next pair
                                        continue
                        except struct.error:
                            pass
                
                pos += 1
            
            print(f"  Found {frame_count} valid energy pairs")
            return gaps[skip:] if len(gaps) > skip else gaps
            
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        return []

def compute_fep_energies(gaps, kT=0.596, alpha=20.0, A=25.0, bins=50, min_pts=10):
    """
    Compute dG* and dG0 using FEP with energy gaps.
    """
    if not gaps:
        return None, None
    
    gaps = np.array(gaps)
    print(f"  Energy gap statistics: min={np.min(gaps):.3f}, max={np.max(gaps):.3f}, mean={np.mean(gaps):.3f}")
    
    # Apply alpha shift to state 2 energies
    gaps_shifted = gaps + alpha
    
    # Bin the energy gaps
    gap_min, gap_max = np.min(gaps_shifted), np.max(gaps_shifted)
    bin_edges = np.linspace(gap_min, gap_max, bins + 1)
    hist, _ = np.histogram(gaps_shifted, bins=bin_edges)
    
    # Filter bins with enough points
    valid_bins = hist >= min_pts
    if not np.any(valid_bins):
        print(f"Error: No bins have >= {min_pts} points")
        return None, None
    
    # Compute probabilities
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    prob = hist / np.sum(hist)
    
    # Free energy profile: -kT * ln(P)
    dG = -kT * np.log(prob, where=(prob > 0))
    dG[~valid_bins] = np.inf
    
    # dG0: Free energy difference
    valid_indices = np.where(np.isfinite(dG))[0]
    if len(valid_indices) < 2:
        print("Error: Insufficient valid bins for dG0 calculation")
        return None, None
    
    dG0 = dG[valid_indices[0]] - dG[valid_indices[-1]]
    
    # dG*: Activation energy
    dG_star = np.max(dG[valid_indices]) - np.min(dG[valid_indices])
    
    # Apply coupling constant A
    if A != 0:
        dG0 -= A
        dG_star -= A
    
    return dG0, dG_star

def main():
    print("=" * 60)
    print("Q FEP Energy Analysis (Qdyn 5.10.27 Format)")
    print("=" * 60)
    
    # Find all .en files
    en_files = sorted(glob.glob("fep_*.en"))
    if not en_files:
        print("ERROR: No fep_*.en files found")
        return
    
    print(f"Found {len(en_files)} energy files")
    
    # Parameters
    kT = 0.596
    skip = 10
    alpha = 20.8
    A = 26.5
    bins = 50
    min_pts = 10
    
    # Process each file with simple parser (more reliable for your format)
    all_gaps = []
    for f in en_files:
        gaps = read_en_file_simple(f, skip=skip)
        if gaps:
            print(f"  Extracted {len(gaps)} gaps from {f}")
            all_gaps.extend(gaps)
        else:
            print(f"  No gaps from {f}")
    
    if not all_gaps:
        print("\nTrying alternative parsing method...")
        all_gaps = []
        for f in en_files:
            gaps = read_en_file_fixed(f, skip=skip)
            if gaps:
                print(f"  Extracted {len(gaps)} gaps from {f}")
                all_gaps.extend(gaps)
    
    if not all_gaps:
        print("\nERROR: No valid energy data extracted")
        print("Files may be corrupted or in unexpected format")
        return
    
    print(f"\nTotal gaps collected: {len(all_gaps)}")
    
    # Compute free energies
    print("\nComputing free energies...")
    dG0, dG_star = compute_fep_energies(all_gaps, kT=kT, alpha=alpha, A=A, bins=bins, min_pts=min_pts)
    
    if dG0 is not None and dG_star is not None:
        print("\n=== FREE ENERGY RESULTS ===")
        print(f"dG0 (reaction free energy): {dG0:.3f} kcal/mol")
        print(f"dG* (activation free energy): {dG_star:.3f} kcal/mol")
        
        # Save gaps to file for verification
        np.savetxt('extracted_gaps.txt', all_gaps, fmt='%.6f')
        print(f"\nEnergy gaps saved to 'extracted_gaps.txt' for verification")
    else:
        print("\nERROR: Failed to compute free energies")

if __name__ == "__main__":
    main()