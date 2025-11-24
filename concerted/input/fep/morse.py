import numpy as np

# Try to import plotting libraries, but don't fail if they're missing
try:
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
    print("Matplotlib found - plots will be generated")
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    print("Matplotlib not available - using text output and ASCII plots")

try:
    from scipy.optimize import curve_fit
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    print("Scipy not available - curve fitting disabled")

class MorseParameterCalculator:
    """Calculate Morse potential parameters for C=N iminium ion"""
    
    def __init__(self):
        # Typical values for C=N iminium ion (these can be adjusted based on your system)
        self.bond_dissociation_energy = 5.5  # eV (approximate for C=N+)
        self.equilibrium_bond_length = 1.25  # Angstroms (typical C=N iminium)
        self.force_constant = 1200  # N/m (approximate)
        
        # Physical constants
        self.planck_constant = 6.626e-34  # J⋅s
        self.speed_of_light = 2.998e10  # cm/s
        self.avogadro = 6.022e23  # mol⁻¹
        self.ev_to_joules = 1.602e-19  # J/eV
        
    def calculate_morse_parameters(self, De=None, re=None, k=None):
        """
        Calculate Morse potential parameters
        
        Parameters:
        De: Dissociation energy (eV)
        re: Equilibrium bond length (Angstroms) 
        k: Force constant (N/m)
        
        Returns:
        dict with Morse parameters
        """
        # Use provided values or defaults
        De = De or self.bond_dissociation_energy
        re = re or self.equilibrium_bond_length
        k = k or self.force_constant
        
        # Convert units
        De_joules = De * self.ev_to_joules  # J
        re_meters = re * 1e-10  # m
        
        # Calculate reduced mass for C=N (approximate atomic masses)
        m_C = 12.011  # amu
        m_N = 14.007  # amu
        reduced_mass = (m_C * m_N) / (m_C + m_N)  # amu
        reduced_mass_kg = reduced_mass * 1.66054e-27  # kg
        
        # Calculate Morse parameter 'a'
        # a = sqrt(k / (2 * De * mu))
        a = np.sqrt(k / (2 * De_joules * reduced_mass_kg))  # m⁻¹
        a_angstrom = a * 1e-10  # Å⁻¹
        
        # Calculate vibrational frequency
        omega_e = np.sqrt(k / reduced_mass_kg) / (2 * np.pi)  # Hz
        omega_e_cm = omega_e / (self.speed_of_light)  # cm⁻¹
        
        # Anharmonicity constant
        xe = (self.planck_constant * omega_e) / (4 * De_joules)
        
        results = {
            'De_eV': De,
            'De_joules': De_joules,
            'De_cm': De_joules / (self.planck_constant * self.speed_of_light * 100),
            're_angstrom': re,
            're_meters': re_meters,
            'a_per_angstrom': a_angstrom,
            'a_per_meter': a,
            'force_constant_N_per_m': k,
            'reduced_mass_amu': reduced_mass,
            'omega_e_cm': omega_e_cm,
            'anharmonicity_xe': xe,
            'vibrational_levels': int(1/(2*xe) - 0.5) if xe > 0 else 0
        }
        
        return results
    
    def morse_potential(self, r, De, re, a):
        """Morse potential function"""
        return De * (1 - np.exp(-a * (r - re)))**2
    
    def create_ascii_plot(self, r, V, params, width=70, height=20):
        """Create ASCII art plot"""
        print("\n" + "="*75)
        print("ASCII PLOT: MORSE POTENTIAL FOR C=N IMINIUM ION")
        print("="*75)
        
        # Find plot range
        V_min = min(-0.5, V.min())
        V_max = max(params['De_eV'] * 1.2, V.max())
        r_min, r_max = r.min(), r.max()
        
        # Create plotting grid
        plot_grid = [[' ' for _ in range(width)] for _ in range(height)]
        
        # Draw axes
        for i in range(height):
            plot_grid[i][0] = '|'
        for j in range(width):
            plot_grid[height-1][j] = '-'
        plot_grid[height-1][0] = '+'
        
        # Plot the curve
        for i in range(len(r)):
            if i % 10 == 0:  # Sample every 10th point
                x_pos = int((r[i] - r_min) / (r_max - r_min) * (width - 2)) + 1
                y_pos = height - 1 - int((V[i] - V_min) / (V_max - V_min) * (height - 1))
                
                if 0 <= x_pos < width and 0 <= y_pos < height:
                    plot_grid[y_pos][x_pos] = '*'
        
        # Mark equilibrium position
        re_pos = int((params['re_angstrom'] - r_min) / (r_max - r_min) * (width - 2)) + 1
        if 0 <= re_pos < width:
            for y in range(height):
                if plot_grid[y][re_pos] == ' ':
                    plot_grid[y][re_pos] = '|'
        
        # Print the plot
        for row in plot_grid:
            print(''.join(row))
        
        print(f"\nX-axis: Bond length from {r_min:.1f} to {r_max:.1f} Å")
        print(f"Y-axis: Energy from {V_min:.1f} to {V_max:.1f} eV")
        print(f"Curve: * = Morse potential")
        print(f"Vertical line: | = Equilibrium bond length (re = {params['re_angstrom']:.2f} Å)")
        print("="*75)
    
    def plot_morse_potential(self, params, r_range=(0.8, 2.5), num_points=500):
        """Plot the Morse potential curve"""
        r = np.linspace(r_range[0], r_range[1], num_points)
        V = self.morse_potential(r, params['De_eV'], params['re_angstrom'], 
                                params['a_per_angstrom'])
        
        if MATPLOTLIB_AVAILABLE:
            plt.figure(figsize=(10, 6))
            plt.plot(r, V, 'b-', linewidth=2, label='Morse Potential')
            plt.axhline(y=0, color='k', linestyle='--', alpha=0.5)
            plt.axvline(x=params['re_angstrom'], color='r', linestyle='--', 
                       alpha=0.7, label=f're = {params["re_angstrom"]:.3f} Å')
            
            # Mark dissociation energy
            plt.axhline(y=params['De_eV'], color='g', linestyle='--', 
                       alpha=0.7, label=f'De = {params["De_eV"]:.2f} eV')
            
            plt.xlabel('Bond Length (Å)')
            plt.ylabel('Potential Energy (eV)')
            plt.title('Morse Potential for C=N Iminium Ion')
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.ylim(-0.5, params['De_eV'] * 1.5)
            plt.show()
        else:
            # Use ASCII plot
            self.create_ascii_plot(r, V, params)
        
        return r, V
    
    def save_data_file(self, r, V, params, filename="morse_cn_data.txt"):
        """Save calculation results to a file"""
        try:
            with open(filename, 'w') as f:
                f.write("# Morse Potential Data for C=N Iminium Ion\n")
                f.write(f"# Dissociation Energy (De): {params['De_eV']:.3f} eV\n")
                f.write(f"# Equilibrium Bond Length (re): {params['re_angstrom']:.3f} A\n")
                f.write(f"# Morse Parameter (a): {params['a_per_angstrom']:.3f} A^-1\n")
                f.write(f"# Force Constant: {params['force_constant_N_per_m']:.0f} N/m\n")
                f.write("# \n")
                f.write("# Bond_Length(A)\tPotential_Energy(eV)\n")
                
                for r_val, v_val in zip(r, V):
                    f.write(f"{r_val:.6f}\t{v_val:.6f}\n")
            
            print(f"\nData saved to: {filename}")
            print("You can plot this data in Excel, LibreOffice, or online plotters!")
            return True
        except Exception as e:
            print(f"Could not save file: {e}")
            return False
        
    def print_results(self, params):
        """Print formatted results"""
        print("\n" + "=" * 65)
        print("MORSE POTENTIAL PARAMETERS FOR C=N IMINIUM ION")
        print("=" * 65)
        print(f"Dissociation Energy (De):")
        print(f"  {params['De_eV']:.3f} eV")
        print(f"  {params['De_joules']:.3e} J")
        print(f"  {params['De_cm']:.1f} cm⁻¹")
        print()
        print(f"Equilibrium Bond Length (re):")
        print(f"  {params['re_angstrom']:.3f} Å")
        print(f"  {params['re_meters']:.3e} m")
        print()
        print(f"Morse Parameter (a):")
        print(f"  {params['a_per_angstrom']:.3f} Å⁻¹")
        print(f"  {params['a_per_meter']:.3e} m⁻¹")
        print()
        print(f"Additional Parameters:")
        print(f"  Force Constant: {params['force_constant_N_per_m']:.0f} N/m")
        print(f"  Reduced Mass: {params['reduced_mass_amu']:.3f} amu")
        print(f"  Vibrational Frequency: {params['omega_e_cm']:.1f} cm⁻¹")
        print(f"  Anharmonicity (xe): {params['anharmonicity_xe']:.4f}")
        print(f"  Estimated Vibrational Levels: {params['vibrational_levels']}")
        print("=" * 65)

def main():
    """Main execution function"""
    print("MORSE PARAMETER CALCULATOR FOR C=N IMINIUM ION")
    print("=" * 50)
    
    # Create calculator
    calc = MorseParameterCalculator()
    
    # Calculate with default parameters
    print("\n1. USING DEFAULT PARAMETERS:")
    params = calc.calculate_morse_parameters()
    calc.print_results(params)
    
    # Create plot and data
    print("\n2. GENERATING PLOT:")
    r, V = calc.plot_morse_potential(params)
    
    # Save data file
    print("\n3. SAVING DATA:")
    calc.save_data_file(r, V, params)
    
    # Example with custom parameters
    print("\n4. CUSTOM PARAMETERS EXAMPLE:")
    print("-" * 40)
    custom_params = calc.calculate_morse_parameters(
        De=6.0,  # eV - stronger bond
        re=1.22,  # Å - slightly shorter
        k=1400   # N/m - stiffer
    )
    calc.print_results(custom_params)
    
    # Show some key values in a table format
    print("\n5. COMPARISON TABLE:")
    print("-" * 60)
    print("Parameter        | Default | Custom  | Units")
    print("-" * 60)
    print(f"De               | {params['De_eV']:7.2f} | {custom_params['De_eV']:7.2f} | eV")
    print(f"re               | {params['re_angstrom']:7.3f} | {custom_params['re_angstrom']:7.3f} | Å")
    print(f"a                | {params['a_per_angstrom']:7.3f} | {custom_params['a_per_angstrom']:7.3f} | Å⁻¹")
    print(f"k                | {params['force_constant_N_per_m']:7.0f} | {custom_params['force_constant_N_per_m']:7.0f} | N/m")
    print(f"ωe               | {params['omega_e_cm']:7.1f} | {custom_params['omega_e_cm']:7.1f} | cm⁻¹")
    print(f"xe               | {params['anharmonicity_xe']:7.4f} | {custom_params['anharmonicity_xe']:7.4f} | -")
    print(f"Vib. levels      | {params['vibrational_levels']:7.0f} | {custom_params['vibrational_levels']:7.0f} | -")
    print("-" * 60)
    
    if not MATPLOTLIB_AVAILABLE:
        print("\n6. TO GET BETTER PLOTS:")
        print("Install matplotlib with: pip install matplotlib")
    
    print("\nCalculation completed successfully!")

# Run the program
if __name__ == "__main__":
    main()