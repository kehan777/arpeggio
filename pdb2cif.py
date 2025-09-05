#!/usr/bin/env python3
"""
Enhanced PDB to CIF converter that includes chemical component information
This addresses the arpeggio requirement for _chem_comp category
python pdb2cif.py xxxx.pdb
python pdb2cif.py .
"""

import sys
import os
from Bio import PDB
from Bio.PDB import MMCIFIO
import tempfile
import subprocess

def get_standard_residue_info():
    """Return standard amino acid and nucleotide information"""
    standard_aa = {
        'ALA': {'type': 'L-PEPTIDE LINKING', 'formula': 'C3 H7 N O2'},
        'ARG': {'type': 'L-PEPTIDE LINKING', 'formula': 'C6 H14 N4 O2'},
        'ASN': {'type': 'L-PEPTIDE LINKING', 'formula': 'C4 H8 N2 O3'},
        'ASP': {'type': 'L-PEPTIDE LINKING', 'formula': 'C4 H7 N O4'},
        'CYS': {'type': 'L-PEPTIDE LINKING', 'formula': 'C3 H7 N O2 S'},
        'GLN': {'type': 'L-PEPTIDE LINKING', 'formula': 'C5 H10 N2 O3'},
        'GLU': {'type': 'L-PEPTIDE LINKING', 'formula': 'C5 H9 N O4'},
        'GLY': {'type': 'L-PEPTIDE LINKING', 'formula': 'C2 H5 N O2'},
        'HIS': {'type': 'L-PEPTIDE LINKING', 'formula': 'C6 H9 N3 O2'},
        'ILE': {'type': 'L-PEPTIDE LINKING', 'formula': 'C6 H13 N O2'},
        'LEU': {'type': 'L-PEPTIDE LINKING', 'formula': 'C6 H13 N O2'},
        'LYS': {'type': 'L-PEPTIDE LINKING', 'formula': 'C6 H14 N2 O2'},
        'MET': {'type': 'L-PEPTIDE LINKING', 'formula': 'C5 H11 N O2 S'},
        'PHE': {'type': 'L-PEPTIDE LINKING', 'formula': 'C9 H11 N O2'},
        'PRO': {'type': 'L-PEPTIDE LINKING', 'formula': 'C5 H9 N O2'},
        'SER': {'type': 'L-PEPTIDE LINKING', 'formula': 'C3 H7 N O3'},
        'THR': {'type': 'L-PEPTIDE LINKING', 'formula': 'C4 H9 N O3'},
        'TRP': {'type': 'L-PEPTIDE LINKING', 'formula': 'C11 H12 N2 O2'},
        'TYR': {'type': 'L-PEPTIDE LINKING', 'formula': 'C9 H11 N O3'},
        'VAL': {'type': 'L-PEPTIDE LINKING', 'formula': 'C5 H11 N O2'},
        # Common water and ions
        'HOH': {'type': 'NON-POLYMER', 'formula': 'H2 O'},
        'WAT': {'type': 'NON-POLYMER', 'formula': 'H2 O'},
        'NA': {'type': 'NON-POLYMER', 'formula': 'Na'},
        'CL': {'type': 'NON-POLYMER', 'formula': 'Cl'},
        'MG': {'type': 'NON-POLYMER', 'formula': 'Mg'},
        'CA': {'type': 'NON-POLYMER', 'formula': 'Ca'},
        'ZN': {'type': 'NON-POLYMER', 'formula': 'Zn'},
        'FE': {'type': 'NON-POLYMER', 'formula': 'Fe'},
    }
    return standard_aa

def extract_residues_from_pdb(pdb_file):
    """Extract unique residue names from PDB file"""
    residues = set()
    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith(('ATOM', 'HETATM')):
                    resname = line[17:20].strip()
                    if resname:
                        residues.add(resname)
    except Exception as e:
        print(f"Warning: Could not extract residues from PDB: {e}")
    return residues

def create_enhanced_cif(pdb_file, output_cif=None):
    """Create CIF file with proper chemical component information"""
    if output_cif is None:
        output_cif = os.path.splitext(pdb_file)[0] + ".cif"
    
    # First, try using BioPython for basic conversion
    try:
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("structure", pdb_file)
        
        # Create temporary CIF file
        temp_cif = tempfile.NamedTemporaryFile(mode='w', suffix='.cif', delete=False)
        temp_cif.close()
        
        io = MMCIFIO()
        io.set_structure(structure)
        io.save(temp_cif.name)
        
        # Read the generated CIF and enhance it
        with open(temp_cif.name, 'r') as f:
            cif_content = f.read()
        
        # Extract residues from original PDB
        residues_in_structure = extract_residues_from_pdb(pdb_file)
        standard_residues = get_standard_residue_info()
        
        # Add chemical component information
        chem_comp_section = "\n# Chemical component information\n"
        chem_comp_section += "loop_\n"
        chem_comp_section += "_chem_comp.id\n"
        chem_comp_section += "_chem_comp.type\n"
        chem_comp_section += "_chem_comp.formula\n"
        chem_comp_section += "_chem_comp.name\n"
        
        for residue in sorted(residues_in_structure):
            if residue in standard_residues:
                info = standard_residues[residue]
                chem_comp_section += f"{residue} '{info['type']}' '{info['formula']}' {residue}\n"
            else:
                # For unknown residues, provide minimal info
                chem_comp_section += f"{residue} 'NON-POLYMER' '?' {residue}\n"
        
        chem_comp_section += "#\n"
        
        # Insert chemical component section after the data block header
        lines = cif_content.split('\n')
        enhanced_lines = []
        data_block_found = False
        
        for line in lines:
            enhanced_lines.append(line)
            if line.startswith('data_') and not data_block_found:
                enhanced_lines.append(chem_comp_section)
                data_block_found = True
        
        # Write enhanced CIF
        with open(output_cif, 'w') as f:
            f.write('\n'.join(enhanced_lines))
        
        # Clean up temp file
        os.unlink(temp_cif.name)
        
        print(f"✅ Enhanced CIF conversion: {pdb_file} -> {output_cif}")
        print(f"   Added chemical component info for {len(residues_in_structure)} residue types")
        return True
        
    except Exception as e:
        print(f"❌ Error in enhanced CIF conversion: {e}")
        return False

def try_gemmi_conversion(pdb_file, output_cif):
    """Try using gemmi for PDB to CIF conversion (if available)"""
    try:
        cmd = ['gemmi', 'convert', pdb_file, output_cif]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print(f"✅ Gemmi conversion successful: {pdb_file} -> {output_cif}")
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False

def pdb_to_cif_enhanced(input_pdb, output_cif=None):
    """Enhanced PDB to CIF conversion with multiple fallback methods"""
    if output_cif is None:
        output_cif = os.path.splitext(input_pdb)[0] + "_enhanced.cif"
    
    print(f"Converting {input_pdb} to enhanced CIF format...")
    
    # Method 1: Try gemmi (most reliable for arpeggio)
    if try_gemmi_conversion(input_pdb, output_cif):
        return output_cif
    
    # Method 2: Enhanced BioPython conversion
    if create_enhanced_cif(input_pdb, output_cif):
        return output_cif
    
    print(f"❌ All conversion methods failed for {input_pdb}")
    return None

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage:")
        print("  python pdb2cif.py file.pdb [output.cif]")
        print("  python pdb2cif.py .  # Convert all PDB files in directory")
        print("")
        print("This converter adds chemical component information required by arpeggio.")
        print("It tries multiple conversion methods for best compatibility.")
        sys.exit(1)

    target = sys.argv[1]

    if os.path.isdir(target):
        # Batch mode
        pdb_files = [f for f in os.listdir(target) if f.lower().endswith(".pdb")]
        if not pdb_files:
            print("❌ No .pdb files found in current directory")
            sys.exit(1)
        
        success_count = 0
        for pdb_file in pdb_files:
            full_path = os.path.join(target, pdb_file)
            result = pdb_to_cif_enhanced(full_path)
            if result:
                success_count += 1
        
        print(f"\n✅ Successfully converted {success_count}/{len(pdb_files)} files")
    else:
        # Single file mode
        output_cif = sys.argv[2] if len(sys.argv) > 2 else None
        result = pdb_to_cif_enhanced(target, output_cif)
        if not result:
            sys.exit(1)