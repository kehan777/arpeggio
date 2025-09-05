#!/usr/bin/env python3
"""
Integrated Arpeggio Analysis and Visualization Tool
Combines pdbe-arpeggio analysis with PyMOL visualization in one command
Now supports automatic PyMOL execution and PSE file generation!

conda activate stab-ddg
cd /mnt/d/conda/jupyter/stab/arpeggio

# Basic usage
python arpeggio_visualizer.py 1crn.cif -s /A/10/

# With automatic PyMOL execution (generates PSE files automatically)
python arpeggio_visualizer.py 1crn.cif -s /A/10/ --auto-pymol

# 多个选择一次处理 (each selection gets separate PyMOL session)
python arpeggio_visualizer.py 1crn.cif -s /A/10/ /B/15/ /C/20/ --auto-pymol

# 不同链和残基 (different chains and residues)
python arpeggio_visualizer.py protein.pdb -s /A/10/ /B/25/ /C/30/ --auto-pymol
"""

import json
import argparse
import os
import sys
import subprocess
import tempfile
import shutil
import time

def run_arpeggio_analysis(structure_file, selection, output_dir="out"):
    """Run pdbe-arpeggio analysis"""
    print(f"Running arpeggio analysis for selection {selection}...")
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Build arpeggio command
    cmd = [
        "pdbe-arpeggio",
        "-s", selection,
        "-o", output_dir,
        structure_file
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print("Arpeggio analysis completed successfully!")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running arpeggio: {e}")
        print(f"stdout: {e.stdout}")
        print(f"stderr: {e.stderr}")
        return False
    except FileNotFoundError:
        print("Error: pdbe-arpeggio not found. Make sure it's installed and in your PATH.")
        return False

def parse_selection(selection_str):
    """Parse selection string like /A/10/ or A:10 into components"""
    if selection_str.startswith('/') and selection_str.endswith('/'):
        # Format: /A/10/
        parts = selection_str.strip('/').split('/')
        if len(parts) >= 2:
            chain = parts[0] if parts[0] else 'A'
            resid = parts[1] if parts[1] else None
            return chain, resid
    elif ':' in selection_str:
        # Format: A:10
        chain, resid = selection_str.split(':')
        return chain, resid
    else:
        # Assume it's just a residue number
        return None, selection_str
    
    return None, None

def generate_enhanced_pymol_script(json_file, structure_file, selection_str, output_script):
    """Generate enhanced PyMOL script with better visualization"""
    
    # Parse selection for filtering
    chain, resid = parse_selection(selection_str)
    filter_selection = f"{chain}:{resid}" if chain and resid else resid
    
    # Enhanced PyMOL configuration with better colors and styling
    pymol_config = {
        "dashcolor": {
            "hbond": {"vdwclash": "purple", "vdw": "purple", "proximal": "purple"},
            "polar": {"vdwclash": "lightblue", "vdw": "lightblue", "proximal": "lightblue"},
            "weakhbond": {"vdwclash": "violet", "vdw": "violet", "proximal": "violet"},
            "weakpolar": {"vdwclash": "paleblue", "vdw": "paleblue", "proximal": "paleblue"},
            "aromatic": {"vdwclash": "red", "vdw": "red", "proximal": "red"},
            "hydrophobic": {"vdwclash": "blue", "vdw": "blue", "proximal": "blue"},
            "carbonyl": {"vdwclash": "orange", "vdw": "orange", "proximal": "orange"},
            "ionic": {"vdwclash": "yellow", "vdw": "yellow", "proximal": "yellow"},
            "metalcomplex": {"covalent": "magenta", "vdwclash": "magenta", "vdw": "magenta"},
            "xbond": {"vdwclash": "cyan", "vdw": "cyan", "proximal": "cyan"},
            "undefined": {"covalent": "white", "vdwclash": "lightcoral", "vdw": "green", "proximal": "grey60"},
            "clash": {"clash": "lightcoral"},
            "covalent": {"covalent": "white"},
            "vdwclash": {"vdwclash": "lightcoral"},
            "vdw": {"vdw": "green"},
            "proximal": {"proximal": "grey60"}
        },
        "dashradius": {
            "hbond": {"vdwclash": 0.03, "vdw": 0.03, "proximal": 0.03},
            "polar": {"vdwclash": 0.025, "vdw": 0.025, "proximal": 0.025},
            "weakhbond": {"vdwclash": 0.02, "vdw": 0.02, "proximal": 0.02},
            "weakpolar": {"vdwclash": 0.02, "vdw": 0.02, "proximal": 0.02},
            "aromatic": {"vdwclash": 0.035, "vdw": 0.035, "proximal": 0.035},
            "hydrophobic": {"vdwclash": 0.025, "vdw": 0.025, "proximal": 0.025},
            "carbonyl": {"vdwclash": 0.025, "vdw": 0.025, "proximal": 0.025},
            "ionic": {"vdwclash": 0.04, "vdw": 0.04, "proximal": 0.04},
            "metalcomplex": {"covalent": 0.04, "vdwclash": 0.04, "vdw": 0.04},
            "xbond": {"vdwclash": 0.025, "vdw": 0.025, "proximal": 0.025},
            "undefined": {"covalent": 0.02, "vdwclash": 0.02, "vdw": 0.02, "proximal": 0.015},
            "clash": {"clash": 0.02},
            "covalent": {"covalent": 0.02},
            "vdwclash": {"vdwclash": 0.02},
            "vdw": {"vdw": 0.02},
            "proximal": {"proximal": 0.015}
        },
        "dashgap": {
            "hbond": {"vdwclash": 0.15, "vdw": 0.15, "proximal": 0.15},
            "polar": {"vdwclash": 0.2, "vdw": 0.2, "proximal": 0.2},
            "weakhbond": {"vdwclash": 0.25, "vdw": 0.25, "proximal": 0.25},
            "weakpolar": {"vdwclash": 0.25, "vdw": 0.25, "proximal": 0.25},
            "aromatic": {"vdwclash": 0.1, "vdw": 0.1, "proximal": 0.1},
            "hydrophobic": {"vdwclash": 0.2, "vdw": 0.2, "proximal": 0.2},
            "carbonyl": {"vdwclash": 0.2, "vdw": 0.2, "proximal": 0.2},
            "ionic": {"vdwclash": 0.1, "vdw": 0.1, "proximal": 0.1},
            "metalcomplex": {"covalent": 0.1, "vdwclash": 0.1, "vdw": 0.1},
            "xbond": {"vdwclash": 0.2, "vdw": 0.2, "proximal": 0.2},
            "undefined": {"covalent": 0.0, "vdwclash": 0.3, "vdw": 0.3, "proximal": 0.8},
            "clash": {"clash": 0.3},
            "covalent": {"covalent": 0.0},
            "vdwclash": {"vdwclash": 0.3},
            "vdw": {"vdw": 0.3},
            "proximal": {"proximal": 0.8}
        },
        "dashlength": {
            "hbond": {"vdwclash": 0.2, "vdw": 0.2, "proximal": 0.2},
            "polar": {"vdwclash": 0.15, "vdw": 0.15, "proximal": 0.15},
            "weakhbond": {"vdwclash": 0.1, "vdw": 0.1, "proximal": 0.1},
            "weakpolar": {"vdwclash": 0.1, "vdw": 0.1, "proximal": 0.1},
            "aromatic": {"vdwclash": 0.25, "vdw": 0.25, "proximal": 0.25},
            "hydrophobic": {"vdwclash": 0.15, "vdw": 0.15, "proximal": 0.15},
            "carbonyl": {"vdwclash": 0.15, "vdw": 0.15, "proximal": 0.15},
            "ionic": {"vdwclash": 0.3, "vdw": 0.3, "proximal": 0.3},
            "metalcomplex": {"covalent": 0.3, "vdwclash": 0.3, "vdw": 0.3},
            "xbond": {"vdwclash": 0.15, "vdw": 0.15, "proximal": 0.15},
            "undefined": {"covalent": 0.5, "vdwclash": 0.1, "vdw": 0.1, "proximal": 0.05},
            "clash": {"clash": 0.1},
            "covalent": {"covalent": 0.5},
            "vdwclash": {"vdwclash": 0.1},
            "vdw": {"vdw": 0.1},
            "proximal": {"proximal": 0.05}
        }
    }
    
    script_lines = []
    
    def add_command(cmd):
        script_lines.append(cmd)
    
    # Enhanced PyMOL setup
    add_command('# Enhanced Arpeggio Visualization Script')
    add_command(f'# Generated for selection: {selection_str}')
    add_command('reinitialize')
    add_command('')
    
    # Enhanced visual settings
    add_command('# Enhanced visual settings')
    add_command('set valence, 1')
    add_command('set stick_rad, 0.12')
    add_command('set line_width, 2')
    add_command('set mesh_width, 0.5')
    add_command('set label_size, 12')
    add_command('set sphere_scale, 0.2')
    add_command('set dash_round_ends, 1')
    add_command('set dash_gap, 0.2')
    add_command('set dash_length, 0.15')
    add_command('')
    
    # Enhanced color definitions
    add_command('# Enhanced color palette')
    add_command("set_color lightblue, (173, 216, 230)")
    add_command("set_color lightcoral, (240, 128, 128)")
    add_command("set_color paleblue, (175, 238, 238)")
    add_command("set_color violet, (238, 130, 238)")
    add_command("set_color gold, (255, 215, 0)")
    add_command("set_color forest, (34, 139, 34)")
    add_command("set_color crimson, (220, 20, 60)")
    add_command('')
    
    # Premium quality settings
    add_command('# Premium quality settings')
    add_command('set line_smooth, 1')
    add_command('set antialias, 3')
    add_command('set cartoon_fancy_helices, 1')
    add_command('set cartoon_smooth_loops, 1')
    add_command('set depth_cue, 1')
    add_command('set specular, 1')
    add_command('set shininess, 50')
    add_command('set surface_quality, 2')
    add_command('set stick_quality, 20')
    add_command('set sphere_quality, 3')
    add_command('set cartoon_sampling, 20')
    add_command('set ribbon_sampling, 15')
    add_command('set transparency_mode, 2')
    add_command('set stick_ball, 1')
    add_command('set stick_ball_ratio, 2.0')
    add_command('set ray_shadows, 1')
    add_command('rebuild')
    add_command('')
    
    # Load structure
    add_command(f'load {structure_file}')
    add_command('select binding_site, None')
    add_command('select target_residue, None')
    add_command('set defer_update, 1')
    add_command('')
    
    # Parse and process contacts
    contacts = parse_json_contacts(json_file, filter_selection)
    interaction_styling = []  # Initialize styling list
    
    if not contacts:
        print("Warning: No contacts found for the specified selection!")
        add_command('# No contacts found for the specified selection')
        add_command('print "No contacts found for the specified selection"')
    else:
        print(f"Found {len(contacts)} contacts for visualization")
        
        # Highlight target residue
        if chain and resid:
            add_command(f'# Highlight target residue {chain}:{resid}')
            add_command(f'select target_residue, chain {chain} and resi {resid}')
            add_command('show sticks, target_residue')
            add_command('color yellow, target_residue')
            add_command('set stick_radius, 0.2, target_residue')
            add_command('')
        
        used_interaction_types = set()
        
        # Process contacts with enhanced grouping
        interaction_counts = {}
        
        for i, contact in enumerate(contacts):
            interaction_type = contact['interaction_type']
            dist_flag = contact['dist_flag']
            
            # Count interactions by type
            if interaction_type not in interaction_counts:
                interaction_counts[interaction_type] = 0
            interaction_counts[interaction_type] += 1
            
            # Handle special distance flags
            if interaction_type in ['clash', 'covalent', 'vdwclash', 'vdw', 'proximal']:
                dist_flag = interaction_type
            
            label = f'{interaction_type}-{dist_flag}'
            
            add_command(f"distance {label}, ({contact['bgn_sel']}), ({contact['end_sel']})")
            add_command(f"show sticks, byres ({contact['bgn_sel']})")
            add_command(f"show sticks, byres ({contact['end_sel']})")
            add_command(f"select binding_site, binding_site + byres ({contact['bgn_sel']})")
            add_command(f"select binding_site, binding_site + byres ({contact['end_sel']})")
            
            used_interaction_types.add((interaction_type, dist_flag))
        
        # Store styling info for later application (after util.cbaw)
        interaction_styling = []
        for interaction_type, flag in used_interaction_types:
            label = f'{interaction_type}-{flag}'
            
            if interaction_type in pymol_config['dashcolor'] and flag in pymol_config['dashcolor'][interaction_type]:
                color = pymol_config['dashcolor'][interaction_type][flag]
                radius = pymol_config['dashradius'][interaction_type][flag]
                gap = pymol_config['dashgap'][interaction_type][flag]
                length = pymol_config['dashlength'][interaction_type][flag]
                
                interaction_styling.append((label, color, radius, gap, length))
        
        # Print interaction summary
        add_command('')
        add_command('# Interaction summary')
        for interaction_type, count in sorted(interaction_counts.items()):
            add_command(f'print "Found {count} {interaction_type} interactions"')
    
    # Enhanced final visualization
    add_command('')
    add_command('# Enhanced final visualization - IMPORTANT: util.cbaw resets colors!')
    add_command('hide labels')
    add_command('util.cbaw')  # This MUST come before interaction styling
    add_command('bg_color white')
    add_command('show cartoon')
    add_command('set cartoon_side_chain_helper, 1')
    add_command('set cartoon_transparency, 0.3')
    add_command('# Keep main chain visible, only hide unnecessary lines')
    #add_command('hide lines, het')  # Only hide lines for heteroatoms
    add_command('hide everything, het')
    add_command('show sticks, het')
    add_command('show spheres, het')
    add_command('color atomic, het')
    
    # Apply interaction styling AFTER util.cbaw to prevent color reset
    if contacts and interaction_styling:
        add_command('')
        add_command('# Apply interaction-specific styling AFTER util.cbaw')
        for label, color, radius, gap, length in interaction_styling:
            add_command(f'color {color}, {label}')
            add_command(f'set dash_radius, {radius}, {label}')
            add_command(f'set dash_gap, {gap}, {label}')
            add_command(f'set dash_length, {length}, {label}')
    
    # Create interaction legend
    add_command('')
    add_command('# Create interaction legend')
    add_command('select legend_area, None')
    
    # Enhanced binding site visualization
    add_command('')
    add_command('# Enhanced binding site visualization')
    add_command('# Surface display disabled per user preference')
    add_command('')
    add_command('# Label interacting residues with amino acid code and position')
    add_command('label binding_site and name CA, oneletter + resi')
    add_command('set label_color, black')
    add_command('set label_size, 14')
    add_command('set label_outline_color, black')
    add_command('set label_position, (0, 0, 2)')  # Position labels above atoms
    
    # Disable less important interactions for clarity
    add_command('disable undefined-proximal')
    add_command('disable proximal-proximal')
    
    # Final settings
    add_command('set defer_update, 0')
    add_command('zoom binding_site')
    add_command('')
    
    # Save session with descriptive name
    base_name = os.path.splitext(os.path.basename(structure_file))[0]
    selection_clean = selection_str.replace('/', '_').replace(':', '_')
    pse_file = f"{base_name}_{selection_clean}.pse"
    add_command(f'save {pse_file}')
    add_command('')
    add_command(f'print "Visualization complete! Session saved as {pse_file}"')
    add_command(f'print "Target selection: {selection_str}"')
    add_command(f'print "Total contacts: {len(contacts) if contacts else 0}"')
    
    # Write script
    with open(output_script, 'w') as f:
        f.write('\n'.join(script_lines))
    
    return pse_file

def parse_json_contacts(json_file, selection=None):
    """Parse arpeggio JSON output and extract contact information"""
    if not os.path.exists(json_file):
        print(f"Warning: JSON file {json_file} not found")
        return []
    
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    # Parse selection if provided
    filter_chain = None
    filter_resid = None
    if selection:
        if ':' in selection:
            filter_chain, filter_resid = selection.split(':')
            filter_resid = int(filter_resid)
        else:
            filter_resid = int(selection)

    contacts = []
    
    for contact in data:
        if contact.get('type') != 'atom-atom':
            continue
            
        atom_bgn = contact.get('bgn', {})
        atom_end = contact.get('end', {})
        
        # Apply selection filter
        if selection:
            bgn_chain = atom_bgn.get('auth_asym_id', 'A')
            bgn_resid = atom_bgn.get('auth_seq_id', 1)
            end_chain = atom_end.get('auth_asym_id', 'A')
            end_resid = atom_end.get('auth_seq_id', 1)
            
            bgn_matches = True
            end_matches = True
            
            if filter_chain:
                bgn_matches = (bgn_chain == filter_chain and bgn_resid == filter_resid)
                end_matches = (end_chain == filter_chain and end_resid == filter_resid)
            else:
                bgn_matches = (bgn_resid == filter_resid)
                end_matches = (end_resid == filter_resid)
            
            if not (bgn_matches or end_matches):
                continue
        
        # Build atom selections
        bgn_sel = f"chain {atom_bgn.get('auth_asym_id', 'A')} and resi {atom_bgn.get('auth_seq_id', 1)} and name {atom_bgn.get('auth_atom_id', 'CA')}"
        end_sel = f"chain {atom_end.get('auth_asym_id', 'A')} and resi {atom_end.get('auth_seq_id', 1)} and name {atom_end.get('auth_atom_id', 'CA')}"
        
        interaction_types = contact.get('contact', [])
        distance = contact.get('distance', 0.0)
        
        # Determine distance flag
        dist_flag = 'proximal'
        for contact_type in interaction_types:
            if contact_type in ['clash']:
                dist_flag = 'clash'
                break
            elif contact_type in ['covalent']:
                dist_flag = 'covalent'
                break
            elif contact_type in ['vdw_clash']:
                dist_flag = 'vdwclash'
                break
            elif contact_type in ['vdw']:
                dist_flag = 'vdw'
                break
            elif contact_type in ['proximal']:
                dist_flag = 'proximal'
                break
        
        # Map contact types
        contact_type_mapping = {
            'clash': 'clash', 'covalent': 'covalent', 'vdw_clash': 'vdwclash',
            'vdw': 'vdw', 'proximal': 'proximal', 'hbond': 'hbond',
            'weak_hbond': 'weakhbond', 'xbond': 'xbond', 'ionic': 'ionic',
            'metal': 'metalcomplex', 'aromatic': 'aromatic', 'hydrophobic': 'hydrophobic',
            'carbonyl': 'carbonyl', 'polar': 'polar', 'weak_polar': 'weakpolar'
        }
        
        for contact_type in interaction_types:
            interaction_type = contact_type_mapping.get(contact_type, 'undefined')
            contacts.append({
                'bgn_sel': bgn_sel,
                'end_sel': end_sel,
                'interaction_type': interaction_type,
                'dist_flag': dist_flag,
                'original_type': contact_type,
                'distance': distance
            })
    
    return contacts

def run_pymol_automation(structure_file, script_file, pse_file):
    """Automatically run PyMOL with the generated script and save PSE file"""
    print(f"Starting PyMOL automation for {script_file}...")
    
    # Create a temporary PyMOL command script
    temp_script = f"temp_pymol_commands_{os.getpid()}.pml"
    
    try:
        with open(temp_script, 'w') as f:
            f.write(f"# Automated PyMOL execution script\n")
            f.write(f"load {structure_file}\n")
            f.write(f"run {script_file}\n")
            f.write(f"save {pse_file}\n")
            f.write("quit\n")
        
        # Run PyMOL with the command script
        cmd = ["pymol", "-c", temp_script]  # -c for command line mode
        
        print(f"Executing: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)  # 5 minute timeout
        
        if result.returncode == 0:
            print(f"✓ PyMOL automation completed successfully!")
            print(f"  - PSE file saved: {pse_file}")
            return True
        else:
            print(f"✗ PyMOL automation failed with return code {result.returncode}")
            if result.stderr:
                print(f"Error output: {result.stderr}")
            return False
            
    except subprocess.TimeoutExpired:
        print("✗ PyMOL automation timed out (5 minutes)")
        return False
    except FileNotFoundError:
        print("✗ PyMOL not found. Make sure PyMOL is installed and in your PATH.")
        print("  You can still manually run: pymol", script_file)
        return False
    except Exception as e:
        print(f"✗ Error during PyMOL automation: {e}")
        return False
    finally:
        # Clean up temporary script
        if os.path.exists(temp_script):
            try:
                os.remove(temp_script)
            except:
                pass

def main():
    parser = argparse.ArgumentParser(description='''
    
    ############################################
    # INTEGRATED ARPEGGIO VISUALIZATION TOOL  #
    ############################################
    
    Combines arpeggio analysis and PyMOL visualization in one command.
    
    Examples:
        # Analyze chain A residue 10
        python arpeggio_visualizer.py 1crn.cif -s /A/10/
        
        # Analyze with automatic PyMOL execution and PSE generation
        python arpeggio_visualizer.py 1crn.cif -s /A/10/ --auto-pymol
        
        # Analyze chain B residue 25 with custom output
        python arpeggio_visualizer.py protein.pdb -s /B/25/ -o results
        
        # Multiple residues with automatic PyMOL (each gets separate PSE file)
        python arpeggio_visualizer.py protein.pdb -s /A/10/ /A/15/ /B/20/ --auto-pymol
    
    ''', formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument('structure_file', help='Path to PDB or CIF structure file')
    parser.add_argument('-s', '--selection', nargs='+', required=True, 
                       help='Selection(s) to analyze (e.g., /A/10/ or A:10)')
    parser.add_argument('-o', '--output', default='out', 
                       help='Output directory (default: out)')
    parser.add_argument('--skip-analysis', action='store_true',
                       help='Skip arpeggio analysis (use existing JSON files)')
    parser.add_argument('--pymol-only', action='store_true',
                       help='Generate PyMOL script only, don\'t run arpeggio')
    parser.add_argument('--auto-pymol', action='store_true',
                       help='Automatically run PyMOL to generate PSE files')
    parser.add_argument('--no-auto-pymol', action='store_true',
                       help='Disable automatic PyMOL execution (default behavior)')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.structure_file):
        print(f"Error: Structure file {args.structure_file} not found")
        sys.exit(1)
    
    # Process each selection
    for selection in args.selection:
        print(f"\n{'='*60}")
        print(f"Processing selection: {selection}")
        print(f"{'='*60}")
        
        # Determine JSON file name
        base_name = os.path.splitext(os.path.basename(args.structure_file))[0]
        json_file = os.path.join(args.output, f"{base_name}.json")
        
        # Run arpeggio analysis if needed
        if not args.skip_analysis and not args.pymol_only:
            success = run_arpeggio_analysis(args.structure_file, selection, args.output)
            if not success:
                print(f"Failed to analyze selection {selection}")
                continue
        
        # Generate PyMOL visualization
        selection_clean = selection.replace('/', '_').replace(':', '_')
        script_name = f"{base_name}_{selection_clean}.pml"
        
        print(f"Generating enhanced PyMOL visualization...")
        pse_file = generate_enhanced_pymol_script(json_file, args.structure_file, selection, script_name)
        
        # Automatic PyMOL execution
        pymol_success = False
        if args.auto_pymol and not args.no_auto_pymol:
            print(f"Running PyMOL automation...")
            pymol_success = run_pymol_automation(args.structure_file, script_name, pse_file)
            if pymol_success:
                print(f"  - PSE file automatically generated: {pse_file}")
        
        print(f"\n✓ Analysis complete for selection {selection}")
        print(f"  - PyMOL script: {script_name}")
        if pymol_success:
            print(f"  - PSE file: {pse_file} (automatically generated)")
        else:
            print(f"  - Expected PSE file: {pse_file}")
            print(f"  - To visualize manually: pymol {script_name}")
    
    print(f"\n{'='*60}")
    print("All selections processed successfully!")
    if args.auto_pymol and not args.no_auto_pymol:
        print("PyMOL sessions have been automatically generated.")
        print("You can open the PSE files directly in PyMOL.")
    else:
        print("Use 'pymol <script_name>' to view the results.")
        print("Or use --auto-pymol flag to automatically generate PSE files.")
    print(f"{'='*60}")

if __name__ == '__main__':
    main()
