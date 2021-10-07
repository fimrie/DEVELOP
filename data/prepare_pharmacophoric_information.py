#!/usr/bin/env/python
"""
Usage:
    prepare_pharmacophoric_information.py [options]

Options:
    -h --help                Show this screen
    --sdf_path FILE          Path to SD file containing conformers of reference molecules
    --input_path FILE        Path to input file containing examples
    --output_path_base FILE  Base path to output files
    --design_task NAME       Options: linker, elaboration
    --verbose                Print progress and updates to terminal
"""

import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../analysis/'))

from rdkit import Chem
from rdkit.Chem import AllChem

import frag_utils

from docopt import docopt

if __name__ == "__main__":
    # Parse args
    args = docopt(__doc__)
    sdf_path = args.get('--sdf_path')
    input_path = args.get('--input_path')
    output_path_base = args.get('--output_path_base')
    verbose = args.get('--verbose')
    design_task = args.get('--design_task')
    if design_task not in ["linker", "elaboration"]:
        print("Invalid choice for design_task. Must be 'linker' or 'elaboration'.")
    
    # Load data
    fragmentations = []
    with open(input_path, 'r') as f:
        for line in f:
            fragmentations.append(line.strip().split())
    if design_task == "elaboration":
        fragmentations = [f+[0,0] for f in fragmentations]
		
    if verbose:
        print("Num entries in input file: \t%d" % len(fragmentations))

    # Calculate pharmacophoric information
    starting_structures_path = output_path_base+'_starting_structures.sdf'
    pharmacophores_path = output_path_base+'_pharmacophores.sdf'
    fragmentations_pharm, fails = frag_utils.create_frags_pharma_sdf_dataset(fragmentations, 
                                                                             sdf_path, dataset="CASF",
                                                                             sdffile=starting_structures_path,
                                                                             sdffile_pharm=pharmacophores_path,
                                                                             prot="", verbose=True)
	
    print("Done")
