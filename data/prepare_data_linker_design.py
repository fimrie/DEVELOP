#!/usr/bin/env/python
"""
Usage:
    prepare_data_linker_design.py [options]

Options:
    -h --help                Show this screen
    --data_path FILE         Path to data file containing fragments and reference molecules
    --dataset_name NAME      Name of dataset (for use in output file naming)
    --output_dir DIR         Path to directory to save data [default: ./]
    --test_mode              To prepare the data for DeLinker in test mode
    --include_pharms         Include pharmacophoric features in structural data
    --add_idx                Add index denoting line number under sdf_idx
"""

import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit.Chem import rdFMCS
import glob
import json
import numpy as np
from utils import bond_dict, dataset_info, need_kekulize, to_graph_mol, graph_to_adj_mat
import utils
import pickle
import random
from docopt import docopt

from align_utils import align_mol_to_frags

dataset = 'zinc'

# Compute pharmacophores
from rdkit.Chem import ChemicalFeatures
from rdkit import RDConfig
fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
factory = ChemicalFeatures.BuildFeatureFactory(fdefName)

pharms = ['Donor', 'Acceptor', 'Aromatic']
def get_pharm_dict(smi):
    feats = factory.GetFeaturesForMol(Chem.MolFromSmiles(smi))
    dic = {'Donor': 0, 'Acceptor': 0, 'Aromatic': 0}
    for feat in feats:
        if feat.GetFamily() in pharms:
            dic[feat.GetFamily()] +=1
    return [dic['Donor'], dic['Acceptor'], dic['Aromatic']]


def read_file(file_path, add_idx=False, calc_pharm_counts=False):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    num_lines = len(lines)
    data = []
    for i, line in enumerate(lines):
        # Parse input data
        toks = line.strip().split(' ')
        if len(toks) == 3:
            smi_frags, abs_dist, angle = toks
            smi_mol = smi_frags
            smi_linker = ''
            if add_idx:
                idx=i
            else:
                idx=-1
        elif len(toks) == 5:
            smi_mol, smi_linker, smi_frags, abs_dist, angle = toks 
            if add_idx:
                idx=i
            else:
                idx=-1
        elif len(toks) == 6:
            smi_mol, smi_linker, smi_frags, abs_dist, angle, idx = toks
        else:
            print("Incorrect input format. Please check the README for useage.")
            exit()
        # Calculate pharmacophoric information
        if calc_pharm_counts == True:
            pharm_count = get_pharm_dict(smi_linker)
        else:
            pharm_count = []
        struct_data = [float(abs_dist), float(angle)]
        struct_data.extend(pharm_count)
        # Add to dataset
        data.append({'smi_mol': smi_mol, 'smi_linker': smi_linker, 
                     'smi_frags': smi_frags,
                     'abs_dist': struct_data,
                     'idx': int(idx)
                     })
        if i % 2000 == 0:
            print('Finished reading: %d / %d' % (i, num_lines), end='\r')
    print('Finished reading: %d / %d' % (num_lines, num_lines))
    return data

def preprocess(raw_data, dataset, name, output_dir='', test=False):
    print('Parsing smiles as graphs.')
    processed_data =[]
    total = len(raw_data)
    for i, (smi_mol, smi_frags, smi_link, abs_dist, idx) in enumerate([(mol['smi_mol'], mol['smi_frags'], 
                                                                         mol['smi_linker'], mol['abs_dist'], mol['idx']) for mol in raw_data]):
        if test:
            smi_mol = smi_frags
            smi_link = ''
        (mol_out, mol_in), nodes_to_keep, exit_points = align_mol_to_frags(smi_mol, smi_link, smi_frags)
        if mol_out == []:
            continue
        nodes_in, edges_in = to_graph_mol(mol_in, dataset)
        nodes_out, edges_out = to_graph_mol(mol_out, dataset)
        if min(len(edges_in), len(edges_out)) <= 0:
            continue
        entry = {
                'graph_in': edges_in,
                'graph_out': edges_out,
                'node_features_in': nodes_in,
                'node_features_out': nodes_out,
                'smiles_out': smi_mol,
                'smiles_in': smi_frags,
                'v_to_keep': nodes_to_keep,
                'exit_points': exit_points,
                'abs_dist': abs_dist
                }
        if idx!=-1:
            entry['sdf_idx']=idx
        processed_data.append(entry)
        # Progress
        if i % 500 == 0:
            print('Processed: %d / %d' % (i, total), end='\r')
    print('Processed: %d / %d' % (total, total))
    print('Saving data')
    # save the dataset
    with open(output_dir + '/molecules_%s.json' % name, 'w') as f:
        json.dump(processed_data, f)
    print('Length raw data: \t%d' % total)
    print('Length processed data: \t%d' % len(processed_data))
          

if __name__ == "__main__":
    # Parse args
    args = docopt(__doc__)
    if args.get('--data_path') and args.get('--dataset_name'):
        data_paths = [args.get('--data_path')]
        names = [args.get('--dataset_name')]
    else:
        data_paths = ['linker_design/data_zinc_train.txt', 'linker_design/data_zinc_valid.txt', 'linker_design/data_zinc_test.txt', 'linker_design/data_casf.txt', 'linker_design/data_pdbbind_2019_refined.txt']
        names = ['zinc_train', 'zinc_valid', 'zinc_test', 'casf_test', 'pdbbind_test']

    output_dir = args.get('--output_dir')
    test_mode = args.get('--test_mode')
    add_idx = args.get('--add_idx')
    calc_pharm_counts = args.get('--include_pharms')

    for data_path, name in zip(data_paths, names):
        print("Preparing: %s" % name)
        raw_data = read_file(data_path, add_idx=add_idx, calc_pharm_counts=calc_pharm_counts)
        preprocess(raw_data, dataset, name, output_dir, test_mode)
