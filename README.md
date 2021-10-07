# DEVELOP - Deep Generative Design with 3D Pharmacophoric Constraints

This repository contains our implementation of [Deep Generative Design with 3D Pharmacophoric Constraints](https://www.biorxiv.org/content/10.1101/2021.04.27.441676v1.full) (DEVELOP). 

If you found DEVELOP useful, please cite our preprint:

```
@article{Imrie2021DEVELOP,
  title = {Deep Generative Design with 3D Pharmacophoric Constraints},
  author = {Imrie, Fergus and Hadfield, Thomas E and Bradley, Anthony R and Deane, Charlotte M},
  journal = {bioRxiv},
  doi = {10.1101/2021.04.27.441676},
  year = {2021},
  URL = {https://www.biorxiv.org/content/10.1101/2021.04.27.441676v1.full},
  eprint = {https://www.biorxiv.org/content/10.1101/2021.04.27.441676v1.full.pdf}
}
```

# Acknowledgements

We thank the authors of [Constrained Graph Variational Autoencoders for Molecule Design](https://papers.nips.cc/paper/8005-constrained-graph-variational-autoencoders-for-molecule-design) for releasing their code. The code in this repository is based on their source code release ([link](https://github.com/microsoft/constrained-graph-variational-autoencoder)). If you find this code useful, please consider citing their work.

We thank the authors of [libmolgrid](https://pubs.acs.org/doi/10.1021/acs.jcim.9b01145), which we use for molecular gridding ([link to project](https://github.com/gnina/libmolgrid)). If you use DEVELOP, please consider citing their work.

# Requirements

This code was tested in Python 3.7 with Tensorflow 1.14 and molgrid 0.2.1. 

A yaml file containing all requirements is provided. This can be readily setup using conda.

```
conda env create -f DEVELOP_env.yml
conda activate DEVELOP_env
```

Note that a GPU is required to run DEVELOP as a result of the molgrid dependency.

# Data Extraction

We use three primary datasets (ZINC, CASF, PDBbind).

To preprocess these datasets for DEVELOP, please go to `data` directory and run `prepare_data_linker_design.py` or `prepare_data_scaffold_elaboration.py`.

E.g.
```
python prepare_data_linker_design.py --output_dir ./linker_design --add_idx --include_pharms
```

# Running DEVELOP

DEVELOP can be used to perform both linker design and scaffold elaboration. The below commands demonstrate how to train and generate molecules for linker design. Scaffold elaboration can be performed analogously. 

We provide two settings to generate molecules with DEVELOP. The first setting generates molecules with the same number of atoms as the reference molecule. The second setting generates elaborations with a specified number of atoms. 

To train and generate molecules using the first setting, use:

```
python DEVELOP.py --dataset zinc --config '{"num_epochs": 10, "epoch_to_generate": 10, "train_file": "data/linker_design/molecules_zinc_train.json", "train_struct_file": "data/linker_design/zinc_train_structs.types", "valid_file": "data/linker_design/molecules_zinc_test.json", "valid_struct_file": "data/linker_design/zinc_test_structs.types", "struct_data_root": "./data/linker_design/"}'
```

To generate molecules with a pretrained model using the first setting, use:

```
python DEVELOP.py --dataset zinc --restore models/linker_design/pretrained_DEVELOP_model.pickle --config '{"generation": true, "number_of_generation_per_valid": 10, "batch_size": 1, "train_file": "data/linker_design/molecules_casf_test.json", "train_struct_file": "data/linker_design/casf_test_structs.types", "valid_file": "data/linker_design/molecules_casf_test.json", "valid_struct_file": "data/linker_design/casf_test_structs.types", "struct_data_root": "./data/linker_design/"}'
```

To generate molecules using the second setting, first prepare the data for using the `--test_mode` option:

```
python prepare_data_linker_design.py --add_idx --include_pharms --test_mode --data_path linker_design/data_casf.txt --dataset_name casf_test_setting_2 --output_dir ./linker_design
```

Then use:

```
python DEVELOP.py --dataset zinc --restore models/linker_design/pretrained_DEVELOP_model.pickle --config '{"generation": true, "number_of_generation_per_valid": 10, "batch_size": 1, "train_file": "data/linker_design/molecules_casf_test_setting_2.json", "train_struct_file": "data/linker_design/casf_test_structs.types", "valid_file": "data/linker_design/molecules_casf_test_setting_2.json", "valid_struct_file": "data/linker_design/casf_test_structs.types", "struct_data_root": "./data/linker_design/", "min_atoms": 5, "max_atoms": 6}'
```

In both cases, the output is of the following format:

```
Input fragments (SMILES) Ground truth molecule/fragments (SMILES) Generated molecule (SMILES)
```

More configurations can be found in `default_params` in `DEVELOP.py`.

# Pretrained Models and Generated Molecules

We provide pretrained models for both linker design and scaffold elaboration:

```
models/linker_design/pretrained_DEVELOP_model.pickle
models/scaffold_elaboration/pretrained_DEVELOP_model.pickle
```

Generated molecules can be obtained upon request.

# Contact (Questions/Bugs/Requests)

Please submit a Github issue or contact the Oxford Protein Informatics Group (OPIG) [deane@stats.ox.ac.uk](mailto:deane@stats.ox.ac.uk).
