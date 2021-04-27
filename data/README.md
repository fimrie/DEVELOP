# Data 

We have provided three primary datasets (ZINC, CASF, PDBbind), for which we have pre-computed fragmentations and structural information.

To preprocess these datasets for DEVELOP, please go to `data` directory and run `prepare_data_linker_design.py` or `prepare_data_scaffold_elaboration.py`.

E.g.
```
python prepare_data_linker_design.py --output_dir ./linker_design --add_idx --include_pharms
```

If you want to use DEVELOP with the second setting (which generates linkers/elaborations with a specified number of atoms), run the prepare\_data script with the following arguments:

```
python prepare_data_linker_design.py --data_path PATH_TO_DATA --dataset_name NAME_OF_DATASET --test_mode
```

`prepare_data_linker_design.py` takes three possible input formats, listed below. `prepare_data_scaffold_elaboration.py` takes the same input format, except omitting the distance and angle information.

```
Fragments (SMILES) Distance (Angstrom) Angle (Radians)
Full molecule (SMILES) Linker (SMILES) Fragments (SMILES) Distance (Angstrom) Angle (Radians)
Full molecule (SMILES) Linker (SMILES) Fragments (SMILES) Distance (Angstrom) Angle (Radians) Idx (Integer)
```

Below are examples for linker design and scaffold elaboration:

Linker design
```
COc1ccccc1CCC(=O)c1cccc(F)c1 O=C(CC[*:2])[*:1] COc1ccccc1[*:2].Fc1cccc([*:1])c1 4.69 2.00
```
Scaffold elaboration:
```
COc1cc(Cl)cc(C(=O)Nc2ccc(Cl)cn2)c1NC(=O)c1scc(CN(C)C2=NCCO2)c1Cl O=C(Nc1ccc(Cl)cn1)[*:1] COc1cc(Cl)cc(c1NC(=O)c1scc(CN(C)C2=NCCO2)c1Cl)[*:1]
```

In addition, we provide .types files (see [molgrid](https://github.com/gnina/libmolgrid) documentation for more details) and preprocessed 3D structures (in .gninatypes format, but .sdf is also compatible with molgrid) for the datasets employed.

Due to the large number of files and GitHub file size contraints, structures for the training dataset datasets must be downloaded from [http://opig.stats.ox.ac.uk/resources](http://opig.stats.ox.ac.uk/resources).
