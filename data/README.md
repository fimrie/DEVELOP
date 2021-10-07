# Data 

## Provided datasets

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

## Prepraing your own dataset

In [Examples] (https://github.com/oxpig/DEVELOP/tree/master/examples/), we provide two example notebooks that demonstrate all required preprocessing and the use of DEVELOP for both linker design and scaffold elaboration.

### Step 1: Compute fragmentations

(Skip this step if you have prepared your own fragmentations).

If you want to simply provide an SD file containing a set of molecules, run `prepare_data_from_sdf.py`.

```
python prepare_data_from_sdf --sdf_path PATH_TO_DATA --output_path PATH_TO_FILE --verbose --design_task TASK
```

`design_task` should either be `linker` or `elaboration`, depending on whether you are using DEVELOP for linker design or scaffold elaboration, respectively.

This will compute fragmentations and structural information, as per the criteria described in our paper. If you do not want to filter the fragmentations using the 2D chemical property filters described in our paper, add the flag `--no_filters` to the above command.

### Step 1b: Compute structural information 

If you are using DEVELOP for linker design (and did not use the above code to compute fragmentations), you must first calculate the structural information by running `calculate_distance_angle.py`. 

All other users should skip this step!

You will need to supply a data file containing a list of fragments and molecules, and an SD file containing a conformation of each molecule.

```
python calculate_distance_angle.py --data_path PATH_TO_FILE --sdf_path PATH_TO_FILE --output_path PATH_TO_FILE --verbose
```

The format of the data file is: 

```
Fragments (SMILES) Full molecule (SMILES)
```

For example:

```
COc1ccccc1[*:2].Fc1cccc([*:1])c1 COc1ccccc1CCC(=O)c1cccc(F)c1
```

### Step 2: Compute pharmacophoric information

To compute the pharmacophoric information for the set of fragmentations, you need to run `prepare_pharmacophoric_information.py`.

Similarly to step 1b above, you will need to supply a data file containing a list of fragments and molecules, and an SD file containing a conformation of each molecule.

```
python prepare_pharmacophoric_information.py --input_path PATH_TO_FILE --sdf_path PATH_TO_FILE --output_path_base BASE_PATH --verbose --design_task TASK
```

Again, `design_task` should either be `linker` or `elaboration`, depending on whether you are using DEVELOP for linker design or scaffold elaboration, respectively.

This will output two SD files containing the 3D representations of the starting substructure(s) and pharmacophoric points required as input to DEVELOP.

### Step 3: Process SD files and create .types file

Finally, we will create a `.types` file and two directories containing the structures as `.gninatypes` files, the custom file format adopted by [molgrid](https://github.com/gnina/libmolgrid). Note that molgrid also accepts other file formats (e.g. .sdf). See the notebooks in `Examples` such an example.

To do this, run the following command, with the case `BASE_PATH` and `TASK` used in Step 2, while `TYPES_NAME` will determine the naming of the generate files (e.g. ZINC, PDBbind).

```
bash create_gninatypes_and_types_files.sh BASE_PATH TYPES_NAME TASK
```