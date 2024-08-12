# Carbone honeycomb modelling

To build honeycomb models from .dat or .pdb files. Allows intercalation with other structures.

## To run the program

1) Create init_data/{structure}/ folder with AAARAS.PDB or ljout.dat files.
E.g. init_data/A1-7_h3/AAARAS.PDB, where structure=A1-7_h3.
2) Run python main.py {action} {structure}.

## Note

1) Available actions you can check with 'python main.py help' command.
2) If we don't have .pbd init file - run as a first action 'convert_init_dat_to_pdb'.
3) To put al into carbon channel set parameters in result_data/{structure}/structure_settings.json JSON file.
