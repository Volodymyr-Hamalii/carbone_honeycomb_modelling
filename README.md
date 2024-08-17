# Carbone honeycomb modelling

To build honeycomb models from .dat or .pdb files. Allows intercalation with other structures.

## To run the program

1. Run
   pip install -r requirements.txt

2. Create init_data/{structure}/ folder with AAARAS.PDB or ljout.dat files.
   E.g. init_data/A1-7_h3/AAARAS.PDB, where structure=A1-7_h3.

3. Run
   python main.py {action} {structure}

## Description

main.py takes the following parameters: action, structure_folder, set:

1. action -- any action from the list above,
2. structure_folder -- structure to process from 'init_data' or 'result_data' folders,
3. set -- parameters custumization (to build bonds, to filter atoms etc.); just put 'set' as third parameter to customize building,
4. optional arguments -- specific parameters if some method needs them.

For example:
python3 main.py show_init_structure A1-7_h3 set

If you don't want to specify some argument -- just provide '\_' on this place.

## Note

1. You have one processed example for A1-7_h3 carbone honeycomb structure.
2. Available actions you can check with 'python main.py help' command.
3. If we don't have .pbd init file (but you have ljout.dat) -- run as a first action 'convert_init_dat_to_pdb'.
4. To put Al into carbon channel set parameters in result_data/{structure}/structure_settings.json JSON file.
5. To work in dev mode (to see error tracebacks and other info) just set in .env DEV_MODE=true.
