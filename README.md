# Carbon honeycomb modelling

To build honeycomb models from .dat or .pdb files. Allows intercalation with other structures.

## To run the program

1. Run
   pip install -r requirements.txt

2. Create init_data/{structure}/ folder with AAARAS.PDB or ljout.dat files.
   E.g. init_data/A1-7_h3/AAARAS.PDB, where structure=A1-7_h3.

3. Run
   python main.py

## Note

- To work in dev mode (to see error tracebacks and other info) just set in .env DEV_MODE=true.
