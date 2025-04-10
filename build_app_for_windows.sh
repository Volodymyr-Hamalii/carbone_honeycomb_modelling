#!/bin/bash

pyinstaller --noconfirm --onefile --windowed run_app.py --hidden-import=customtkinter --hidden-import=MDAnalysis.lib.formats.cython_util
rm -r build

# Rename run_app.exe to carbon_honeycomb_manager.exe
mv dist/run_app.exe dist/carbon_honeycomb_manager.exe

# Zip init_data and result_data folders and move to dist .zip files
zip -r dist/init_data.zip init_data
zip -r dist/result_data.zip result_data

mv dist/init_data.zip dist/init_data.zip
mv dist/result_data.zip dist/result_data.zip

echo "App built successfully to 'dist' folder."
