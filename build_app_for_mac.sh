#!/bin/bash

# Use PyInstaller to build the application
pyinstaller --noconfirm --onefile --windowed run_app.py --hidden-import=customtkinter --hidden-import=MDAnalysis.lib.formats.cython_util

# Remove the build directory
rm -rf build

mkdir carbon_honeycomb_modelling

# Rename run_app to carbon_honeycomb_modelling
# cp -r dist/* carbon_honeycomb_modelling/
cp dist/run_app carbon_honeycomb_modelling/run_app

rm -rf dist

# Copy init_data and result_data folders and move to dist .zip files
cp -r init_data carbon_honeycomb_modelling/init_data
cp -r result_data carbon_honeycomb_modelling/result_data

# Confirm the app was built successfully
echo "App built successfully to 'carbon_honeycomb_modelling' folder."
