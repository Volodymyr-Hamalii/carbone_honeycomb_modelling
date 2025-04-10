#!/bin/bash

# Use PyInstaller to build the application
pyinstaller --noconfirm --onefile --windowed run_app.py --hidden-import=customtkinter --hidden-import=MDAnalysis.lib.formats.cython_util

# Remove the build directory
rm -rf build

# Rename run_app to carbon_honeycomb_manager
mv dist/run_app dist/carbon_honeycomb_manager

# Zip init_data and result_data folders and move to dist .zip files
zip -r dist/init_data.zip init_data
zip -r dist/result_data.zip result_data

# Confirm the app was built successfully
echo "App built successfully to 'dist' folder."
