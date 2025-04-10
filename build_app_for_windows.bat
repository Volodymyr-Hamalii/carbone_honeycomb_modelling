:: build_app_for_windows.bat

@echo off

:: Use PyInstaller to build the application
pyinstaller --noconfirm --onefile --windowed run_app.py --hidden-import=customtkinter --hidden-import=MDAnalysis.lib.formats.cython_util

:: Remove the build directory
rmdir /s /q build

:: Rename run_app.exe to carbon_honeycomb_manager.exe
move dist\run_app.exe dist\carbon_honeycomb_manager.exe

:: Zip init_data and result_data folders and move to dist .zip files
powershell -Command "Compress-Archive -Path init_data -DestinationPath dist\init_data.zip"
powershell -Command "Compress-Archive -Path result_data -DestinationPath dist\result_data.zip"

:: Confirm the app was built successfully
echo App built successfully to 'dist' folder.
