@echo off
rem Executa a ferramenta a partir da pasta deste arquivo (caminhos relativos)
cd /d "%~dp0"
python main.py
if errorlevel 1 pause
