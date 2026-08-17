@echo off
title Aaron's Archive - DNA Analysis Engine
echo.
echo  =====================================================
echo    Aaron's Archive - DNA Analysis Engine v3.2
echo  =====================================================
echo.
python --version >nul 2>&1
if errorlevel 1 (
    echo  ERROR: Python not found. Download from python.org
    pause & exit /b 1
)
echo  Installing dependencies...
pip install flask waitress --quiet 2>nul
echo.
echo  Starting production server...
echo  Open browser at: http://localhost:4040
echo  Press Ctrl+C to stop.
echo.
python -c "from waitress import serve; from app import app; serve(app, host='0.0.0.0', port=4040, threads=4)"
pause
