@echo off
setlocal
cd /d "%~dp0"

echo ============================================
echo  RTD Reactor Calculator - starting...
echo ============================================

docker compose up -d --build
if errorlevel 1 (
    echo.
    echo Failed to start. Is Docker Desktop running?
    pause
    exit /b 1
)

echo.
echo Waiting for the app to come up...
timeout /t 3 /nobreak >nul

start http://localhost:8000

echo.
echo RTD Reactor Calculator is running at http://localhost:8000
echo.
echo To stop it later, run: docker compose down
echo (or double-click RTD_cal_stop.bat)
echo.
pause
