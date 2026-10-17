@echo off
setlocal
cd /d "%~dp0"

echo Stopping RTD Reactor Calculator...
docker compose down

echo.
pause
