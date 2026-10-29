@echo off
echo ====================================================================
echo   PUSHING LyE-Running-Biomechanics TO GITHUB
echo ====================================================================
echo.

cd /d "%~dp0"

echo 1. Initializing Git repository...
git init
if %errorlevel% neq 0 (
    echo Git is not installed or not in PATH. Please install Git.
    pause
    exit /b
)

echo.
echo 2. Adding all files...
git add .

echo.
echo 3. Creating initial commit...
git commit -m "Release v1.0: Validated, Tunable Lyapunov Exponent Framework for Running Biomechanics"

echo.
echo 4. Setting main branch...
git branch -M main

echo.
set /p REPO_URL="Enter your GitHub repository URL (e.g. https://github.com/YourUsername/LyE-Running-Biomechanics.git): "

if "%REPO_URL%"=="" (
    echo No URL entered. Skipping push. You can push manually later.
    pause
    exit /b
)

git remote remove origin 2>nul
git remote add origin %REPO_URL%

echo.
echo 5. Pushing to GitHub...
git push -u origin main

echo.
echo ====================================================================
echo   SUCCESSFULLY PUSHED TO GITHUB!
echo ====================================================================
echo Now you can copy your repository URL: %REPO_URL%
echo and paste it into docs/LINKEDIN_POST.md before posting on LinkedIn!
echo.
pause
