@echo off

build\Release\qd_test_c.exe
IF %ERRORLEVEL% NEQ 0 pause

build\Release\qd_test_huge.exe
IF %ERRORLEVEL% NEQ 0 pause

build\Release\qd_test_more.exe
IF %ERRORLEVEL% NEQ 0 pause

build\Release\qd_test_pslq.exe
IF %ERRORLEVEL% NEQ 0 pause

build\Release\qd_test_qd.exe
IF %ERRORLEVEL% NEQ 0 pause

build\Release\qd_test_quadt.exe
IF %ERRORLEVEL% NEQ 0 pause

build\Release\qd_bench.exe
IF %ERRORLEVEL% NEQ 0 pause

echo All passed
pause
