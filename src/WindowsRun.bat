where /q pyw
IF ERRORLEVEL 1 (
    ECHO 'pyw' is missing, try 'pythonw'
    where /q pythonw
    IF ERRORLEVEL 1 (
        ECHO 'pythonw' also missing. Abort.
    ) ELSE (
       start pythonw magpie.py
    )
) ELSE (
    start pyw magpie.py
)



