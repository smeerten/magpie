rmdir /S/Q dist
rmdir /S/Q build
call buildMagpie.bat
copy magpie\LICENSE dist\LICENSE
copy magpie\README.md dist\README.md
copy magpie\Manual.pdf dist\Manual.pdf
xcopy /E magpie\samples dist\samples\
xcopy /E magpie\pulseSeqs dist\pulseSeqs\
xcopy /E magpie\shapefiles dist\shapefiles\
copy gpl.rtf dist
copy modern_magpie.nsi dist
"C:\Program Files (x86)\NSIS\makensis.exe" dist\modern_magpie.nsi