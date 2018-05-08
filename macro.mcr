#!MC 1410
$!READDATASET ' "tecplotfiles/out-tecplot-CC-0.2000.plt"  "tecplotfiles/out-tecplot-F1-0.2000.plt"  "tecplotfiles/out-tecplot-F2-0.2000.plt" '
READDATAOPTION = NEW
RESETSTYLE = YES
VARLOADMODE = BYNAME
ASSIGNSTRANDIDS = YES
VARNAMELIST = '"x" "y" "f" "fdrop" "u.x" "u.y" "omega" "pressure"'
$!WRITEDATASET  "tecplotfiles/out-bin-tecplot-CC-0.2000.plt"
INCLUDETEXT = NO
INCLUDEGEOM = NO
INCLUDEDATASHARELINKAGE = YES
BINARY = YES
USEPOINTFORMAT = NO
PRECISION = 9
TECPLOTVERSIONTOWRITE = TECPLOTCURRENT