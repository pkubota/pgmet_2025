dset ^ZonalWind.bin
*OPTIONS  big_endian          
undef  -99999.0
title GPCP MODEL 1.0  degree
xdef   93 linear   -58   0.25
ydef   61 linear   -35   0.25
tdef   120  LINEAR  00Z24MAR2004  1hr
zdef    37 levels 1000 975 950 925 900 875 850 825 800 775 
                   750 700 650 600 550 500 450 400 350 300 
                   250 225 200 175 150 125 100  70  50  30 
                    20  10   7   5   3   2   1             
VARS 1
uvel    37   99    U component of wind         (m/sec )    
ENDVARS
