dset ^ModelWD.bin
*OPTIONS  big_endian          
undef  -99999.0
title GPCP MODEL 1.0  degree
xdef   93 linear   -58   0.25
ydef   61 linear   -35   0.25
tdef   96  LINEAR  0024mar2004  1hr
zdef    37 levels 1000 975 950 925 900 875 850 825 800 775 
                   750 700 650 600 550 500 450 400 350 300 
                   250 225 200 175 150 125 100  70  50  30 
                    20  10   7   5   3   2   1             
VARS 5
thetas    0   99    Temperatura potencial na superficie
ug_sfc    0   99    Velocidade geostrofica zonal na superficie
vg_sfc    0   99    Velocidade geostrofica meridional na superficie
ug_LND    0   99    Velocidade geostrofica zonal no level no divegence
vg_LND    0   99    Velocidade geostrofica meridional no level no divegence
ENDVARS
