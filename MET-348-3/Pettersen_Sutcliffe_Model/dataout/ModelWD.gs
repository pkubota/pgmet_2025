'reinit'
'open ModelSC.ctl'
'open ModelWD.ctl'
'set display color white'
'clear'
xmin=0.5
xmax=10.0
ymin=1.0
ymax=7.0
'set parea 'xmin' 'xmax' 'ymin' 'ymax''
'set gxout shaded'

'run rgbset.gs'
'set gxout shaded'
it=1
while(it<=96)
'set t 'it''
'set lat -34 -21'
'set lon -57 -36'

'set clevs   -0.6    -0.5 -0.4  -0.3    0.3   0.4   0.5  0.6'
'set ccols 48    46  44  42     0      22   24   26   28'
'aa= terma.1+termb.1+termc.1+termd.1'
*'aa= terma.1+termc.1+termd.1'
*'aa= terma.1'
*'aa= termb.1'
*'aa= termc.1'
*'aad= termd.1'
'd aa*1e12'
*'d aad*1e14'
'cbarn_local.gs'

'set gxout stream'
'd ug_sfc.2;vg_sfc.2;mag(ug_LND.2,vg_LND.2)'

'!sleep 0.5'
'clear'
it=it+1

endwhile
