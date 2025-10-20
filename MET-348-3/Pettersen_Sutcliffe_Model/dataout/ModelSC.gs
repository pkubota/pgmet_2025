'reinit'
'open ModelA.ctl'
'set display color white'
'clear'
xmin=1.0
xmax=10.5
ymin=1.0
ymax=7.0
'set parea 'xmin' 'xmax' 'ymin' 'ymax''
'set gxout shaded'

'run rgbset.gs'
'set gxout shaded'
it=1
while(it<=20)
'set t 'it''
'set lat -34 -21'
'set lon -57 -36'

'set clevs   -0.6    -0.5 -0.4  -0.3    0.3   0.4   0.5  0.6'
'set ccols 48    46  44  42     0      22   24   26   28'
'd terma*1e11'
'cbarn_local.gs'
'!sleep 0.5'
'clear'
it=it+1

endwhile
