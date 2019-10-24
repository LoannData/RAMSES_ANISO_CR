pro fait_plot,printps=printps

rd_1D_Temp,d1,file='run_1Temp_nocond.log',/one
rd_1D_Temp,d1c,file='run_1Temp.log',/one
rd_1D_Temp,d2,file='run_2Temp.log'

n1=n_elements(d1)-1L
n1c=n_elements(d1c)-1L
n2=n_elements(d2)-1L

o1=*d1[n1]
o1c=*d1c[n1c]
o2=*d2[n2]

print,o1.t,o1c.t,o2.t,' Myr'

;  o1=*d1[10]
;  o1c=*d1c[10]
;  o2=*d2[10]

ud=1.66d-24
ut=3.154464d13
ul=3.08d21

ut2=1.66d-24/1.38d-16*ul^2d0/ut^2d0
up=ud*ul^2d0/ut^2d0*1d-1
mug=0.5882353d0
mue=1.1363636d0
mui=1.2195118d0
xoff=0.2d0

if keyword_set(printps)then begin
    set_plot,'ps'
    device,filename='spitzer.eps',/encaps,ysize=30,xsize=40,/color
    !p.font=1
    th=5.
    !p.charthick=th
    !p.charsize=2.5
    !p.thick=th
endif else begin
    set_plot,'x'
    window,0,xs=1536,ys=1024
    th=1
endelse
!p.multi=[0,2,2]
plot,o1.x/1d3,o1.p/o1.d*ut2*mug,xtitle='x (Mpc)',ytitle='T (K)',xth=th,yth=th,th=th,yr=[3d7,1d8],/ys
tek_color
oplot,o1c.x/1d3+xoff,o1c.p/o1c.d*ut2*mug,color=2

oplot,o2.x/1d3,(o2.p+o2.pnt)/o2.d*ut2*mug,color=4
oplot,o2.x/1d3,(o2.pnt)/o2.d*ut2*mue,color=4,lines=2
oplot,o2.x/1d3,(o2.p)/o2.d*ut2*mui,color=4,lines=1

plot,o1.x/1d3,o1.d,xtitle='x (Mpc)',ytitle='n (cm!u-3!n)',xth=th,yth=th,th=th,yr=[9d-5,2.5d-4],/ys
oplot,o1c.x/1d3+xoff,o1c.d,color=2

oplot,o2.x/1d3,o2.d,color=4

plot,o1.x/1d3,o1.p*up,xtitle='x (Mpc)',ytitle='P (Pa)',xth=th,yth=th,th=th
oplot,o1c.x/1d3+xoff,o1c.p*up,color=2

oplot,o2.x/1d3,(o2.p+o2.pnt)*up,color=4
oplot,o2.x/1d3,o2.pnt*up,color=4,lines=2
oplot,o2.x/1d3,o2.p*up,color=4,lines=1

plot,o1.x/1d3,o1.u*ul/ut/1d5,xtitle='x (Mpc)',ytitle='u (km/s)',xth=th,yth=th,th=th
oplot,o1c.x/1d3+xoff,o1c.u*ul/ut/1d5,color=2

oplot,o2.x/1d3,o2.u*ul/ut/1d5,color=4

if keyword_set(printps)then begin
    device,/close
    !p.font=-1
    !p.charthick=1.
    !p.thick=1.
endif

end
