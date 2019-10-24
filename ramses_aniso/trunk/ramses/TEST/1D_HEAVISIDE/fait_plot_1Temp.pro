pro fait_plot
  
  ;rd_1d,file='res_1Temp.log',out    
  rd_1d_temp,file='res1d_cr.log',out    
  ana=out

nout=7

er_tot = fltarr(nout)
time   =  fltarr(nout)
tout   =  fltarr(nout)
tout(0)=1
tout(1)=2
tout(2)=nout-1
xx=findgen(100)/200.0

for iout=0,nout-1 do begin
   dat = *out[iout]

   er_tot(iout)=total(dat.pnt*0.5^dat.l)
  
   time(iout)=dat.t
endfor


set_plot,'x'
window,0
;ana=fltarr(100)
set_plot,'ps'
device,filename='heaviside.eps',/encaps,yoffset=5.3,ysize=17,xsize=20,/color
!p.font=0
!p.charthick=2.5
!p.charsize=1.5
!p.thick=2.5

!p.position=[0.15,0.4,0.9,0.95]
iplot=0
for iout=0,2 do begin
   iiout=(tout(iout))
   ana = 0.6d0+0.2d0*erf((xx-0.25)/(4.0d0*time(iiout))^0.5)
   ana = 6.d0+2d0*erf((xx-0.25)/(4.0d0*time(iiout))^0.5)


   dat=*out[iiout]   

   if(iout eq 0)then begin
      plot,xx,ana,xr=[0,1],yr=[2,10],/noerase,xtickn=[' ',' ',' ',' ',' ',' '],/ynozero,ystyle=8,ytitle=textoidl("P")
      ;; plot,xx,ana,xr=[0,1],yr=[0.2,1],/noerase,xtickn=[' ',' ',' ',' ',' ',' '],/ynozero,ystyle=8,ytitle=textoidl("P")
      oplot,1.0d0-xx,ana
      oplot,dat.x,dat.pnt,psym=1
      oplot,dat.x,8/5.*(dat.l-3)+2,line=2
      ;; oplot,dat.x,0.8/3.*(dat.l-4)+0.2,line=2
print,dat.l

   endif else begin
      oplot,xx,ana,color=2*iout
   oplot,1.0d0-xx,ana,color=2*iout

   oplot,dat.x,dat.pnt,psym=2*iout+2,color=2*iout
   
;   oplot,dat.x,0.8/3.*(dat.l-4)+0.2,line=2,color=2*iout
   oplot,dat.x,8./5.*(dat.l-3)+2,line=2,color=2*iout

print,dat.l
endelse
print,dat.t
plots,[0.48,0.52],[8,8]
plots,[0.48,0.52],[8,8],color=2
plots,[0.48,0.52],[7.95,7.95],color=4
;; oplot,out_dir.x,1.2*(out_dir.l-5)*1.25,line=1

endfor
axis,yaxis=1,ytitle='!17 AMR level',yr=[0,5]

!p.position=[0.15,0.1,0.9,0.4]
iplot=1
for iout=0,2 do begin
   iiout=(tout(iout))

   dat=*out[iiout]

   for icell=0,dat.nc-1 do begin
      if(dat.x(icell) le 0.5)then begin
         dat.d(icell)=6d0+2d0*erf((dat.x(icell)-0.25)/(4.0d0*dat.t)^0.5)
      endif else begin
         dat.d(icell)=6d0+2d0*erf((0.75-dat.x(icell))/(4.0d0*dat.t)^0.5)
      endelse
   endfor

      if(iout eq 0)then plot,dat.x,alog10(abs(dat.pnt-dat.d)/dat.d),xtitle='x',ytitle='log(relative error)',/noerase,yr=[-4.8,-1.1],/ys; ,ytickn=['10!u-6',' ','10!u-4',,' ','10!u-2',' ',' '] ;,ytickn=[' ','10!u-6','10!u-4','10!u-2',' ']
   if(iout ne 0)then  oplot,dat.x,alog10(abs(dat.pnt-dat.d)/dat.d),color=iout*2
print,abs(dat.pnt-dat.d)/dat.d
print,dat.pnt
print,dat.d
print,'   '
endfor

device,/close

;; iplot=0

;; plot,xx,alog10(ana2),line=2,ytitle=textoidl("log(E_r)"),/ynozero,ystyle=8,yr=[0,6]


;; device,filename='dirac.eps',/encaps,yoffset=5.3,ysize=17,xsize=20,/color
;; out_ana=*ana[2]
;; out_dir=*dir[8]
;; out_rob=*rob[8]
;; out_vn=*vn[8]

;; xx=findgen(100)/100.
;; chi=1.d10;/7.5d0
;; x0=0.5
;; E0=1.d5
;; tt=out_ana.t
;; anac = E0/(2.0d0*(chi*3.1415*tt)^.5)
;; ana2 = 1.+anac*exp(-(xx-x0)^2./(4.0d0*tt*chi));^0.5)

;; out_dir.p=1.+anac*exp(-(out_dir.x-x0)^2./(4.0d0*tt*chi));^0.5)
;; out_rob.p=1.+anac*exp(-(out_rob.x-x0)^2./(4.0d0*tt*chi));^0.5)
;; out_vn.p=1.+anac*exp(-(out_vn.x-x0)^2./(4.0d0*tt*chi));^0.5)

;; out_dir.p=abs(out_dir.p-out_dir.er)/out_dir.p
;; out_rob.p=abs(out_rob.p-out_rob.er)/out_rob.p
;; out_vn.p=abs(out_vn.p-out_vn.er)/out_vn.p


;; ;oplot,xx,alog10(ana2),psym=6

;; ;; plot,out_ana.x,alog10(out_ana.er),line=1,ytitle=textoidl("log(E_r)"),xtitle='x',/ynozero,ystyle=8

;; !p.position=[0.15,0.4,0.9,0.95]
;; iplot=0

;; plot,xx,alog10(ana2),line=2,ytitle=textoidl("log(E_r)"),/ynozero,ystyle=8,yr=[0,6],/noerase,xtickn=[' ',' ',' ',' ',' ',' ']

;; oplot,out_dir.x,alog10(out_dir.er),color=0
;; oplot,out_rob.x,alog10(out_rob.er),color=2
;; oplot,out_vn.x,alog10(out_vn.er),color=4

;; axis,yaxis=1,ytitle='!17 AMR level',yr=[0,4]
;; oplot,out_dir.x,1.2*(out_dir.l-5)*1.25,line=1

;; ;xyouts,0.05,0.3,'!17(a)'

;; t1=textoidl("Analytic")
;; t2=textoidl("\alpha=1")
;; t3=textoidl("\alpha=0.5")
;; t4=textoidl("\alpha=0")
;; t5=textoidl("AMR level")
;; legend,[t5,t1,t2,t3,t4],/fill,lines=[1,2,0,0,0],color=[0,0,0,2,4],box=0,spacing=1., /top, /right, thick=[2.5,2.5,2.5,2.5,2.5],linsize=0.3

;; !p.position=[0.15,0.1,0.9,0.4]
;; iplot=1
;; plot_io,out_dir.x,out_dir.p,xtitle='x',ytitle='Relative error',/noerase,yr=[1.d-5,10.],ytickn=['10!u-5',' ','10!u-3',' ','10!u-1',' ','10'];,ytickn=[' ','10!u-6','10!u-4','10!u-2',' ']
;; oplot,out_rob.x,out_rob.p,color=2,line=0
;; oplot,out_vn.x,out_vn.p,color=4

;; device,/close
;; print,out_dir.l

;; print,time_ana,er_ana
;; print,time(*,0),er_tot(*,0)
;; print,time(*,1),er_tot(*,1)
;; print,time(*,2),er_tot(*,2)
end
