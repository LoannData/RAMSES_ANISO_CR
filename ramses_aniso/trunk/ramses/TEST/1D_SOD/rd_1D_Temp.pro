;+
; NAME:
;       RD_1D_TEMP
;
; PURPOSE:
;       This procedure reads 1D profiles from a RAMSES LOG ASCII file.
;
; CATEGORY:
;       Input/Output.
;
; CALLING SEQUENCE:
;       RD_1D, Data, FILE=file, NMAX=nmax, LMAX=lmax
;
; OPTIONAL INPUTS:
;       FILE:   if set, input the scalar string containing the name of
;       the file to be read. Otherwise, a PICKFILE widget is launched. 
;
;       NMAX:  if set, the maximum number of time steps to be
;       read from the file. Default: 100.
;
;       LAX:   if set, the maximum number of levels to be read from
;       the file. Default: 10.
;
; OUTPUTS:
;       Data: A pointer array of structure containing the hydro
;       variables of the RAMSES 1D run for each output time.
;
; EXAMPLE:
;       To extract profiles from a RAMSES LOG ASCII file, type:
;
;               RD_1D, data, file='run12.log'
;
;       To get the mesh structure for the 3rd output time, type:
;
;               out3 = *data[2]
;
;       The 1D mesh structure is organized as follows:
;
;               out3.t    : time                  (DOUBLE)
;               out3.nc   : number of cells       (INTEGER)
;               out3.np   : number of particles   (INTEGER)
;               out3.x    : cell center position  (DOUBLE ARRAY)
;               out3.l    : level of refinement   (DOUBLE ARRAY)
;               out3.d    : mass density          (DOUBLE ARRAY)
;               out3.u    : velocity              (DOUBLE ARRAY)
;               out3.p    : pressure              (DOUBLE ARRAY)
;               out3.xp   : particle positions    (DOUBLE ARRAY)
;               out3.vp   : particle velocities   (DOUBLE ARRAY)
;               out3.temp : electron temperature  (DOUBLE ARRAY)
;
; MODIFICATION HISTORY:
;       Written by:     Romain Teyssier, 01/01/2000.
;                       e-mail: Romain.Teyssier@cea.fr
;       Fevrier, 2001:  Comments and header added by Romain Teyssier.
;-
;###################################################
;###################################################
;###################################################
pro rd_1d_Temp,data,file=file,part=part,nmax=nmax

if not keyword_set(file) then file=pickfile(/READ)
if not keyword_set(nmax) then nmax=31000L
print,'Reading ',file
print,nmax
output=PTRARR(nmax)

openr,1,file
st=''
i=0
on_ioerror, cycle
while not eof(1) do begin
    ipos=-1
    while ipos eq -1 and not eof(1) do begin
        readf,1,st
        ipos=strpos(st,'Output')
    endwhile
    if not eof(1) then begin
        ncell=fix(strmid(st,ipos+6,10))

        if keyword_set(part) then begin
            readf,1,st
            npart=MAX([fix(strmid(st,ipos+6,10)),1])
        endif
        ipos=-1
        ; Define profile structure
        if keyword_set(part) then begin
            time={t:0.0d0, nc:ncell, np:npart, l:intarr(ncell) $
                  , x:dblarr(ncell), d: dblarr(ncell), u:dblarr(ncell) $
                  , p:dblarr(ncell),temp:dblarr(ncell), xp:dblarr(npart), vp:dblarr(npart) }
        endif else begin
            time={t:0.0d0, nc:ncell, l:intarr(ncell) $
                  , x:dblarr(ncell), d: dblarr(ncell), u:dblarr(ncell) $
                  , p:dblarr(ncell), temp:dblarr(ncell), pnt:dblarr(ncell) }
        endelse
 readf,1,st
 readf,1,st


        for j=0,ncell-1 do begin

            ;; readf,1,ll,rr,dd,uu,vv,ww,pp,aa,bb,cc;,ee
            ;; ppnt=pp
            readf,1,ll,rr,dd,uu,vv,ww,ppnt,pp;,aa,bb,cc;,ee
;print,ll,rr,dd,uu,pp,ipos
            time.l(j)=ll
            time.x(j)=rr
            time.d(j)=dd
            time.p(j)=pp
            time.u(j)=uu
            time.pnt(j)=ppnt
;            time.temp(j)=tt         
        endfor
        if keyword_set(part) then begin
            for j=0,npart-1 do begin
                readf,1,xx,vv
                time.xp(j)=xx
                time.vp(j)=vv
            endfor
        endif
        ipos=-1
        while ipos eq -1 and not eof(1) do begin
            readf,1,st
            ipos=strpos(st,'Fine step')
        endwhile
        if(not eof(1))then begin
            t=double(strmid(st,ipos+20,21))
            time.t=t
            pc=ptr_new(time)
            output(i)=pc
            i=i+1
        endif
        cycle:
    endif
endwhile
close,1

if(i le 0) then begin
    print,'No valid time read'
    return
endif

ntime=i
print,i
data=output(0:ntime-1)

end
;#########################################################
;#########################################################
;#########################################################
