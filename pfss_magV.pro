pro pfss_magV,date = date,$ ; date of the event yyyymmdd
              hh = hh, $    ; time of pfss (00,06,12 or 18)
              rss  = rss    ; solar surface height
;    
@pfss_data_block  ;  sets pfss common block variables
if not keyword_set(date) then date = '20210101'
if not keyword_set(rss) then rss = 2.5
if not keyword_set(hh) then hh = '00'
YYYY= strmid(date,0,4)
MM =  strmid(date,4,2)
DD =  strmid(date,6,2)
print, date
pfss_restore,pfss_time2file(YYYY+'-'+MM+'-'+DD+'_'+hh+':04:00',/ssw_cat,/url)  ;  gets file


;sfield=pfss_surffield_restore(pfss_time2file(YYYY+'-'+MM+'-'+DD+'_'+hh+':04:00',/ssw_cat,/url,/surffield))  ;  for all users
nlat0=192  ;  number of latitudinal gridpoints in magnetogram
;pfss_mag_create,mag,0,nlat0,file=sfield;,/quiet
;pfss_get_potl_coeffs,mag,rtop=rss



pfss_potl_field,rss


cth=cos(theta)
monopole=(mean_dtheta(total((Br(*,*,0)),1),cth)/nlon)(0)
r2inv=rebin(reform(1./rix^2,1,1,nr),nlon,nlat,nr)
Br=Br-float(monopole*r2inv)  ;  removes monopole from coronal field

save,rix,theta,phi,filename=getenv("HOME")+'/Documents/Events/'+date+'/coordB.sav'
save,Br,Bth,Bph,filename=getenv("HOME")+'/Documents/Events/'+date+'/magvB.sav'

print,'saved'
end