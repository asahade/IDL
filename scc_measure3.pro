;+
; $Id: scc_measure3.pro, 2024/04/24 sahade Exp $
;
; Project     :	STEREO - SECCHI
;
; Name        :	SCC_MEASURE3
;
; Purpose     :	Measure 3D coordinates from three STEREO images
;
; Category    :	STEREO, SECCHI, Coordinates
;
; Explanation :	Widget application to allow a user to select with the cursor a
;               feature appearing in three images.  The selected coordinates are
;               used to locate the point in three-dimensional space.
;
;               The user first selects a point on one image.  The program then
;               displays a line on the other two images representing the
;               line-of-sight from the first image.  The user then selects the
;               point along this line with the same feature, in any of the other images.
;               The 3D coordinates are then calculated and displayed within the widget
;               as (Earth-based) Stonyhurst heliographic longitude and
;               latitude, along with radial distance in solar radii.
;
;               In essence, this is a demonstration procedure showing how to
;               derive 3D information from multi-viewpoint data.
;               Tested with:  EUVI A & B, EIT, AIA and EUI
;
; Syntax      :	SCC_MEASURE3, IMAGE_AHEAD, IMAGE_BEHIND, IMAGE_OTHER, INDEX_AHEAD, INDEX_BEHIND, INDEX_OTHER
; 
; 
; or
; 
;               SCC_MEASURE, FILE_AHEAD, FILE_BEHIND, FILE_OTHER
;               Only works if Other can be read by sccreadfits

;
; Inputs      :	IMAGE_AHEAD  = Image for the Ahead observatory
;               IMAGE_BEHIND = Image for the Behind observatory
;               IMAGE_OTHER  = Image for the Other observatory
;               INDEX_AHEAD  = FITS index structure for Ahead
;               INDEX_BEHIND = FITS index structure for Behind
;               INDEX_OTHER  = FITS index structure for Other
;
;               Prefered method
;               
;               Not tested:
; Opt. Inputs : FILE_AHEAD  = Filename for the Ahead observatory.
;               FILE_BEHIND = Filename for the Behind observatory
;               FILE_OTHER  = Filename for other observatory
;
; ;               
;
; Outputs     :	None.
;
; Opt. Outputs:	None.
;
; Keywords    :	WSIZE   = Window size to use for displaying the images.  The
;                         default is 512.
;
;               OUTFILE = Name of file to store results in.  Use the "Store"
;                         button to write each selected measurement into this
;                         file.  The stored points will also be displayed on
;                         the images.  Use the "Clear stored" button to clear
;                         these points out from memory.  At the same time, a
;                         blank line will be written to the output file to mark
;                         the separation between collections of points.
;
;               APPEND  = Append to existing output file.
;
;               FORCESAVE = If set, saves last results in given output
;                           file upon exit, even if manual 'Save' has
;                           or has not been used.  Has no effect if no
;                           output file was defined.
;
;               CROP = If set, restricts chosing of points to the
;                      inner region so that a crop to half size
;                      can later be done (outside of scc_measure)
;
;               NO_BLOCK = Call XMANAGER with /NO_BLOCK
;
;   	    	/DEBUG	print diagnostic information
;
;
; Calls       :	SCCREADFITS, FITSHEAD2WCS, EXPTV, SCC_MEASURE3_EVENT,
;               SCC_MEASURE3_SELECT_1ST, SCC_MEASURE3_SELECT_2ND, BOOST_ARRAY,
;               SCC_MEASURE3_REDISPLAY, SSC_MEASURE3_REPLOT, SSC_MEASURE3_CLEANUP,
;               SETWINDOW, GET_TV_SCALE, WCS_GET_COORD, CONVERT_SUNSPICE_COORD
;
; Common      :	SCC_MEASURE3 = Internal common block used by the widget program.
;
; Restrictions:	None.
;
; Side effects:	None.
;
;
; Adapted from SCC_MEASURE (William Thompson)
; 
; Written     :	Version 1, 26-Jan-2023, Abril Sahade

;-
;
;==============================================================================
;
pro scc_measure3_redisplay, window=k_window
;
;  Procure to redisplay the windows after various operations such as changing
;  the color table.
;
common scc_measure3, ahead, h_ahead, behind, h_behind, eother, h_eother, $
  wcs_ahead, wcs_behind, wcs_eother, win_left, win_right, win_other, $
  image_left, image_right, image_other, lon_wid, lat_wid, rsun_wid, $
  zoom_wid, unzoom_wid, maxz, pix_paraml, pix_paramr, pix_paramo, $
  pix_left, pix_right, pix_other, heeq_left, heeq_right, heeq_other, $
  in_progress, win_last, lzoom, rzoom, ozoom, lfrz_wid, rfrz_wid, $
  ofrz_wid, color, color1, color2, lmin, lmax, rmin, rmax, omin, omax, $
  lmin_wid, lmax_wid, rmin_wid, rmax_wid, omin_wid, omax_wid, $
  subimage_left, subimage_right, subimage_other, $
  origin_left, origin_right, origin_other, outlun, out_wid, $
  clear_wid, oplot_wid, n_stored, store_left, store_right, store_other, $
  pixpair, pixforce, roi, debugon

;
if n_elements(k_window) eq 1 then window = k_window else window = -1
;
;  If the left (or no) window was selected, then redisplay the left window.
;
if (window eq win_left) or (window eq -1) then begin
    setwindow, win_left
    get_tv_scale, sx, sy, mx, my, jx, jy
    tv, image_left, jx, jy
    if n_stored gt 0 then $
      plots, store_left[0,*], store_left[1,*], psym=1, color=color
endif
;
;  If the right (or no) window was selected, then redisplay the right window.
;
if (window eq win_right) or (window eq -1) then begin
    setwindow, win_right
    get_tv_scale, sx, sy, mx, my, jx, jy
    tv, image_right, jx, jy
    if n_stored gt 0 then $
      plots, store_right[0,*], store_right[1,*], psym=1, color=color
endif
;
;;  If the other (or no) window was selected, then redisplay the right window.
;
if (window eq win_other) or (window eq -1) then begin
  setwindow, win_other
  get_tv_scale, sx, sy, mx, my, jx, jy
  tv, image_other, jx, jy
  if n_stored gt 0 then $
    plots, store_other[0,*], store_other[1,*], psym=1, color=color
endif
;
;  If no particular window was selected, and some points have been selected,
;  then recreate the plots on top of the image.
;
if (window eq -1) and (win_last ne '') then begin
;
;  If the point selection is in progress, then plot the point on the 1st
;  window, and the line on the 2nd and other window.
;
    if in_progress then begin
        if win_last eq 'LEFT' then begin
            setwindow, win_left
            plots, pix_left[0],  pix_left[1],  psym=1, symsize=3, color=color
            setwindow, win_right
            plots, pix_paramr[0,*], pix_paramr[1,*], color=color
            x = [0,wcs_behind.naxis[1]-1]
            paramr = poly_fit(pix_paramr[0,*], pix_paramr[1,*], 1)
            plots, x, poly(x, paramr), line=1, color=color
            setwindow, win_other
            plots, pix_paramo[0,*], pix_paramo[1,*], color=color
            x2 = [0,wcs_eother.naxis[1]-1]
            paramo = poly_fit(pix_paramo[0,*], pix_paramo[1,*], 1)
            plots, x2, poly(x2, paramo), line=1, color=color
        end else if win_last eq 'RIGHT' then begin
            setwindow, win_right
            plots, pix_right[0],  pix_right[1],  psym=1, symsize=3, color=color
            setwindow, win_left
            plots, pix_paraml[0,*], pix_paraml[1,*], color=color
            x = [0,wcs_ahead.naxis[1]-1]
            paraml = poly_fit(pix_paraml[0,*], pix_paraml[1,*], 1)
            plots, x, poly(x, paraml), line=1, color=color
            setwindow, win_other
            plots, pix_paramo[0,*], pix_paramo[1,*], color=color
            x2 = [0,wcs_eother.naxis[1]-1]
            paramo = poly_fit(pix_paramo[0,*], pix_paramo[1,*], 1)
            plots, x2, poly(x2, paramo), line=1, color=color
        end else begin
          setwindow, win_other
          plots, pix_other[0],  pix_other[1],  psym=1, symsize=3, color=color
          setwindow, win_left
          plots, pix_paraml[0,*], pix_paraml[1,*], color=color
          x = [0,wcs_ahead.naxis[1]-1]
          paraml = poly_fit(pix_paraml[0,*], pix_paraml[1,*], 1)
          plots, x, poly(x, paraml), line=1, color=color
          setwindow, win_right
          plots, pix_paramr[0,*], pix_paramr[1,*], color=color
          x2 = [0,wcs_behind.naxis[1]-1]
          paramr = poly_fit(pix_paramr[0,*], pix_paramr[1,*], 1)
          plots, x2, poly(x2, paramr), line=1, color=color
        endelse
;
;  Otherwise, plot both points.
;
    end else begin
        setwindow, win_left
        plots, pix_left[0],  pix_left[1],  psym=1, symsize=3, color=color
        setwindow, win_right
        plots, pix_right[0], pix_right[1], psym=1, symsize=3, color=color
        setwindow, win_other
        ;x = [0,wcs_eother.naxis[1]-1]
        ;plots, pix_param2[0,*], pix_param2[1,*], color=color
        ;plots, x, poly(x, param2), line=1, color=color
        ;plots, pix_param3[0,*], pix_param3[1,*], color=color1
        ;plots, x, poly(x, param3), line=1, color=color1
        plots, pix_other[0], pix_other[1], psym=1, symsize=3, color=color ; 
    endelse
endif
;
end
;
;==============================================================================
;
pro scc_measure3_oplot_from_file ;without modifications from 2img version
;
;  Procedure to overplot traces already saved in the output file.
;
common scc_measure3
;
;  Read in the file.
;
infile = (fstat(outlun)).name
openr, inlun, infile, /get_lun
delvarx, table, iseg
line = 'string'
noffset0 = 0
values = dblarr(7)
while not eof(inlun) do begin
    readf, inlun, line
;
;  Only use lines containing 7 numbers.
;
    on_ioerror, finish_segment
    reads, line, values
    boost_array, table, values
    goto, next_line
;
;  Otherwise, end the segment.
;
finish_segment:
    noffset = n_elements(table) / 7
    if noffset ne noffset0 then begin
        boost_array, iseg, noffset
        noffset0 = noffset
    endif
;
next_line:
;
endwhile
noffset = n_elements(table) / 7
if noffset ne noffset0 then boost_array, iseg, noffset
free_lun, inlun
;
;  Overplot the points in the left and right windows
;
for j=0,n_elements(iseg)-1 do begin
    if j eq 0 then i0=0 else i0 = iseg[j-1]
    i1 = iseg[j] - 1
    setwindow, win_left
    plots, table[5,i0:i1], table[6,i0:i1], psym=0, color=color
    setwindow, win_right
    plots, table[3,i0:i1], table[4,i0:i1], psym=0, color=color
endfor
;
end
;
;==============================================================================
;
pro scc_measure3_replot, left=left, right=right, other=other
;
;  Procedure to replot the images, e.g. after changing the image range.
;
common scc_measure3
;
;  Unless only the right image was selected, then plot the left image.  Get the
;  current min and max values.
;
if not keyword_set(right) then begin
  if not keyword_set(other) then begin
    setwindow, win_left
    widget_control, lmin_wid, get_value=lmin1
    widget_control, lmax_wid, get_value=lmax1
    if lmax1 gt lmin1 then begin
        lmin = lmin1
        lmax = lmax1
    endif
    widget_control, lmin_wid, set_value=lmin
    widget_control, lmax_wid, set_value=lmax
    exptv, subimage_left, /data, /nobox, /noexact, origin=origin_left, $
      min=lmin, max=lmax, bscaled=image_left
  endif
endif
;
;  Unless only the left image was selected, then plot the right image.  Get the
;  current min and max values.
;
if not keyword_set(left) then begin
  if not keyword_set(other) then begin
    setwindow, win_right
    widget_control, rmin_wid, get_value=rmin1
    widget_control, rmax_wid, get_value=rmax1
    if rmax1 gt rmin1 then begin
        rmin = rmin1
        rmax = rmax1
    endif
    widget_control, rmin_wid, set_value=rmin
    widget_control, rmax_wid, set_value=rmax
    exptv, subimage_right, /data, /nobox, /noexact, origin=origin_right, $
      min=rmin, max=rmax, bscaled=image_right
  endif
endif
;  Copied for other
;  
;
if not keyword_set(left) then begin
  if not keyword_set(right) then begin
    setwindow, win_other
    widget_control, omin_wid, get_value=omin1
    widget_control, omax_wid, get_value=omax1
    if omax1 gt omin1 then begin
      omin = omin1
      omax = omax1
    endif
    widget_control, omin_wid, set_value=omin
    widget_control, omax_wid, set_value=omax
    exptv, subimage_other, /data, /nobox, /noexact, origin=origin_other, $
      min=omin, max=omax, bscaled=image_other
  endif
endif
;
;  Redisplay any overplots.
;
scc_measure3_redisplay
;
end
;
;==============================================================================
;
pro scc_measure3_select_1st, ev, win_1st, win_2nd, win_3rd,$
                            image_1st, image_2nd, image_3rd, $
                            wcs_1st, wcs_2nd, wcs_3rd, sc_1st, sc_2nd, sc_3rd, $
                            heeq_1st, pix_param_2nd, pix_param_3rd, pix_1st
;
;  Procedure to select the first point.
;
common scc_measure3
;
;  Start by converting from device pixels to pixels within the image.  Overplot
;  the selected point.
;
scc_measure3_redisplay, window=win_1st
setwindow, win_1st
pix = convert_coord(ev.x, ev.y, /device, /to_data)
pix_1st = pix[0:1]

if roi[9] then begin            ;ROI active
    if (pix_1st[0] lt roi[0] or pix_1st[0] gt roi[1] or $
        pix_1st[1] lt roi[2] or pix_1st[1] gt roi[3]) then begin
        ; out of our region of interest, force an early routine
        roi[8]=0                ; failed
        return
    end else begin
        roi[8]=1                ; success
    endelse
endif

if (sc_1st eq 'ahead')  then pixpair[0:1]=pix_1st 
if (sc_1st eq 'behind') then pixpair[2:3]=pix_1st
if (sc_1st eq 'eother') then pixpair[4:5]=pix_1st
if (pixforce ne 0) then pixforce=2; mark as 'needs to be saved'

plots, pix_1st[0], pix_1st[1], psym=1, symsize=3, color=color
;
;  Convert from image pixel into helioprojective-cartesian coordinates, in
;  radians.
;
coord = wcs_get_coord(wcs_1st, pix_1st)
conv = !dpi / 180.d0
case wcs_1st.cunit[0] of
    'arcmin': conv = conv / 60.d0
    'arcsec': conv = conv / 3600.d0
    'mas':    conv = conv / 3600.d3
    'rad':    conv = 1.d0
    else:     conv = conv
endcase
coord = coord * conv
;
;  Calculate the equivalent heliocentric coordinates for distances D within
;  +/- maxz of Dsun.
;
dsun = wcs_1st.position.dsun_obs
d = dsun + maxz * [-1,1]
cosy = cos(coord[1])
x = d * cosy * sin(coord[0])
y = d * sin(coord[1])
z = dsun - d * cosy * cos(coord[0])
;
;  Determine the spacecraft parameter to pass to convert_sunspice_coord.
;
spacecraft = sc_1st
test = execute('header = h_'+spacecraft)
obsrvtry = 'Earth'
if datatype(header) eq 'STC' then begin
    if tag_exist(header, 'OBSRVTRY') then begin
        obsrvtry = header.obsrvtry
    end else if tag_exist(header, 'TELESCOP') then begin
        obsrvtry = header.telescop
    endif
end else begin
    temp = fxpar(header, 'OBSRVTRY', count=count)
    if count gt 0 then obsrvtry = temp else begin
        temp = fxpar(header, 'TELESCOP', count=count)
        if count gt 0 then obsrvtry = temp
    endelse
endelse
spacecraft = parse_sunspice_name(obsrvtry, /earth_default)
if wcs_1st.position.soho and (not wcs_1st.position.pos_assumed) then $
  spacecraft = 'SOHO'

IF debugon THEN print,'Selected ',coord,' ',wcs_1st.cunit[0]
if debugon THEN help,dsun,spacecraft
;
;  Convert to HEEQ coordinates, with rearranging into HGRTN format as an
;  intermediate state.
;
coord = transpose([[z],[x],[y]])
convert_sunspice_coord, wcs_1st.time.observ_date, coord, 'HGRTN', 'HEEQ', $
  spacecraft=spacecraft
heeq_1st = coord
;
;  Determine the spacecraft parameter to pass to convert_sunspice_coord.
;
spacecraft = sc_2nd
test = execute('header = h_'+spacecraft)
obsrvtry = 'Earth'
if datatype(header) eq 'STC' then begin
  if tag_exist(header, 'OBSRVTRY') then begin
    obsrvtry = header.obsrvtry
  end else if tag_exist(header, 'TELESCOP') then begin
    obsrvtry = header.telescop
  endif
end else begin
  temp = fxpar(header, 'OBSRVTRY', count=count)
  if count gt 0 then obsrvtry = temp else begin
    temp = fxpar(header, 'TELESCOP', count=count)
    if count gt 0 then obsrvtry = temp
  endelse
endelse
spacecraft = parse_sunspice_name(obsrvtry, /earth_default)
if wcs_2nd.position.soho and (not wcs_2nd.position.pos_assumed) then $
  spacecraft = 'SOHO'
;
;  Switch to the other window, and convert from HEEQ to heliocentric-cartesian
;  x,y,z values for the other perspective.
;
scc_measure3_redisplay, window=win_2nd
coord2=coord
convert_sunspice_coord, wcs_2nd.time.observ_date, coord2, 'HEEQ', 'HGRTN', $
  spacecraft=spacecraft
x = reform(coord2[1,*])
y = reform(coord2[2,*])
z = reform(coord2[0,*])
;
;  Convert from heliocentric-cartesian to helioprojective-cartesian.  Put into
;  the target units.
;
dsun = wcs_2nd.position.dsun_obs
d = sqrt(x^2 + y^2 + (dsun-z)^2)
coord2 = transpose([[atan(x,dsun-z)], [asin(y/d)]])
conv = 180.d0 / !dpi
case wcs_2nd.cunit[0] of
  'arcmin': conv = conv * 60.d0
  'arcsec': conv = conv * 3600.d0
  'mas':    conv = conv * 3600.d3
  'rad':    conv = 1.d0
  else:     conv = conv
endcase
coord2 = coord2 * conv
IF debugon THEN print,'Translating to ',coord2,' ',wcs_2nd.cunit[0]
IF debugon THEN help,dsun,spacecraft
;
;  Convert from helioprojective-cartesian to image pixel coordinates, and
;  overplot the points.
;
pix_param_2nd = wcs_get_pixel(wcs_2nd, coord2)
plots, pix_param_2nd[0,*], pix_param_2nd[1,*], color=color
;
;  Extrapolate over the full X-range of the image.
;
param_2nd = poly_fit(pix_param_2nd[0,*], pix_param_2nd[1,*], 1)
x = [0,wcs_2nd.naxis[1]-1]
plots, x, poly(x, param_2nd), line=1, color=color
;
;  Determine the spacecraft parameter to pass to convert_sunspice_coord.
;
spacecraft = sc_3rd
test = execute('header = h_'+spacecraft)
obsrvtry = 'Earth'
if datatype(header) eq 'STC' then begin
  if tag_exist(header, 'OBSRVTRY') then begin
    obsrvtry = header.obsrvtry
  end else if tag_exist(header, 'TELESCOP') then begin
    obsrvtry = header.telescop
  endif
end else begin
  temp = fxpar(header, 'OBSRVTRY', count=count)
  if count gt 0 then obsrvtry = temp else begin
    temp = fxpar(header, 'TELESCOP', count=count)
    if count gt 0 then obsrvtry = temp
  endelse
endelse
spacecraft = parse_sunspice_name(obsrvtry, /earth_default)
if wcs_3rd.position.soho and (not wcs_3rd.position.pos_assumed) then $
  spacecraft = 'SOHO'
;
;  Switch to the other window, and convert from HEEQ to heliocentric-cartesian
;  x,y,z values for the other perspective.
;
scc_measure3_redisplay, window=win_3rd
coord3=coord
convert_sunspice_coord, wcs_3rd.time.observ_date, coord3, 'HEEQ', 'HGRTN', $
  spacecraft=spacecraft
x = reform(coord3[1,*])
y = reform(coord3[2,*])
z = reform(coord3[0,*])
;
;  Convert from heliocentric-cartesian to helioprojective-cartesian.  Put into
;  the target units.
;
dsun = wcs_3rd.position.dsun_obs
d = sqrt(x^2 + y^2 + (dsun-z)^2)
coord3 = transpose([[atan(x,dsun-z)], [asin(y/d)]])
conv = 180.d0 / !dpi
case wcs_3rd.cunit[0] of
  'arcmin': conv = conv * 60.d0
  'arcsec': conv = conv * 3600.d0
  'mas':    conv = conv * 3600.d3
  'rad':    conv = 1.d0
  else:     conv = conv
endcase
coord3 = coord3 * conv
IF debugon THEN print,'Translating to ',coord3,' ',wcs_3rd.cunit[0]
IF debugon THEN help,dsun,spacecraft
;
;  Convert from helioprojective-cartesian to image pixel coordinates, and
;  overplot the points.
;
pix_param_3rd = wcs_get_pixel(wcs_3rd, coord3)
plots, pix_param_3rd[0,*], pix_param_3rd[1,*], color=color
;
;  Extrapolate over the full X-range of the image.
;
param_3rd = poly_fit(pix_param_3rd[0,*], pix_param_3rd[1,*], 1)
x3 = [0,wcs_3rd.naxis[1]-1]
plots, x3, poly(x3, param_3rd), line=1, color=color
;
;  Deactivate the zoom button, and clear out the text widgets.
;
widget_control, zoom_wid, sensitive=0
widget_control, out_wid,  sensitive=0
widget_control, lon_wid, set_value=''
widget_control, lat_wid, set_value=''
widget_control, rsun_wid, set_value=''
;
end
;
;==============================================================================
;
pro scc_measure3_select_2nd, ev, win_2nd, win_3rd, image_2nd, image_3rd,$
                             wcs_2nd, wcs_3rd, sc_2nd, sc_3rd, $
                             heeq_1st, heeq_2nd, pix_param_2nd, pix_param_3rd, $
                             pix_2nd, pix_3rd
;
;  Procedure to select the 2nd point.
;
common scc_measure3
;
;  Start by converting from device pixels to pixels within the image.
;
scc_measure3_redisplay, window=win_2nd
pix = convert_coord(ev.x, ev.y, /device, /to_data)
pix = pix[0:1]

;  Find the intersection between the linear fit from the other window
;  (PARAM_1ST) and the orthogonal line represented by the current selection.
;  Overplot the selection.
;
param_2nd = poly_fit(pix_param_2nd[0,*], pix_param_2nd[1,*], 1)
a0 = param_2nd[0]
a1 = param_2nd[1]
x = (pix[0] + a1*(pix[1]-a0)) / (1 + a1^2)
pix_2nd = [x, a0 + a1*x]
;
; Now check our region of interest again
;
if roi[9] then begin            ;ROI active
    if (pix_2nd[0] lt roi[4] or pix_2nd[0] gt roi[5] or $
        pix_2nd[1] lt roi[6] or pix_2nd[1] gt roi[7]) then begin
        ; out of our region of interest, force an early routine
        roi[8]=0                ; failed
        return
    end else begin
        roi[8]=1                ; success
    end
endif
if (sc_2nd eq 'ahead')  then pixpair[0:1]=pix_2nd
if (sc_2nd eq 'behind') then pixpair[2:3]=pix_2nd
if (sc_2nd eq 'eother') then pixpair[4:5]=pix_2nd
if (pixforce ne 0) then pixforce=2; mark as 'needs to be saved'

plots, pix_2nd[0], pix_2nd[1], psym=1, symsize=3, color=color
;
;  Convert from image pixel into helioprojective-cartesian coordinates, in
;  radians.  Store the angles in the common block
;
coord = wcs_get_coord(wcs_2nd, pix_2nd)
conv = !dpi / 180.d0
case wcs_2nd.cunit[0] of
    'arcmin': conv = conv / 60.d0
    'arcsec': conv = conv / 3600.d0
    'mas':    conv = conv / 3600.d3
    'rad':    conv = 1.d0
    else:     conv = conv
endcase
coord = coord * conv
;
;  Calculate the equivalent heliocentric coordinates for distances D within
;  +/- maxz of Dsun.
;
dsun = wcs_2nd.position.dsun_obs
d = dsun + maxz * [-1,1]
cosy = cos(coord[1])
x = d * cosy * sin(coord[0])
y = d * sin(coord[1])
z = dsun - d * cosy * cos(coord[0])
;
;  Determine the spacecraft parameter to pass to convert_sunspice_coord.
;
spacecraft = sc_2nd
test = execute('header = h_'+spacecraft)
obsrvtry = 'Earth'
if datatype(header) eq 'STC' then begin
    if tag_exist(header, 'OBSRVTRY') then begin
        obsrvtry = header.obsrvtry
    end else if tag_exist(header, 'TELESCOP') then begin
        obsrvtry = header.telescop
    endif
end else begin
    temp = fxpar(header, 'OBSRVTRY', count=count)
    if count gt 0 then obsrvtry = temp else begin
        temp = fxpar(header, 'TELESCOP', count=count)
        if count gt 0 then obsrvtry = temp
    endelse
endelse
spacecraft = parse_sunspice_name(obsrvtry, /earth_default)
if wcs_2nd.position.soho and (not wcs_2nd.position.pos_assumed) then $
  spacecraft = 'SOHO'
;
;  Convert to HEEQ coordinates, with rearranging into HGRTN format as an
;  intermediate state.  Store the coordinates in the common block
;
coord = transpose([[z],[x],[y]])
convert_sunspice_coord, wcs_2nd.time.observ_date, coord, 'HGRTN', $
  'HEEQ', spacecraft=spacecraft
heeq_2nd = coord
;
spacecraft = sc_3rd
test = execute('header = h_'+spacecraft)
obsrvtry = 'Earth'
if datatype(header) eq 'STC' then begin
  if tag_exist(header, 'OBSRVTRY') then begin
    obsrvtry = header.obsrvtry
  end else if tag_exist(header, 'TELESCOP') then begin
    obsrvtry = header.telescop
  endif
end else begin
  temp = fxpar(header, 'OBSRVTRY', count=count)
  if count gt 0 then obsrvtry = temp else begin
    temp = fxpar(header, 'TELESCOP', count=count)
    if count gt 0 then obsrvtry = temp
  endelse
endelse
spacecraft = parse_sunspice_name(obsrvtry, /earth_default)
if wcs_3rd.position.soho and (not wcs_3rd.position.pos_assumed) then $
  spacecraft = 'SOHO'
;
;  Switch to the other window, and convert from HEEQ to heliocentric-cartesian
;  x,y,z values for the other perspective.
; Plot the other line
scc_measure3_redisplay, window=win_3rd
coord3=coord
convert_sunspice_coord, wcs_3rd.time.observ_date, coord3, 'HEEQ', 'HGRTN', $
  spacecraft=spacecraft
x = reform(coord3[1,*])
y = reform(coord3[2,*])
z = reform(coord3[0,*])
;
;  Convert from heliocentric-cartesian to helioprojective-cartesian.  Put into
;  the target units.
;
dsun = wcs_3rd.position.dsun_obs
d = sqrt(x^2 + y^2 + (dsun-z)^2)
coord3 = transpose([[atan(x,dsun-z)], [asin(y/d)]])
conv = 180.d0 / !dpi
case wcs_3rd.cunit[0] of
  'arcmin': conv = conv * 60.d0
  'arcsec': conv = conv * 3600.d0
  'mas':    conv = conv * 3600.d3
  'rad':    conv = 1.d0
  else:     conv = conv
endcase
coord3 = coord3 * conv
IF debugon THEN print,'Translating to ',coord3,' ',wcs_3rd.cunit[0]
IF debugon THEN help,dsun,spacecraft
;
;  Convert from helioprojective-cartesian to image pixel coordinates, and
;  overplot the points.
;
plots, pix_param_3rd[0,*], pix_param_3rd[1,*], color=color
pix_param4 = wcs_get_pixel(wcs_3rd, coord3)
plots, pix_param4[0,*], pix_param4[1,*], color=color1
;
;  Extrapolate over the full X-range of the image.
;
param_4 = poly_fit(pix_param4[0,*], pix_param4[1,*], 1)
x3 = [0,wcs_3rd.naxis[1]-1]
plots, x3, poly(x3, param_4), line=1, color=color1
;
;  Based on the HEEQ coordinates from the left and right images, find the
;  intersection on the equatorial (x-y) plane.
;
p1st = poly_fit(heeq_1st[0,*], heeq_1st[1,*], 1)
p2nd = poly_fit(heeq_2nd[0,*], heeq_2nd[1,*], 1)
x = (p1st[0] - p2nd[0]) / (p2nd[1] - p1st[1])
y = (poly(x,p1st) + poly(x,p2nd)) / 2
;
;  Using the same point, find the Z position.
;
p1st = poly_fit(heeq_1st[0,*], heeq_1st[2,*], 1)
p2nd = poly_fit(heeq_2nd[0,*], heeq_2nd[2,*], 1)
z = (poly(x,p1st) + poly(x,p2nd)) / 2
;
;  Store 3rd image coordinates of the point
coordhq = transpose([[x],[y],[z]])
convert_sunspice_coord, wcs_3rd.time.observ_date, coordhq, 'HEEQ', 'HGRTN', $
  spacecraft=spacecraft
x3 = reform(coordhq[1,*])
y3 = reform(coordhq[2,*])
z3 = reform(coordhq[0,*])
;
;  Convert from heliocentric-cartesian to helioprojective-cartesian.  Put into
;  the target units.
;
dsun = wcs_3rd.position.dsun_obs
d = sqrt(x3^2 + y3^2 + (dsun-z3)^2)
coordhq = transpose([[atan(x3,dsun-z3)], [asin(y3/d)]])
conv = 180.d0 / !dpi
case wcs_3rd.cunit[0] of
  'arcmin': conv = conv * 60.d0
  'arcsec': conv = conv * 3600.d0
  'mas':    conv = conv * 3600.d3
  'rad':    conv = 1.d0
  else:     conv = conv
endcase
coordhq= coordhq * conv
;
;  Convert from helioprojective-cartesian to image pixel coordinates, and
;  overplot the points.
;
pix_3rd = wcs_get_pixel(wcs_3rd, coordhq)

; Now check our region of interest again
;
if roi[9] then begin            ;ROI active
  if (pix_3rd[0] lt roi[10] or pix_3rd[0] gt roi[11] or $
      pix_3rd[1] lt roi[12] or pix_3rd[1] gt roi[13]) then begin
    ; out of our region of interest, force an early routine
    roi[8]=0                ; failed
    return
  end else begin
    roi[8]=1                ; success
  end
endif
if (sc_3rd eq 'ahead')  then pixpair[0:1]=pix_3rd
if (sc_3rd eq 'behind') then pixpair[2:3]=pix_3rd
if (sc_3rd eq 'eother') then pixpair[4:5]=pix_3rd
plots, pix_3rd[0], pix_3rd[1], psym=2, symsize=3, color=color2
;
;  Populate the widgets.
;
rad = sqrt(x^2 + y^2 + z^2)
lon = atan(y, x) * 180 / !dpi
widget_control, lon_wid, set_value=ntrim(float(lon))
lat = asin(z / rad) * 180 / !dpi
widget_control, lat_wid, set_value=ntrim(float(lat))
rad = rad / 6.95508e8
widget_control, rsun_wid, set_value=ntrim(float(rad))
;
;  Activate the zoom and store buttons.
;
widget_control, zoom_wid, /sensitive
if outlun gt 0 then widget_control, out_wid, /sensitive
;
end
;
;==============================================================================
;
pro scc_measure3_cleanup, id ;without modifications from 2img version
;
;  Cleanup procedure, to make sure that the output file is properly closed.
;
common scc_measure3
if outlun gt 0 then free_lun, outlun
end
;
;==============================================================================
;
pro store_scc_measure3
  common scc_measure3

  widget_control, lon_wid, get_value=lon
  widget_control, lat_wid, get_value=lat
  widget_control, rsun_wid, get_value=rad
  printf, outlun, lon, lat, rad, pixpair, format='(9(F11.5,1x))'
;        printf, outlun, lon, lat, rad, format='(3F5.5)'
  flush, outlun
  widget_control, out_wid, sensitive=0
  widget_control, clear_wid, /sensitive
  widget_control, oplot_wid, /sensitive
  n_stored = n_stored + 1
  if n_stored eq 1 then begin
      store_left = pix_left
      store_right = pix_right
      store_other = pix_other
  end else begin
      boost_array, store_left,  pix_left
      boost_array, store_right, pix_right
      boost_array, store_other, pix_other
  endelse

end

;==============================================================================
;
pro scc_measure3_event, ev
;
;  Widget event handler.
;
common scc_measure3
;
widget_control, ev.id, get_uvalue=uvalue
case uvalue of
    'EXIT': begin
        if (pixforce eq 2) then store_scc_measure3
        widget_control, /destroy, ev.top
    end
;
    'XLOADCT': xloadct, group=ev.top, updatecallback="scc_measure3_redisplay"
;
;  The plot color was changed.
;
    'COLOR': begin
        widget_control, ev.id, get_value=color
        color = 0 > color < (!d.table_size - 1)
        widget_control, ev.id, set_value=color
        scc_measure3_replot
    end
;
;  The left window was selected.  Start by converting from device pixels to
;  pixels within the image.  Overplot the selected point.
;
    'LEFT': if ev.press gt 0 then begin
        if win_last eq 'LEFT' then in_progress = 0
        if in_progress then begin
          if win_last eq 'RIGHT' then begin
            scc_measure3_select_2nd, ev, win_left, win_other, $
               image_left, image_other, wcs_behind, wcs_eother, $
              'behind','eother', heeq_right, heeq_left, $
               pix_paraml, pix_paramo, pix_left, pix_other
            if roi[8] or (roi[9] eq 0) then in_progress = 0
          end else begin
            scc_measure3_select_2nd, ev, win_left, win_right, $
              image_left, image_right, wcs_behind, wcs_ahead, $
              'behind','ahead', heeq_other, heeq_left, $
              pix_paraml, pix_paramr, pix_left, pix_right
            if roi[8] or (roi[9] eq 0) then in_progress = 0
          endelse
        end else begin
            scc_measure3_select_1st, ev, win_left, win_right, win_other, $
               image_left, image_right, image_other, wcs_behind, $
               wcs_ahead, wcs_eother, 'behind', 'ahead', 'eother', $
               heeq_left, pix_paramr, pix_paramo, pix_left
            if roi[8] or (roi[9] eq 0) then in_progress = 1
        endelse
        if roi[8] or (roi[9] eq 0) then win_last = 'LEFT'
    endif
; 
;  The right window was selected.  Start by converting from device pixels to
;  pixels within the image.
;
    'RIGHT': if ev.press gt 0 then begin
        if win_last eq 'RIGHT' then in_progress = 0
        if in_progress then begin
          if win_last eq 'LEFT' then begin
            scc_measure3_select_2nd, ev, win_right, win_other, $
               image_right, image_other, wcs_ahead, wcs_eother, $
               'ahead','eother', heeq_left, heeq_right, $
               pix_paramr, pix_paramo, pix_right, pix_other
            if roi[8] or (roi[9] eq 0) then in_progress = 0
          end else begin
            scc_measure3_select_2nd, ev, win_right, win_left, $
              image_right, image_left, wcs_ahead, wcs_behind, $
              'ahead','behind', heeq_other, heeq_right, $
              pix_paramr, pix_paraml, pix_right, pix_left
            if roi[8] or (roi[9] eq 0) then in_progress = 0
          endelse
        end else begin
            scc_measure3_select_1st, ev, win_right, win_left, win_other, $
              image_right, image_left, image_other, wcs_ahead,$
              wcs_behind, wcs_eother, 'ahead', 'behind', 'eother', $
              heeq_right, pix_paraml, pix_paramo, pix_right
            if roi[8] or (roi[9] eq 0) then in_progress = 1
        endelse
        if roi[8] or (roi[9] eq 0) then win_last = 'RIGHT'
    endif
;
;  The right window was selected.  Start by converting from device pixels to
;  pixels within the image.
;
'OTHER': if ev.press gt 0 then begin
  if win_last eq 'OTHER' then in_progress = 0
  if in_progress then begin
    if win_last eq 'RIGHT' then begin
      scc_measure3_select_2nd, ev, win_other, win_left, $
        image_other, image_left, wcs_eother, wcs_behind, $
        'eother','behind', heeq_right, heeq_other, $
        pix_paramo, pix_paraml, pix_other, pix_left
        if roi[8] or (roi[9] eq 0) then in_progress = 0
    end else begin
      scc_measure3_select_2nd, ev, win_other, win_right, $
        image_other, image_right, wcs_eother, wcs_ahead, $
        'eother','ahead', heeq_left, heeq_other, $
        pix_paramo, pix_paramr, pix_other, pix_right
      if roi[8] or (roi[9] eq 0) then in_progress = 0
    endelse  
  end else begin
    scc_measure3_select_1st, ev, win_other, win_right, win_left, $
      image_other, image_right, image_left, wcs_eother,$
      wcs_ahead, wcs_behind, 'eother', 'ahead', 'behind', $
      heeq_other, pix_paramr, pix_paraml, pix_other
    if roi[8] or (roi[9] eq 0) then in_progress = 1
  endelse
  if roi[8] or (roi[9] eq 0) then win_last = 'OTHER'
endif
;
;  Modify the image ranges.
;
    'LMIN': scc_measure3_replot, /left
    'LMAX': scc_measure3_replot, /left
    'RMIN': scc_measure3_replot, /right
    'RMAX': scc_measure3_replot, /right
    'OMIN': scc_measure3_replot, /other
    'OMAX': scc_measure3_replot, /other
;
;  Store any selected points.
;
    'STORE': begin
        store_scc_measure3 
        if (pixforce ne 0) then pixforce=1; mark as 'saved'
    end
;
;  Clear the previously stored points from memory.
;
    'CLEAR': begin
        n_stored = 0
        printf, outlun
        flush, outlun
        widget_control, clear_wid, sensitive=0
        scc_measure3_redisplay
    end
;
;  Overplot previously stored points from the output file.
;
    'OPLOT': scc_measure3_oplot_from_file
;
;  Zoom in on the previously selected points.
;
    'ZOOM': begin
        zoomed = 0
        widget_control, lfrz_wid, get_value=lfrz
        if (not lfrz) and (min(wcs_behind.naxis/(2*lzoom)) gt 4) then begin
            zoomed = 1
            lzoom = lzoom * 2
            ;IF h_behind.telescop EQ 'SOLO' THEN lzoom = lzoom * 1.5
            naxis = wcs_behind.naxis
            nn = naxis / lzoom
            origin_left = floor(pix_left - nn/2) > 0
            if origin_left[0]+nn[0] gt naxis[0] then $
              origin_left[0] = naxis[0] - nn[0]
            if origin_left[1]+nn[1] gt naxis[1] then $
              origin_left[1] = naxis[1] - nn[1]
            ix = origin_left[0]
            iy = origin_left[1]
            subimage_left = behind[ix:ix+nn[0]-1, iy:iy+nn[1]-1]
            setwindow, win_left
            exptv, subimage_left, /data, /nobox, /noexact, origin=origin_left, $
                   min=lmin, max=lmax, bscaled=image_left
            plots, pix_left[0], pix_left[1], psym=1, symsize=3, color=color1
        endif
;
        widget_control, rfrz_wid, get_value=rfrz
        if (not rfrz) and (min(wcs_ahead.naxis/(2*rzoom)) gt 4) then begin
            zoomed = 1
            rzoom = rzoom * 2
            ;IF h_ahead.telescop EQ 'SOLO' THEN rzoom = rzoom * 1.5
            naxis = wcs_ahead.naxis
            nn = naxis / rzoom
            origin_right = floor(pix_right - nn/2) > 0
            if origin_right[0]+nn[0] gt naxis[0] then $
              origin_right[0] = naxis[0] - nn[0]
            if origin_right[1]+nn[1] gt naxis[1] then $
              origin_right[1] = naxis[1] - nn[1]
            ix = origin_right[0]
            iy = origin_right[1]
            subimage_right = ahead[ix:ix+nn[0]-1, iy:iy+nn[1]-1]
            setwindow, win_right
            exptv, subimage_right, /data, /nobox, /noexact, origin=origin_right, $
                   min=rmin, max=rmax, bscaled=image_right
            plots, pix_right[0], pix_right[1], psym=1, symsize=3, color=color1
        endif
; 
        widget_control, ofrz_wid, get_value=ofrz
        if (not ofrz) and (min(wcs_eother.naxis/(2*ozoom)) gt 4) then begin
          zoomed = 1
          ozoom = ozoom * 2
          IF h_eother.telescop EQ 'SOLO' THEN ozoom = ozoom * 2.
          naxis = wcs_eother.naxis
          nn = naxis / ozoom
          origin_other = floor(pix_other - nn/2) > 0
          if origin_other[0]+nn[0] gt naxis[0] then $
            origin_other[0] = naxis[0] - nn[0]
          if origin_other[1]+nn[1] gt naxis[1] then $
            origin_other[1] = naxis[1] - nn[1]
          ix = origin_other[0]
          iy = origin_other[1]
          subimage_other = eother[ix:ix+nn[0]-1, iy:iy+nn[1]-1]
          setwindow, win_other
          exptv, subimage_other, /data, /nobox, /noexact, origin=origin_other, $
            min=omin, max=omax, bscaled=image_other
          plots, pix_other[0], pix_other[1], psym=1, symsize=3, color=color1
        endif
;
;
;  Redisplay any overplots.
;
        if zoomed then begin
            widget_control, unzoom_wid, /sensitive
            scc_measure3_redisplay
        endif
    end
;
;  Zoom out from the previously zoomed image.
;
    'UNZOOM': begin
        unzoomed = 0
        widget_control, lfrz_wid, get_value=lfrz
        if (not lfrz) and (lzoom gt 1) then begin
            unzoomed = 1
            lzoom = lzoom / 2
            if lzoom eq 1 then begin
                setwindow, win_left
                subimage_left = behind
                origin_left = [0,0]
                exptv, behind, /data, /nobox, /noexact, min=lmin, max=lmax, $
                       bscaled=image_left
                if not in_progress then plots, pix_left[0], pix_left[1], psym=1, $
                                               symsize=3, color=color1
            end else begin
                naxis = wcs_behind.naxis
                nn = naxis / lzoom
                origin_left = floor(pix_left - nn/2) > 0
                if origin_left[0]+nn[0] gt naxis[0] then $
                  origin_left[0] = naxis[0] - nn[0]
                if origin_left[1]+nn[1] gt naxis[1] then $
                  origin_left[1] = naxis[1] - nn[1]
                ix = origin_left[0]
                iy = origin_left[1]
                subimage_left = behind[ix:ix+nn[0]-1, iy:iy+nn[1]-1]
                setwindow, win_left
                exptv, subimage_left, /data, /nobox, /noexact, min=lmin, $
                       max=lmax, origin=origin_left, bscaled=image_left
                if not in_progress then plots, pix_left[0], pix_left[1], psym=1, $
                                               symsize=3, color=color1
            endelse
        endif
;
        widget_control, rfrz_wid, get_value=rfrz
        if (not rfrz) and (rzoom gt 1) then begin
            unzoomed = 1
            rzoom = rzoom / 2
            if rzoom eq 1 then begin
                setwindow, win_right
                subimage_right = ahead
                origin_right = [0,0]
                exptv, ahead, /data, /nobox, /noexact, min=rmin, max=rmax, $
                       bscaled=image_right
                if not in_progress then plots, pix_right[0], pix_right[1], $
                                               psym=1, symsize=3, color=color1
            end else begin
                naxis = wcs_ahead.naxis
                nn = naxis / rzoom
                origin_right = floor(pix_right - nn/2) > 0
                if origin_right[0]+nn[0] gt naxis[0] then $
                  origin_right[0] = naxis[0] - nn[0]
                if origin_right[1]+nn[1] gt naxis[1] then $
                  origin_right[1] = naxis[1] - nn[1]
                ix = origin_right[0]
                iy = origin_right[1]
                subimage_right = ahead[ix:ix+nn[0]-1, iy:iy+nn[1]-1]
                setwindow, win_right
                exptv, subimage_right, /data, /nobox, /noexact, min=rmin, $
                       max=rmax, origin=origin_right, bscaled=image_right
                if not in_progress then plots, pix_right[0], pix_right[1], $
                                               psym=1, symsize=3, color=color1
            endelse
        endif
 ;
        widget_control, ofrz_wid, get_value=ofrz
        if (not ofrz) and (ozoom gt 1) then begin
          unzoomed = 1
          ozoom = ozoom / 2
          if ozoom eq 1 then begin
            setwindow, win_other
            subimage_other = eother
            origin_other = [0,0]
            exptv, eother, /data, /nobox, /noexact, min=omin, max=omax, $
              bscaled=image_other
            if not in_progress then plots, pix_other[0], pix_other[1], $
              psym=1, symsize=3, color=color
          end else begin
            naxis = wcs_ahead.naxis
            nn = naxis / ozoom
            origin_other = floor(pix_other - nn/2) > 0
            if origin_other[0]+nn[0] gt naxis[0] then $
              origin_other[0] = naxis[0] - nn[0]
            if origin_other[1]+nn[1] gt naxis[1] then $
              origin_other[1] = naxis[1] - nn[1]
            ix = origin_other[0]
            iy = origin_other[1]
            subimage_other = eother[ix:ix+nn[0]-1, iy:iy+nn[1]-1]
            setwindow, win_other
            exptv, subimage_other, /data, /nobox, /noexact, min=omin, $
              max=omax, origin=origin_other, bscaled=image_other
            if not in_progress then plots, pix_other[0], pix_other[1], $
              psym=1, symsize=3, color=color1
          endelse
        endif
;
        

        if (lzoom eq 1) and (rzoom eq 1) then widget_control, unzoom_wid, sensitive=0
;
;  Redisplay any overplots.
;
        if unzoomed then scc_measure3_redisplay
    end
;
;  Dummy events for the LFRZ and RFRZ checkboxes.
;
    'LFRZ':
    'RFRZ':
    'OFRZ':
;
    else: message, /continue, 'Unrecognized event ' + uvalue
endcase
;
end
;
;==============================================================================
;
pro scc_measure3, file_behind, file_ahead, file_eother, index_behind, index_ahead, index_eother, $
                 wsize=k_wsize, color=c0, color1=c1, color2=c2, $
                 outfile=outfile, secchiprep=secchiprep, $
                 append=append, forcesave=forcesave, crop=crop, $
                 no_block=no_block, debug=debug, _extra=_extra
; ahead, behind, eother \\ // left, right, other
;  Main procedure to set up the widget.
;
common scc_measure3
;
IF keyword_set(DEBUG) THEN debugon=1 ELSE debugon=0

if (n_elements(forcesave) ne 0) then pixforce=1 else pixforce=0
pixpair=fltarr(6)
if (xregistered("scc_measure3") NE 0) then return
;
if datatype(file_ahead,1) eq 'String' then begin
    if keyword_set(secchiprep) then begin
        secchi_prep, file_ahead,  h_ahead,  ahead,  _extra=_extra
        secchi_prep, file_behind, h_behind, behind, _extra=_extra
        secchi_prep, file_eother, h_eother, eother, _extra=_extra
    end else begin
        ahead  = sccreadfits(file_ahead,  h_ahead)
        behind = sccreadfits(file_behind, h_behind)
        eother = sccreadfits(file_eother, h_eother)
    endelse
    wcs_ahead = fitshead2wcs(h_ahead,KEYWORD_NULL=0)
    wcs_behind = fitshead2wcs(h_behind,KEYWORD_NULL=0)
    wcs_eother = fitshead2wcs(h_eother,KEYWORD_NULL=0)
end else begin
    ahead = file_ahead
    h_ahead = index_ahead
    wcs_ahead = fitshead2wcs(index_ahead)
    behind = file_behind
    h_behind = index_behind
    wcs_behind = fitshead2wcs(index_behind)
    eother = file_eother
    h_eother = index_eother
    wcs_eother = fitshead2wcs(index_eother)
endelse
;
; Get correct label for windows
;
IF datatype(h_ahead) NE 'STC' THEN h_ahead=fitshead2struct(h_ahead)
IF datatype(h_behind) NE 'STC' THEN h_behind=fitshead2struct(h_behind)
IF datatype(h_eother) NE 'STC' THEN h_eother=fitshead2struct(h_eother)
IF tag_exist(h_ahead, 'telescop') THEN BEGIN
    IF h_ahead.telescop EQ 'STEREO' THEN BEGIN
        IF h_ahead.obsrvtry EQ 'STEREO_A' THEN scright='Ahead' ELSE $
        IF h_ahead.obsrvtry EQ 'STEREO_B' THEN scright='Behind'
    ENDIF ELSE scright=h_ahead.telescop
ENDIF ELSE scright = 'First'
IF tag_exist(h_behind, 'telescop') THEN BEGIN
    IF h_behind.telescop EQ 'STEREO' THEN BEGIN
        IF h_behind.obsrvtry EQ 'STEREO_A' THEN scleft='Ahead' ELSE $
        IF h_behind.obsrvtry EQ 'STEREO_B' THEN scleft='Behind'
    ENDIF ELSE scleft=h_behind.telescop
ENDIF ELSE scleft = 'Second'
IF tag_exist(h_eother, 'telescop') THEN BEGIN
  IF h_eother.telescop EQ 'STEREO' THEN BEGIN
    IF h_eother.obsrvtry EQ 'STEREO_A' THEN scother='Ahead' ELSE $
      IF h_eother.obsrvtry EQ 'STEREO_B' THEN scother='Behind'
  ENDIF ELSE scother=h_eother.telescop
ENDIF ELSE scother = 'Third'
;
; LASCO level-0.5 has CUNIT undefined
;
IF scright EQ 'SOHO' THEN wcs_ahead.cunit =['arcsec','arcsec']
IF scleft  EQ 'SOHO' THEN wcs_behind.cunit=['arcsec','arcsec']
IF scother EQ 'SOHO' THEN wcs_eother.cunit =['arcsec','arcsec']
;
; If setting up for cropping, store the allowable region of interest,
; otherwise defaultq region is 'the entire image'
;
roi =fltarr(14)
if (n_elements(crop) ne 0) then begin
    roi[0]=wcs_ahead.naxis[0]/4
    roi[1]=3*wcs_ahead.naxis[0]/4
    roi[2]=wcs_ahead.naxis[1]/4
    roi[3]=3*wcs_ahead.naxis[1]/4
    roi[4]=wcs_behind.naxis[0]/4
    roi[5]=3*wcs_behind.naxis[0]/4
    roi[6]=wcs_behind.naxis[1]/4
    roi[7]=3*wcs_behind.naxis[1]/4
    roi[8]=1; status flag
    roi[9]=1; ROI active
    roi[10]=wcs_eother.naxis[0]/4
    roi[11]=3*wcs_eother.naxis[0]/4
    roi[12]=wcs_eother.naxis[1]/4
    roi[13]=3*wcs_eother.naxis[1]/4

    ; also mark our ROI within the data
    ahead[roi[0]:roi[1],roi[2]]=max(ahead)
    ahead[roi[0]:roi[1],roi[3]]=max(ahead)
    ahead[roi[0],roi[2]:roi[3]]=max(ahead)
    ahead[roi[1],roi[2]:roi[3]]=max(ahead)

    behind[roi[4]:roi[5],roi[6]]=max(behind)
    behind[roi[4]:roi[5],roi[7]]=max(behind)
    behind[roi[4],roi[6]:roi[7]]=max(behind)
    behind[roi[5],roi[6]:roi[7]]=max(behind)
    
    other[roi[10]:roi[11],roi[12]]=max(eother)
    other[roi[10]:roi[11],roi[13]]=max(eother)
    other[roi[10]:roi[12],roi[13]]=max(eother)
    other[roi[11],roi[12]:roi[13]]=max(eother)

endif



;
;  Set up the main widget base.
;
main = widget_base(title="3 viewpoint 3D coordinate measuring tool", /column)
;
;  Set up the two graphics windows.
;
if n_elements(k_wsize) eq 1 then wsize=k_wsize else wsize=512
if n_elements(c0) eq 1 then color= c0 else color =254b ;white
if n_elements(c1) eq 1 then color1=c1 else color1=250b ;blue
if n_elements(c2) eq 1 then color2=c2 else color2=252b ;magenta
show = widget_base(main, /row)
lcolumn = widget_base(show, /column, /frame)
dummy = widget_label(lcolumn, value=scleft)
left  = widget_draw(lcolumn, xsize=wsize, ysize=wsize, retain=2, $
                    uvalue="LEFT", /button_events)
dummy = widget_base(lcolumn, /row)
lmin_wid = cw_field(dummy, title='Image minimum:', xsize=15, /floating, $
                uvalue='LMIN', /return_events)
lmax_wid = cw_field(dummy, title='maximum:', xsize=15, /floating, $
                uvalue='LMAX', /return_events)
lfrz_wid = cw_bgroup(dummy, 'Freeze size', uvalue='LFRZ', /nonexclusive)
;
rcolumn = widget_base(show, /column, /frame)
dummy = widget_label(rcolumn, value=scright)
right = widget_draw(rcolumn, xsize=wsize, ysize=wsize, retain=2, $
                    uvalue="RIGHT", /button_events)
dummy = widget_base(rcolumn, /row)
rmin_wid = cw_field(dummy, title='Image minimum:', xsize=15, /floating, $
                uvalue='RMIN', /return_events)
rmax_wid = cw_field(dummy, title='maximum:', xsize=15, /floating, $
                uvalue='RMAX', /return_events)
rfrz_wid = cw_bgroup(dummy, 'Freeze size', uvalue='RFRZ', /nonexclusive)
;
ocolumn = widget_base(show, /column, /frame)
dummy = widget_label(ocolumn, value=scother)
other = widget_draw(ocolumn, xsize=wsize, ysize=wsize, retain=2, $
  uvalue="OTHER", /button_events)
dummy = widget_base(ocolumn, /row)
omin_wid = cw_field(dummy, title='Image minimum:', xsize=15, /floating, $
  uvalue='OMIN', /return_events)
omax_wid = cw_field(dummy, title='maximum:', xsize=15, /floating, $
  uvalue='OMAX', /return_events)
ofrz_wid = cw_bgroup(dummy, 'Freeze size', uvalue='OFRZ', /nonexclusive)
;
;  Define a row for the output, and define fields for the heliographic
;  longitude, latitude, and radial distance.
;
output = widget_base(main, /row)
dummy = widget_label(output, value='Heliographic')
lon_wid  = cw_field(output, title='Longitude:',  xsize=10, /noedit)
lat_wid  = cw_field(output, title='Latitude:',  xsize=10, /noedit)
rsun_wid = cw_field(output, title='Solar radii:', xsize=10, /noedit)
out_wid = widget_button(output, value='Store', uvalue='STORE')
clear_wid = widget_button(output, value='Clear stored', uvalue='CLEAR')
oplot_wid = widget_button(output, value='Overplot from file', uvalue='OPLOT')
;
;  Define a row for buttons.  Place the exit button a little off from the other
;  buttons.
;
control = widget_base(main, /row)
dummy = widget_button(control, uvalue="XLOADCT", value='Adjust color table')
zoom_wid = widget_button(control, uvalue="ZOOM", value='Zoom in')
unzoom_wid = widget_button(control, uvalue="UNZOOM", value='Zoom out')
color = !d.table_size - 1
;  set colors, corresponds to ANA color table 47
color_wid = widget_slider(control, uvalue="COLOR", value=color, minimum=0, $
                          maximum=!d.table_size-1, title='Plot color')
dummy = widget_label(control, value="          ")
dummy = widget_button(control, uvalue="EXIT", Value="Exit", /frame)
;
;  Realize the widget, and draw the images.
;
widget_control, main, /realize
widget_control, left, get_value=win_left
setwindow, win_left
test = sigrange(behind)
IF scleft EQ 'SOLO' THEN test = hist_equal(alog(behind>1),per=0.1)
lmin = min(test, max=lmax)
exptv, test, /data, /nobox, /noexact, bscaled=image_left
subimage_left = behind
origin_left = [0,0]
widget_control, lmin_wid, set_value=lmin
widget_control, lmax_wid, set_value=lmax
;
widget_control, right, get_value=win_right
setwindow, win_right
test = sigrange(ahead)
IF scright EQ 'SOLO' THEN test = hist_equal(alog(ahead>1),per=0.1)
rmin = min(test, max=rmax)
exptv, test, /data, /nobox, /noexact, bscaled=image_right
subimage_right = ahead
origin_right = [0,0]
widget_control, rmin_wid, set_value=rmin
widget_control, rmax_wid, set_value=rmax
;
widget_control, other, get_value=win_other
setwindow, win_other
test = sigrange(eother)
IF scother EQ 'SOLO' THEN test = hist_equal(alog(eother>1),per=0.1)
omin = min(test, max=omax)
exptv, test, /data, /nobox, /noexact, bscaled=image_other
subimage_other = eother
origin_other = [0,0]
widget_control, omin_wid, set_value=omin
widget_control, omax_wid, set_value=omax
;
setwindow, win_left             ;Makes sure that plot parameters are stored.
;
;  Deactivate the zoom, unzoom, and store widgets.
;
widget_control, zoom_wid, sensitive=0
widget_control, unzoom_wid, sensitive=0
widget_control, out_wid, sensitive=0
widget_control, clear_wid, sensitive=0
if not keyword_set(append) then widget_control, oplot_wid, sensitive=0
;
;  Determine the maximum scale of the two images, in meters.  Make sure that
;  it's at least 3 solar radii
;
scale0 = max(wcs_ahead.naxis*wcs_ahead.cdelt)
conv = !dpi / 180.d0
case wcs_ahead.cunit[0] of
    'arcmin': conv = conv / 60.d0
    'arcsec': conv = conv / 3600.d0
    'mas':    conv = conv / 3600.d3
    'rad':    conv = 1.d0
    else:     conv = conv
endcase
scale0 = scale0 * conv * wcs_ahead.position.dsun_obs
;
scale1 = max(wcs_behind.naxis*wcs_behind.cdelt)
conv = !dpi / 180.d0
case wcs_behind.cunit[0] of
    'arcmin': conv = conv / 60.d0
    'arcsec': conv = conv / 3600.d0
    'mas':    conv = conv / 3600.d3
    'rad':    conv = 1.d0
    else:     conv = conv
endcase
scale1 = scale1 * conv * wcs_behind.position.dsun_obs
;
scale2 = max(wcs_eother.naxis*wcs_eother.cdelt)
conv = !dpi / 180.d0
case wcs_eother.cunit[0] of
  'arcmin': conv = conv / 60.d0
  'arcsec': conv = conv / 3600.d0
  'mas':    conv = conv / 3600.d3
  'rad':    conv = 1.d0
  else:     conv = conv
endcase
scale2 = scale2 * conv * wcs_eother.position.dsun_obs
maxz = scale0 > scale1 > scale2 > 2.1e9
;
;  If the OUTFILE parameter was passed, then open a file for writing.
;
if n_elements(outfile) eq 1 then $
  openw, outlun, outfile, /get_lun, append=append else $
  outlun = -1
;
;  Initialize some control parameters, and set the widget going.
;
in_progress = 0
win_last = ''
lzoom = 1
rzoom = 1
ozoom = 1
n_stored = 0
;
xmanager, 'scc_measure3', main, cleanup='scc_measure3_cleanup', no_block=no_block
end