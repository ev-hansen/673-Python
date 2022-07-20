;+
;PROCEDURE:	calc_hope_pressure
;PURPOSE:
;  load a .cdf ECT-HOPE level 3 file and calculate particle pressure 
;
;INPUTS:
;  probe:	'a' or 'b'
;  factor:      Factor to be multiplied to the particle fluxes
;
;KEYWORDS:
;  time_range:	2 element vector specifying the time range.
;
;CREATED BY:	Cristian Ferradas, 2021-06-01 
;UPDATED BY:	Cristian Ferradas, 2021-10-20, added option to use relativistic energies
;-
;----------------------------------------------------------------------

PRO calc_hope_pressure, probe=probe, $
                        factor=factor, $
                        low_energy=low_energy, $
                        up_energy=up_energy, $
                        pot_corr=pot_corr, $
                        relat=relat, $
                        swindow=swindow, $
                        time_range=trange ; optional

  COMMON get_error, get_err_no, get_err_msg, default_verbose

;----------------------------------------------------------------------
; Manually enter timespan
;----------------------------------------------------------------------
 IF 0 EQ 1 THEN BEGIN

   time_s_str = '2013-03-17/00:00:00'  ;'2013-02-08/10:00:00'
   time_e_str = '2013-03-18/00:00:00'  ;'2013-02-08/18:00:00'

   time_s = time_double(time_s_str)  ;-4.5*3600
   time_e = time_double(time_e_str)  ;time_s+9.0*3600

   time = time_s

   timespan, time, time_e-time_s, /SECONDS
   probe = 'a'           ;'a','b'
   factor = 1.0
   low_energy = [1e2,1e2]   ; [i+, e-] in eV
   up_energy = [5.5e4,5.5e4]    ; [i+, e-] in eV
   pot_corr = 1        ; Use sc potential corrected fluxes? 0:No, 1:Yes
   relat = 1           ; Use relativistic energy equation? 0:No, 1:Yes
   swindow = 60        ; Smoothing window in sec; swindow=0: No smoothing

 ENDIF

;----------------------------------------------------------------------

 ;tplot_options, 'title', 'Van Allen Probes-'+strupcase(probe)

 IF NOT keyword_set(trange) THEN get_timespan, tr ELSE tr=time_double(trange)
 IF n_elements(tr) EQ 1 THEN tr = [tr, tr+86399d0]

 ts = time_struct(tr(0))
 te = time_struct(tr(1))

;----------------------------------------------------------------------
; Calculate number of days
;----------------------------------------------------------------------
 nds_flt = (tr[1]-tr[0])/60.0/60.0/24.0
     nds = ceil(nds_flt)

;**---------
  dify = te.year-ts.year
  difm = te.month-ts.month
  difd = te.doy-ts.doy

  IF (difd LT 0) AND (dify LT 0) THEN stop, 'Wrong timespan!'
  IF dify EQ 0 THEN BEGIN
    ndys = difd+1
    IF (te.hour EQ 0) AND (te.min EQ 0) AND (te.sec EQ 0) THEN ndys = difd
  ENDIF ELSE BEGIN
    ndys = nds
  ENDELSE
;**---------

 print, '# of days = ', ndys

;----------------------------------------------------------------------
; Get data for specified days
;----------------------------------------------------------------------
; Initialize system variables for HOPE
  rbsp_hope_init, probe=probe

FOR idy = 0, ndys-1 DO BEGIN                       ; each day

  it = tr[0] + idy*24.0*60.0*60.0                  ; read in the daily files one by one
  t  = time_struct(it)

  yr = t.year
  mn = t.month
  dt = t.date

  yr_str = string(yr,'(I4.4)')
  mn_str = string(mn,'(I2.2)')
  dt_str = string(dt,'(I2.2)')

  ymd_str = yr_str + mn_str + dt_str

  dir = '~/data/rbsp/rbsp'+probe+'/ect/hope/level3/PA/rel04/'+yr_str+'/'    ; public files
  nm_hd   = 'rbsp'+probe+'_rel04_ect-hope-*-*3_'           ; public files
  fln_tmp = nm_hd+ymd_str+'_v*.cdf'

; Check if cdf file exists in local directory, if not then download it
  rp = !rbsp_hope.remote_data_dir+'rbsp'+probe+'/hope/level3/pitchangle/'+yr_str+'/'
  rf = fln_tmp
  lp = !rbsp_hope.local_data_dir+'rbsp'+probe+'/ect/hope/level3/PA/rel04/'+yr_str+'/'

  dir_fln = file_search(dir+fln_tmp, fold_case=1, count=nfl)
  IF nfl EQ 0 THEN BEGIN
    files = spd_download(remote_path=rp, remote_file=rf, local_path=lp, /LAST_VERSION, /NO_UPDATE) ;local_file=rf
    dir_fln = file_search(dir+fln_tmp, fold_case=1, count=nfl)
  ENDIF

  IF nfl EQ 0 THEN goto, nofile

  fln_r_wdir = dir_fln[n_elements(dir_fln)-1]

  ; print the info of the .cdf file
  ; print_cdf_info, fln_r_wdir
  print, fln_r_wdir

;----------------------------------------------------------------------
; Read data in the .cdf file
;----------------------------------------------------------------------
  dd = get_cdf_data(file=fln_r_wdir)

; Time
  t_ion     = reform(dd.epoch_ion.data)
  t_ele     = reform(dd.epoch_ele.data)
  t_ion_str = time_double(t_ion, /EPOCH)
  t_ele_str = time_double(t_ele, /EPOCH)

; Energy
; FLOAT     Array[72, 3468]
  e_data_ion = transpose(dd.hope_energy_ion.data)   ; UNITS           STRING    'eV'
  e_data_ele = transpose(dd.hope_energy_ele.data)   ; UNITS           STRING    'eV'

; Ion mode of operation
; Indicates instrument mode for ion data.
; Each mode has a different set of energy channels.
; See corresponding HOPE_ENERGY_ION record.
; (0) Apogee mode: normal operation.
; (1) Perigee mode: minimum energy channel raised during pass through perigee.
  mode_ion = reform(dd.mode_ion.data)

; Electron mode of operation
; Indicates instrument mode for electron data.
; Each mode has a different set of energy channels.
; See corresponding HOPE_ENERGY_Ion record.
; (0) Apogee mode: normal operation.
; (1) Perigee mode: minimum energy channel raised during pass through perigee.
; (2) burst mode: subset of energies sampled at rapid cadence.
  mode_ele = reform(dd.mode_ele.data)

; Variables to save
  var_arr = ['Pressure_H','Pressure_HE','Pressure_O','Pressure_E']   ;,'COUNTS_P','COUNTS_HE','COUNTS_O','COUNTS_E']
  nm_arr  = ['H', 'He', 'O', 'e']   ;, 'H', 'He', 'O', 'e']
  nvar    = n_elements(var_arr)

;;----------------------------------------------------------------------
FOR ivar = 0, nvar-1 DO BEGIN                                         ; for each species

; FPDU:
; Proton flux as a function of pitch angle (11) as computed from magnetic field direction crossed with spacecraft spin axis and energy (72).  
; Note that data have been expanded to the highest possible resolution, even though before 9/2013 the transmitted data are collapsed 
; (36 energies, 4-8-16-8-4 spin angles) to fit the telemetry allocation.
; UNITS           STRING    's!E-1!Ncm!E-2!Nster!E-1!NkeV!E-1!N'
; FILLVAL         FLOAT      -1.00000e+31
  IF ivar EQ 0 THEN ddata = dd.fpdu.data          ; H+ unidirectional flux  ; FLOAT     = Array[72, 11, 3468]
  IF ivar EQ 1 THEN ddata = dd.fhedu.data         ; He+ unidirectional flux ; FLOAT     = Array[72, 11, 3468]
  IF ivar EQ 2 THEN ddata = dd.fodu.data          ; O+ unidirectional flux  ; FLOAT     = Array[72, 11, 3468]
  IF ivar EQ 3 THEN ddata = dd.fedu.data          ; e- unidirectional flux  ; FLOAT     = Array[72, 11, 3468]
  ;IF ivar EQ 4 THEN ddata = dd.COUNTS_P.data          ; H+ unidirectional counts
  ;IF ivar EQ 5 THEN ddata = dd.COUNTS_HE.data         ; He+ unidirectional counts
  ;IF ivar EQ 6 THEN ddata = dd.COUNTS_O.data          ; O+ unidirectional counts
  ;IF ivar EQ 7 THEN ddata = dd.COUNTS_E.data          ; e- unidirectional counts

; Mark fill values (including the filling value:  -1.00000e+31 [?]) as NaN
  nan_sub = where(ddata EQ -1.00000e+31, nn)
  IF nn GT 0 THEN ddata[nan_sub] = !Values.F_NaN

  daty_avg_int = transpose(ddata)

; Multiply by a factor
  IF ivar EQ 0 THEN daty_avg_int = daty_avg_int*factor
  IF ivar EQ 1 THEN daty_avg_int = daty_avg_int*factor
  IF ivar EQ 2 THEN daty_avg_int = daty_avg_int*factor
  ;IF ivar EQ 3 THEN daty_avg_int = daty_avg_int*factor

;----------------------------------------------------------------------
; Create variables with the ion/electron data. If timespan is longer than
; one day then combine the data for several days
;----------------------------------------------------------------------
   IF idy EQ 0 THEN BEGIN

    if ivar eq 0 then dat_t_ion = t_ion_str
    if ivar eq 0 then daty_e_ion = e_data_ion
    if ivar eq 0 then dat_mode_ion = mode_ion
    if ivar eq 0 then daty_avg_int_H1  = daty_avg_int
    if ivar eq 1 then daty_avg_int_He1 = daty_avg_int
    if ivar eq 2 then daty_avg_int_O1  = daty_avg_int
    if ivar eq 3 then dat_t_ele = t_ele_str
    if ivar eq 3 then daty_e_ele = e_data_ele
    if ivar eq 3 then dat_mode_ele = mode_ele
    if ivar eq 3 then daty_avg_int_e1  = daty_avg_int
    ;if ivar eq 4 then daty_counts_H1 = daty_avg_int
    ;if ivar eq 5 then daty_counts_He1 = daty_avg_int
    ;if ivar eq 6 then daty_counts_O1 = daty_avg_int
    ;if ivar eq 7 then daty_counts_e1 = daty_avg_int

   ENDIF ELSE BEGIN

    if ivar eq 0 then dat_t_ion = [dat_t_ion, t_ion_str]
    if ivar eq 0 then daty_e_ion = [daty_e_ion, e_data_ion]
    if ivar eq 0 then dat_mode_ion = [dat_mode_ion, mode_ion]
    if ivar eq 0 then daty_avg_int_H1  = [daty_avg_int_H1,  daty_avg_int]
    if ivar eq 1 then daty_avg_int_He1 = [daty_avg_int_He1, daty_avg_int]
    if ivar eq 2 then daty_avg_int_O1  = [daty_avg_int_O1,  daty_avg_int]
    if ivar eq 3 then dat_t_ele = [dat_t_ele, t_ele_str]
    if ivar eq 3 then daty_e_ele = [daty_e_ele, e_data_ele]
    if ivar eq 3 then dat_mode_ele = [dat_mode_ele, mode_ele]
    if ivar eq 3 then daty_avg_int_e1 = [daty_avg_int_e1, daty_avg_int]
    ;if ivar eq 4 then daty_counts_H1 = [daty_counts_H1, daty_avg_int]
    ;if ivar eq 5 then daty_counts_He1 = [daty_counts_He1, daty_avg_int]
    ;if ivar eq 6 then daty_counts_O1 = [daty_counts_O1, daty_avg_int]
    ;if ivar eq 7 then daty_counts_e1 = [daty_counts_e1, daty_avg_int]

   ENDELSE

ENDFOR                                             ; each species
;----------------------------------------------------------------------

  nofile:  IF nfl EQ 0 THEN print, 'File missing: ', dir+fln_tmp

ENDFOR                                             ; each day
;----------------------------------------------------------------------

; Limit data arrays to time interval requested
  get_timespan, tt
  itime_ion = where(dat_t_ion GE tt[0] AND dat_t_ion LE tt[1], c_itime_ion)
  itime_ele = where(dat_t_ele GE tt[0] AND dat_t_ele LE tt[1], c_itime_ele)
  IF c_itime_ion LE 1 THEN BEGIN
    get_err_no = 1
    get_err_msg = 'Less than 2 ion data points found for time interval'
    MESSAGE, get_err_msg, /CONTINUE
    RETURN
  ENDIF
  IF c_itime_ele LE 1 THEN BEGIN
    get_err_no = 1
    get_err_msg = 'Less than 2 ele data points found for time interval'
    MESSAGE, get_err_msg, /CONTINUE
    RETURN
  ENDIF

  dat_t_ion        = dat_t_ion[itime_ion]
  daty_e_ion       = daty_e_ion[itime_ion,*]
  dat_mode_ion     = dat_mode_ion[itime_ion]
  daty_avg_int_H1  = daty_avg_int_H1[itime_ion,*,*]
  daty_avg_int_He1 = daty_avg_int_He1[itime_ion,*,*]
  daty_avg_int_O1  = daty_avg_int_O1[itime_ion,*,*]
  ;daty_counts_H1   = daty_counts_H1[itime_ion,*,*]
  ;daty_counts_He1  = daty_counts_He1[itime_ion,*,*]
  ;daty_counts_O1   = daty_counts_O1[itime_ion,*,*]

  dat_t_ele       = dat_t_ele[itime_ele]
  daty_e_ele      = daty_e_ele[itime_ele,*]
  dat_mode_ele    = dat_mode_ele[itime_ele]
  daty_avg_int_e1 = daty_avg_int_e1[itime_ele,*,*]
  ;daty_counts_e1  = daty_counts_e1[itime_ele,*,*]

; Limit data arrays to energy range
  id = 0
  m0 = dat_mode_ion[id]
  WHILE m0 NE 0 DO BEGIN
    id = id+1
    m0 = dat_mode_ion[id]
  ENDWHILE
  en_ind_ion = where(daty_e_ion[id,*] GE low_energy[0] AND daty_e_ion[id,*] LE up_energy[0], nechi)
  en_ind_ele = where(daty_e_ele[id,*] GE low_energy[1] AND daty_e_ele[id,*] LE up_energy[1], neche)

; First and last energy bins considered
  eion0 = daty_e_ion[id,en_ind_ion[0]]
  eion1 = daty_e_ion[id,en_ind_ion[nechi-1]]
  eele0 = daty_e_ele[id,en_ind_ele[0]]
  eele1 = daty_e_ele[id,en_ind_ele[neche-1]]

;----------------------------------------------------------------------
; Delete data on perigee mode and electron burst mode
;----------------------------------------------------------------------
; Delete ion data on perigee mode
  indi = where(dat_mode_ion EQ 1)
  daty_e_ion[indi,*] = !Values.F_NAN
  daty_avg_int_H1[indi,*,*] = !Values.F_NAN
  daty_avg_int_He1[indi,*,*] = !Values.F_NAN
  daty_avg_int_O1[indi,*,*] = !Values.F_NAN

; Delete electron data on burst mode
  mode_e = float(dat_mode_ele)
  inde = where(dat_mode_ele EQ 2)
  daty_e_ele[inde,*] = !Values.F_NAN
  daty_avg_int_e1[inde,*,*] = !Values.F_NAN
  mode_e[inde] = !Values.F_NAN

  FOR ie=0, 71 DO BEGIN
    daty_e_ele[*,ie] = interpol(daty_e_ele[*,ie], dat_t_ele, dat_t_ele, /NAN)
    FOR ip=0, 10 DO BEGIN
      ifin = where(finite(daty_avg_int_e1[*,ip,ie]) EQ 1, cfin)
      daty_avg_int_e1[*,ip,ie] = interpol(daty_avg_int_e1[*,ip,ie], dat_t_ele, dat_t_ele, /NAN)
      daty_avg_int_e1[0:ifin[0],ip,ie] = !Values.F_NAN
    ENDFOR
  ENDFOR
  ;dat_mode_ele_int = interpol(mode_e, dat_t_ele, dat_t_ele, /NAN)

; Delete electron data on perigee mode
  inde = where(dat_mode_ele EQ 1)
  daty_e_ele[inde,*] = !Values.F_NAN
  daty_avg_int_e1[inde,*,*] = !Values.F_NAN

  dat_mode_ion_sm     = dat_mode_ion
  dat_mode_ele_sm     = dat_mode_ele
  ;dat_mode_ele_int_sm = dat_mode_ele_int

; Assign fluxes to new variables
  daty_h  = daty_avg_int_H1
  daty_he = daty_avg_int_He1
  daty_o  = daty_avg_int_O1
  daty_e  = daty_avg_int_e1

;----------------------------------------------------------------------
; Correct fluxes in each energy channel with the measured spacecraft potential
;----------------------------------------------------------------------
  IF pot_corr EQ 1 THEN BEGIN

; Get meaured spacecraft potential from EFW
  ;read_cdf_efw_l3, probe=probe, smooth=0
  ;read_txt_efw_scpot, probe=probe
  tplot_restore, f='~/data/rbsp/rbsp'+probe+'/efw/rbsp'+probe+'_scpot.tplot'
  get_data, 'rbsp'+probe+'_scpot', data=scpot        ; negative potential
  pot = scpot.y*(-1)                                 ; potential

  scpot_ion = interpol(pot, scpot.x, dat_t_ion)
  scpot_ele = interpol(pot, scpot.x, dat_t_ele)

; Add the spacecraft potential to each HOPE energy channel
  emod_ion = fltarr(c_itime_ion, 72)
  emod_ele = fltarr(c_itime_ele, 72)

  FOR it=0, c_itime_ion-1 DO emod_ion[it,*] = daty_e_ion[it,*]+scpot_ion[it]
  FOR it=0, c_itime_ele-1 DO emod_ele[it,*] = daty_e_ele[it,*]-scpot_ele[it]

; Interpolate fluxes back to the original HOPE energy channels
  daty_h_corr  = fltarr(c_itime_ion, 11, 72)
  daty_he_corr = fltarr(c_itime_ion, 11, 72)
  daty_o_corr  = fltarr(c_itime_ion, 11, 72)
  daty_e_corr  = fltarr(c_itime_ele, 11, 72)
; Ions
  FOR it=0, c_itime_ion-1 DO BEGIN
    FOR ip=0, 10 DO BEGIN
      daty_h_corr[it,ip,*]  = interpol(daty_avg_int_H1[it,ip,*], emod_ion[it,*], daty_e_ion[it,*])
      daty_he_corr[it,ip,*] = interpol(daty_avg_int_He1[it,ip,*], emod_ion[it,*], daty_e_ion[it,*])
      daty_o_corr[it,ip,*]  = interpol(daty_avg_int_O1[it,ip,*], emod_ion[it,*], daty_e_ion[it,*])
    ; Do not extrapolate in energy
      ipos = where(daty_avg_int_H1[it,ip,*] GT 0, cpos)
      IF cpos GT 0 THEN BEGIN
        iex = where(daty_e_ion[it,*] LT emod_ion[it,ipos[0]] OR daty_e_ion[it,*] GT emod_ion[it,ipos[cpos-1]], cex)
        IF cex GT 0 THEN daty_h_corr[it,ip,iex] = 0.0 ;!Values.F_NAN
      ENDIF
      ipos = where(daty_avg_int_He1[it,ip,*] GT 0, cpos)
      IF cpos GT 0 THEN BEGIN
        iex = where(daty_e_ion[it,*] LT emod_ion[it,ipos[0]] OR daty_e_ion[it,*] GT emod_ion[it,ipos[cpos-1]], cex)
        IF cex GT 0 THEN daty_he_corr[it,ip,iex] = 0.0 ;!Values.F_NAN
      ENDIF
      ipos = where(daty_avg_int_O1[it,ip,*] GT 0, cpos)
      IF cpos GT 0 THEN BEGIN
        iex = where(daty_e_ion[it,*] LT emod_ion[it,ipos[0]] OR daty_e_ion[it,*] GT emod_ion[it,ipos[cpos-1]], cex)
        IF cex GT 0 THEN daty_o_corr[it,ip,iex] = 0.0 ;!Values.F_NAN
      ENDIF
    ENDFOR
  ENDFOR
; Electrons
  FOR it=0, c_itime_ele-1 DO BEGIN
    FOR ip=0, 10 DO BEGIN
      daty_e_corr[it,ip,*] = interpol(daty_avg_int_e1[it,ip,*], emod_ele[it,*], daty_e_ele[it,*])
    ; Do not extrapolate in energy
      ipos = where(daty_avg_int_e1[it,ip,*] GT 0, cpos)
      IF cpos GT 0 THEN BEGIN
        iex = where(daty_e_ele[it,*] LT emod_ele[it,ipos[0]] OR daty_e_ele[it,*] GT emod_ele[it,ipos[cpos-1]], cex)
        IF cex GT 0 THEN daty_e_corr[it,ip,iex] = 0.0 ;!Values.F_NAN
      ENDIF
    ENDFOR
  ENDFOR

  daty_h  = daty_h_corr
  daty_he = daty_he_corr
  daty_o  = daty_o_corr
  daty_e  = daty_e_corr

  ENDIF

;----------------------------------------------------------------------
; Calculate pressure
;----------------------------------------------------------------------
; Use cgs units (to match the differential flux units: 1/s-cm2-sr-keV)
; Speed of light, c
  c = 2.9979e10                  ; cm/s

; Mass
  mass_ele = 9.11e-28            ; grams
  mass_pro = 1.67e-24            ; grams
  mass_hel = 4*mass_pro          ; grams
  mass_oxy = 16*mass_pro         ; grams

; Energy
  ev_to_erg  = 1.60219e-12              ; erg/eV
  en_ion_erg = daty_e_ion*ev_to_erg     ; erg
  en_ele_erg = daty_e_ele*ev_to_erg     ; erg

  en_ion_kev = daty_e_ion/1e3    ; keV
  en_ele_kev = daty_e_ele/1e3    ; keV

  del_ion_en = fltarr(c_itime_ion,72)
  del_ele_en = fltarr(c_itime_ele,72)

  FOR it=0, c_itime_ion-1 DO BEGIN
    FOR ie=1, 70 DO del_ion_en[it,ie] = (en_ion_kev[it,ie+1]-en_ion_kev[it,ie-1])/2.0
    del_ion_en[it,0]  = en_ion_kev[it,1]-en_ion_kev[it,0]
    del_ion_en[it,71] = en_ion_kev[it,71]-en_ion_kev[it,70]
  ENDFOR
  FOR it=0, c_itime_ele-1 DO BEGIN
    FOR ie=1, 70 DO del_ele_en[it,ie] = (en_ele_kev[it,ie+1]-en_ele_kev[it,ie-1])/2.0
    del_ele_en[it,0]  = en_ele_kev[it,1]-en_ele_kev[it,0]
    del_ele_en[it,71] = en_ele_kev[it,71]-en_ele_kev[it,70]
  ENDFOR

; Gamma
  gam_pro = 1 + en_ion_erg/(mass_pro*c^2)
  gam_hel = 1 + en_ion_erg/(mass_hel*c^2)
  gam_oxy = 1 + en_ion_erg/(mass_oxy*c^2)
  gam_ele = 1 + en_ele_erg/(mass_ele*c^2)

; Pitch angle
  pa_arr = dd.PITCH_ANGLE.data*!pi/180.0          ; rad
  npa    = n_elements(pa_arr)
  del_pa = fltarr(npa)
  FOR ii=1, npa-2 DO del_pa[ii] = (pa_arr[ii+1]-pa_arr[ii-1])/2.0
  del_pa[0]     = pa_arr[1]-pa_arr[0]
  del_pa[npa-1] = pa_arr[npa-1]-pa_arr[npa-2]

; Create variables
  hperp  = fltarr(c_itime_ion,nechi)
  hpara  = fltarr(c_itime_ion,nechi)
  heperp = fltarr(c_itime_ion,nechi)
  hepara = fltarr(c_itime_ion,nechi)
  operp  = fltarr(c_itime_ion,nechi)
  opara  = fltarr(c_itime_ion,nechi)
  eperp  = fltarr(c_itime_ele,neche)
  epara  = fltarr(c_itime_ele,neche)

  p_para_h  = dblarr(c_itime_ion)
  p_perp_h  = dblarr(c_itime_ion)
  p_para_he = dblarr(c_itime_ion)
  p_perp_he = dblarr(c_itime_ion)
  p_para_o  = dblarr(c_itime_ion)
  p_perp_o  = dblarr(c_itime_ion)
  p_para_e  = dblarr(c_itime_ele)
  p_perp_e  = dblarr(c_itime_ele)

; Pressure calculation
  ang_perp = 0.5*sin(pa_arr)^3
  ang_para = sin(pa_arr)*cos(pa_arr)^2

; ions
  ;FOR it=0, c_itime_ion-1 DO BEGIN
  ;  FOR ie=en_ind_ion[0], en_ind_ion[nechi-1] DO BEGIN
  ;    FOR ip=0, npa-1 DO BEGIN
  ;      del_ion = del_ion_en[it,ie]*del_pa[ip]
  ;      ang_perp = 0.5*sin(pa_arr[ip])^3
  ;      ang_para = sin(pa_arr[ip])*cos(pa_arr[ip])^2
  ;      ;IF finite(daty_avg_int_O1[it,ip,ie]) EQ 0 THEN daty_avg_int_O1[it,ip,ie] = 0.0
  ;      p_perp_h[it]  = p_perp_h[it] + 2*!pi*sqrt(2*mass_pro*en_ion_erg[it,ie])*daty_avg_int_H1[it,ip,ie]*ang_perp*del_ion
  ;      p_para_h[it]  = p_para_h[it] + 2*!pi*sqrt(2*mass_pro*en_ion_erg[it,ie])*daty_avg_int_H1[it,ip,ie]*ang_para*del_ion
  ;      p_perp_he[it] = p_perp_he[it] + 2*!pi*sqrt(2*mass_hel*en_ion_erg[it,ie])*daty_avg_int_He1[it,ip,ie]*ang_perp*del_ion
  ;      p_para_he[it] = p_para_he[it] + 2*!pi*sqrt(2*mass_hel*en_ion_erg[it,ie])*daty_avg_int_He1[it,ip,ie]*ang_para*del_ion
  ;      p_perp_o[it]  = p_perp_o[it] + 2*!pi*sqrt(2*mass_oxy*en_ion_erg[it,ie])*daty_avg_int_O1[it,ip,ie]*ang_perp*del_ion
  ;      p_para_o[it]  = p_para_o[it] + 2*!pi*sqrt(2*mass_oxy*en_ion_erg[it,ie])*daty_avg_int_O1[it,ip,ie]*ang_para*del_ion
  ;    ENDFOR
  ;  ENDFOR
  ;ENDFOR

  FOR it=0, c_itime_ion-1 DO BEGIN
    FOR ie=en_ind_ion[0], en_ind_ion[nechi-1] DO BEGIN
      hperp[it,ie-en_ind_ion[0]]  = total(daty_h[it,*,ie]*ang_perp*del_pa, /NAN)
      hpara[it,ie-en_ind_ion[0]]  = total(daty_h[it,*,ie]*ang_para*del_pa, /NAN)
      heperp[it,ie-en_ind_ion[0]] = total(daty_he[it,*,ie]*ang_perp*del_pa, /NAN)
      hepara[it,ie-en_ind_ion[0]] = total(daty_he[it,*,ie]*ang_para*del_pa, /NAN)
      operp[it,ie-en_ind_ion[0]]  = total(daty_o[it,*,ie]*ang_perp*del_pa, /NAN)
      opara[it,ie-en_ind_ion[0]]  = total(daty_o[it,*,ie]*ang_para*del_pa, /NAN)
    ENDFOR
  ENDFOR

  FOR it=0, c_itime_ion-1 DO BEGIN
    IF relat EQ 0 THEN BEGIN
      p_perp_h[it]  = 2*!pi*total(sqrt(2*mass_pro*en_ion_erg[it,en_ind_ion])*hperp[it,*]*del_ion_en[it,en_ind_ion], /NAN);, /NAN
      p_para_h[it]  = 2*!pi*total(sqrt(2*mass_pro*en_ion_erg[it,en_ind_ion])*hpara[it,*]*del_ion_en[it,en_ind_ion], /NAN)
      p_perp_he[it] = 2*!pi*total(sqrt(2*mass_hel*en_ion_erg[it,en_ind_ion])*heperp[it,*]*del_ion_en[it,en_ind_ion], /NAN)
      p_para_he[it] = 2*!pi*total(sqrt(2*mass_hel*en_ion_erg[it,en_ind_ion])*hepara[it,*]*del_ion_en[it,en_ind_ion], /NAN)
      p_perp_o[it]  = 2*!pi*total(sqrt(2*mass_oxy*en_ion_erg[it,en_ind_ion])*operp[it,*]*del_ion_en[it,en_ind_ion], /NAN)
      p_para_o[it]  = 2*!pi*total(sqrt(2*mass_oxy*en_ion_erg[it,en_ind_ion])*opara[it,*]*del_ion_en[it,en_ind_ion], /NAN)
    ENDIF
    IF relat EQ 1 THEN BEGIN
      p_perp_h[it]  = 2*!pi*total(mass_pro*c*sqrt(gam_pro[it,en_ind_ion]^2-1)*hperp[it,*]*del_ion_en[it,en_ind_ion], /NAN)
      p_para_h[it]  = 2*!pi*total(mass_pro*c*sqrt(gam_pro[it,en_ind_ion]^2-1)*hpara[it,*]*del_ion_en[it,en_ind_ion], /NAN)
      p_perp_he[it] = 2*!pi*total(mass_hel*c*sqrt(gam_hel[it,en_ind_ion]^2-1)*heperp[it,*]*del_ion_en[it,en_ind_ion], /NAN)
      p_para_he[it] = 2*!pi*total(mass_hel*c*sqrt(gam_hel[it,en_ind_ion]^2-1)*hepara[it,*]*del_ion_en[it,en_ind_ion], /NAN)
      p_perp_o[it]  = 2*!pi*total(mass_oxy*c*sqrt(gam_oxy[it,en_ind_ion]^2-1)*operp[it,*]*del_ion_en[it,en_ind_ion], /NAN)
      p_para_o[it]  = 2*!pi*total(mass_oxy*c*sqrt(gam_oxy[it,en_ind_ion]^2-1)*opara[it,*]*del_ion_en[it,en_ind_ion], /NAN)
    ENDIF
  ENDFOR

; electrons
  ;FOR it=0, c_itime_ele-1 DO BEGIN
  ;  FOR ie=en_ind_ele[0], en_ind_ele[neche-1] DO BEGIN
  ;    FOR ip=0, npa-1 DO BEGIN
  ;      del_ele = del_ele_en[it,ie]*del_pa[ip]
  ;      ang_perp = 0.5*sin(pa_arr[ip])^3
  ;      ang_para = sin(pa_arr[ip])*cos(pa_arr[ip])^2
  ;      p_perp_e[it] = p_perp_e[it] + 2*!pi*sqrt(2*mass_ele*en_ele_erg[it,ie])*daty_avg_int_e1[it,ip,ie]*ang_perp*del_ele
  ;      p_para_e[it] = p_para_e[it] + 2*!pi*sqrt(2*mass_ele*en_ele_erg[it,ie])*daty_avg_int_e1[it,ip,ie]*ang_para*del_ele
  ;    ENDFOR
  ;  ENDFOR
  ;ENDFOR

  FOR it=0, c_itime_ele-1 DO BEGIN
    FOR ie=en_ind_ele[0], en_ind_ele[neche-1] DO BEGIN
      eperp[it,ie-en_ind_ele[0]] = total(daty_e[it,*,ie]*ang_perp*del_pa, /NAN)
      epara[it,ie-en_ind_ele[0]] = total(daty_e[it,*,ie]*ang_para*del_pa, /NAN)
    ENDFOR
  ENDFOR

  FOR it=0, c_itime_ele-1 DO BEGIN
    IF relat EQ 0 THEN BEGIN
      p_perp_e[it] = 2*!pi*total(sqrt(2*mass_ele*en_ele_erg[it,en_ind_ele])*eperp[it,*]*del_ele_en[it,en_ind_ele], /NAN);, /NAN
      p_para_e[it] = 2*!pi*total(sqrt(2*mass_ele*en_ele_erg[it,en_ind_ele])*epara[it,*]*del_ele_en[it,en_ind_ele], /NAN)
    ENDIF
    IF relat EQ 1 THEN BEGIN
      p_perp_e[it] = 2*!pi*total(mass_ele*c*sqrt(gam_ele[it,en_ind_ele]^2-1)*eperp[it,*]*del_ele_en[it,en_ind_ele], /NAN)
      p_para_e[it] = 2*!pi*total(mass_ele*c*sqrt(gam_ele[it,en_ind_ele]^2-1)*epara[it,*]*del_ele_en[it,en_ind_ele], /NAN)
    ENDIF
  ENDFOR

; Convert units
  p_perp_h  = p_perp_h/10.0*1e9     ; Ba (dyne/cm^2) --> nPa
  p_para_h  = p_para_h/10.0*1e9     ; Ba (dyne/cm^2) --> nPa
  p_perp_he = p_perp_he/10.0*1e9    ; Ba (dyne/cm^2) --> nPa
  p_para_he = p_para_he/10.0*1e9    ; Ba (dyne/cm^2) --> nPa
  p_perp_o  = p_perp_o/10.0*1e9     ; Ba (dyne/cm^2) --> nPa
  p_para_o  = p_para_o/10.0*1e9     ; Ba (dyne/cm^2) --> nPa
  p_perp_e  = p_perp_e/10.0*1e9     ; Ba (dyne/cm^2) --> nPa
  p_para_e  = p_para_e/10.0*1e9     ; Ba (dyne/cm^2) --> nPa

  p_perp_h_sm  = p_perp_h
  p_para_h_sm  = p_para_h
  p_perp_he_sm = p_perp_he
  p_para_he_sm = p_para_he
  p_perp_o_sm  = p_perp_o
  p_para_o_sm  = p_para_o
  p_perp_e_sm  = p_perp_e
  p_para_e_sm  = p_para_e

;----------------------------------------------------------------------
; Smooth data
;----------------------------------------------------------------------
  IF swindow GT 0 THEN BEGIN
    dt = tt[1]-tt[0]
    ti = dindgen(dt/swindow)*swindow + tt[0]

    dat_t_ion_nogap = dat_t_ion
    dat_t_ele_nogap = dat_t_ele
    p_perp_h_nogap  = p_perp_h
    p_para_h_nogap  = p_para_h
    p_perp_he_nogap = p_perp_he
    p_para_he_nogap = p_para_he
    p_perp_o_nogap  = p_perp_o
    p_para_o_nogap  = p_para_o
    p_perp_e_nogap  = p_perp_e
    p_para_e_nogap  = p_para_e

  ; Look for data gaps to fill them with NANs
  ; ions
    tdif = ts_diff(dat_t_ion, 1)*(-1)
    indgap = where(tdif GT 1e2, ngap)
    IF ngap GT 0 THEN BEGIN
      dat_t_ion_nogap = dat_t_ion[0:indgap[0]]
      p_perp_h_nogap  = p_perp_h[0:indgap[0]]
      p_para_h_nogap  = p_para_h[0:indgap[0]]
      p_perp_he_nogap = p_perp_he[0:indgap[0]]
      p_para_he_nogap = p_para_he[0:indgap[0]]
      p_perp_o_nogap  = p_perp_o[0:indgap[0]]
      p_para_o_nogap  = p_para_o[0:indgap[0]]
      FOR ig=0, ngap-1 DO BEGIN
        idum = indgap[ig]
        lgap = tdif[idum]
        nbin = floor(lgap/22.0)          ; default time cadence of 22 seconds
        tgap = (dindgen(nbin-1)+1)*22 + dat_t_ion[idum]
        pgap = make_array(nbin-1, /FLOAT, value=!VALUES.F_NAN)
        IF ig LT ngap-1 THEN BEGIN
          dat_t_ion_nogap = [dat_t_ion_nogap, tgap, dat_t_ion[idum+1:indgap[ig+1]]]
          p_perp_h_nogap  = [p_perp_h_nogap, pgap, p_perp_h[idum+1:indgap[ig+1]]]
          p_para_h_nogap  = [p_para_h_nogap, pgap, p_para_h[idum+1:indgap[ig+1]]]
          p_perp_he_nogap = [p_perp_he_nogap, pgap, p_perp_he[idum+1:indgap[ig+1]]]
          p_para_he_nogap = [p_para_he_nogap, pgap, p_para_he[idum+1:indgap[ig+1]]]
          p_perp_o_nogap  = [p_perp_o_nogap, pgap, p_perp_o[idum+1:indgap[ig+1]]]
          p_para_o_nogap  = [p_para_o_nogap, pgap, p_para_o[idum+1:indgap[ig+1]]]
        ENDIF ELSE BEGIN
          dat_t_ion_nogap = [dat_t_ion_nogap, tgap, dat_t_ion[idum+1:-1]]
          p_perp_h_nogap  = [p_perp_h_nogap, pgap, p_perp_h[idum+1:-1]]
          p_para_h_nogap  = [p_para_h_nogap, pgap, p_para_h[idum+1:-1]]
          p_perp_he_nogap = [p_perp_he_nogap, pgap, p_perp_he[idum+1:-1]]
          p_para_he_nogap = [p_para_he_nogap, pgap, p_para_he[idum+1:-1]]
          p_perp_o_nogap  = [p_perp_o_nogap, pgap, p_perp_o[idum+1:-1]]
          p_para_o_nogap  = [p_para_o_nogap, pgap, p_para_o[idum+1:-1]]
        ENDELSE
      ENDFOR
    ENDIF
  ; electrons
    tdif = ts_diff(dat_t_ele, 1)*(-1)
    indgap = where(tdif GT 1e2, ngap)
    IF ngap GT 0 THEN BEGIN
      dat_t_ele_nogap = dat_t_ele[0:indgap[0]]
      p_perp_e_nogap  = p_perp_e[0:indgap[0]]
      p_para_e_nogap  = p_para_e[0:indgap[0]]
      FOR ig=0, ngap-1 DO BEGIN
        idum = indgap[ig]
        lgap = tdif[idum]
        nbin = floor(lgap/22.0)          ; default time cadence of 22 seconds
        tgap = (dindgen(nbin-1)+1)*22 + dat_t_ele[idum]
        pgap = make_array(nbin-1, /FLOAT, value=!VALUES.F_NAN)
        IF ig LT ngap-1 THEN BEGIN
          dat_t_ele_nogap = [dat_t_ele_nogap, tgap, dat_t_ele[idum+1:indgap[ig+1]]]
          p_perp_e_nogap  = [p_perp_e_nogap, pgap, p_perp_e[idum+1:indgap[ig+1]]]
          p_para_e_nogap  = [p_para_e_nogap, pgap, p_para_e[idum+1:indgap[ig+1]]]
        ENDIF ELSE BEGIN
          dat_t_ele_nogap = [dat_t_ele_nogap, tgap, dat_t_ele[idum+1:-1]]
          p_perp_e_nogap  = [p_perp_e_nogap, pgap, p_perp_e[idum+1:-1]]
          p_para_e_nogap  = [p_para_e_nogap, pgap, p_para_e[idum+1:-1]]
        ENDELSE
      ENDFOR
    ENDIF

  ; Look for a data gap at the end to fill it with NANs
  ; ions
    dtend = tt[1]-dat_t_ion_nogap[-1]
    IF dtend GT 1e2 THEN BEGIN
      nbin_end = floor(dtend/22.0)          ; default time cadence of 22 seconds
      tgap_end = (dindgen(nbin_end-1)+1)*22 + dat_t_ion_nogap[-1]
      pgap_end = make_array(nbin_end-1, /FLOAT, value=!VALUES.F_NAN)
      dat_t_ion_nogap = [dat_t_ion_nogap, tgap_end]
      p_perp_h_nogap  = [p_perp_h_nogap, pgap_end]
      p_para_h_nogap  = [p_para_h_nogap, pgap_end]
      p_perp_he_nogap = [p_perp_he_nogap, pgap_end]
      p_para_he_nogap = [p_para_he_nogap, pgap_end]
      p_perp_o_nogap  = [p_perp_o_nogap, pgap_end]
      p_para_o_nogap  = [p_para_o_nogap, pgap_end]
    ENDIF
  ; electrons
    dtend = tt[1]-dat_t_ele_nogap[-1]
    IF dtend GT 1e2 THEN BEGIN
      nbin_end = floor(dtend/22.0)          ; default time cadence of 22 seconds
      tgap_end = (dindgen(nbin_end-1)+1)*22 + dat_t_ele_nogap[-1]
      pgap_end = make_array(nbin_end-1, /FLOAT, value=!VALUES.F_NAN)
      dat_t_ele_nogap = [dat_t_ele_nogap, tgap_end]
      p_perp_e_nogap  = [p_perp_e_nogap, pgap_end]
      p_para_e_nogap  = [p_para_e_nogap, pgap_end]
    ENDIF

  ; Regrid time with swindow cadence
    p_perp_h_sm  = interpol(p_perp_h_nogap, dat_t_ion_nogap, ti);, /NAN
    p_para_h_sm  = interpol(p_para_h_nogap, dat_t_ion_nogap, ti)
    p_perp_he_sm = interpol(p_perp_he_nogap, dat_t_ion_nogap, ti)
    p_para_he_sm = interpol(p_para_he_nogap, dat_t_ion_nogap, ti)
    p_perp_o_sm  = interpol(p_perp_o_nogap, dat_t_ion_nogap, ti)
    p_para_o_sm  = interpol(p_para_o_nogap, dat_t_ion_nogap, ti)
    p_perp_e_sm  = interpol(p_perp_e_nogap, dat_t_ele_nogap, ti)
    p_para_e_sm  = interpol(p_para_e_nogap, dat_t_ele_nogap, ti)

    dat_mode_ion_sm = interpol(dat_mode_ion, dat_t_ion, ti)
    dat_mode_ele_sm = interpol(dat_mode_ele, dat_t_ele, ti)
    dat_mode_ele_sm = interpol(dat_mode_ele, dat_t_ele, ti)

  ; Do not extrapolate
    iex = where(ti LT dat_t_ion_nogap[0] OR ti GT dat_t_ion_nogap[-1], cex)
    IF cex GT 0 THEN BEGIN
      p_perp_h_sm[iex]  = !Values.F_NAN
      p_para_h_sm[iex]  = !Values.F_NAN
      p_perp_he_sm[iex] = !Values.F_NAN
      p_para_he_sm[iex] = !Values.F_NAN
      p_perp_o_sm[iex]  = !Values.F_NAN
      p_para_o_sm[iex]  = !Values.F_NAN
    ENDIF
    iex = where(ti LT dat_t_ele_nogap[0] OR ti GT dat_t_ele_nogap[-1], cex)
    IF cex GT 0 THEN BEGIN
      p_perp_e_sm[iex] = !Values.F_NAN
      p_para_e_sm[iex] = !Values.F_NAN
    ENDIF
  ENDIF

;----------------------------------------------------------------------
; Replace zeros with NaNs
;----------------------------------------------------------------------
  inan = where(p_perp_h_sm EQ 0.0, cnan)
  p_perp_h_sm[inan] = !Values.F_NAN
  inan = where(p_para_h_sm EQ 0.0, cnan)
  p_para_h_sm[inan] = !Values.F_NAN
  inan = where(p_perp_he_sm EQ 0.0, cnan)
  p_perp_he_sm[inan] = !Values.F_NAN
  inan = where(p_para_he_sm EQ 0.0, cnan)
  p_para_he_sm[inan] = !Values.F_NAN
  inan = where(p_perp_o_sm EQ 0.0, cnan)
  p_perp_o_sm[inan] = !Values.F_NAN
  inan = where(p_para_o_sm EQ 0.0, cnan)
  p_para_o_sm[inan] = !Values.F_NAN

  inan = where(p_perp_e_sm EQ 0.0, cnan)
  p_perp_e_sm[inan] = !Values.F_NAN
  inan = where(p_para_e_sm EQ 0.0, cnan)
  p_para_e_sm[inan] = !Values.F_NAN

;----------------------------------------------------------------------
; Average pressures
;----------------------------------------------------------------------
  p_avg_h  = (2.0 * p_perp_h_sm + p_para_h_sm) / 3.0
  p_avg_he = (2.0 * p_perp_he_sm + p_para_he_sm) / 3.0
  p_avg_o  = (2.0 * p_perp_o_sm + p_para_o_sm) / 3.0
  p_avg_e  = (2.0 * p_perp_e_sm + p_para_e_sm) / 3.0

  IF swindow EQ 0 THEN BEGIN
    p_avg_e_int = interpol(p_avg_e, dat_t_ele, dat_t_ion, /NAN)
    p_avg_t = p_avg_h + p_avg_he + p_avg_o + p_avg_e_int
    p_species = [[p_avg_t],[p_avg_h],[p_avg_he],[p_avg_o],[p_avg_e_int]]
  ENDIF ELSE BEGIN
    p_avg_t = p_avg_h + p_avg_he + p_avg_o + p_avg_e
    p_species = [[p_avg_t],[p_avg_h],[p_avg_he],[p_avg_o],[p_avg_e]]
  ENDELSE

; Create vector arrays with Pperp, Ppara, Pavg, and Pavg for all species
  p_h  = [[p_perp_h_sm],[p_para_h_sm],[p_avg_h]]
  p_he = [[p_perp_he_sm],[p_para_he_sm],[p_avg_he]]
  p_o  = [[p_perp_o_sm],[p_para_o_sm],[p_avg_o]]
  p_e  = [[p_perp_e_sm],[p_para_e_sm],[p_avg_e]]

;----------------------------------------------------------------------
; Create TPLOT variables for the ion and electron data
;----------------------------------------------------------------------
  erange_ion = string(low_energy[0]/1e3, '(F0.1)')+'-'+string(up_energy[0]/1e3, '(F0.1)')
  erange_ele = string(low_energy[1]/1e3, '(F0.1)')+'-'+string(up_energy[1]/1e3, '(F0.1)')

  estr_ion = string(eion0, '(I0.1)')+'eV-'+string(eion1/1e3, '(I0.1)')+'keV'
  estr_ele = string(eele0, '(I0.1)')+'eV-'+string(eele1/1e3, '(I0.1)')+'keV'

  IF relat EQ 0 THEN relstr = '_nonrel' ELSE relstr = '_rel'

  IF swindow EQ 0 THEN BEGIN
    tion = dat_t_ion
    tele = dat_t_ele
  ENDIF ELSE BEGIN
    tion = ti
    tele = ti
  ENDELSE

FOR ivar = 0, nvar-1 DO BEGIN

  nmh = 'rbsp'+probe+'_hope_'
  var = nmh+var_arr[ivar]+relstr

  IF ivar EQ 0 THEN store_data, var+'_'+erange_ion+'_perp', dat={x:tion, y:p_perp_h_sm}, dlim={ytitle:nm_arr[ivar]+'!U+!N Pperp (nPa)!C'+estr_ion}
  IF ivar EQ 0 THEN store_data, var+'_'+erange_ion+'_para', dat={x:tion, y:p_para_h_sm}, dlim={ytitle:nm_arr[ivar]+'!U+!N Ppara (nPa)!C'+estr_ion}
  IF ivar EQ 0 THEN store_data, var+'_'+erange_ion+'_avg', dat={x:tion, y:p_avg_h}, dlim={ytitle:nm_arr[ivar]+'!U+!N Pavg (nPa)!C'+estr_ion}
  IF ivar EQ 0 THEN store_data, var+'_'+erange_ion, dat={x:tion, y:p_h}, dlim={ytitle:nm_arr[ivar]+'!U+!N P (nPa)!C'+estr_ion}

  IF ivar EQ 1 THEN store_data, var+'_'+erange_ion+'_perp', dat={x:tion, y:p_perp_he_sm}, dlim={ytitle:nm_arr[ivar]+'!U+!N Pperp (nPa)!C'+estr_ion}
  IF ivar EQ 1 THEN store_data, var+'_'+erange_ion+'_para', dat={x:tion, y:p_para_he_sm}, dlim={ytitle:nm_arr[ivar]+'!U+!N Ppara (nPa)!C'+estr_ion}
  IF ivar EQ 1 THEN store_data, var+'_'+erange_ion+'_avg', dat={x:tion, y:p_avg_he}, dlim={ytitle:nm_arr[ivar]+'!U+!N Pavg (nPa)!C'+estr_ion}
  IF ivar EQ 1 THEN store_data, var+'_'+erange_ion, dat={x:tion, y:p_he}, dlim={ytitle:nm_arr[ivar]+'!U+!N P (nPa)!C'+estr_ion}

  IF ivar EQ 2 THEN store_data, var+'_'+erange_ion+'_perp', dat={x:tion, y:p_perp_o_sm}, dlim={ytitle:nm_arr[ivar]+'!U+!N Pperp (nPa)!C'+estr_ion}
  IF ivar EQ 2 THEN store_data, var+'_'+erange_ion+'_para', dat={x:tion, y:p_para_o_sm}, dlim={ytitle:nm_arr[ivar]+'!U+!N Ppara (nPa)!C'+estr_ion}
  IF ivar EQ 2 THEN store_data, var+'_'+erange_ion+'_avg', dat={x:tion, y:p_avg_o}, dlim={ytitle:nm_arr[ivar]+'!U+!N Pavg (nPa)!C'+estr_ion}
  IF ivar EQ 2 THEN store_data, var+'_'+erange_ion, dat={x:tion, y:p_o}, dlim={ytitle:nm_arr[ivar]+'!U+!N P (nPa)!C'+estr_ion}

  IF ivar EQ 3 THEN store_data, var+'_'+erange_ele+'_perp', dat={x:tele, y:p_perp_e_sm}, dlim={ytitle:nm_arr[ivar]+'!U-!N Pperp (nPa)!C'+estr_ele}
  IF ivar EQ 3 THEN store_data, var+'_'+erange_ele+'_para', dat={x:tele, y:p_para_e_sm}, dlim={ytitle:nm_arr[ivar]+'!U-!N Ppara (nPa)!C'+estr_ele}
  IF ivar EQ 3 THEN store_data, var+'_'+erange_ele+'_avg', dat={x:tele, y:p_avg_e}, dlim={ytitle:nm_arr[ivar]+'!U-!N Pavg (nPa)!C'+estr_ele}
  IF ivar EQ 3 THEN store_data, var+'_'+erange_ele, dat={x:tele, y:p_e}, dlim={ytitle:nm_arr[ivar]+'!U-!N P (nPa)!C'+estr_ele}
  ;IF ivar EQ 4 THEN store_data, var, dat={x:dat_t_ion, y:daty_counts_H1,   v:daty_e_ion}, dlim=lim
  ;IF ivar EQ 5 THEN store_data, var, dat={x:dat_t_ion, y:daty_counts_He1,  v:daty_e_ion}, dlim=lim
  ;IF ivar EQ 6 THEN store_data, var, dat={x:dat_t_ion, y:daty_counts_O1,   v:daty_e_ion}, dlim=lim
  ;IF ivar EQ 7 THEN store_data, var, dat={x:dat_t_ele, y:daty_counts_e1,   v:daty_e_ele}, dlim=lim

  options, var+'*', 'x_no_interp', 1
  options, var+'*', 'y_no_interp', 1

  options, var+'_'+erange_ion, 'colors', [6,2,0]
  options, var+'_'+erange_ion, 'labels', ['Pperp','Ppara','Pavg']
  options, var+'_'+erange_ion, 'labflag', -1

  options, var+'_'+erange_ele, 'colors', [6,2,0]
  options, var+'_'+erange_ele, 'labels', ['Pperp','Ppara','Pavg']
  options, var+'_'+erange_ele, 'labflag', -1

ENDFOR                                             ; each var

  store_data, 'rbsp'+probe+'_hope_Pressure_All_species'+relstr+'_'+erange_ion, dat={x:tion, y:p_species}, dlim={ytitle:'Pressure (nPa)!C'+estr_ion}
  store_data, 'rbsp'+probe+'_MODE_ION', dat={x:tion, y:dat_mode_ion_sm}, dlim=lim
  store_data, 'rbsp'+probe+'_MODE_ELE', dat={x:tele, y:dat_mode_ele_sm}, dlim=lim
  ;store_data, 'rbsp'+probe+'_MODE_ELE_int', dat={x:tele, y:dat_mode_ele_int_sm}, dlim=lim

  options, 'rbsp'+probe+'_hope_Pressure_All_species'+relstr+'_'+erange_ion, 'colors', [0,2,3,6,4]
  options, 'rbsp'+probe+'_hope_Pressure_All_species'+relstr+'_'+erange_ion, 'labels', ['Tot','H!U+','He!U+','O!U+','e!U-']
  options, 'rbsp'+probe+'_hope_Pressure_All_species'+relstr+'_'+erange_ion, 'labflag', -1

;tplot,'*'


END
