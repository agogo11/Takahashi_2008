PRO yzreprod_fig4_takahashi2008

  ;.RESET_SESSION
  COMMON share, tau, tau0, alphabl, alphafa,tau_lcl
  path = 'G:\MyIDLpackage\RadCooling_Takahashi2008\'

  ;Reference parameters
  sigma = 5.67e-8 ; stefan-bolzmann constant Wm-2K-4
  RdCp    =   0.286;      % -             Ratio of Rd/Cp
  g=9.80665d ; gravity const.
  BOLTZMAN = 1.380658e-23;            % BOLTZMAN,  Boltzman constant (J/K)
  AVOGADRO = .602214199e24;           % AVOGADRO, Avogadro constant (1/mol)
  MD       = 28.9644e-3;              % MD, molar mass dry air (kg/mol)
  MV  = 18.0153e-3;                   % MD, molae mass water vapor (kg/mol)
  r_v = (AVOGADRO)*(BOLTZMAN)/(MV); % r_v, gas constant for water vapor (J/(kg-K))
  r_d = (AVOGADRO)*(BOLTZMAN)/(MD); % r_d, gas constant for dry air (J/(kg-K))
  cp  = 7./2.*(r_d);                   % cp, specific heat of air (J/(kg-K))

  ;simulation-specific parameters
  Ps = 1000. ; hPa
  thetas = 300.; K surface air (potential) temeprature
  n = 4 ; exponent that relates p to tau
  alphabl = 4.*R_d/(cp*n)


  ;simulation specific arrays
  gamma_theta_array = [6.5, 6.5, 5.]/1000.  ;K/m
  p_lcl_array = [900., 950., 900.]  ;hPa
  nrun = N_ELEMENTS(gamma_theta_array)

  ;basic-state column tau of the atmosphere
  tau0_array = arrgen(1, 100., 1.)
  ntau0 = N_ELEMENTS(tau0_array)


  ;initialize plot/output variables to zero FOR new simulation
  epsilon_lcl_plot       = FLTARR(nrun, ntau0);
  epsilon_sfc_plot       = FLTARR(nrun, ntau0);
  epsilon_toa_plot       = FLTARR(nrun, ntau0);

  ;; Loop over values tau0
  FOR irun = 0, nrun - 1 DO BEGIN
    gamma_theta = gamma_theta_array[irun]
    p_lcl = p_lcl_array[irun]
    alphafa = 4.*R_d*gamma_theta/(g*n)

    FOR itau0 = 0, ntau0 - 1 DO BEGIN
      tau0 = tau0_array[itau0]
      tau_LCL = tau0*(P_lcl/Ps)^4.

      ;1. cal radiative fluxes at surface
      tau = tau0

      ;cal upwelling fluxes
      IF tau GT tau_lcl THEN BEGIN
        inteR_atmos = QROMB('f_u_bl',tau,tau0)
        U = sigma*(thetas^4.)*(EXP(-(tau0 - tau))) + sigma*(thetas^4.)*inteR_atmos
      ENDIF ELSE BEGIN
        inteR_atmos = QROMB('f_u_fa',tau,tau_lcl)
        inteR_atmos_bl = QROMB('f_u_bl',tau_lcl,tau0)
        U = sigma*(thetas^4.)*(EXP(-(tau0 - tau))) + sigma*(thetas^4.)*(inteR_atmos + inteR_atmos_bl)
      ENDELSE

      ;cal downwelling fluxes
      IF tau GT tau_lcl THEN BEGIN
        inteR_atmos = QROMB('f_d_bl',tau_lcl, tau)
        inteR_atmos_fa = QROMB('f_d_fa',0, tau_lcl)
        D = sigma*(thetas^4.)*(inteR_atmos + inteR_atmos_fa)
      ENDIF ELSE BEGIN
        inteR_atmos = QROMB('f_d_fa',0., tau);
        D = sigma*(thetas^4.)*inteR_atmos
      ENDELSE

      epsilon_sfc_plot[irun, itau0] = (U-D)/(sigma*(thetas^4.))

      ;2. cal radiative fluxes at LCL
      tau = tau_lcl

      ;cal upwelling fluxes
      IF tau GT tau_lcl THEN BEGIN
        inteR_atmos = QROMB('f_u_bl',tau,tau0)
        U = sigma*(thetas^4.)*(EXP(-(tau0 - tau))) + sigma*(thetas^4.)*inteR_atmos
      ENDIF ELSE BEGIN
        inteR_atmos = QROMB('f_u_fa',tau,tau_lcl)
        inteR_atmos_bl = QROMB('f_u_bl',tau_lcl,tau0)
        U = sigma*(thetas^4.)*(EXP(-(tau0 - tau))) + sigma*(thetas^4.)*(inteR_atmos + inteR_atmos_bl)
      ENDELSE

      ;cal downwelling fluxes
      IF tau GT tau_lcl THEN BEGIN
        inteR_atmos = QROMB('f_d_bl',tau_lcl, tau)
        inteR_atmos_fa = QROMB('f_d_fa',0, tau_lcl)
        D = sigma*(thetas^4.)*(inteR_atmos + inteR_atmos_fa)
      ENDIF ELSE BEGIN
        inteR_atmos = QROMB('f_d_fa',0., tau);
        D = sigma*(thetas^4.)*inteR_atmos
      ENDELSE

      epsilon_lcl_plot[irun, itau0] = (U-D)/(sigma*(thetas^4.))

      ;3. cal radiative fluxes at LCL
      tau = 0.

      ;cal upwelling fluxes
      IF tau GT tau_lcl THEN BEGIN
        inteR_atmos = QROMB('f_u_bl',tau,tau0)
        U = sigma*(thetas^4.)*(EXP(-(tau0 - tau))) + sigma*(thetas^4.)*inteR_atmos
      ENDIF ELSE BEGIN
        inteR_atmos = QROMB('f_u_fa',tau,tau_lcl)
        inteR_atmos_bl = QROMB('f_u_bl',tau_lcl,tau0)
        U = sigma*(thetas^4.)*(EXP(-(tau0 - tau))) + sigma*(thetas^4.)*(inteR_atmos + inteR_atmos_bl)
      ENDELSE

      ;cal downwelling fluxes
      IF tau GT tau_lcl THEN BEGIN
        inteR_atmos = QROMB('f_d_bl',tau_lcl, tau)
        inteR_atmos_fa = QROMB('f_d_fa',0, tau_lcl)
        D = sigma*(thetas^4.)*(inteR_atmos + inteR_atmos_fa)
      ENDIF ELSE BEGIN
        inteR_atmos = QROMB('f_d_fa',0., tau);
        D = sigma*(thetas^4.)*inteR_atmos
      ENDELSE

      epsilon_toa_plot[irun, itau0] = (U-D)/(sigma*(thetas^4.))
    ENDFOR

  ENDFOR

  ;plot
  cgps_open, Path + 'Fig4_Takahashi2008.eps',/CMYK,$
    /nomatch,Font = 1,/Times, xsize =10/2.54, ysize = 10/2.54, fontsize =14

  pos = cglayout([1,1], OXMargin=[5,1], OYMargin=[4,1])

  mycharsize = 1.
  myxrange = [1., 100.]
  myyrange = [0., 0.9]

  myxtitle = cggreek('tau') + '0'
  myytitle = cggreek('epsilon')

  mycolors = ['TG7','TG4','TG1']
  mylinestyle = [0,1,2]

  location_lg = [5., 0.85]
  text_lg = [cggreek('gamma') + '=6.5 K/km, P$\subLCL$ = 900 mb',$
    cggreek('gamma') + '=6.5 K/km, P$\subLCL$ = 950 mb',$
    cggreek('gamma') + '=5.0 K/km, P$\subLCL$ = 900 mb']

  linestyle_lg = [0,1,2]

  ;Fig 1
  cgplot,tau0_array, epsilon_lcl_plot[0,*],/noerase,/nodata, /xlog,$
    xtitle = myxtitle, ytitle = myytitle,$
    xrange = myxrange,yrange = myyrange,$
    ytickinterval = 0.1,$
    charsize = mycharsize, POSITION = pos[*,0]

  FOR irun = 0, nrun - 1 DO BEGIN
    cgoplot,tau0_array, epsilon_sfc_plot[irun,*], thick = 2., color = mycolors[0], linestyle = mylinestyle[irun]
    cgoplot,tau0_array, epsilon_lcl_plot[irun,*], thick = 2., color = mycolors[1], linestyle = mylinestyle[irun]
    cgoplot,tau0_array, epsilon_toa_plot[irun,*], thick = 2., color = mycolors[2], linestyle = mylinestyle[irun]
  ENDFOR

  cglegend, title = text_lg, linestyle = linestyle_lg, Location = location_lg, /data,$
    vspace = 0.8, charsize = 0.7*mycharsize

  location_lg = [20., 0.7]

  text_lg = [cggreek('epsilon') + '$\subsfc$',$
    cggreek('epsilon') + '$\subLCL$',$
    cggreek('epsilon') + '$\subTOA$']

  cglegend, title = text_lg, color = mycolors, Location = location_lg, /data,$
    vspace = 0.8, length = 0., charsize = 0.9*mycharsize, tcolors = mycolors

  cgps_close, /png, density = 1000

  PRINT, ' '

END
