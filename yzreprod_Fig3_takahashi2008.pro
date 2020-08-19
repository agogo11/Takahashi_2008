PRO yzreprod_Fig3_takahashi2008

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
  P_lcl = 900. ; hPa
  Ps = 1000. ; hPa
  gamma_theta = 6.5/1000. ; K/m
  thetas = 300.; K surface air (potential) temeprature
  n = 4 ; exponent that relates p to tau
  alphabl = 4.*R_d/(cp*n)
  alphafa = 4.*R_d*gamma_theta/(g*n)

  ;simulation specific arrays
  ;basic-state column tau of the atmosphere
  tau0_array = [100., 25., 6., 2., 0.5]
  ntau0 = N_ELEMENTS(tau0_array)

  tau_LCL_array = tau0_array*(P_lcl/Ps)^4.

  p_array = arrgen(ps - 0.5, 0.5, 10.)
  np = N_ELEMENTS(p_array)

  ;initialize plot/output variables to zero FOR new simulation
  ntau = np
  p_tau_eq_1_plot       = FLTARR(ntau0);
  epsilon_tau_eq_1_plot       = FLTARR(ntau0);
  epsilon_lcl_plot       = FLTARR(ntau0);

  epsilon_plot = FLTARR(ntau0, ntau)
  epsilon_plot1 = FLTARR(ntau0, ntau)
  tau_plot = FLTARR(ntau0, ntau)
  p_plot = FLTARR(ntau0, ntau)

  ;; Loop over values tau0
  FOR itau0 = 0, ntau0 - 1 DO BEGIN
    tau0 = tau0_array[itau0]
    tau_LCL = tau_LCL_array[itau0]

    ;cal pressure with tau = 1
    p_tau_eq_1_plot[itau0] = ps*((1./tau0)^(1./n))

    FOR ip = 0, np - 1 DO BEGIN
      tau = tau0*((p_array[ip]/ps)^(n))
      PRINT, ' '
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

      epsilon_plot[itau0, ip] = (U-D)/(sigma*(thetas^4.))
      tau_plot[itau0, ip] = tau
    ENDFOR

    ;find the value at LCL
    ind_p900 = yzclosest(p_array, P_lcl)
    epsilon_lcl_plot[itau0] = epsilon_plot[itau0, ind_p900]
  ENDFOR

  ;plot
  cgps_open, Path + 'Fig3_Takahashi2008.eps',/CMYK,$
    /nomatch,Font = 1,/Times, xsize =10/2.54, ysize = 10/2.54, fontsize =14

  pos = cglayout([1,1], OXMargin=[5,1], OYMargin=[4,1])

  mycharsize = 1.
  myxrange = [0., 1.]
  myyrange = [1000., 0.]

  myxtitle = cggreek('epsilon')
  myytitle = 'p (hPa)'

  mycolors = ['TG7','TG5','TG4','TG3','TG1']
  ;mycolors = ['RYB7','RYB5','RYB4','RYB3','RYB1']
  ;mycolors = ['YGB7','YGB5','YGB4','YGB3','YGB1']

  location_lg = [0.02, 50]
  text_lg = [cggreek('tau') + '0=100',cggreek('tau') + '0=25',cggreek('tau') + '0=6',cggreek('tau') + '0=2',cggreek('tau') + '0=0.5']
  linestyle_lg = MAKE_ARRAY(ntau0, value = 0)

  ;Fig 1
  cgplot,epsilon_plot[0, *], p_array,/noerase,/nodata, $
    xtitle = myxtitle, ytitle = myytitle,$
    xrange = myxrange,yrange = myyrange,$
    xtickinterval = 0.2, ytickinterval = 100,$
    charsize = mycharsize, POSITION = pos[*,0]

  FOR itau0 = 0, ntau0 - 1 DO BEGIN
    cgoplot,epsilon_plot[itau0, *], p_array, thick = 2., color = mycolors[itau0]
    cgoplot, epsilon_lcl_plot[itau0], p_lcl, psym = cgsymcat('opencircle'), color = mycolors[itau0], symsize = 0.7
  ENDFOR

  cglegend, title = text_lg, linestyle = linestyle_lg, color = mycolors, Location = location_lg, /data,$
    vspace = 0.8, charsize = 0.7*mycharsize, tcolors = mycolors

  cgps_close, /png, density = 1000

  PRINT, ' '

END

FUNCTION f_u_bl, x
  COMMON share

  RETURN,  ((x/tau0)^alphabl)*(EXP(-(x - tau)))
END

FUNCTION f_u_fa, x
  COMMON share

  RETURN,  (((tau_lcl/tau0)^alphabl)*((x/tau_lcl)^alphafa))*(EXP(-(x - tau)))
END

FUNCTION f_d_fa, x
  COMMON share
  RETURN,  (((tau_lcl/tau0)^alphabl)*((x/tau_lcl)^alphafa))*(EXP(-(tau - x)))
END

FUNCTION f_d_bl, x
  COMMON share
  RETURN,  ((x/tau0)^alphabl)*(EXP(-(tau - x)))
END

