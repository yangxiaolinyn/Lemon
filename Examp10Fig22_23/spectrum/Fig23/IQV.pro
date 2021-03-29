      pro IQV

      oldn=!D.name & set_plot,'ps'
 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      ns = 201
      pi1=3.141592653589793 / 2.
      phi = fltarr(ns)
      for i=0, ns-1 do begin
          phi(i) =  2.*!Pi / 200. * i ;(ns - i-1) * 90.0 / (ns - 1)
      endfor 
      ;print, theta, theta2
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      ;rfiles, infile = './IQUV_1.txt',I0, I1, I2, I3, ns
      ;rdfile, infile = './IQUV352.txt', mu = 3, times = 4, IQUV_1, ns
      rdfile, infile = './dataTe=3.520E-01tau=5.000E-02.dat', n_iquv = 4, $
                       mu = 2, num_scatter = 7, IQUV_1, n
      ;rdfile, infile = './IQUV_1.txt', mu = 3, times = 4, IQUV_1, ns
      ;rfiles, infile = './I_Theta=05_2.txt',I10, I11, I12, I13, ns
      ;rfiles, infile = './I_Theta=10_1.txt',I20, I21, I22, I23, ns
      ;rfiles, infile = './I_Theta=10_2.txt',I30, I31, I32, I33, ns
      ;rfile3, infile = './I_Theta=05_1.txt',I0, I1, I2, ns
      ;print, I0(*)
      ;I0(0, 200:*) = 0.
      ;I0(1, 200:*) = 0.
      ;I0(2, 200:*) = 0.
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      ;readfile, infile = './I12345_mu50.txt', Ij0, Ij1, Ij2, Ij3, Ij4, Ij5, Ij6, ns
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ;rfilequv, infile = './IQUVpmu11.txt', I11, Q11, U11, V11, d11, ns
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ;rfilequv, infile = './IQUV_E.txt', I0, I1, I2, I3, I4, I5, I6, ns
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      pi1=3.141592653589793 / 2.
      lg10E = fltarr(ns)
      nu_2 = 1.
      nu_1 = -5.
      for i=0, ns-1 do begin
          lg10E(i) = ( nu_2 - nu_1 ) / (ns-1) * i + nu_1
      endfor
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      ;readfile, infile = './Iquv60.txt', Iquv_I60, Iquv_q60, Iquv_u60, Iquv_v60, nf
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      ;readfile, infile = './Iquv80.txt', Iquv_I80, Iquv_q80, Iquv_u80, Iquv_v80, nf
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
 

      ratio=1.6/1.0
      l=16 & xxss=l*(ratio) & yyss=l
      ;xoff=(LL-xxss)/2.,yoff=(3*LL/2.-yyss)/2.,
      !p.font = 0
      device,filename='./fig352_3.ps',xsize=xxss,ysize=yyss,bits_per_pixel=8,$
      /color,xoff=(2-xxss)/2.0,yoff=(2-yyss)/2.,$
      set_font='Times-Roman';, /tt_font

      loadct,30
      RRR=bytscl(findgen(256))
      GGG=bytscl(findgen(256))
      BBB=bytscl(findgen(256))
      RRR[243:255]=[0,255,0  ,0  ,0  ,255,255,200,0  ,0  ,200,200,0  ]
      GGG[243:255]=[0,0  ,255,0  ,255,0  ,255,0  ,200,0  ,200,0  ,200]
      BBB[243:255]=[0,0  ,0  ,255,255,255,0  ,0  ,0  ,200,0  ,200,200]
      ;245=green,243=black,244=red,251=yellow,246=blue,253=cyan,254=magenta,255=white
      TVLCT,RRR,GGG,BBB
      blue = 246
      red = 244
      green = 245
      black = 243
      yellow = 247
      cyan = 253
      magenta = 254
 
      xgap1=0.1
      xgap2=0.01
      xgap =0.005

      ;xlen2= ( 1. - xgap1 - xgap2 - xgap3 ) / 3.
      dyy = 0.35
      y0 = 5
      rat = 1.0
      theta0 = acos(0.8)
      n_1 =  2 

      m = 1
      ygap_up = 0.01
      ygap_d = 0.1
      ygap = 0.01
      ylen = (1. - ygap_up - ygap_d ) / 1

      xlow = -5
      xup = 1
      ylow = .5
      yup =    1.2

      h_gauss = 0.04
      pdel = 1.0

      xlen = (1. - xgap1 - xgap2 - xgap * (m-1.) ) / 1.

      xgs = fltarr(3)
      ygs = fltarr(3)
      xgs(0) = 0.08
      xgs(1) = 0.01
      xgs(2) = 0.08
      ygs(0) = 0.02
      ygs(1) = 0.1
      ygs(2) = 0.005
      n_x = 2
      n_y = 2
      ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      panel_set, m_x = n_x, m_y = n_y, xgaps = xgs, ygaps = ygs, $
                   xl_arr, xr_arr, yu_arr, yd_arr
      ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      for j=0, n_x-1 do begin    
        for i=0, n_y-1 do begin    
 
          ;xlen  = ( 1. - xgap1 - xgap2 ) 
          ;x1 = xgap1
          ;x2 = x1 + xlen  
          IF (i EQ 0) and (j EQ 0) THEN BEGIN    
              ylow = -7
              yup  = 1
          endif   
          IF (i EQ 0) and (j EQ 1) THEN BEGIN    
              ylow = -7
              yup  = 1
          endif   
          IF (i EQ 1) and (j EQ 0) THEN BEGIN    
              n_ytick = 4
              ylow = -0.05
              yup  =0.15
          endif    
          IF (i EQ 1) and (j EQ 1) THEN BEGIN    
              ylow = -0.05
              yup  =0.15
              n_ytick = 4
          endif    
   
          posup=[xl_arr[j], yd_arr[n_y-1-i] , xr_arr[j], yu_arr[n_y-1-i] ] 
          plot,[xlow,xup],[ylow,yup],pos=[posup],/noerase,/nodata,/device,$
             xrange=[xlow,xup],yrange=[ylow,yup],/ynozero,/normal,xstyle=4+1,ystyle=4+1
   
          IF (i EQ 0) and (j EQ 0) THEN BEGIN 
              ;ylow = -5
              ;yup  = 1
              axis,yaxis=0,ytitle=textoidl('Log(\nuF_\nu)'),$
                        yticks=(yup-ylow),yminor=5, $
                yrange=[ylow,yup],ystyle=1, charsize=1.1
              axis,yaxis=1,ytickname=replicate(' ',12),yticks=(yup-ylow),yminor=5
              axis,xaxis=1,xticks=6,xminor=5,xtickname=replicate(' ',15) 
              axis,xaxis=0,xticks=6,xminor=5,xtickname=replicate(' ',15) 
              ;axis,xaxis=0,xticks=6,xminor=5,xrange=[xlow,xup],xstyle=1,$ 
              ;charsize=1.3,xtitle=textoidl('log10( h\nu / m_ec^2 )')
              ;,xthick=tickth,color=colorsfont=-1,$;,xtickname=replicate(' ',6),$;,
  
              ns = 0
              i_mu = 0
              max_I =  max( 10.^2 ) 
              max_I2 =  max( 2. ) 
              i_thi = 3
              oplot, lg10E, alog10( IQUV_1[0, i_mu, 6, *] / max( IQUV_1[0, i_mu, 6, *] ) ), $
                        thick=i_thi, color=blue, linestyle=6;,psym=-4   

              oplot, lg10E, alog10( IQUV_1[0, i_mu, 0, *] / max( IQUV_1[0, i_mu, 6, *] ) ), $
                        thick=i_thi, color=red, linestyle=6;,psym=-4 

              oplot, lg10E, alog10( IQUV_1[0, i_mu, 1, *] / max( IQUV_1[0, i_mu, 6, *] ) ), $
                        thick=i_thi, color=green, linestyle=1;,psym=-4  

              oplot, lg10E, alog10( IQUV_1[0, i_mu, 2, *] / max( IQUV_1[0, i_mu, 6, *] ) ), $
                        thick=i_thi, color=cyan, linestyle=2;,psym=-4 

              oplot, lg10E, alog10( IQUV_1[0, i_mu, 3, *] / max( IQUV_1[0, i_mu, 6, *] ) ), $
                        thick=i_thi, color=248, linestyle=3;,psym=-4 

              oplot, lg10E, alog10( IQUV_1[0, i_mu, 4, *] / max( IQUV_1[0, i_mu, 6, *] ) ), $
                        thick=i_thi, color=247, linestyle=4;,psym=-4 

              ;KED2,  lg10E, IQUV_1[0, 0, *],  I_kde, 0.00011
              ;oplot, lg10E, I_kde, thick=4, color=blue, linestyle=6;,psym=-4  
   
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             n_arr = 2
             charac_arr = sindgen(n_arr);" S(1,1,0,0), S(1,0,1,0), S(1,0,0,1) "
             charac_arr(0) = textoidl('\mu_{obs} = 0.11   \tau = 0.05')
             charac_arr(1) = textoidl('T_{bb} = 10 eV    T_e = 352 keV')
             ;charac_arr(0) = textoidl(' \mu = 0.1')
             ;charac_arr(3) = textoidl(' S_{in} = (1, 0, 0, 1)')
             ;charac_arr(3) = textoidl('S_{in}=(1,1/2,1/2,1/\sqrt{2})')
             dy = 1.
             y_bias = dy * 0.2
             cha_size = 1.2
             x0 = -3
             y0 = 0.4
             length = 0.
             colors = indgen(n_arr)
             ;colors(2) = white
             colors(1) = green
             colors(0) = green
             ;colors(3) = green
             line_thick = 4
          ;curves, dy = dy, charac_arr = charac_arr, n_arr = n_arr, cha_size = cha_size, x0 = x0, $
          ;  y0 = y0, length = length, colors = colors, line_thick = line_thick, y_bias = y_bias

          xyout, dy = dy, charac_arr = charac_arr, n_arr = n_arr, cha_size = cha_size, x0 = x0, $
              y0 = y0, colors = colors, line_thick = line_thick, y_bias = y_bias 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
 
   
          endif 

          IF (i EQ 0) and (j EQ 1) THEN BEGIN  
              axis,yaxis=0,ytitle=textoidl('Log(\nuF_\nu) '),$
                        yticks=(yup-ylow),yminor=5, $
                yrange=[ylow,yup],ystyle=1, charsize=1.1
              axis,yaxis=1,ytickname=replicate(' ',12),yticks=(yup-ylow),yminor=5
              axis,xaxis=1,xticks=6,xminor=5,xtickname=replicate(' ',15) 
              axis,xaxis=0,xticks=6,xminor=5,xtickname=replicate(' ',15) 
              ;axis,xaxis=0,xticks=6,xminor=5,xrange=[xlow,xup],xstyle=1,$ 
              ;charsize=1.3,xtitle=textoidl('log10( h\nu / m_ec^2 )')
 
 
            IF 1 eq 1 THEN BEGIN   
 
              ns = 0
              i_mu = 1
              max_I =  max( 10.^2 ) 
              max_I2 =  max( 2. ) 
              i_thi = 3
             
              oplot, lg10E, alog10( IQUV_1[0, i_mu, 6, *] / max( IQUV_1[0, i_mu, 6, *] ) ), $
                        thick=i_thi, color=blue, linestyle=6;,psym=-4  

              oplot, lg10E, alog10( IQUV_1[0, i_mu, 0, *] / max( IQUV_1[0, i_mu, 6, *] ) ), $
                        thick=i_thi, color=red, linestyle=6;,psym=-4 

              oplot, lg10E, alog10( IQUV_1[0, i_mu, 1, *] / max( IQUV_1[0, i_mu, 6, *] ) ), $
                        thick=i_thi, color=green, linestyle=1;,psym=-4  

              oplot, lg10E, alog10( IQUV_1[0, i_mu, 2, *] / max( IQUV_1[0, i_mu, 6, *] ) ), $
                        thick=i_thi, color=cyan, linestyle=2;,psym=-4 

              oplot, lg10E, alog10( IQUV_1[0, i_mu, 3, *] / max( IQUV_1[0, i_mu, 6, *] ) ), $
                        thick=i_thi, color=248, linestyle=3;,psym=-4 

              oplot, lg10E, alog10( IQUV_1[0, i_mu, 4, *] / max( IQUV_1[0, i_mu, 6, *] ) ), $
                        thick=i_thi, color=247, linestyle=4;,psym=-4 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             n_arr = 1
             charac_arr = sindgen(n_arr);" S(1,1,0,0), S(1,0,1,0), S(1,0,0,1) "
             charac_arr(0) = textoidl('\mu_{obs} = 0.50 ')
             ;charac_arr(1) = textoidl('T_{bb} = 10 eV    T_e = 352 keV')
             ;charac_arr(0) = textoidl(' \mu = 0.1')
             ;charac_arr(3) = textoidl(' S_{in} = (1, 0, 0, 1)')
             ;charac_arr(3) = textoidl('S_{in}=(1,1/2,1/2,1/\sqrt{2})')
             dy = 1. 
             y_bias = dy * 0.2
             cha_size = 1.2
             x0 = -2
             y0 = -0.5
             length = 0.
             colors = indgen(n_arr)
             ;colors(2) = white
             ;colors(1) = white
             colors(0) = red
             ;colors(3) = green
             line_thick = 4
          ;curves, dy = dy, charac_arr = charac_arr, n_arr = n_arr, cha_size = cha_size, x0 = x0, $
          ;  y0 = y0, length = length, colors = colors, line_thick = line_thick, y_bias = y_bias

          xyout, dy = dy, charac_arr = charac_arr, n_arr = n_arr, cha_size = cha_size, x0 = x0, $
              y0 = y0, colors = colors, line_thick = line_thick, y_bias = y_bias
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            endif

            IF 1 eq 0 THEN BEGIN  
              KED2, lg10E, sqrt( IQUV_1[2, 1, *]^2 + IQUV_1[2, 2, *]^2 ) / IQUV_1[2, 0, *], $
                           Y_array_kde, 0.001
               max3 = 1.; max(Y_array_kde)
              oplot, lg10E,  Y_array_kde / max3, thick=3,$
                             color=green, linestyle=6;,psym=-4 
               ;print, sqrt( IQUV_1[0, 1, *]^2 + IQUV_1[0, 2, *]^2 )


              KED2, lg10E, sqrt( IQUV_1[1, 1, *]^2 + IQUV_1[1, 2, *]^2 ) / IQUV_1[1, 0, *], $
                           Y_array_kde, 0.001
               max2 = max(Y_array_kde)
              oplot, lg10E,  Y_array_kde / max3 , thick=3,$
                             color=red, linestyle=6;,psym=-4 

              KED2, lg10E, sqrt( IQUV_1[0, 1, *]^2 + IQUV_1[0, 2, *]^2 ) / IQUV_1[0, 0, *], $
                           Y_array_kde, 0.001
               max2 = max(Y_array_kde)
              oplot, lg10E,  Y_array_kde / max3, thick=3,$
                             color=blue, linestyle=6;,psym=-4   
           endif 
            IF 1 eq 0 THEN BEGIN  
              KED2, lg10E, IQUV_1[2, 1, *], $ ;/ IQUV_1[0, 0, *], $
                           Y_array_kde, 0.001
               max3 = max(Y_array_kde)
              oplot, lg10E,  Y_array_kde / max3 * 1, thick=3,$
                             color=green, linestyle=6;,psym=-4 
               ;print, sqrt( IQUV_1[0, 1, *]^2 + IQUV_1[0, 2, *]^2 )


              KED2, lg10E, IQUV_1[1, 1, *] / IQUV_1[1, 0, *], $
                           Y_array_kde, 0.05
               max2 = 1;max(Y_array_kde)
              oplot, lg10E,  Y_array_kde / max3 * 1, thick=3,$
                             color=red, linestyle=6;,psym=-4 

              KED2, lg10E, IQUV_1[0, 1, *] / IQUV_1[2, 0, *], $
                           Y_array_kde, 0.05
               max2 = max(Y_array_kde)
              oplot, lg10E,  Y_array_kde / max3 * 1, thick=3,$
                             color=blue, linestyle=6;,psym=-4   
           endif 
     
          endif 

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          IF (i EQ 1) and (j EQ 0) THEN BEGIN  
              axis,yaxis=0,ytitle=textoidl(' Polarization Degree \delta '),$
                        yticks=n_ytick,yminor=5, $
                yrange=[ylow,yup],ystyle=1, charsize=1.1
              axis,yaxis=1,ytickname=replicate(' ',12),yticks=n_ytick,yminor=5
              axis,xaxis=1,xticks=6,xminor=5,xtickname=replicate(' ',15) 
              ;axis,xaxis=0,xticks=6,xminor=2,xtickname=replicate(' ',15) 
              axis,xaxis=0,xticks=6,xminor=5,xrange=[xlow,xup],xstyle=1,$ 
              charsize=1.3,xtitle=textoidl('Log( h\nu / m_ec^2 )')

 
              i_thick = 3
              oplot, lg10E, IQUV_1[1, 0, 6, *]  /  IQUV_1[0, 0, 6, *] , $
                        thick=i_thick, color=blue, linestyle=6;,psym=-4  
 
              ;i_thick = 3
              ;oplot, lg10E, IQUV_1[2, 2, 6, *]  /  IQUV_1[0, 2, 6, *] , $
              ;          thick=i_thick, color=blue, linestyle=6;,psym=-4  
 

              ;KED2, lg10E,  IQUV_1[0, 3, *], Y_array_kde, 0.001
              ;oplot, lg10E,  Y_array_kde / max(abs(Y_array_kde))*max(abs( IQUV_1[0, 3, *] )), thick=4,$
              ;               color=blue, linestyle=6;,psym=-4 
 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             n_arr = 3
             charac_arr = sindgen(n_arr);" S(1,1,0,0), S(1,0,1,0), S(1,0,0,1) "
             charac_arr(2) = textoidl(' \mu = 0.9')
             charac_arr(1) = textoidl(' \mu = 0.5')
             charac_arr(0) = textoidl(' \mu = 0.1')
             ;charac_arr(3) = textoidl(' S_{in} = (1, 0, 0, 1)')
             ;charac_arr(3) = textoidl('S_{in}=(1,1/2,1/2,1/\sqrt{2})')
             dy = 0.8
             y_bias = dy * 0.2
             cha_size = 1.3
             x0 = -2
             y0 = 7
             length = 1.2
             colors = indgen(n_arr)
             colors(2) = blue
             colors(1) = red
             colors(0) = green
             ;colors(3) = green
             line_thick = 4
          curves, dy = dy, charac_arr = charac_arr, n_arr = n_arr, cha_size = cha_size, x0 = x0, $
            y0 = y0, length = length, colors = colors, line_thick = line_thick, y_bias = y_bias

            ;xyouts, 1.4, y0-5*dy,  textoidl('S_{in} = (1, 1 / 2, 1 / 2, 1 / 2^{1/2})'), charsize=1.3

            x_length = length * 4.
            y_length = dy * (n_arr + 1 )
            xll = x0 - dy
            yll = y0 + dy - y_length
            ;rectangle, xll = xll, yll = yll, x_length = x_length, y_length = y_length  
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          endif   
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          IF (i EQ 1) and (j EQ 1) THEN BEGIN  
              axis,yaxis=0,ytitle=textoidl(' Polarization Degree \delta '),$
                        yticks=n_ytick,yminor=5, $
                yrange=[ylow,yup],ystyle=1, charsize=1.1
              axis,yaxis=1,ytickname=replicate(' ',12),yticks=n_ytick,yminor=5
              axis,xaxis=1,xticks=6,xminor=5,xtickname=replicate(' ',15) 
              ;axis,xaxis=0,xticks=6,xminor=2,xtickname=replicate(' ',15) 
              axis,xaxis=0,xticks=6,xminor=5,xrange=[xlow,xup],xstyle=1,$ 
              charsize=1.3,xtitle=textoidl('Log( h\nu / m_ec^2 )')
 
              i_thick = 3
              oplot, lg10E, IQUV_1[1, 1, 6, *] /  IQUV_1[0, 1, 6, *] , $
                        thick=i_thick, color=blue, linestyle=6;,psym=-4  

              oplot, lg10E, IQUV_1[1, 1, 6, *] /  max( abs( IQUV_1[1, 1, 6, *] ) ) * 0.1 , $
                        thick=i_thick, color=red, linestyle=6;,psym=-4  

          endif   
        endfor 
      endfor 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      device,/close 

      ;free_lun,lunAo10
      set_plot,oldn

      end

PRO rfiles, infile = infile, I0, I1, I2, I3, n
;Demonstrate reading from a file 
;infile = !Bowman + ’data/table.txt’ ;Input file name
n = FILE_LINES(infile) 
      I0=fltarr(n)
      I1=fltarr(n)
      I2=fltarr(n)
      I3=fltarr(n) 
;x = FLTARR(n)  
 print, ' data number= ', n  
 
 OPENR, iunit, infile, /GET_LUN ;Open input file
 Point_lun, iunit, 0
   
      FOR i = 0, n-1 DO BEGIN
         ;ReadF,lunAo10, x1, x2, x3, x4;, FORMAT = "(4F22.17)"
         READF, iunit, x0, x1, x2, x3;, FORMAT = "(4F22.17)"
         I0[i] = x0
         I1[i] = x1
         I2[i] = x2
         I3[i] = x3 
         ;print, x3, x4, x5, x6
      ENDFOR 
FREE_LUN, iunit ;Close input file
;FOR i = 0, n−1 DO PRINT, x[i], logx[i], $ ;Print values to terminal
;    FORMAT = "(2F12.5)"
END

 
PRO rfile3, infile = infile, I0, I1, I2, n
;Demonstrate reading from a file 
;infile = !Bowman + ’data/table.txt’ ;Input file name
n = FILE_LINES(infile) 
      I0=fltarr(n)
      I1=fltarr(n)
      I2=fltarr(n)  
;x = FLTARR(n)  
 print, ' data number= ', n  
 
 OPENR, iunit, infile, /GET_LUN ;Open input file
 Point_lun, iunit, 0
   
      FOR i = 0, n-1 DO BEGIN
         ;ReadF,lunAo10, x1, x2, x3, x4;, FORMAT = "(4F22.17)"
         READF, iunit, x0, x1, x2;, FORMAT = "(4F22.17)"
         I0[i] = x0
         I1[i] = x1
         I2[i] = x2 
         ;print, x3, x4, x5, x6
      ENDFOR 
FREE_LUN, iunit ;Close input file
;FOR i = 0, n−1 DO PRINT, x[i], logx[i], $ ;Print values to terminal
;    FORMAT = "(2F12.5)"
END  


PRO rdfile, infile = infile, n_iquv = n_iquv, mu = mu, num_scatter = num_scatter, I1_6, n
;Demonstrate reading from a file
;infile = !Bowman + ’data/table.txt’ ;Input file name
n = FILE_LINES(infile) / num_scatter
      I1_6=fltarr(n_iquv, mu, num_scatter, n)
     ;xf=fltarr( times ) 
;x = FLTARR(n)  
 print, ' data number= ', n  
 
 OPENR, iunit, infile, /GET_LUN ;Open input file
 Point_lun, iunit, 0
  
      FOR i = 0, num_scatter-1 DO BEGIN
        FOR j = 0, n-1 DO BEGIN
         ;ReadF,lunAo10, x1, x2, x3, x4;, FORMAT = "(4F22.17)"
         ;if times EQ 7 then begin
             ;READF, iunit, x0[0], x0[1], x0[2], x0[3], x0[4], x0[5], x0[6];, FORMAT = "(4F22.17)"
             READF, iunit, x0, x1, x2, x3, x4, x5, x6, x7;, x8, x9, x10, x11;, FORMAT = "(4F22.17)"
             I1_6[0, 0, i, j] = x0
             I1_6[1, 0, i, j] = x1
             I1_6[2, 0, i, j] = x2
             I1_6[3, 0, i, j] = x3

             I1_6[0, 1, i, j] = x4
             I1_6[1, 1, i, j] = x5
             I1_6[2, 1, i, j] = x6
             I1_6[3, 1, i, j] = x7

            ; I1_6[0, 2, i, j] = x8
            ; I1_6[1, 2, i, j] = x9
            ; I1_6[2, 2, i, j] = x10
            ; I1_6[3, 2, i, j] = x11
             ;print, 'tss', x6, times
         ;endif
         ;if times EQ 4 then begin
             ;READF, iunit, x0[0], x0[1], x0[2], x0[3], x0[4];, FORMAT = "(4F22.17)"
             ;READF, iunit, x0, x1, x2, x3;, FORMAT = "(4F22.17)"
             ;print, 'tss', x0[6], times
             ;I1_6[i, 0, j] = x0
             ;I1_6[i, 1, j] = x1
             ;I1_6[i, 2, j] = x2
             ;I1_6[i, 3, j] = x3 
         ;endif 
        ENDFOR
      ENDFOR
FREE_LUN, iunit ;Close input file
;FOR i = 0, n−1 DO PRINT, x[i], logx[i], $ ;Print values to terminal
;    FORMAT = "(2F12.5)"
END




  
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PRO panel_set, m_x = m_x, m_y = m_y, xgaps = xgaps, ygaps = ygaps, $
                   xl_arr, xr_arr, yu_arr, yd_arr  
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      xlen = ( 1. - xgaps[0] - xgaps[1] - xgaps[2] * (m_x - 1.) ) / m_x
      ylen = ( 1. - ygaps[0] - ygaps[1] - ygaps[2] * (m_y - 1.) ) / m_y
      xl_arr = fltarr(m_x)*0.
      xr_arr = fltarr(m_x)*0.

      yu_arr = fltarr(m_y)*0.
      yd_arr = fltarr(m_y)*0. 

      FOR i = 0, m_x-1 DO BEGIN 
          xl_arr[i] = xgaps[0] + ( xlen + xgaps[2] ) * i
          xr_arr[i] = xgaps[0] + ( xlen + xgaps[2] ) * i + xlen
      endfor

      FOR i = 0, m_y-1 DO BEGIN 
          yd_arr[i] = ygaps[1] + ( ylen + ygaps[2] ) * i
          yu_arr[i] = ygaps[1] + ( ylen + ygaps[2] ) * i + ylen
      endfor
   
END




PRO KED2, X_array, Y_array, Y_array_kde, h_gauss


      ;h_gauss = 0.02
      ;nums = size(Y_array)
      ;ns = ;nums[3]
      ns = n_elements(Y_array)
      print, ns  

      ;Imu11_0_kde = fltarr(ns)
      ;for i=0, ns-1 do begin
      ;    for j=0, ns-1 do begin
      ;        ;x_j = lg10E(j)IQmu11
      ;        Imu11_0_kde(i) = Imu11_0_kde(i) + Imu11_0(j) * exp( - (lg10E(i) - lg10E(j))^2 / 2./h_gauss )
      ;    endfor
      ;endfor 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      Y_array_kde = fltarr(ns)
      Y_array_kde[*] = 0.
      for i=0, ns-1 do begin
          for j=0, ns-1 do begin 
            ;Dx_ij = abs( X_array(i) - X_array(j) ) 
            Y_array_kde(i) = Y_array_kde(i) + Y_array(j) * exp( - ( X_array(i) - X_array(j) )^2 / 2./h_gauss ) 
          endfor
      endfor   
END


PRO KED3, X_array, Y_array, h_gauss, n_end = n_end, Y_array_kde
 
      ns = n_elements(Y_array)
      print, ns   
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      Y_array_kde = fltarr(ns)
      Y_array_kde[*] = 0.
      for i=0, n_end do begin
          for j=0, n_end do begin 
            ;Dx_ij = abs( X_array(i) - X_array(j) ) 
            Y_array_kde(i) = Y_array_kde(i) + Y_array(j) * exp( - ( X_array(i) - X_array(j) )^2 / 2./h_gauss ) 
          endfor
      endfor   
END


;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PRO curves, dy = dy, charac_arr = charac_arr, n_arr = n_arr, cha_size = cha_size, x0 = x0, $
              y0 = y0, length = length, colors = colors, line_thick = line_thick, y_bias = y_bias
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    FOR i = 0, n_arr-1 DO BEGIN  
        oplot, [x0, x0 + length], [y0-dy*i, y0-dy*i], thick=line_thick, color=colors(i), linestyle=6;,psym=-4  
        xyouts, x0 + length*1.01, y0-dy*i-y_bias, charac_arr(i), charsize=cha_size
    ENDFOR 
END



;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PRO rectangle, xll = xll, yll = yll, x_length = x_length, y_length = y_length
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
            oplot, [xll, xll], [yll, yll+y_length], thick=4, color=black, linestyle=6;,psym=-4  
        oplot, [xll+x_length, xll+x_length], [yll, yll+y_length], thick=4, color=black, linestyle=6;,psym=-4  
            oplot, [xll, xll+x_length], [yll, yll], thick=4, color=black, linestyle=6;,psym=-4  
        oplot, [xll, xll+x_length], [yll+y_length, yll+y_length], thick=4, color=black, linestyle=6;,psym=-4  
   
END



;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PRO xyout, dy = dy, charac_arr = charac_arr, n_arr = n_arr, cha_size = cha_size, x0 = x0, $
              y0 = y0, colors = colors, line_thick = line_thick, y_bias = y_bias
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    FOR i = 0, n_arr-1 DO BEGIN  
      ;oplot, [x0, x0 + length], [y0-dy*i, y0-dy*i], thick=line_thick, color=colors(i), linestyle=6;,psym=-4  
      xyouts, x0, y0-dy*i-y_bias, charac_arr(i), charsize=cha_size
    ENDFOR 
END


