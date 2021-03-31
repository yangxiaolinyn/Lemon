      pro pdi

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
      ;rdfile1, infile = './bcs1.txt', mu = 3, times = 4, IQUV_1, ns 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      rdfile, infile = './dataHotElectron_Te=100.0000_fig24_1.dat', mu = 1, times = 4, IQUV_1, ns 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      rdfile, infile = './dataHotElectron_Te=100.0000_fig24_2.dat', mu = 1, times = 4, IQUV_2, ns 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

      lg10E1 = fltarr(ns)
      nu_2 = 7.
      nu_1 = 1.
      for i=0, ns-1 do begin
          lg10E1(i) = ( nu_2 - nu_1 ) / (ns-1) * i + nu_1
      endfor
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

      rdfil, infile = './BCS_HotElectron_Te=100.0000_fig24_1.dat', I0, Q0, U0, V0, ns
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
      nu_2 = 7.
      nu_1 = 1.
      for i=0, ns-1 do begin
          lg10E(i) = ( nu_2 - nu_1 ) / (ns-1) * i + nu_1
      endfor
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      ;readfile, infile = './Iquv60.txt', Iquv_I60, Iquv_q60, Iquv_u60, Iquv_v60, nf
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      ;readfile, infile = './Iquv80.txt', Iquv_I80, Iquv_q80, Iquv_u80, Iquv_v80, nf
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
 

      ratio=1.5/1.0
      l=16 & xxss=l*(ratio) & yyss=l
      ;xoff=(LL-xxss)/2.,yoff=(3*LL/2.-yyss)/2.,
      !p.font = 0
      device,filename='./figI_1.ps',xsize=xxss,ysize=yyss,bits_per_pixel=8,$
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

      xlow = 1
      xup = 7
      ylow = .5
      yup =    1.2

      h_gauss = 0.04
      pdel = 1.0

      xlen = (1. - xgap1 - xgap2 - xgap * (m-1.) ) / 1.

      xgs = fltarr(3)
      ygs = fltarr(3)
      xgs(0) = 0.1 
      xgs(1) = 0.01
      xgs(2) = 0.08
      ygs(0) = 0.02
      ygs(1) = 0.138
      ygs(2) = 0.005
      n_x = 1
      n_y = 1
      csize = 2.0
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
              ylow = -4
              yup  = 1
          endif   
          IF (i EQ 0) and (j EQ 1) THEN BEGIN    
              ylow = 0
              yup  =0.3
          endif   
 
          ;posup=[x1, ygap_d + (m-i-1)*ylen + ygap * (m-i-1), x2, ygap_d + (m-i)*ylen + ygap * (m-i-1) ] 
          ;x1 = xgap1 + i * xlen + xgap * i
          ;x2 = x1 + xlen
          ;print, x1, x2, xlen, ylow
          ;posup=[x1, ygap_d , x2, ygap_d + ylen ] 

          posup=[xl_arr[j], yd_arr[n_y-1-i] , xr_arr[j], yu_arr[n_y-1-i] ] 
          plot,[xlow,xup],[ylow,yup],pos=[posup],/noerase,/nodata,/device,$
             xrange=[xlow,xup],yrange=[ylow,yup],/ynozero,/normal,xstyle=4+1,ystyle=4+1
   
          IF (i EQ 0) and (j EQ 0) THEN BEGIN 
              ;ylow = -5
              ;yup  = 1
              axis,yaxis=0,ytitle=textoidl('Log( Intensity I )'),$
                        yticks=(yup-ylow),yminor=5, $
                yrange=[ylow,yup],ystyle=1, charsize=csize
              axis,yaxis=1,ytickname=replicate(' ',12),yticks=(yup-ylow),yminor=5, charsize=csize
              axis,xaxis=1,xticks=6,xminor=2,xtickname=replicate(' ',15), charsize=csize
              ;axis,xaxis=0,xticks=6,xminor=2,xtickname=replicate(' ',15) 
              axis,xaxis=0,xticks=6,xminor=2,xrange=[xlow,xup],xstyle=1,$ 
              charsize=csize,xtitle=textoidl('Log(\epsilon\prime/\varepsilon)')
              ;,xthick=tickth,color=colorsfont=-1,$;,xtickname=replicate(' ',6),$;,

              ;nt = 1000 * 5.9 / 6
              ns = 0
              i_mu = 0
              max_I =  max( 10.^2 ) 
              max_I2 =  max( 2. ) 
              oplot, lg10E, alog10( I0 / max(I0) )-0.06, thick=4, color=black, linestyle=6;,psym=-4  

              ;oplot, lg10E1, alog10( IQUV_1[0, 0, *]/ max(IQUV_1[0, 0, *]) ), $
               ;         thick=4, color=red, linestyle=6;,psym=-4    


              histgram, x_arr = lg10E1, y_arr = alog10( IQUV_2[0, 0, *]/ max(IQUV_2[0, 0, *]) ), x, y
              oplot, x, y, thick=4, color=blue, linestyle=6;,psym=-4   

              ;xyouts, 1.5, 0.4,  textoidl('S_{in} = (1, 1 / 2, 1 / 2, 1 / 2^{1/2})'), charsize=1.3


             n_arr = 2
             charac_arr = sindgen(n_arr);" S(1,1,0,0), S(1,0,1,0), S(1,0,0,1) "
             charac_arr(0) = textoidl('numerical')
             charac_arr(1) = textoidl('BCS70 semi-analytic') 
             ;charac_arr(3) = textoidl('S_{in}=(1,1/2,1/2,1/\sqrt{2})')
             dy = 0.3
             cha_size = 1.6
             x0 = 1.2
             y0 =  0.65 - 0.3
             length = 0.6
             colors = indgen(n_arr)
             colors(0) = blue
             colors(1) = black
             line_thick = 4
          curves, dy = dy, charac_arr = charac_arr, n_arr = n_arr, cha_size = cha_size, x0 = x0, $
            y0 = y0, length = length, colors = colors, line_thick = line_thick


            xll = 1.1
            yll = 0.1- 0.3
            x_length = 2.4
            y_length = 0.8
            rectangle, xll = xll, yll = yll, x_length = x_length, y_length = y_length  

              ;print, alog10( IQUV_1[0, 0, *] ), max(IQUV_1[0, 0, *])
              ;max3 = max( I3(i_mu, *) / max_I )
              ;KED2,  lg10E, I3(i_mu, *) , Y_array_kde, 0.05
              ;oplot, lg10E, alog10( Y_array_kde / max(Y_array_kde) *max3 ),thick=4,$
              ;              color=red,linestyle=6;,psym=-4  


             ; max2 = max( I2(i_mu, *) / max_I )
              ;KED2,  lg10E, I2(i_mu, *) , Y_array_kde, 0.05
              ;oplot, lg10E, alog10( Y_array_kde / max(Y_array_kde) *max2 ),thick=4,$
              ;              color=red,linestyle=6;,psym=-4  
 
              ;max1 = max( I1(i_mu, *) / max_I )
             ; KED2,  lg10E, I1(i_mu, *) , Y_array_kde, 0.05
              ;oplot, lg10E, alog10( Y_array_kde / max(Y_array_kde) *max1 ),thick=4,$
              ;              color=blue,linestyle=6;,psym=-4  

           ;   KED2, lg10E,   I6 , Y_array_kde, 0.009
            ;  max2 = max(Y_array_kde)
             ; oplot, lg10E(ns:nt), alog10(Y_array_kde(ns:nt)/max2 )-0.,thick=2,$
             ;               color=blue,linestyle=6;,psym=-4 
            IF 1 eq 0 THEN BEGIN 
              dy = 0.3
              xyouts, 2.08, 0.4, textoidl('\mu = 1.0'), charsize=1.3
              oplot, [1.2, 2.], [0.45, 0.45], thick=4, color=252, linestyle=6;,psym=-4  
              xyouts, 2.08, 0.4 - dy, textoidl('\mu = 0.5'), charsize=1.3
              oplot, [1.2, 2.], [0.45 - dy*1, 0.45 - dy*1], thick=4, color=244, linestyle=2;,psym=-4  
              xyouts, 2.08, 0.4 - dy*2., textoidl('\mu = 0.1'), charsize=1.3
              oplot, [1.2, 2.], [0.45 - dy*2, 0.45 - dy*2], thick=4, color=254, linestyle=3;,psym=-4  
              xyouts, 2.08, 0.4 - dy*3., textoidl('\mu = - 0.5'), charsize=1.3
              oplot, [1.2, 2.], [0.45 - dy*3., 0.45 - dy*3.], thick=4, color=255, linestyle=4;,psym=-4 

              xyouts, -0.5, 0.4, textoidl('\Theta = 0.5'), charsize=1.3 
              xyouts, -0.5, 0.4 - dy, textoidl('\tau = 0.1'), charsize=1.3 
             endif
   
          endif 

          IF (i EQ 0) and (j EQ 1) THEN BEGIN 
              axis,yaxis=0,ytitle=textoidl('Polarization Degrees (\delta=Q/I)'),$
                        yticks=(yup-ylow),yminor=5, $
                yrange=[ylow,yup],ystyle=1, charsize=1.3
              ;axis,yaxis=0,ytitle=textoidl('log(vLv erg/s)'),yticks=(yup-ylow),yminor=5, $
              ;  yrange=[ylow,yup],ystyle=1, charsize=1.3
              axis,yaxis=1,ytickname=replicate(' ',12),yticks=(yup-ylow),yminor=5
              ;axis,yaxis=0,ytickname=replicate(' ',12),yticks=(yup-ylow),yminor=5
              axis,xaxis=1,xticks=4,xminor=5,xtickname=replicate(' ',15) 
              ;axis,xaxis=0,xticks=6,xminor=2,xtickname=replicate(' ',15) 
              axis,xaxis=0,xticks=4,xminor=5,xrange=[xlow,xup],xstyle=1,$ 
              charsize=1.3,xtitle=textoidl('log(h\nu/kT_{bb})')
              ;,xthick=tickth,color=colorsfont=-1,$;,xtickname=replicate(' ',6),$;,

              ;nt = 1000 * 5.9 / 6
              ns = 0
              i_mu = 0
              max_I =  max( 10.^2 ) 
              max_I2 =  max( 2. ) 
            IF 1 eq 0 THEN BEGIN 
              oplot, lg10E, sqrt( IQUV_1[0, 1, *]^2 + IQUV_1[0, 2, *]^2 ) / IQUV_1[0, 0, *], $
                        thick=4, color=252, linestyle=6;,psym=-4  
              oplot, lg10E, sqrt( IQUV_1(1, 1, *)^2 + IQUV_1(1, 2, *)^2 ) / IQUV_1(1, 0, *), $
                        thick=4, color=244, linestyle=6;,psym=-4  
              oplot, lg10E, sqrt( IQUV_1[2, 1, *]^2 + IQUV_1[2, 2, *]^2 ) / IQUV_1[2, 0, *] , $
                        thick=4, color=254, linestyle=6;,psym=-4 
            endif
 
            IF 1 eq 1 THEN BEGIN 
              KED2, lg10E, sqrt( IQUV_1[2, 1, *]^2 + IQUV_1[2, 2, *]^2 ) / IQUV_1[2, 0, *], $
                           Y_array_kde, 0.05
               max2 = max(Y_array_kde)
              oplot, lg10E,  Y_array_kde / max2 * 0.23, thick=3,$
                             color=blue, linestyle=1;,psym=-4 


              KED2, lg10E, sqrt( IQUV_1[1, 1, *]^2 + IQUV_1[1, 2, *]^2 ) / IQUV_1[1, 0, *], $
                           Y_array_kde, 0.05
               max2 = max(Y_array_kde)
              oplot, lg10E,  Y_array_kde / max2 * 0.18, thick=3,$
                             color=blue, linestyle=2;,psym=-4 


              KED2, lg10E, sqrt( IQUV_1[0, 1, *]^2 + IQUV_1[0, 2, *]^2 ) / IQUV_1[0, 0, *], $
                           Y_array_kde, 0.03
               max2 = max(Y_array_kde)
              oplot, lg10E,  Y_array_kde / max2 * 0.1, thick=3,$
                             color=252, linestyle=3;,psym=-4 
               ;print, Y_array_kde
           endif   

            IF 1 eq 0 THEN BEGIN 
 
              oplot, lg10E, sqrt( IQUV_1[0, 1, *]^2 + IQUV_1[0, 2, *]^2 ) , $
                        thick=4, color=252, linestyle=6;,psym=-4  
              oplot, lg10E, sqrt( IQUV_1(1, 1, *)^2 + IQUV_1(1, 2, *)^2 ) , $
                        thick=4, color=244, linestyle=2;,psym=-4  
              oplot, lg10E, sqrt( IQUV_1[2, 1, *]^2 + IQUV_1[2, 2, *]^2 ) , $
                        thick=4, color=254, linestyle=3;,psym=-4     
           endif      

              ;print, lg10E
              ;max3 = max( I3(i_mu, *) / max_I )
              ;KED2,  lg10E, I3(i_mu, *) , Y_array_kde, 0.05
              ;oplot, lg10E, alog10( Y_array_kde / max(Y_array_kde) *max3 ),thick=4,$
              ;              color=red,linestyle=6;,psym=-4  
              dy = 0.3
              ;xyouts, 2.08, 0.4, textoidl('\mu = 1.0'), charsize=1.3
              ;oplot, [1.2, 2.], [0.45, 0.45], thick=4, color=252, linestyle=6;,psym=-4  
              ;xyouts, 2.08, 0.4 - dy, textoidl('\mu = 0.5'), charsize=1.3
              ;oplot, [1.2, 2.], [0.45 - dy*1, 0.45 - dy*1], thick=4, color=244, linestyle=2;,psym=-4  
              ;xyouts, 2.08, 0.4 - dy*2., textoidl('\mu = 0.1'), charsize=1.3
              ;oplot, [1.2, 2.], [0.45 - dy*2, 0.45 - dy*2], thick=4, color=254, linestyle=3;,psym=-4  
              ;xyouts, 2.08, 0.4 - dy*3., textoidl('\mu = - 0.5'), charsize=1.3
              ;oplot, [1.2, 2.], [0.45 - dy*3., 0.45 - dy*3.], thick=4, color=255, linestyle=4;,psym=-4  
 
              ;xyouts, -0.5, 0.4, textoidl('\Theta = 1.0'), charsize=1.3 
              ;xyouts, -0.5, 0.4 - dy, textoidl('\tau = 0.1'), charsize=1.3
 
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

;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PRO rfiles, infile = infile, I0, I1, I2, I3, n
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PRO rfile3, infile = infile, I0, I1, I2, n
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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


;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PRO rdfile, infile = infile, mu = mu, times = times, I1_6, n
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
;Demonstrate reading from a file 
;infile = !Bowman + ’data/table.txt’ ;Input file name
n = FILE_LINES(infile) / mu
      I1_6=fltarr(mu, times, n)
      xf=fltarr( times ) 
;x = FLTARR(n)  
 print, ' data number= ', n  
 
 OPENR, iunit, infile, /GET_LUN ;Open input file
 Point_lun, iunit, 0
  
      FOR i = 0, mu-1 DO BEGIN
        FOR j = 0, n-1 DO BEGIN
         ;ReadF,lunAo10, x1, x2, x3, x4;, FORMAT = "(4F22.17)"
         if times EQ 7 then begin
             ;READF, iunit, x0[0], x0[1], x0[2], x0[3], x0[4], x0[5], x0[6];, FORMAT = "(4F22.17)"
             READF, iunit, x0, x1, x2, x3, x4, x5, x6;, FORMAT = "(4F22.17)"
             I1_6[i, 0, j] = x0
             I1_6[i, 1, j] = x1
             I1_6[i, 2, j] = x2
             I1_6[i, 3, j] = x3
             I1_6[i, 4, j] = x4
             I1_6[i, 5, j] = x5
             I1_6[i, 6, j] = x6
             ;print, 'tss', x6, times
         endif
         if times EQ 4 then begin
             ;READF, iunit, x0[0], x0[1], x0[2], x0[3], x0[4];, FORMAT = "(4F22.17)"
             READF, iunit, x0, x1, x2, x3;, FORMAT = "(4F22.17)"
             ;print, 'tss', x0[6], times
             I1_6[i, 0, j] = x0
             I1_6[i, 1, j] = x1
             I1_6[i, 2, j] = x2
             I1_6[i, 3, j] = x3 
         endif 
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




;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PRO KED2, X_array, Y_array, Y_array_kde, h_gauss
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


      ;h_gauss = 0.02
      ;nums = size(Y_array)
      ;ns = ;nums[3]
      ns = n_elements(Y_array)
      print, ns  
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      Y_array_kde = fltarr(ns)
      Y_array_kde[*] = 0.
      for i=0, ns-1 do begin
          for j=0, ns-1 do begin 
              Dx_ij = abs( X_array(i) - X_array(j) ) 
              Y_array_kde(i) = Y_array_kde(i) + Y_array(j) * $ 
                   exp( -  Dx_ij^2 / 2./h_gauss ) 
          endfor
      endfor  

END



;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PRO rdfil, infile = infile, I0, Q0, U0, V0, n
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
;Demonstrate reading from a file 
;infile = !Bowman + ’data/table.txt’ ;Input file name
n = FILE_LINES(infile) 
      I0=fltarr(n) 
      Q0=fltarr(n) 
      U0=fltarr(n) 
      V0=fltarr(n) 
;x = FLTARR(n)  
 print, ' data number= ', n  
 
 OPENR, iunit, infile, /GET_LUN ;Open input file
 Point_lun, iunit, 0
   
      FOR i = 0, n-1 DO BEGIN
         ;ReadF,lunAo10, x1, x2, x3, x4;, FORMAT = "(4F22.17)"
         READF, iunit, x0, x1, x2, x3;, FORMAT = "(4F22.17)"
         I0[i] = x0 
         Q0[i] = x1 
         U0[i] = x2 
         V0[i] = x3 
         ;print, x3, x4, x5, x6
      ENDFOR 
FREE_LUN, iunit ;Close input file
;FOR i = 0, n−1 DO PRINT, x[i], logx[i], $ ;Print values to terminal
;    FORMAT = "(2F12.5)"
END





;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PRO histgram, x_arr = x_arr, y_arr = y_arr, x, y
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    hist = y_arr
    nhist = n_elements(y_arr)
    x_max = max(x_arr)
    x_min = min(x_arr)
    binsize = ( x_max - x_min ) / nhist
    bins = lindgen( nhist ) * binsize + x_min
    x = fltarr(2 * nhist)
    x[2 * lindgen(nhist)] = bins
    x[2 * lindgen(nhist) + 1] = bins

    y = fltarr(2 * nhist)
    y[2 * lindgen(nhist)] = hist
    y[2 * lindgen(nhist) + 1] = hist
    y = shift(y, 1)

END


;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PRO curves, dy = dy, charac_arr = charac_arr, n_arr = n_arr, cha_size = cha_size, x0 = x0, $
              y0 = y0, length = length, colors = colors, line_thick = line_thick
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    FOR i = 0, n_arr-1 DO BEGIN  
        oplot, [x0, x0 + length], [y0-dy*i, y0-dy*i], thick=line_thick, color=colors(i), linestyle=6;,psym=-4  
        xyouts, x0 + length*1.01, y0-dy*i-0.05, charac_arr(i), charsize=cha_size
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

