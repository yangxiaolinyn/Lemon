      pro pd6

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
      ;rfiles, infile = './I_Theta=05_1.txt',I0, I1, I2, I3, ns
      ;rfiles, infile = './I_Theta=05_2.txt',I10, I11, I12, I13, ns
      ;rfiles, infile = './I_Theta=10_1.txt',I20, I21, I22, I23, ns
      ;rfiles, infile = './I_Theta=10_2.txt',I30, I31, I32, I33, ns
      rfiles, infile = './dataTe=2.555E-01tau=1.000E-01.dat',I0, I1, I2, I3, ns
      rfiles, infile = './dataTe=5.110E-01tau=1.000E-01.dat',I20, I21, I22, I23, ns
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
      nu_2 = 5.
      nu_1 = -1.
      for i=0, ns-1 do begin
          lg10E(i) = ( nu_2 - nu_1 ) / (ns-1) * i + nu_1
      endfor
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      ;readfile, infile = './Iquv60.txt', Iquv_I60, Iquv_q60, Iquv_u60, Iquv_v60, nf
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      ;readfile, infile = './Iquv80.txt', Iquv_I80, Iquv_q80, Iquv_u80, Iquv_v80, nf
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
 

      ratio=2.4/1.0
      l=16 & xxss=l*(ratio) & yyss=l
      ;xoff=(LL-xxss)/2.,yoff=(3*LL/2.-yyss)/2.,
      !p.font = 0
      device,filename='./fig_20.ps',xsize=xxss,ysize=yyss,bits_per_pixel=8,$
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

      xlow = -1
      xup = 3
      ylow = .5
      yup =    1.2

      h_gauss = 0.04
      pdel = 1.0

      xlen = (1. - xgap1 - xgap2 - xgap * (m-1.) ) / 1.

      xgs = fltarr(3)
      ygs = fltarr(3)
      xgs(0) = 0.06
      xgs(1) = 0.01
      xgs(2) = 0.08
      ygs(0) = 0.02
      ygs(1) = 0.1
      ygs(2) = 0.005
      n_x = 2
      n_y = 1
      ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      panel_set, m_x = n_x, m_y = n_y, xgaps = xgs, ygaps = ygs, $
                   xl_arr, xr_arr, yu_arr, yd_arr
      ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      for j=0, n_x-1 do begin    
        for i=0, n_y-1 do begin    
 
          ;xlen  = ( 1. - xgap1 - xgap2 ) 
          ;x1 = xgap1
          ;x2 = x1 + xlen  
          IF i EQ 0 THEN BEGIN      
              ylow = -5
              yup  = 1
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
              axis,yaxis=0,ytitle=textoidl('log(\nuF_\nu)   (Photon Flux in Arbitrary Units)'),$
                        yticks=(yup-ylow),yminor=5, $
                yrange=[ylow,yup],ystyle=1, charsize=1.3
              axis,yaxis=1,ytickname=replicate(' ',12),yticks=(yup-ylow),yminor=5
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
              oplot, lg10E, alog10(I0 / max(I1) ), $;- 2.03, $
                        thick=4, color=252, linestyle=6;,psym=-4  
              oplot, lg10E, alog10(I1 / max(I0) )- 0.1 , $
                        thick=4, color=244, linestyle=2;,psym=-4  
              oplot, lg10E, alog10(I2 / max(I0) )- 0.8 , $
                        thick=4, color=254, linestyle=3;,psym=-4  
              oplot, lg10E, alog10(I3 / max(I0) )- 0.1 , $
                        thick=4, color=255, linestyle=4;,psym=-4    

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

          IF (i EQ 0) and (j EQ 1) THEN BEGIN
              axis,yaxis=0,ytitle=textoidl('log(\nuF_\nu)   (Photon Flux in Arbitrary Units)'),$
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
              oplot, lg10E, alog10(I20(*) / max(I21) ) - 0.08 , $;
                        thick=4, color=252, linestyle=6;,psym=-4  
              oplot, lg10E, alog10(I21(*) / max(I21) ) - 0.35 , $
                        thick=4, color=244, linestyle=2;,psym=-4  
              oplot, lg10E, alog10(I22(*) / max(I21) ) - 1.05 , $
                        thick=4, color=254, linestyle=3;,psym=-4  
              oplot, lg10E, alog10(I23(*) / max(I21) ) - 0.35 , $
                        thick=4, color=255, linestyle=4;,psym=-4    

              ;max3 = max( I3(i_mu, *) / max_I )
              ;KED2,  lg10E, I3(i_mu, *) , Y_array_kde, 0.05
              ;oplot, lg10E, alog10( Y_array_kde / max(Y_array_kde) *max3 ),thick=4,$
              ;              color=red,linestyle=6;,psym=-4  
              dy = 0.3
              xyouts, 2.08, 0.4, textoidl('\mu = 1.0'), charsize=1.3
              oplot, [1.2, 2.], [0.45, 0.45], thick=4, color=252, linestyle=6;,psym=-4  
              xyouts, 2.08, 0.4 - dy, textoidl('\mu = 0.5'), charsize=1.3
              oplot, [1.2, 2.], [0.45 - dy*1, 0.45 - dy*1], thick=4, color=244, linestyle=2;,psym=-4  
              xyouts, 2.08, 0.4 - dy*2., textoidl('\mu = 0.1'), charsize=1.3
              oplot, [1.2, 2.], [0.45 - dy*2, 0.45 - dy*2], thick=4, color=254, linestyle=3;,psym=-4  
              xyouts, 2.08, 0.4 - dy*3., textoidl('\mu = - 0.5'), charsize=1.3
              oplot, [1.2, 2.], [0.45 - dy*3., 0.45 - dy*3.], thick=4, color=255, linestyle=4;,psym=-4  
 
              xyouts, -0.5, 0.4, textoidl('\Theta = 1.0'), charsize=1.3 
              xyouts, -0.5, 0.4 - dy, textoidl('\tau = 0.1'), charsize=1.3
 
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
