      pro pd5

      oldn=!D.name & set_plot,'ps'
 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      ns = 201
      pi1=3.141592653589793 / 2.
      fr = fltarr(ns)
      ri = fltarr(ns)
      for i=0, ns-1 do begin
          r =  (1. - 1./3.) / (ns - 1) * i + 1./3. 
          ri(i) = 511. * r
          fr(i) = 1. / r + r + (2. * r - 1)^2 / r^2 -1.
      endfor 
      ;print, theta, theta2
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      ;rfilequv, infile = './Imu1.txt', m1, m2, m3, m4, m5, m6, m7, m8, ns
      rfilequv, infile = './HxmImu1.dat', hm1, hm2, hm3, hm4, hm5, hm6, hm7, hm8, hns
      rfilequv, infile = './Imu1.dat', m1, m2, m3, m4, m5, m6, m7, m8, ns
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      pi1=3.141592653589793 / 2.
      nvs = fltarr(ns)
      nu_2 = 600.
      nu_1 = 0.
      for i=0, ns-1 do begin
          nvs(i) = ( nu_2 - nu_1 ) / (ns-1) * i + nu_1
      endfor
      print, ns, ns, hns
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      ;readfile, infile = './Iquv60.txt', Iquv_I60, Iquv_q60, Iquv_u60, Iquv_v60, nf
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      ;readfile, infile = './Iquv80.txt', Iquv_I80, Iquv_q80, Iquv_u80, Iquv_v80, nf
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
 

      ratio=1./0.618
      l=16 & xxss=l*(ratio) & yyss=l
      ;xoff=(LL-xxss)/2.,yoff=(3*LL/2.-yyss)/2.,
      !p.font = 0
      device,filename='./Hxm_12.ps',xsize=xxss,ysize=yyss,bits_per_pixel=8,$
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

      xlen=0.6
      ylen=0.6
      xlen2=0.3
      ylen2=0.3
      xgap1=0.105
      xgap2=0.03  

      ;xlen2= ( 1. - xgap1 - xgap2 - xgap3 ) / 3.
      dyy = 0.35
      y0 = 5
      rat = 1.0
      theta0 = acos(0.8)
      n_1 =  2 

      m = 1
      ygap_up = 0.02
      ygap_d = 0.125
      ygap = 0.01
      ylen = (1. - ygap_up - ygap_d ) / m

      xlow = 150
      xup = 550
      ylow = .5
      yup =    1.2

      h_gauss = 0.04
      pdel = 0.02
      for i=0,m-1 do begin    
 
          xlen  = ( 1. - xgap1 - xgap2 ) 
          x1 = xgap1
          x2 = x1 + xlen  
          IF i EQ 0 THEN BEGIN      
              ylow = -1.6
              yup  = 0.2
          endif  
          IF i EQ 2 THEN BEGIN      
              ylow = -5
              yup  = -0
          endif 
          IF i EQ 1 THEN BEGIN      
              ylow = -0.
              yup  = 0.1
          endif 
          IF i EQ 3 THEN BEGIN      
              ylow = 0
              yup  = 0.1
          endif   

          posup=[x1, ygap_d + (m-i-1)*ylen + ygap * (m-i-1), x2, ygap_d + (m-i)*ylen + ygap * (m-i-1) ] 
          plot,[xlow,xup],[ylow,yup],pos=[posup],/noerase,/nodata,/device,$
             xrange=[xlow,xup],yrange=[ylow,yup],/ynozero,/normal,xstyle=4+1,ystyle=4+1
 
         ;axis,xaxis=0,xticks=6,xminor=4,xrange=[xlow,xup],xstyle=1,$;font=-1,$;,xtickname=replicate(' ',6),$;,
          ;charsize=1,xtitle=textoidl('\theta');,xthick=tickth,color=colors
         nt = 1000 * 5.99 / 6

          IF i EQ 0 THEN BEGIN
              axis,yaxis=0,ytitle=textoidl('g_1(E, \mu_2) (Photon/keV ster.)'),yticks=9,yminor=2, $
                yrange=[ylow,yup],ystyle=1, charsize=1.8 
              axis,yaxis=1,ytickname=replicate(' ',20),yticks=9, yminor=2
              axis,xaxis=1,xticks=4,xminor=2,xtickname=replicate(' ',15) 
              ;axis,xaxis=0,xticks=6,xminor=2,xtickname=replicate(' ',15) 
              axis,xaxis=0,xticks=4,xminor=2,xrange=[xlow,xup],xstyle=1,$
                      ;font=-1,$;,xtickname=replicate(' ',6),$;,
               charsize=1.8,xtitle=textoidl('h\nu (keV)');,xthick=tickth,color=colors

              ;nt = 1000 * 5.9 / 6
              ns = 0
              max_I =  max( m1 )
              np = where(m1 > 0.)
              ;np = 150
              s1 = 1.; alog10( m1(np) / max_I )
              print, np(0) ;max_I , min(m1), alog10( m1 / m1(np) )
              ;It = I0 + I1 + I2 + I3 + I4 + I5
              ;It =  I2 + I3 + I4 + I5

               tk = 1
               plot_curve_hxm, nvs, hm1, tk, 255, 6, 2 
               plot_curve_hxm, nvs, hm2, tk, 254, 6, 2 
               plot_curve_hxm, nvs, hm3, tk, 253, 6, 2  
               plot_curve_hxm, nvs, hm4, tk, 252, 6, 2  
               plot_curve_hxm, nvs, hm5, tk, 251, 6, 2  
               plot_curve_hxm, nvs, hm6, tk, 250, 6, 2  
               plot_curve_hxm, nvs, hm7, tk, red, 6, 2 
               plot_curve_hxm, nvs, hm8, tk, 248, 6, 8  
                
               tk = 5
               plot_curve, nvs, m1, hm1, tk, black, 1, 2 
               plot_curve, nvs, m2, hm2, tk, black, 1, 2  
               plot_curve, nvs, m3, hm3, tk, black, 1, 2  
               plot_curve, nvs, m4, hm4, tk, black, 1, 2  
               plot_curve, nvs, m5, hm5, tk, black, 1, 2  
               plot_curve, nvs, m6, hm6, tk, black, 1, 2 
               plot_curve, nvs, m7, hm7, tk, black, 1, 2  
               plot_curve, nvs, m8, hm8, tk, black, 1, 8
  
              x_s = 360.0
              yup = 0.11
              dy = 0.09
              csz = 1.5
              xyouts, x_s, yup - dy*0, textoidl('\mu_2=0.05') , charsize=csz
              xyouts, x_s, yup - dy*1, textoidl('\mu_2=0.10') , charsize=csz
              xyouts, x_s, yup - dy*2, textoidl('\mu_2=0.15') , charsize=csz
              xyouts, x_s, yup - dy*3, textoidl('\mu_2=0.25') , charsize=csz

              xyouts, x_s, yup - dy*4, textoidl('\mu_2=0.50') , charsize=csz
              xyouts, x_s, yup - dy*5, textoidl('\mu_2=0.75') , charsize=csz
              xyouts, x_s, yup - dy*6, textoidl('\mu_2=0.95') , charsize=csz
              xyouts, x_s, yup - dy*7, textoidl('\mu_2=1.00') , charsize=csz
              xp1 =  300
              xp2 = 355
              yp2 = yup + 0.01
              yp1 = yp2
              dys = dy
              tk = 2
              oplot, [xp1, xp2], [yp1, yp2] , thick=tk, color=255, linestyle=6;,psym=-4
              oplot, [xp1, xp2], [yp1-dys, yp2-dys] , thick=tk, color=254, linestyle=6;,psym=-4
              oplot, [xp1, xp2], [yp1-2 *dys, yp2-2 *dys] , thick=tk, color=253, linestyle=6;,psym=-4
              oplot, [xp1, xp2], [yp1-3 *dys, yp2-3 *dys] , thick=tk, color=252, linestyle=6;,psym=-4 
              oplot, [xp1, xp2], [yp1-4 *dys, yp2-4 *dys] , thick=tk, color=251, linestyle=6;,psym=-4
              oplot, [xp1, xp2], [yp1-5 *dys, yp2-5 *dys] , thick=tk, color=250, linestyle=6;,psym=-4
              oplot, [xp1, xp2], [yp1-6 *dys, yp2-6 *dys] , thick=tk, color=red, linestyle=6;,psym=-4
              oplot, [xp1, xp2], [yp1-7 *dys, yp2-7 *dys] , thick=tk, color=248, linestyle=6;,psym=-4  
          endif  
 
          IF i EQ 1 THEN BEGIN
              axis,yaxis=0,ytitle=textoidl('\delta'),yticks=(yup-ylow)/pdel,yminor=2,yrange=[ylow,yup],ystyle=1
              axis,yaxis=1,ytickname=replicate(' ',12),yticks=(yup-ylow)/pdel,yminor=2
              axis,xaxis=1,xticks=6,xminor=2,xtickname=replicate(' ',15) 
              axis,xaxis=0,xticks=6,xminor=2,xtickname=replicate(' ',15) 
              ;axis,xaxis=0,xticks=6,xminor=4,xrange=[xlow,xup],xstyle=1,$
                      ;font=-1,$;,xtickname=replicate(' ',6),$;,
              ;charsize=1,xtitle=textoidl('log(E/m_ec^2)');,xthick=tickth,color=colors

              ;nt = 1000 * 5.9 / 6 sqrt( Q11(ns:nt)^2+U11(ns:nt)^2 )
              ns = 1
              ;oplot, lg10E(ns:nt),  sqrt( Q11(ns:nt)^2+U11(ns:nt)^2 ) /  abs( I11(ns:nt) ) , $
              ;          thick=4, color=blue, linestyle=6;,psym=-4  

              ;KED2, lg10E, sqrt( Q11^2+U11^2 ), Y_array_kde, 0.001
;
              ;oplot, lg10E(ns:nt),  Y_array_kde(ns:nt) / abs(I11) , $
              ;         thick=4, color=red, linestyle=6;,psym=-4  
          endif 
 
          IF i EQ 2 THEN BEGIN
              axis,yaxis=0,ytitle=textoidl('log(vLv)'),yticks=(yup-ylow),yminor=2,yrange=[ylow,yup],ystyle=1
              axis,yaxis=1,ytickname=replicate(' ',12),yticks=(yup-ylow),yminor=2
              axis,xaxis=1,xticks=6,xminor=2,xtickname=replicate(' ',15) 
              axis,xaxis=0,xticks=6,xminor=2,xtickname=replicate(' ',15) 
              ;axis,xaxis=0,xticks=6,xminor=4,xrange=[xlow,xup],xstyle=1,$
                      ;font=-1,$;,xtickname=replicate(' ',6),$;,
              ;charsize=1,xtitle=textoidl('log(E/m_ec^2)');,xthick=tickth,color=colors

              ;nt = 1000 * 5 / 6
              ns = 1
              max_I2 = max( Ij6(ns:nt) )
              ;oplot, lg10E(ns:nt), alog10(Ij6(ns:nt) /  max_I2 ) , $
              ;          thick=2, color=blue, linestyle=6;,psym=-4     
          endif  
          IF i EQ 3 THEN BEGIN
              axis,yaxis=0,ytitle=textoidl('\delta'),yticks=(yup-ylow)/pdel,yminor=2,yrange=[ylow,yup],ystyle=1
              axis,yaxis=1,ytickname=replicate(' ',12),yticks= (yup-ylow)/pdel,yminor=2
              axis,xaxis=1,xticks=6,xminor=2,xtickname=replicate(' ',15) 
              axis,xaxis=0,xticks=6,xminor=2,xrange=[xlow,xup],xstyle=1,$
                      ;font=-1,$;,xtickname=replicate(' ',6),$;,
              charsize=1,xtitle=textoidl('log(E/m_ec^2)');,xthick=tickth,color=colors

              ;nt = 1000 * 5 / 6
              ;ns = 1 
              ;oplot, lg10E(ns:nt), sqrt( Q50(ns:nt)^2+U50(ns:nt)^2 ) /  I50(ns:nt)  , $
             ;           thick=4, color=blue, linestyle=6;,psym=-4   
          endif 
      endfor 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      device,/close 

      ;free_lun,lunAo10
      set_plot,oldn

      end

;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PRO rfilequv, infile = infile, m1, m2, m3, m4, m5, m6, m7, m8, n
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
;Demonstrate reading from a file 
;infile = !Bowman + ’data/table.txt’ ;Input file name
n = FILE_LINES(infile)-1
      m1=fltarr(n)
      m2=fltarr(n)
      m3=fltarr(n)
      m4=fltarr(n)
      m5=fltarr(n)
      m6=fltarr(n)
      m7=fltarr(n)
      m8=fltarr(n) 
;x = FLTARR(n)  
 print, ' data number= ', n  
 
 OPENR, iunit, infile, /GET_LUN ;Open input file
 Point_lun, iunit, 0
  
      FOR j = 0, n-1 DO BEGIN
         ;ReadF,lunAo10, x1, x2, x3, x4;, FORMAT = "(4F22.17)"
         READF, iunit, x0, x1, x2, x3, x4, x5, x6, x7;, FORMAT = "(4F22.17)"
         m1[j] = x0
         m2[j] = x1
         m3[j] = x2
         m4[j] = x3
         m5[j] = x4
         m6[j] = x5
         m7[j] = x6
         m8[j] = x7
         ;print, x3, x4, x5, x6
      ENDFOR
FREE_LUN, iunit ;Close input file
;FOR i = 0, n−1 DO PRINT, x[i], logx[i], $ ;Print values to terminal
;    FORMAT = "(2F12.5)"
END
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PRO plot_curve, xarr, yarr, yarrh, tk, colc, ls, num
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nh = where(yarrh > 0.) 
    np = where(yarr > 0.) 
    nc = (max(np) + min(np))/2
    ty = yarrh(nc) / yarrh(nh(0))  
    if num eq 8 then begin  
        yarr(max(np)) = 1.e-30 
    endif else begin
        yarr(max(np)+1) = 1.e-30 
    endelse
    oplot, xarr, alog10( yarr / yarr(nc)*ty ), thick=5, color=colc, linestyle=ls;,psym=-4 
    ;print, '**************************************************************************'
    ;print, ty, nc
 
END
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PRO plot_curve_hxm, xarr, yarr, tk, colc, ls, num
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    tk = 5
    np = where(yarr > 0.)    
    if num eq 8 then begin   
        yarr(max(np)+1) = 1.e-30 
    endif else begin
        yarr(max(np)+1) = 1.e-30 
    endelse
    oplot, xarr, alog10( yarr / yarr(np(0)) ), thick=2, color=colc, linestyle=ls;,psym=-4  
    ;print, alog10( yarr / yarr(np(0)) )

END
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

