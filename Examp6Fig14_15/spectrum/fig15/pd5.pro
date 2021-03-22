      pro pd5

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
      readfile, infile = './ChandraIQUV_phi_mu0=-0.5000.dat', m=4, Ichd, Qchd, Uchd, Vchd, nf
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      ;readfile, infile = './Iquv30.txt', Iquv_I30, Iquv_q30, Iquv_u30, Iquv_v30, nf
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      ;readfile, infile = './Iquv60.txt', Iquv_I60, Iquv_q60, Iquv_u60, Iquv_v60, nf
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      ;readfile, infile = './Iquv80.txt', Iquv_I80, Iquv_q80, Iquv_u80, Iquv_v80, nf
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      readfile, infile = './IQUVphi_mu0=-0.5000.dat', m=4, Imc, Qmc, Umc, Vmc, nf 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      ;readfile, infile = './IQUV30_mu0=08_nozeroQUV.txt', IQUV_2I30, IQUV_2Q30, IQUV_2U30, IQUV_2V30, nf 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      ;readfile, infile = './IQUV60_mu0=08_nozeroQUV.txt', IQUV_2I60, IQUV_2Q60, IQUV_2U60, IQUV_2V60, nf 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      ;readfile, infile = './IQUV80_mu0=08_nozeroQUV.txt', IQUV_2I80, IQUV_2Q80, IQUV_2U80, IQUV_2V80, nf 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      Iquv_phi = fltarr(nf)
      for i=0, nf-1 do begin
          Iquv_phi( i ) = 2.*!Pi / (nf-1) * i
      endfor 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

      ratio=1./1.6
      l=16 & xxss=l*(ratio) & yyss=l
      ;xoff=(LL-xxss)/2.,yoff=(3*LL/2.-yyss)/2.,
      !p.font = 0
      device,filename='./fig14_1.ps',xsize=xxss,ysize=yyss,bits_per_pixel=8,$
      /color,xoff=(2-xxss)/2.0,yoff=(2-yyss)/2.,$
      set_font='Times-Roman';, /tt_font

      loadct,30
      RRR=bytscl(findgen(256))
      GGG=bytscl(findgen(256))
      BBB=bytscl(findgen(256))
      RRR[243:255]=[0,255,0  ,0  ,0  ,255,255,200,0  ,0  ,200,200,0]
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
      xgap1=0.18
      xgap2=0.02  

      ;xlen2= ( 1. - xgap1 - xgap2 - xgap3 ) / 3.
      dyy = 0.35
      y0 = 5
      rat = 1.0
      theta0 = acos(0.8)
      n_1 =  2 

      m = 4
      ygap_up = 0.01
      ygap_d = 0.07
      ygap = 0.01
      ylen = (1. - ygap_up - ygap_d - 3. * ygap) / m

      xlow = 0
      xup = 2*!Pi
      ylow = .5
      yup =    1.2

      h_gauss = 0.04
      
      for i=0,3 do begin    
 
          xlen  = ( 1. - xgap1 - xgap2 ) 
          x1 = xgap1
          x2 = x1 + xlen  
          IF i EQ 0 THEN BEGIN      
              ylow = 0.38
              yup  = 0.68
          endif  
          IF i EQ 1 THEN BEGIN      
              ylow = -.1 
              yup =  .22
          endif  
          IF i EQ 2 THEN BEGIN      
              ylow = -.35
              yup =   .35
          endif  
          IF i EQ 3 THEN BEGIN      
              ylow = -0.8
              yup =   0.8
          endif  

          posup=[x1, ygap_d + (m-i-1)*ylen + ygap * (m-i-1), x2, ygap_d + (m-i)*ylen + ygap * (m-i-1) ] 
          plot,[xlow,xup],[ylow,yup],pos=[posup],/noerase,/nodata,/device,$
             xrange=[xlow,xup],yrange=[ylow,yup],/ynozero,/normal,xstyle=4+1,ystyle=4+1
 
         ;axis,xaxis=0,xticks=6,xminor=4,xrange=[xlow,xup],xstyle=1,$;font=-1,$;,xtickname=replicate(' ',6),$;,
          ;charsize=1,xtitle=textoidl('\theta');,xthick=tickth,color=colors


          ;print, max(abs(Iquv_2q10))/max(abs(Iquv_2I10)), $
          ;     max(abs(IQUV_2Q30))/ max(abs(Iquv_2I30)), max(abs(Iquv_q10))
          IF i EQ 0 THEN BEGIN
              axis,yaxis=0,ytitle=textoidl('I'),yticks=5,yminor=2,yrange=[ylow,yup],ystyle=1
              axis,yaxis=1,ytickname=replicate(' ',12),yticks=5,yminor=2
              axis,xaxis=1,xticks=4,xminor=4,xtickname=replicate(' ',15)
              axis,xaxis=0,xticks=4,xminor=4,xtickname=replicate(' ',15)
 
              IQUV_PLOT, phi, Iquv_phi, Ichd[*, 0], Ichd[*, 1], Ichd[*, 2], Ichd[*, 3], $
              Imc[*, 0], Imc[*, 1], Imc[*, 2], Imc[*, 3], blue, red, magenta, green, black 


              ;max_I1 =  max( abs(Iquv_I10) ) * 1.01  
              ;oplot, phi,      Iquv_I10 , thick=2, color=c5, linestyle=6;,psym=-4
  
              ;max_2I1 = max(abs(IQUV_2I10))
              ;oplot, Iquv_phi, IQUV_2I10 / max_2I1 * max_I1, thick=3, color=c1, linestyle=6;,psym=-4
          endif  

          IF i EQ 1 THEN BEGIN
              axis,yaxis=0,ytitle=textoidl('Q'),yticks=4,yminor=2,yrange=[ylow,yup],ystyle=1
              axis,yaxis=1,ytickname=replicate(' ',12),yticks=4,yminor=2
              axis,xaxis=1,xticks=5,xminor=4,xtickname=replicate(' ',15)
              axis,xaxis=0,xticks=5,xminor=4,xtickname=replicate(' ',15)

              IQUV_PLOT, phi, Iquv_phi, Qchd[*, 0], Qchd[*, 1], Qchd[*, 2], Qchd[*, 3], $
              Qmc[*, 0], Qmc[*, 1], Qmc[*, 2], Qmc[*, 3], blue, red, magenta, green, black 
 

              ;max_I1 =  max( abs(Iquv_q10) ) * 1.01  
              ;oplot, phi,      Iquv_q10 , thick=2, color=c5, linestyle=6;,psym=-4
  
              ;max_2I1 = max(abs(IQUV_2Q10))
              ;oplot, Iquv_phi, IQUV_2Q10 / max_2I1 * max_I1 , thick=3, color=c1, linestyle=6;,psym=-4
 
              ;xyouts, !Pi / 2.2, yup * 0.8, textoidl('\mu=0.05')  
              ;xyouts, !Pi / 2.2, yup * 0.6, textoidl('\mu=0.3')  
              ;xyouts, !Pi / 2.2, yup * 0.4, textoidl('\mu=0.6')  
              ;xyouts, !Pi / 2.2, yup * 0.2, textoidl('\mu=0.85') 
              xp1 = !Pi / 9.5 
              xp2 = !Pi / 2.5 
              yp2 = yup * 0.8
              yp1 = yp2
              dys = yup * 0.2
              ;oplot, [xp1, xp2], [yp1, yp2] , thick=3, color=blue, linestyle=1;,psym=-4
              ;oplot, [xp1, xp2], [yp1-dys, yp2-dys] , thick=3, color=red, linestyle=2;,psym=-4
              ;oplot, [xp1, xp2], [yp1-2 *dys, yp2-2 *dys] , thick=3, color=magenta, linestyle=3;,psym=-4
              ;oplot, [xp1, xp2], [yp1-3 *dys, yp2-3 *dys] , thick=3, color=green, linestyle=4;,psym=-4  
          endif  
          IF i EQ 2 THEN BEGIN
              axis,yaxis=0,ytitle=textoidl('U'),yticks=7,yminor=2,yrange=[ylow,yup],ystyle=1
              axis,yaxis=1,ytickname=replicate(' ',12),yticks=7,yminor=2
              axis,xaxis=1,xticks=5,xminor=4,xtickname=replicate(' ',15)
              axis,xaxis=0,xticks=5,xminor=4,xtickname=replicate(' ',15)

              IQUV_PLOT, phi, Iquv_phi, Uchd[*, 0], Uchd[*, 1], Uchd[*, 2], Uchd[*, 3], $
              Umc[*, 0], Umc[*, 1], Umc[*, 2], Umc[*, 3], blue, red, magenta, green, black 
  
              ;max_I1 =  max( abs(Iquv_u10) ) * 1.01  
              ;oplot, phi,      Iquv_u10 , thick=2, color=black, linestyle=6;,psym=-4
  
              ;max_2I1 = max(abs(IQUV_2U10))
              ;oplot, Iquv_phi, IQUV_2U10 / max_2I1 * max_I1, thick=3, color=red, linestyle=1;,psym=-4
          endif  
          IF i EQ 3 THEN BEGIN
              axis,yaxis=0,ytitle=textoidl('V'),yticks=4,yminor=2,yrange=[ylow,yup],ystyle=1
              axis,yaxis=1,ytickname=replicate(' ',12),yticks=4,yminor=2
              axis,xaxis=1,xticks=4,xminor=4,xtickname=replicate(' ',15)
              labels = ['0', textoidl('\pi/2'), textoidl('\pi'), textoidl('3\pi/4'), textoidl('2\pi  ')]
              axis,xaxis=0,xticks=4,xminor=4,xrange=[xlow,xup],xstyle=1,$
                        xtickname=labels,$; xtickname=replicate( textoidl('\pi'),5),$;,
                      charsize=1,xtitle=textoidl('\varphi');,xthick=tickth,color=colors

              IQUV_PLOT, phi, Iquv_phi, Vchd[*, 0], Vchd[*, 1], Vchd[*, 2], Vchd[*, 3], $
              Vmc[*, 0], Vmc[*, 1], Vmc[*, 2], Vmc[*, 3], blue, red, magenta, green, black  

              dxx = !pi / 1.95
              dyy = 0.02
              xyouts, !Pi / 2.2 + dxx, yup * 0.8 - dyy, textoidl('\mu=0.05')  
              xyouts, !Pi / 2.2 + dxx, yup * 0.6 - dyy, textoidl('\mu=0.30')  
              xyouts, !Pi / 2.2 + dxx, yup * 0.4 - dyy, textoidl('\mu=0.60')  
              xyouts, !Pi / 2.2 + dxx, yup * 0.2 - dyy, textoidl('\mu=0.85') 
              xp1 = !Pi / 9.5  + dxx
              xp2 = !Pi / 2.5  + dxx
              yp2 = yup * 0.8
              yp1 = yp2
              dys = yup * 0.2
              oplot, [xp1, xp2], [yp1, yp2] , thick=3, color=blue, linestyle=1;,psym=-4
              oplot, [xp1, xp2], [yp1-dys, yp2-dys] , thick=3, color=red, linestyle=1;,psym=-4
              oplot, [xp1, xp2], [yp1-2 *dys, yp2-2 *dys] , thick=3, color=magenta, linestyle=1;,psym=-4
              oplot, [xp1, xp2], [yp1-3 *dys, yp2-3 *dys] , thick=3, color=green, linestyle=1;,psym=-4 
          endif   
      endfor 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      device,/close 

      ;free_lun,lunAo10
      set_plot,oldn

      end


PRO READFILE, infile = infile, m=m, I_r, I_l, I_u, I_v, n
;Demonstrate reading from a file 
;infile = !Bowman + ’data/table.txt’ ;Input file name
n = FILE_LINES(infile) / m
      I_l=fltarr(n, m)
      I_r=fltarr(n, m)
      I_u=fltarr(n, m)
      I_v=fltarr(n, m)
;x = FLTARR(n)  
 print, ' data number= ', n  
 
 OPENR, iunit, infile, /GET_LUN ;Open input file
 Point_lun, iunit, 0
  
  FOR j = 0, m-1 DO BEGIN
      FOR i = 0, n-1 DO BEGIN
         ;ReadF,lunAo10, x1, x2, x3, x4;, FORMAT = "(4F22.17)"
         READF, iunit, xr, xl, xu, xv;, FORMAT = "(4F22.17)"
         I_r[i, j] = xr
         I_l[i, j] = xl
         I_u[i, j] = xu
         I_v[i, j] = xv
      ENDFOR
  ENDFOR
FREE_LUN, iunit ;Close input file
;FOR i = 0, n−1 DO PRINT, x[i], logx[i], $ ;Print values to terminal
;    FORMAT = "(2F12.5)"
END





PRO IQUV_PLOT, phi, Iquv_phi, Iquv_I10, Iquv_I30, Iquv_I60, Iquv_I80, $
              IQUV_2I10, IQUV_2I30, IQUV_2I60, IQUV_2I80, c1, c2, c3, c4, c5

;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          max_I1 =  max( abs(Iquv_I10) ) * 1.0
          oplot, phi,      Iquv_I10 , thick=1, color=c5, linestyle=6;,psym=-4
 
          ;KED, Iquv_phi, IQUV_2I10, Y_array_kde, h_gauss
          max_2I1 = max(abs(IQUV_2I10))
          oplot, Iquv_phi, IQUV_2I10 / max_2I1 * max_I1, thick=3, color=c1, linestyle=1;,psym=-4 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          max_I1 =  max(abs( Iquv_I30 )) * 1.0
          oplot, phi,      Iquv_I30 , thick=1, color=c5, linestyle=6;,psym=-4  

          ;KED, Iquv_phi, IQUV_2I30, Y_array_kde, h_gauss 
          max_2I1 = max(abs(IQUV_2I30)) 
          oplot, Iquv_phi, IQUV_2I30 / max_2I1 * max_I1, thick=3, color=c2, linestyle=1;,psym=-4 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
          max_I1 =  max( abs(Iquv_I60) ) * 1.0
          oplot, phi,      Iquv_I60 , thick=1, color=c5, linestyle=6;,psym=-4   

          ;KED, Iquv_phi, IQUV_2I60, Y_array_kde, h_gauss 
          max_2I1 = max(abs(IQUV_2I60)) 
          oplot, Iquv_phi, IQUV_2I60 / max_2I1 * max_I1, thick=3, color=c3, linestyle=1;,psym=-4 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
          max_I1 =  max(abs( Iquv_I80) ) * 1.0
          oplot, phi,      Iquv_I80 , thick=1, color=c5, linestyle=6;,psym=-4   

          ;KED, Iquv_phi, IQUV_2I80, Y_array_kde, h_gauss 
          max_2I1 = max(abs(IQUV_2I80)) 
          oplot, Iquv_phi, IQUV_2I80 / max_2I1 * max_I1, thick=3, color=c4, linestyle=1; ,psym=-4 
  
;FOR i = 0, n−1 DO PRINT, x[i], logx[i], $ ;Print values to terminal
;    FORMAT = "(2F12.5)"
END
   
