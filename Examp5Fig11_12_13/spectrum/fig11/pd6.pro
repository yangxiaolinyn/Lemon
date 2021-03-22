      pro pd6

      oldn=!D.name & set_plot,'ps'
 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      pi1=3.141592653589793 / 2.   
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      readfile, infile = './Iquv.dat', Irluv_I, Irluv_q, Irluv_u, Irluv_v, nf

      Irlu_theta0_180 = fltarr(nf)
      for i=0, nf/2-1 do begin
          Irlu_theta0_180(i) =  (100.0 - i) * 90.0 / (nf/2.-1.)
      endfor
      for i=0, nf/2-1 do begin
          Irlu_theta0_180(i + nf/2) =  - i * 90.0 / (nf/2.-1.)
      endfor
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      readfile, infile = './Iquv90.dat', Irluv_I90, Irluv_q90, Irluv_u90, Irluv_v90, nf

      Irlu_theta_90 = fltarr(nf)
      for i=0,100 do begin
          Irlu_theta_90(i) =  (100.0 - i) * 90.0 / (nf-1)
      endfor  
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      readfile, infile = './IQUVphi=90_mu0=-0.2000.dat', IQUV_I90, IQUV_Q90, IQUV_U90, IQUV_V90, nf 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

      nf1 = nf
      theta90 = fltarr(nf)
      for i=1, nf-1 do begin
          theta90(i) = acos(1./(nf-1.) * i ) / !Pi * 180.
      endfor
      theta90(0) = acos(0.001) / !Pi * 180.  
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      readfile, infile = './IQUVphi=0180_mu0=-0.2000.dat', I_0_180, Q_0_180, U_0_180, V_0_180, nf 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

      Irlu_theta0180 = fltarr(nf)
      for i=1, nf/2-1 do begin
          Irlu_theta0180(i) =  acos(1./(nf/2.-1.) * i ) / !Pi * 180.;(nf/2-1. - i) * 90.0 / (nf/2.-1.)
      endfor
      for i=0, nf/2-2 do begin
          Irlu_theta0180(i + nf/2) = -acos(1./(nf/2.-1.) * (nf/2-1-i) ) / !Pi * 180.;- i * 90.0 / (nf/2.-1.)
      endfor
      Irlu_theta0180(0) = acos(0.001) / !Pi * 180.
      Irlu_theta0180(nf-1) = -acos(0.001) / !Pi * 180. 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
     ;print, Irlu_theta
 

      ratio=1./0.85
      l=16 & xxss=l*(ratio) & yyss=l
      ;xoff=(LL-xxss)/2.,yoff=(3*LL/2.-yyss)/2.,
      !p.font = 0
      device,filename='./reflec_no0QUV_12.ps',xsize=xxss,ysize=yyss,bits_per_pixel=8,$
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
      yellow = 251
      cyan = 253
      magenta = 254

      xlen=0.6
      ylen=0.6
      xlen2=0.3
      ylen2=0.3
      xgap1=0.12
      xgap2=0.02
      xgap3 =0.11
      xlow = -90.
      xup = 90.0
      ylow = - 0.15
      yup =    0.55

      csize = 1.3

      xlen2= ( 1. - xgap1 - xgap2 - xgap3 ) / 3.
      xlen  = xlen2 * 2.
      x1 = xgap1
      x2 = x1 + xlen
      posup=[x1, (1.-ylen)/2., x2, (1.+ylen)/2.] 
      plot,[xlow,xup],[ylow,yup],pos=[posup],/noerase,/nodata,/device,$
         xrange=[xlow,xup],yrange=[ylow,yup],/ynozero,/normal,xstyle=4+1,ystyle=4+1
 
      axis,xaxis=1,xticks=6,xminor=4,xtickname=replicate(' ',15),charsize=csize
      axis,xaxis=0,xticks=6,xminor=4,xrange=[xlow,xup],xstyle=1,$;font=-1,$;,xtickname=replicate(' ',6),$;,
      charsize=csize,xtitle=textoidl('\theta');,xthick=tickth,color=colors

      axis,yaxis=0,ytitle=textoidl('Intensity'),yticks=7,yminor=2,yrange=[ylow,yup],ystyle=1,charsize=csize
      axis,yaxis=1,ytickname=replicate(' ',12),yticks=7,yminor=2,charsize=csize
      dyy = 0.35
      y0 = 5
      rat = 1.0
      theta0 = acos(0.8)
      n_1 =  2
      max_I = max( abs(I_0_180))
      max_I2 = max( abs(Irluv_I) ) * 0.990
      num_scale_v = max(abs(Irluv_v)) / max(abs(V_0_180))
      print, max_I, max_I2 
      for i=0,0 do begin    
          oplot, Irlu_theta0_180,  Irluv_I  ,  thick=1, color=black, linestyle=6;,psym=-4   
          oplot, Irlu_theta0_180,  -Irluv_q ,  thick=1, color=black, linestyle=6;,psym=-4  
          oplot, Irlu_theta0_180,  Irluv_u, thick=1, color=black, linestyle=6;,psym=-4  
          oplot, Irlu_theta0_180,  Irluv_v*0.1+0.3, thick=1, color=black, linestyle=6;,psym=-4  

 
          oplot, Irlu_theta0180,  I_0_180/max_I*max_I2, thick=5, color=black, linestyle=1;,psym=-4  
          oplot, Irlu_theta0180,  -Q_0_180/max_I*max_I2, thick=5, color=black, linestyle=1;,psym=-4  
          oplot, Irlu_theta0180,  U_0_180/max_I*max_I2, thick=5, color=black, linestyle=1;,psym=-4  
          oplot, Irlu_theta0180,  V_0_180*num_scale_v*0.1+0.3, thick=5, color=black, linestyle=1;,psym=-4  

          ;xyouts, 25., 0.9  , 'I'
          ;xyouts, 25., 0.21 , 'Q'
          ;xyouts,  80., 0.2 , 'U'
          ;xyouts,  -45., -0.4 , 'V'

          xyouts, Irlu_theta0_180[nf/2], Irluv_I[nf/2]+0.015  , 'I',charsize=csize
          xyouts, Irlu_theta0_180[nf/4], -Irluv_q[nf/4]+0.015 , 'Q',charsize=csize
          xyouts, Irlu_theta0_180[nf/4], Irluv_u[nf/4]+0.015 , 'U',charsize=csize
          xyouts, Irlu_theta0_180[nf/2], Irluv_v[nf/2]*0.1+0.3+0.015 , 'V',charsize=csize
  
  
      endfor 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      x3 = x2 + xgap2
      x4 = x3 + xlen2

      xlow = 0.
      xup = 90.
      ylow = -0.05
      yup =    0.4

      posup=[x3, (1.-ylen)/2., x4, (1.+ylen)/2.] 
      plot,[xlow,xup],[ylow,yup],pos=[posup],/noerase,/nodata,/device,$
         xrange=[xlow,xup],yrange=[ylow,yup],/ynozero,/normal,xstyle=4+1,ystyle=4+1
 
      axis,xaxis=1,xticks=3,xminor=2,xtickname=replicate(' ',15),charsize=csize
      axis,xaxis=0,xticks=3,xminor=2,xrange=[xlow,xup],xstyle=1,$;font=-1,$;,xtickname=replicate(' ',6),$;,
      charsize=csize,xtitle=textoidl('\theta');,xthick=tickth,color=colors

      ;axis,yaxis=0,ytitle=textoidl('I_r+I_l'),yticks=6,yminor=2,yrange=[ylow,yup],ystyle=1
      axis,yaxis=0,ytickname=replicate(' ',12),yticks=5,yminor=2,charsize=csize
      ;axis,yaxis=1,ytickname=replicate(' ',12),yticks=7,yminor=2
      axis,yaxis=1,ytitle=textoidl('Intensity'),yticks=5,yminor=2,yrange=[ylow,yup],ystyle=1,charsize=csize
      ;max_I = max( abs(IQUV_I90) )
      ;max_I2 = max( abs(Irluv_I90) ) * 1.;025
      scal_1 = max( abs(Irluv_I90) ) / max( abs(IQUV_I90) )*0.99 
      scal_v1 = max(abs(Irluv_v90)) / max(abs(IQUV_V90))
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
           oplot, Irlu_theta_90(*),   Irluv_I90 , thick=1, color=black, linestyle=6;,psym=-4    
           oplot, Irlu_theta_90(*),   Irluv_q90 , thick=1, color=black, linestyle=6;,psym=-4
           oplot, Irlu_theta_90(*),  -Irluv_u90,  thick=1, color=black, linestyle=6;,psym=-4 
           oplot, Irlu_theta_90(*),  -Irluv_v90,  thick=1, color=black, linestyle=6;,psym=-4 


           oplot, theta90(*),   IQUV_I90(*) * scal_1,  thick=5, color=black, linestyle=1;,psym=-4 
           oplot, theta90(*),   IQUV_Q90(*) * scal_1,  thick=5, color=black, linestyle=1;,psym=-4 
           oplot, theta90(*),  -IQUV_U90(*) * scal_1,  thick=5, color=black, linestyle=1;,psym=-4 
           oplot, theta90(*),  -IQUV_V90(*) * scal_v1, thick=5, color=black, linestyle=1;,psym=-4 

          ;xyouts, 60., 0.98 , 'I'
          ;xyouts, 75., 0.25 , 'Q'
          ;xyouts, 15., 0.18 , 'V'
          ;x;youts,  68., 0.4 , 'U'
          xyouts, theta90(nf1/4), Irluv_I90(nf1/4)+0.006 , 'I',charsize=csize
          xyouts, theta90(nf1/4), Irluv_q90(nf1/4)-0.04 , 'Q',charsize=csize
          xyouts, 15, -Irluv_u90(75) , 'V',charsize=csize
          xyouts, 15, -Irluv_v90(75)+0.006 , 'U',charsize=csize
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ;print, Ir + Il / max( Ir + Il )
      device,/close 

      ;free_lun,lunAo10
      set_plot,oldn

      end


PRO READFILE, infile = infile, I_r, I_l, I_u, I_v, n
;Demonstrate reading from a file 
;infile = !Bowman + ’data/table.txt’ ;Input file name
n = FILE_LINES(infile)
      I_l=fltarr(n)
      I_r=fltarr(n)
      I_u=fltarr(n)
      I_v=fltarr(n)
;x = FLTARR(n)  
 print, ' data number= ', n  
 
 OPENR, iunit, infile, /GET_LUN ;Open input file
 Point_lun, iunit, 0
  
      FOR i = 0, n-1 DO BEGIN
         ;ReadF,lunAo10, x1, x2, x3, x4;, FORMAT = "(4F22.17)"
         READF, iunit, xr, xl, xu, xv;, FORMAT = "(4F22.17)"
         I_r[i] = xr
         I_l[i] = xl
         I_u[i] = xu
         I_v[i] = xv
      ENDFOR
FREE_LUN, iunit ;Close input file
;FOR i = 0, n−1 DO PRINT, x[i], logx[i], $ ;Print values to terminal
;    FORMAT = "(2F12.5)"
END
   
