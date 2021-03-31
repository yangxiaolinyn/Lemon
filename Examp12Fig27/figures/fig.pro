      pro fig

      oldn=!D.name & set_plot,'ps'
 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      pi1=3.141592653589793 / 2. 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      ;readfile, infile = './I_Q_V.txt', Q_s, U_s, V_s, nf
      ;readfile, infile = './I_Q_V4.txt', I_ifm, Q_ifm, U_ifm, V_ifm, nf
      readfile, infile = './MC_QUV.dat', I_i1, Q_i1, U_i1, V_i1, nf

      s_i1 = fltarr(nf+1)
      for i=0, nf do begin
          s_i1(i) = (3.20 / nf) * i
      endfor  

      readfile3, infile = './QUV.dat', Q_quv, U_quv, V_quv, nf2
 
      s_quv = fltarr(nf2+1)
      for i=0, nf2 do begin
          s_quv(i) = (3.20 / nf2) * i
      endfor  

;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

      readfile2, infile = './I_Q.dat', I_iq, Q_iq, nf2
      readfile, infile = './MC_IQ.dat', I_i2, Q_i2, U_i2, V_i2, nf

      s_iq = fltarr(nf2+1)
      for i=0, nf2 do begin
          s_iq(i) = (10.0 / nf2) * i
      endfor  
      s_i2 = fltarr(nf+1)
      for i=0, nf do begin
          s_i2(i) = (10.0 / nf) * i
      endfor  
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
     ;print, Irlu_theta
 

      ratio=1./0.4
      l=16 & xxss=l*(ratio) & yyss=l
      ;xoff=(LL-xxss)/2.,yoff=(3*LL/2.-yyss)/2.,
      !p.font = 0
      device,filename='./figs_1.ps',xsize=xxss,ysize=yyss,bits_per_pixel=8,$
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

      ylen=0.9
      ;xlen2=0.3
      ylen2=0.3
      xgap1=0.11
      xgap2=0.02
      xgap3 =0.01
      xlen=1 - xgap3 - xgap1
      ylow = -0.4
      yup =  1.2

      cs = 2.0
      xlen2= ( 1. - xgap1 - xgap2 - xgap3 ) / 3.
      xlen  = xlen2 * 2.
      xlen=1 - xgap3 - xgap1
      x1 = xgap1
      x2 = x1 + xlen
 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      xgs = fltarr(3)
      ygs = fltarr(3)
      xgs(0) = 0.08
      xgs(1) = 0.01
      xgs(2) = 0.06
      ygs(0) = 0.02
      ygs(1) = 0.14
      ygs(2) = 0.005
      n_x = 2
      n_y = 1
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      panel_set, m_x = n_x, m_y = n_y, xgaps = xgs, ygaps = ygs, $
                   xl_arr, xr_arr, yu_arr, yd_arr

      
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      for j=0, n_x-1 do begin    
        for i=0, n_y-1 do begin 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          IF (i EQ 0) and (j EQ 0) THEN BEGIN  
              ylow = -0.4
              yup =  1.2
              xlow = 0.
              xup = 10.0
          endif   
          IF (i EQ 0) and (j EQ 1) THEN BEGIN    
              ylow = -0.8
              yup  = 1.4
              xlow = 0.
              xup = 3.0
          endif  

;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          posup=[xl_arr[j], yd_arr[n_y-1-i] , xr_arr[j], yu_arr[n_y-1-i] ] 
          plot,[xlow,xup],[ylow,yup],pos=[posup],/noerase,/nodata,/device,$
             xrange=[xlow,xup],yrange=[ylow,yup],/ynozero,/normal,xstyle=4+1,ystyle=4+1
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
       IF (i EQ 0) and (j EQ 0) THEN BEGIN 
 
          axis,xaxis=1,xticks=10,xminor=5,xtickname=replicate(' ',15), charsize=cs
          axis,xaxis=0,xticks=10,xminor=5,xrange=[xlow,xup],xstyle=1,$
                 ;font=-1,$;,xtickname=replicate(' ',6),$;,
          charsize=cs,xtitle=textoidl('s');,xthick=tickth,color=colors

          axis, yaxis=0, ytitle=textoidl('Intensity'), yticks=8, yminor=4, $
                     yrange=[ylow,yup],ystyle=1, charsize=cs
          axis, yaxis=1, ytickname=replicate(' ',12), yticks=8, yminor=4, charsize=cs
   
          oplot, s_iq,  I_iq/max( abs( I_iq ) ), thick=1, color=blue, linestyle=6;,psym=-4  
          oplot, s_iq,  Q_iq/max( abs( I_iq ) ), thick=1, color=red, linestyle=6;,psym=-4  
          ;oplot, s_i,  U_ifm/I_ifm, thick=2, color=red, linestyle=6;,psym=-4 
          ;oplot, s_i,  V_ifm/I_ifm, thick=2, color=green, linestyle=6;,psym=-4 

          oplot, s_i2,  I_i2/max( abs( I_i2 ) ), thick=5, color=black, linestyle=1;,psym=-4   
          oplot, s_i2,  Q_i2/max( abs( I_i2 ) ) , thick=5, color=black, linestyle=1;,psym=-4  
          ;oplot, s_i,  U_i/I_i , thick=2, color=red, linestyle=1;,psym=-4 
          ;oplot, s_i,  V_i/I_i , thick=2, color=green, linestyle=1;,psym=-4   
  
          ;print, Q_i/max( abs( Q_i ) )
          ;print, U_i/max( abs( U_i ) )
          xyouts, 5., 1.0  , 'I', charsize=cs
          xyouts, 5., -0.23 , 'Q', charsize=cs
          ;xyouts,  80., 0.2 , 'U'
          ;xyouts,  -45., -0.4 , 'V'
  
      endif 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
       IF (i EQ 0) and (j EQ 1) THEN BEGIN 
 
          axis,xaxis=1,xticks=3,xminor=5,xtickname=replicate(' ',15), charsize=cs
          axis,xaxis=0,xticks=3,xminor=5,xrange=[xlow,xup],xstyle=1,$
                 ;font=-1,$;,xtickname=replicate(' ',6),$;,
          charsize=cs,xtitle=textoidl('s');,xthick=tickth,color=colors

          axis, yaxis=0, ytitle=textoidl(''), yticks=11, yminor=4, $
                     yrange=[ylow,yup],ystyle=1, charsize=cs
          axis, yaxis=1, ytickname=replicate(' ',12), yticks=11, yminor=4, charsize=cs
    
          oplot, s_quv,  Q_quv/max( abs( Q_quv ) ), thick=1, color=blue, linestyle=6;,psym=-4 
          oplot, s_quv,  U_quv/max( abs( Q_quv ) ), thick=1, color=red, linestyle=6;,psym=-4 
          oplot, s_quv,  V_quv/max( abs( Q_quv ) ), thick=1, color=green, linestyle=6;,psym=-4  
          ;oplot, s_i,  U_ifm/I_ifm, thick=2, color=red, linestyle=6;,psym=-4 
          ;oplot, s_i,  V_ifm/I_ifm, thick=2, color=green, linestyle=6;,psym=-4 
   
          oplot, s_i1,  Q_i1/max( abs( Q_i1 ) ) , thick=5, color=black, linestyle=1;,psym=-4  
          oplot, s_i1,  U_i1/max( abs( Q_i1 ) ) , thick=5, color=black, linestyle=1;,psym=-4  
          oplot, s_i1,  V_i1/max( abs( Q_i1 ) ) , thick=5, color=black, linestyle=1;,psym=-4  
          ;oplot, s_i,  U_i/I_i , thick=2, color=red, linestyle=1;,psym=-4 
          ;oplot, s_i,  V_i/I_i , thick=2, color=green, linestyle=1;,psym=-4   
  
          ;print, Q_i/max( abs( Q_i ) )
          ;print, U_i/max( abs( U_i ) )
          ;xyouts, 5., 1.0  , 'I', charsize=cs
          ;xyouts, 5., -0.23 , 'Q', charsize=cs
          ;xyouts,  80., 0.2 , 'U'
          ;xyouts,  -45., -0.4 , 'V'
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
             dy = 0.13 
             dy2 = 0.05
             dx = 0.6
             cha_size = 1.2
             x0 = 0.2
             y0 = 1.2
             oplot, [x0, x0+dx], [y0, y0] , thick=1, color=black, linestyle=6;,psym=-4 
             oplot, [x0, x0+dx], [y0-dy, y0-dy] , thick=1, color=black, linestyle=6;,psym=-4 
             oplot, [x0, x0+dx], [y0-dy*2.0, y0-dy*2.0] , thick=1, color=black, linestyle=6;,psym=-4 

             oplot, [x0, x0+dx], [y0, y0] , thick=5, color=blue, linestyle=1;,psym=-4 
             oplot, [x0, x0+dx], [y0-dy, y0-dy] , thick=5, color=red, linestyle=1;,psym=-4 
             oplot, [x0, x0+dx], [y0-dy*2.0, y0-dy*2.0] , thick=5, color=green, linestyle=1;,psym=-4  

             xyouts, x0+dx*1.1, y0-dy2  , 'Q', charsize=cs
             xyouts, x0+dx*1.1, y0-dy-dy2 , 'U', charsize=cs
             xyouts, x0+dx*1.1, y0-dy*2.0-dy2 , 'V', charsize=cs
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
      endif 
  
        endfor 
      endfor 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      device,/close 

      ;free_lun,lunAo10
      set_plot,oldn

      end


PRO READFILE, infile = infile, I_s, Q_s, U_s, V_s, n
;Demonstrate reading from a file 
;infile = !Bowman + ’data/table.txt’ ;Input file name
n = FILE_LINES(infile)
      I_s=fltarr(n)
      Q_s=fltarr(n)
      U_s=fltarr(n)
      V_s=fltarr(n)
      ;I_u=fltarr(n)
      ;I_v=fltarr(n)
;x = FLTARR(n)  
 print, ' data number= ', n  
 
 OPENR, iunit, infile, /GET_LUN ;Open input file
 Point_lun, iunit, 0
  
      FOR i = 0, n-1 DO BEGIN
         ;ReadF,lunAo10, x1, x2, x3, x4;, FORMAT = "(4F22.17)"
         READF, iunit, x0, x1, x2, x3;, xu, xv;, FORMAT = "(4F22.17)"
         I_s[i] = x0
         Q_s[i] = x1
         U_s[i] = x2
         V_s[i] = x3
         ;I_v[i] = xv
      ENDFOR
FREE_LUN, iunit ;Close input file
;FOR i = 0, n−1 DO PRINT, x[i], logx[i], $ ;Print values to terminal
;    FORMAT = "(2F12.5)"
END

 

PRO READFILE2, infile = infile, I_s, Q_s, n
;Demonstrate reading from a file 
;infile = !Bowman + ’data/table.txt’ ;Input file name
n = FILE_LINES(infile)
      I_s=fltarr(n)
      Q_s=fltarr(n)
 print, ' data number= ', n  
 
 OPENR, iunit, infile, /GET_LUN ;Open input file
 Point_lun, iunit, 0
  
      FOR i = 0, n-1 DO BEGIN
         ;ReadF,lunAo10, x1, x2, x3, x4;, FORMAT = "(4F22.17)"
         READF, iunit, x0, x1;, xu, xv;, FORMAT = "(4F22.17)"
         I_s[i] = x0
         Q_s[i] = x1
         ;I_v[i] = xv
      ENDFOR
FREE_LUN, iunit ;Close input file
;FOR i = 0, n−1 DO PRINT, x[i], logx[i], $ ;Print values to terminal
;    FORMAT = "(2F12.5)"
END
   


PRO READFILE3, infile = infile, Q_s, U_s, V_s, n
;Demonstrate reading from a file 
;infile = !Bowman + ’data/table.txt’ ;Input file name
n = FILE_LINES(infile)
      Q_s=fltarr(n)
      U_s=fltarr(n)
      V_s=fltarr(n)
      ;I_u=fltarr(n)
      ;I_v=fltarr(n)
;x = FLTARR(n)  
 print, ' file row number= ', n  
 
 OPENR, iunit, infile, /GET_LUN ;Open input file
 Point_lun, iunit, 0
  
      FOR i = 0, n-1 DO BEGIN
         ;ReadF,lunAo10, x1, x2, x3, x4;, FORMAT = "(4F22.17)"
         READF, iunit, x1, x2, x3;, xu, xv;, FORMAT = "(4F22.17)"
         Q_s[i] = x1
         U_s[i] = x2
         V_s[i] = x3
         ;I_v[i] = xv
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


;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PRO xyout, dy = dy, charac_arr = charac_arr, n_arr = n_arr, cha_size = cha_size, x0 = x0, $
              y0 = y0, colors = colors, line_thick = line_thick, y_bias = y_bias
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    FOR i = 0, n_arr-1 DO BEGIN  
      ;oplot, [x0, x0 + length], [y0-dy*i, y0-dy*i], thick=line_thick, color=colors(i), linestyle=6;,psym=-4  
      xyouts, x0, y0-dy*i-y_bias, charac_arr(i), charsize=cha_size
    ENDFOR 
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



