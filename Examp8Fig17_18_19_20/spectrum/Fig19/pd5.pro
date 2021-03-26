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
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  I_Theta=05_1.txt
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      ;readfile, infile = './IQUV_E2.txt',I0, I1, I2, I3, I4, I5, I6, ns
      ;readfile, infile = './Fig3_1.txt',I0, I1, I2, I3, I4, I5, I6, ns
      readfile1, infile = './dataTe=9.607E-02tau=6.000E-02.dat', I_mu_6, ns
      ;print, I0(0, *)
     ; I0(0, 200:*) = 0.
     ; I0(1, 200:*) = 0.
     ; I0(2, 200:*) = 0.
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
 

      ratio=3./0.95
      l=16 & xxss=l*(ratio) & yyss=l
      ;xoff=(LL-xxss)/2.,yoff=(3*LL/2.-yyss)/2.,
      !p.font = 0
      device,filename='./fig_12.ps',xsize=xxss,ysize=yyss,bits_per_pixel=8,$
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
 
      xgap1=0.05
      xgap2=0.005
      xgap =0.01

      ;xlen2= ( 1. - xgap1 - xgap2 - xgap3 ) / 3.
      dyy = 0.35
      y0 = 5
      rat = 1.0
      theta0 = acos(0.8)
      n_1 =  2 

      m = 3
      ygap_up = 0.015
      ygap_d = 0.145
      ygap = 0.01
      ylen = (1. - ygap_up - ygap_d ) / 1

      xlow = -1
      xup = 5
      ylow = .5
      yup =    1.2

      h_gauss = 0.04
      pdel = 1.0

      xlen = (1. - xgap1 - xgap2 - xgap * (m-1.) ) / 3.
      for i=0,m-1 do begin    
 
          ;xlen  = ( 1. - xgap1 - xgap2 ) 
          ;x1 = xgap1
          ;x2 = x1 + xlen  
          csize = 2. 
          IF i EQ 0 THEN BEGIN      
              ylow = -8
              yup  = -0
          endif  
          IF i EQ 2 THEN BEGIN      
              ylow = -8.
              yup  = -0
          endif 
          IF i EQ 1 THEN BEGIN      
              ylow = -8
              yup  = 0
          endif 
          IF i EQ 3 THEN BEGIN      
              ylow = -8.
              yup  = 0
          endif   
 
          ;posup=[x1, ygap_d + (m-i-1)*ylen + ygap * (m-i-1), x2, ygap_d + (m-i)*ylen + ygap * (m-i-1) ] 
          x1 = xgap1 + i * xlen + xgap * i
          x2 = x1 + xlen
          print, x1, x2, xlen, ylow
          posup=[x1, ygap_d , x2, ygap_d + ylen ] 
          plot,[xlow,xup],[ylow,yup],pos=[posup],/noerase,/nodata,/device,$
             xrange=[xlow,xup],yrange=[ylow,yup],/ynozero,/normal,xstyle=4+1,ystyle=4+1
 
         ;axis,xaxis=0,xticks=6,xminor=4,xrange=[xlow,xup],xstyle=1,$;font=-1,$;,xtickname=replicate(' ',6),$;,
          ;charsize=1,xtitle=textoidl('\theta');,xthick=tickth,color=colors
         nt = 1000 * 5.99 / 6

          IF i EQ 0 THEN BEGIN
              axis,yaxis=0,ytitle=textoidl('log(vLv erg/s)'),yticks=(yup-ylow),yminor=2, $
                yrange=[ylow,yup],ystyle=1, charsize= csize
              axis,yaxis=1,ytickname=replicate(' ',12),yticks=(yup-ylow),yminor=2
              axis,xaxis=1,xticks=6,xminor=5,xtickname=replicate(' ',15) 
              ;axis,xaxis=0,xticks=6,xminor=2,xtickname=replicate(' ',15) 
              axis,xaxis=0,xticks=6,xminor=5,xrange=[xlow,xup],xstyle=1,$ 
              charsize=csize,xtitle=textoidl('log(h\nu/kT_{bb})')
              ;,xthick=tickth,color=colorsfont=-1,$;,xtickname=replicate(' ',6),$;,

              ;nt = 1000 * 5.9 / 6
              ns = 0
              i_mu = 0
              tk1 = 6.5
              max_I =  max( I_mu_6(i_mu, 6, *) ) 
              oplot, lg10E, alog10( I_mu_6(i_mu, 0, *) / max_I ) , $
                        thick=tk1, color=blue, linestyle=6;,psym=-4  
              oplot, lg10E, alog10( I_mu_6(i_mu, 1, *) / max_I ) , $
                        thick=tk1, color=green, linestyle=6;,psym=-4  
              oplot, lg10E, alog10( I_mu_6(i_mu, 2, *) / max_I ) , $
                        thick=tk1, color=cyan, linestyle=6;,psym=-4  
              oplot, lg10E, alog10( I_mu_6(i_mu, 3, *) / max_I ) , $
                        thick=tk1, color=255, linestyle=6;,psym=-4  
              oplot, lg10E, alog10( I_mu_6(i_mu, 4, *) / max_I ) , $
                        thick=tk1, color=254, linestyle=6;,psym=-4 
              oplot, lg10E, alog10( I_mu_6(i_mu, 5, *) / max_I ) , $
                        thick=tk1, color=red, linestyle=6;,psym=-4    
              oplot, lg10E, alog10( I_mu_6(i_mu, 6, *) / max_I ) , $
                        thick=tk1, color=black, linestyle=6;,psym=-4  


              dy = 0.5
              xyouts, 2.5, -0.9, textoidl('\Theta = 0.188'), charsize=1.8
              xyouts, 2.5, -0.9 - dy, textoidl('\tau = 0.06'), charsize=1.8
              xyouts, 2.5, -0.9 - 2*dy, textoidl('\mu = 0.1'), charsize=1.8
   
          endif  
 
          IF i EQ 1 THEN BEGIN
              ;axis,yaxis=0,ytitle=textoidl('log(vLv)'),yticks=(yup-ylow)/pdel,yminor=2,$
              ;              yrange=[ylow,yup],ystyle=1
              axis,yaxis=0,ytickname=replicate(' ',12),yticks=(yup-ylow)/pdel,yminor=2
              axis,yaxis=1,ytickname=replicate(' ',12),yticks=(yup-ylow)/pdel,yminor=2
              axis,xaxis=1,xticks=6,xminor=5,xtickname=replicate(' ',15) 
              ;axis,xaxis=0,xticks=6,xminor=2,xtickname=replicate(' ',15) 
              axis,xaxis=0,xticks=6,xminor=5,xrange=[xlow,xup],xstyle=1,$
                      ;font=-1,$;,xtickname=replicate(' ',6),$;,
              charsize=csize,xtitle=textoidl('log(h\nu/kT_{bb})');,xthick=tickth,color=colors

              i_mu = 1 
              max_I =  max( I_mu_6(i_mu, 6, *) ) 
              tk1 = 6.5
              oplot, lg10E, alog10( I_mu_6(i_mu, 0, *) / max_I ) , $
                        thick=tk1, color=blue, linestyle=6;,psym=-4  
              oplot, lg10E, alog10( I_mu_6(i_mu, 1, *) / max_I ) , $
                        thick=tk1, color=green, linestyle=6;,psym=-4  
              oplot, lg10E, alog10( I_mu_6(i_mu, 2, *) / max_I ) , $
                        thick=tk1, color=cyan, linestyle=6;,psym=-4  
              oplot, lg10E, alog10( I_mu_6(i_mu, 3, *) / max_I ) , $
                        thick=tk1, color=255, linestyle=6;,psym=-4  
              oplot, lg10E, alog10( I_mu_6(i_mu, 4, *) / max_I ) , $
                        thick=tk1, color=254, linestyle=6;,psym=-4 
              oplot, lg10E, alog10( I_mu_6(i_mu, 5, *) / max_I ) , $
                        thick=tk1, color=red, linestyle=6;,psym=-4     
              oplot, lg10E, alog10( I_mu_6(i_mu, 6, *) / max_I ) , $
                        thick=tk1, color=black, linestyle=6;,psym=-4 

              dy = 0.5
              ;xyouts, 2.5, -0.9, textoidl('\Theta = 0.188'), charsize=1.8
              ;xyouts, 2.5, -0.9 - dy, textoidl('\tau = 0.06'), charsize=1.8
              xyouts, 2.5, -0.9 - 2*dy, textoidl('\mu = 0.5'), charsize=1.8
               
              ;oplot, lg10E, alog10(I6(0, *) / max_I  ) , $
              ;          thick=4, color=black, linestyle=6;,psym=-4 

              ;oplot, lg10E, alog10(I6(2, *) / max_I  ) , $
              ;          thick=4, color=black, linestyle=6;,psym=-4 
          endif 
 
          IF i EQ 2 THEN BEGIN
              ;axis,yaxis=0,ytitle=textoidl('log(vLv)'),yticks=(yup-ylow),yminor=2,yrange=[ylow,yup],ystyle=1
              axis,yaxis=1,ytickname=replicate(' ',12),yticks=(yup-ylow),yminor=2
              axis,yaxis=0,ytickname=replicate(' ',12),yticks=(yup-ylow)/pdel,yminor=2
              axis,xaxis=1,xticks=6,xminor=5,xtickname=replicate(' ',15) 
              ;axis,xaxis=0,xticks=6,xminor=2,xtickname=replicate(' ',15) 
              axis,xaxis=0,xticks=6,xminor=5,xrange=[xlow,xup],xstyle=1,$ 
              charsize=csize,xtitle=textoidl('log(h\nu/kT_{bb})')
              ;axis,xaxis=0,xticks=6,xminor=4,xrange=[xlow,xup],xstyle=1,$
                      ;font=-1,$;,xtickname=replicate(' ',6),$;,
              ;charsize=1,xtitle=textoidl('log(E/m_ec^2)');,xthick=tickth,color=colors

              print, ylow
              ;nt = 1000 * 5 / 6
              i_mu = 2   
              tk1 = 6.5
              max_I =  max( I_mu_6(i_mu, 6, *) ) 
              oplot, lg10E, alog10( I_mu_6(i_mu, 0, *) / max_I ) , $
                        thick=tk1, color=blue, linestyle=6;,psym=-4  
              oplot, lg10E, alog10( I_mu_6(i_mu, 1, *) / max_I ) , $
                        thick=tk1, color=green, linestyle=6;,psym=-4  
              oplot, lg10E, alog10( I_mu_6(i_mu, 2, *) / max_I ) , $
                        thick=tk1, color=cyan, linestyle=6;,psym=-4  
              oplot, lg10E, alog10( I_mu_6(i_mu, 3, *) / max_I ) , $
                        thick=tk1, color=255, linestyle=6;,psym=-4  
              oplot, lg10E, alog10( I_mu_6(i_mu, 4, *) / max_I ) , $
                        thick=tk1, color=254, linestyle=6;,psym=-4 
              oplot, lg10E, alog10( I_mu_6(i_mu, 5, *) / max_I ) , $
                        thick=tk1, color=red, linestyle=6;,psym=-4      
              oplot, lg10E, alog10( I_mu_6(i_mu, 6, *) / max_I ) , $
                        thick=tk1, color=black, linestyle=6;,psym=-4 
              ;print, I_mu_6(i_mu, 0, *)
              dy = 0.5
              ;xyouts, 2.5, -0.9, textoidl('\Theta = 0.188'), charsize=1.8
              ;xyouts, 2.5, -0.9 - dy, textoidl('\tau = 0.06'), charsize=1.8
              xyouts, 2.5, -0.9 - 2*dy, textoidl('\mu = 0.9'), charsize=1.8
          endif  
          IF i EQ 3 THEN BEGIN
              axis,yaxis=0,ytitle=textoidl('\delta'),yticks=(yup-ylow)/pdel,yminor=2,yrange=[ylow,yup],ystyle=1
              axis,yaxis=1,ytickname=replicate(' ',12),yticks= (yup-ylow)/pdel,yminor=2
              axis,xaxis=1,xticks=6,xminor=2,xtickname=replicate(' ',15) 
              axis,xaxis=0,xticks=6,xminor=2,xrange=[xlow,xup],xstyle=1,$
                      ;font=-1,$;,xtickname=replicate(' ',6),$;,
              charsize=1,xtitle=textoidl('log(E/m_ec^2)');,xthick=tickth,color=colors

              ;nt = 1000 * 5 / 6
              ns = 1 
              oplot, lg10E(ns:nt), sqrt( Q50(ns:nt)^2+U50(ns:nt)^2 ) /  I50(ns:nt)  , $
                        thick=4, color=blue, linestyle=6;,psym=-4   
          endif 
      endfor 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      device,/close 

      ;free_lun,lunAo10
      set_plot,oldn

      end


PRO READFILE, infile = infile, I0, I1, I2, I3, I4, I5, I6, n
;Demonstrate reading from a file 
;infile = !Bowman + ’data/table.txt’ ;Input file name
n = FILE_LINES(infile) / 3
      I0=fltarr(3, n)
      I1=fltarr(3,n)
      I2=fltarr(3,n)
      I3=fltarr(3,n)
      I4=fltarr(3,n)
      I5=fltarr(3,n)
      I6=fltarr(3,n) 
;x = FLTARR(n)  
 print, ' data number= ', n  
 
 OPENR, iunit, infile, /GET_LUN ;Open input file
 Point_lun, iunit, 0
  
    FOR j = 0, 2 DO BEGIN
      FOR i = 0, n-1 DO BEGIN
         ;ReadF,lunAo10, x1, x2, x3, x4;, FORMAT = "(4F22.17)"
         READF, iunit, x0, x1, x2, x3, x4, x5, x6;, FORMAT = "(4F22.17)"
         I0[j, i] = x0
         I1[j, i] = x1
         I2[j, i] = x2
         I3[j, i] = x3
         I4[j, i] = x4
         I5[j, i] = x5
         I6[j, i] = x6
         ;print, x3, x4, x5, x6
      ENDFOR
    ENDFOR
FREE_LUN, iunit ;Close input file
;FOR i = 0, n−1 DO PRINT, x[i], logx[i], $ ;Print values to terminal
;    FORMAT = "(2F12.5)"
END
   



PRO READFILE1, infile = infile, I_mu_6, n
    ;Demonstrate reading from a file 
    ;infile = !Bowman + ’data/table.txt’ ;Input file name
    n = FILE_LINES(infile) / 3
    I_mu_6 = fltarr( 3, 7, n )

      ;I0=fltarr(3, n)
      ;I1=fltarr(3,n)
      ;I2=fltarr(3,n)
      ;I3=fltarr(3,n)
      ;I4=fltarr(3,n)
      ;I5=fltarr(3,n)
      ;I6=fltarr(3,n) 
;x = FLTARR(n)  
 print, ' data number= ', n  
 
 OPENR, iunit, infile, /GET_LUN ;Open input file
 Point_lun, iunit, 0
  
    FOR j = 0, 2 DO BEGIN
      FOR i = 0, n-1 DO BEGIN
         ;ReadF,lunAo10, x1, x2, x3, x4;, FORMAT = "(4F22.17)"
         READF, iunit, x0, x1, x2, x3, x4, x5, x6;, FORMAT = "(4F22.17)"
         I_mu_6[j, 0, i] = x0
         I_mu_6[j, 1, i] = x1
         I_mu_6[j, 2, i] = x2
         I_mu_6[j, 3, i] = x3
         I_mu_6[j, 4, i] = x4
         I_mu_6[j, 5, i] = x5
         I_mu_6[j, 6, i] = x6
         ;print, x3, x4, x5, x6
      ENDFOR
    ENDFOR
FREE_LUN, iunit ;Close input file
;FOR i = 0, n−1 DO PRINT, x[i], logx[i], $ ;Print values to terminal
;    FORMAT = "(2F12.5)"
END
