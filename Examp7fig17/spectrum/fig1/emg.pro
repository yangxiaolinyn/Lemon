      pro emg

      oldn=!D.name & set_plot,'ps'
 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      rdfile, infile = './IQUVphi_mu0=0.1000.dat', mu1, nf
      rdfile, infile = './IQUVphi_mu0=0.2000.dat', mu2, nf
      rdfile, infile = './IQUVphi_mu0=0.4000.dat', mu4, nf
      rdfile, infile = './IQUVphi_mu0=0.5000.dat', mu5, nf
      rdfile, infile = './IQUVphi_mu0=0.8000.dat', mu8, nf
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      ns = nf / 2
      theta = fltarr(nf)
      for i=0, ns-1 do begin
          theta(i) = acos(1./(ns-1) * i) / !Pi * 180. 
          theta(2*ns-1-i) = - acos(1./(ns-1) * i ) / !Pi * 180.
      endfor
      theta2 = fltarr(ns)
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      ;readfile, infile = './Iquv60.txt', Iquv_I60, Iquv_q60, Iquv_u60, Iquv_v60, nf
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      ;readfile, infile = './Iquv80.txt', Iquv_I80, Iquv_q80, Iquv_u80, Iquv_v80, nf
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
 

      ratio=4./ 1.4
      l=16 & xxss=l*(ratio) & yyss=l
      ;xoff=(LL-xxss)/2.,yoff=(3*LL/2.-yyss)/2.,
      !p.font = 0
      device,filename='./emg_8.ps',xsize=xxss,ysize=yyss,bits_per_pixel=8,$
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
 
      xgap1=0.035
      xgap2=0.01
      xgap =0.03 

      ;xlen2= ( 1. - xgap1 - xgap2 - xgap3 ) / 3.
      dyy = 0.35
      y0 = 5
      rat = 1.0
      theta0 = acos(0.8)
      n_1 =  2 

      m = 4
      ygap_up = 0.03
      ygap_d = 0.1
      ygap = 0.01
      ylen = (1. - ygap_up - ygap_d ) / 1

      xlow = -90
      xup = 90 

      h_gauss = 0.04
      pdel = 1.0
      Total_Num = 2.*10.^5
      y00 = 0.2

      xlen = (1. - xgap1 - xgap2 - xgap * (m-1.) ) / m
      for i=0,m-1 do begin    
 
          ;xlen  = ( 1. - xgap1 - xgap2 ) 
          ;x1 = xgap1
          ;x2 = x1 + xlen  
          IF i EQ 0 THEN BEGIN
              ylow = 0
              yup  = 20
          endif  
          IF i EQ 1 THEN BEGIN 
              ylow = -2
              yup  = 5
          endif  
          IF i EQ 2 THEN BEGIN 
              ylow = -1
              yup  = 4
          endif  
          IF i EQ 3 THEN BEGIN 
              ylow = -1
              yup  = 9
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
              y00 = 0.4
              axis,yaxis=0,ytitle=textoidl('I'),yticks=5,yminor=4, $
                yrange=[ylow,yup],ystyle=1, charsize=1.3
              axis,yaxis=1,ytickname=replicate(' ',30),yticks=(yup-ylow),yminor=2
              axis,xaxis=1,xticks=6,xminor=2,xtickname=replicate(' ',15) 
              ;axis,xaxis=0,xticks=6,xminor=2,xtickname=replicate(' ',15) 
              axis,xaxis=0,xticks=6,xminor=4,xrange=[xlow,xup],xstyle=1,$ 
              charsize=1.3,xtitle=textoidl('\theta')
              ;,xthick=tickth,color=colorsfont=-1,$;,xtickname=replicate(' ',6),$;,
              oplot,  theta, mu1(0, *)/ Total_Num + 8 * y00 , $
                                thick=3, color=black, linestyle=6;,psym=-4 

              ;Total_Num = max( mu2(0, *) )
              oplot,  theta, mu2(0, *)/ Total_Num + 6 * y00 , $
                                thick=3, color=black, linestyle=6;,psym=-4 

              ;Total_Num = max( mu4(0, *) )
              oplot,  theta, mu4(0, *)/ Total_Num + 4 * y00 , $
                                thick=3, color=black, linestyle=6;,psym=-4 

              ;Total_Num = max( mu5(0, *) )
              oplot,  theta, mu5(0, *)/ Total_Num + 2 * y00 , $
                                thick=3, color=black, linestyle=6;,psym=-4 

              ;Total_Num = max( mu8(0, *) )
              oplot,  theta, mu8(0, *)/ Total_Num + 0 * y00 , $
                                thick=3, color=black, linestyle=6;,psym=-4  
          endif  
 
          IF i EQ 1 THEN BEGIN
              axis,yaxis=0,ytitle=textoidl('Q'),yticks=(yup-ylow)/pdel,yminor=2,$
                            yrange=[ylow,yup],ystyle=1
              y00 = 0.4
              ;axis,yaxis=0,ytickname=replicate(' ',30),yticks=(yup-ylow)/pdel,yminor=2
              axis,yaxis=1,ytickname=replicate(' ',30),yticks=(yup-ylow)/pdel,yminor=2
              axis,xaxis=1,xticks=6,xminor=2,xtickname=replicate(' ',15) 
              ;axis,xaxis=0,xticks=6,xminor=2,xtickname=replicate(' ',15) 
              axis,xaxis=0,xticks=6,xminor=4,xrange=[xlow,xup],xstyle=1,$
                      ;font=-1,$;,xtickname=replicate(' ',6),$;,
              charsize=1.3,xtitle=textoidl('\theta');,xthick=tickth,color=colors
 
              oplot,  theta, mu1(1, *)/ Total_Num + 8 * y00 , $
                                thick=3, color=black, linestyle=6;,psym=-4  
              ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              pointfind, mus = mu1(1, *), xarr = theta, x_0
              pointdraw, x0 = x_0, yc = 8 * y00, col = black
              ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
              oplot,  theta, mu2(1, *)/ Total_Num + 6 * y00 , $
                                thick=3, color=black, linestyle=6;,psym=-4 
              ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              pointfind, mus = mu2(1, *), xarr = theta, x_0
              pointdraw, x0 = x_0, yc = 6 * y00, col = black
              ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

              ;Total_Num = max( mu4(0, *) )
              oplot,  theta, mu4(1, *)/ Total_Num + 4 * y00 , $
                                thick=3, color=black, linestyle=6;,psym=-4 
              ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              pointfind, mus = mu4(1, *), xarr = theta, x_0
              pointdraw, x0 = x_0, yc = 4 * y00, col = black
              ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

              ;Total_Num = max( mu5(0, *) )
              oplot,  theta, mu5(1, *)/ Total_Num + 2 * y00 , $
                                thick=3, color=black, linestyle=6;,psym=-4 
              ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              pointfind, mus = mu5(1, *), xarr = theta, x_0
              pointdraw, x0 = x_0, yc = 2 * y00, col = black
              ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

              ;Total_Num = max( mu8(0, *) )
              oplot,  theta, mu8(1, *)/ Total_Num + 0 * y00 , $
                                thick=3, color=black, linestyle=6;,psym=-4 
              ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              pointfind, mus = mu8(1, *), xarr = theta, x_0
              pointdraw, x0 = x_0, yc = 0 * y00, col = black
              ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          endif 
          IF i EQ 2 THEN BEGIN
              Total_Num = max( abs( mu1(2, *) ) )
              mi = 2
              ;labels=[-90, -60, -30, 0, 30, 60, '  ']
              axis,yaxis=0,ytitle=textoidl('U'),yticks=(yup-ylow)/pdel,yminor=2,$
                            yrange=[ylow,yup],ystyle=1
              ;axis,yaxis=0,ytickname=replicate(' ',30),yticks=(yup-ylow)/pdel,yminor=2
              axis,yaxis=1,ytickname=replicate(' ',30),yticks=(yup-ylow)/pdel,yminor=2
              axis,xaxis=1,xticks=6,xminor=2,xtickname=replicate(' ',15) 
              ;axis,xaxis=0,xticks=6,xminor=2,xtickname=replicate(' ',15) 
              axis,xaxis=0,xticks=6,xminor=4,xrange=[xlow,xup],xstyle=1,$
                      ;font=-1,$;,xtickname=replicate(' ',6),$;,
              charsize=1.3,xtitle=textoidl('\theta'), xtickname=labels;,xthick=tickth,color=colors

              oplot,  theta, mu1(mi, *)/ Total_Num + 8 * y00 , $
                                thick=3, color=black, linestyle=6;,psym=-4 

              ;Total_Num = max( mu2(mi, *) )
              oplot,  theta, mu2(mi, *)/ Total_Num + 6 * y00 , $
                                thick=3, color=black, linestyle=6;,psym=-4 

              ;Total_Num = max( mu4(mi, *) )
              oplot,  theta, mu4(mi, *)/ Total_Num + 4 * y00 , $
                                thick=3, color=black, linestyle=6;,psym=-4 

              ;Total_Num = max( mu5(mi, *) )
              oplot,  theta, mu5(mi, *)/ Total_Num + 2 * y00 , $
                                thick=3, color=black, linestyle=6;,psym=-4 

              ;Total_Num = max( mu8(mi, *) )
              oplot,  theta, mu8(mi, *)/ Total_Num + 0 * y00 , $
                                thick=3, color=black, linestyle=6;,psym=-4  
          endif   
          IF i EQ 3 THEN BEGIN
              Total_Num = 2.*10.^5
              ;Total_Num = max(mu1(2, *))
              mi = 3
              axis,yaxis=0,ytitle=textoidl('V'),yticks=(yup-ylow)/pdel,yminor=2,$
                            yrange=[ylow,yup],ystyle=1
              ;axis,yaxis=0,ytickname=replicate(' ',30),yticks=(yup-ylow)/pdel,yminor=2
              axis,yaxis=1,ytickname=replicate(' ',30),yticks=(yup-ylow)/pdel,yminor=2
              axis,xaxis=1,xticks=6,xminor=2,xtickname=replicate(' ',15) 
              ;axis,xaxis=0,xticks=6,xminor=2,xtickname=replicate(' ',15) 
              axis,xaxis=0,xticks=6,xminor=4,xrange=[xlow,xup],xstyle=1,$
                      ;font=-1,$;,xtickname=replicate(' ',6),$;,
              charsize=1.3,xtitle=textoidl('\theta');,xthick=tickth,color=colors

              oplot,  theta, mu1(mi, *)/ Total_Num + 12 * y00 , $
                                thick=3, color=black, linestyle=6;,psym=-4 

              ;Total_Num = max( mu2(0, *) )
              oplot,  theta, mu2(mi, *)/ Total_Num + 9 * y00 , $
                                thick=3, color=black, linestyle=6;,psym=-4 

              ;Total_Num = max( mu4(0, *) )
              oplot,  theta, mu4(mi, *)/ Total_Num + 6 * y00 , $
                                thick=3, color=black, linestyle=6;,psym=-4 

              ;Total_Num = max( mu5(0, *) )
              oplot,  theta, mu5(mi, *)/ Total_Num + 3 * y00 , $
                                thick=3, color=black, linestyle=6;,psym=-4 

              ;Total_Num = max( mu8(0, *) )
              oplot,  theta, mu8(mi, *)/ Total_Num + 0 * y00 , $
                                thick=3, color=black, linestyle=6;,psym=-4  
          endif      
      endfor 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      device,/close 

      ;free_lun,lunAo10
      set_plot,oldn

      end

PRO pointfind, mus = mus, xarr = xarr, x0
 
      x0=fltarr(5)*0.
      n = n_elements(mus) 
      m = 0
      FOR i = 0, n-2 DO BEGIN
          if(mus(i) * mus(i+1) < 0.)then begin
              x0(m) = ( xarr(i) + xarr(i+1) ) / 2.
              m = m + 1
          endif
      endfor   
END

PRO pointdraw, x0 = x0, yc = yc, col = col
    
      FOR i = 0, 4 DO BEGIN
          if(x0(i) ne 0.)then begin
              xy = Convert_Coord(x0(i), yc, /Data, /To_Device) 
              Polyfill, CIRCLE(xy(0), xy(1), 90.), /Fill, Color=col, /Device 
          endif
      endfor   
END

FUNCTION CIRCLE, xcenter, ycenter, radius
points = (2 * !PI / 99.0) * FIndGen(100)
x = xcenter + radius * Cos(points)
y = ycenter + radius * Sin(points)
RETURN, Transpose([[x],[y]])
END

   

PRO rdfile, infile = infile, I0, n
;Demonstrate reading from a file 
;infile = !Bowman + ’data/table.txt’ ;Input file name
 
n = FILE_LINES(infile) 
      I0=fltarr(4, n)
      ;I1=fltarr(4, n)
      ;I2=fltarr(4, n)
      ;I3=fltarr(4, n) 
;x = FLTARR(n)  
 print, ' data number= ', n  
 
 OPENR, iunit, infile, /GET_LUN ;Open input file
 Point_lun, iunit, 0
  
    ;FOR j = 0, m-1 DO BEGIN
      FOR i = 0, n-1 DO BEGIN
         ;ReadF,lunAo10, x1, x2, x3, x4;, FORMAT = "(4F22.17)"
         READF, iunit, x0, x1, x2, x3;, FORMAT = "(4F22.17)"
         I0[0, i] = x0
         I0[1, i] = x1
         I0[2, i] = x2
         I0[3, i] = x3
         ;I1[j, i] = x1
         ;I2[j, i] = x2
         ;I3[j, i] = x3 
         ;print, x3, x4, x5, x6
      ENDFOR
    ;ENDFOR
FREE_LUN, iunit ;Close input file
;FOR i = 0, n−1 DO PRINT, x[i], logx[i], $ ;Print values to terminal
;    FORMAT = "(2F12.5)"
END



