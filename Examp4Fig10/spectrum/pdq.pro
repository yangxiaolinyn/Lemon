      pro pdq

      oldn=!D.name & set_plot,'ps'

      ns = 21
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      READFILE, infile = './polarizationDegree_ChanD.txt', lv, ns
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      READFILE, infile = './ChandraI.txt', ChandI, ns 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      vv = fltarr(ns)
      for i=0, ns-1 do begin
          vv(i) = 0. + 1./(ns - 1.0) * i
      endfor
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      ;ns = 101  
 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      READFILE, infile = './EIQ_I_tau=0.2000.dat', IQ_I02, ns
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      READFILE, infile = './EIQ_Q_tau=0.2000.dat', IQ_Q02, ns
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      READFILE, infile = './EIQ_I_tau=0.5000.dat', IQ_I05, ns
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      READFILE, infile = './EIQ_Q_tau=0.5000.dat', IQ_Q05, ns
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      READFILE, infile = './EIQ_I_tau=1.0000.dat', IQ_I1, ns
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      READFILE, infile = './EIQ_Q_tau=1.0000.dat', IQ_Q1, ns
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      READFILE, infile = './EIQ_I_tau=2.0000.dat', IQ_I2, ns
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      READFILE, infile = './EIQ_Q_tau=2.0000.dat', IQ_Q2, ns
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      READFILE, infile = './EIQ_I_tau=5.0000.dat', IQ_I5, ns
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      READFILE, infile = './EIQ_Q_tau=5.0000.dat', IQ_Q5, ns
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      cosq1 = fltarr(ns)
      for i=0, ns-1 do begin
          cosq1(i) = ( i ) * 1. / (ns - 1) 
      endfor
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      ratio=1./0.79
      l=16 & xxss=l*(ratio) & yyss=l
      ;xoff=(LL-xxss)/2.,yoff=(3*LL/2.-yyss)/2.,
      !p.font = 0
      device,filename='./pd_Chand_q_1.ps',xsize=xxss,ysize=yyss,bits_per_pixel=8,$
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

      xlen=0.7
      ylen=0.7
      xlow = 0.
      xup = 1.0
      ylow = - 0.
      yup = .18
      cs = 1.5

      posup=[(1.-xlen)/2.,(1.-ylen)/2.,(1.+xlen)/2.,(1.+ylen)/2.] 
      plot,[xlow,xup],[ylow,yup],pos=[posup],/noerase,/nodata,/device,$
         xrange=[xlow,xup],yrange=[ylow,yup],/ynozero,/normal,xstyle=4+1,ystyle=4+1
 
      axis,xaxis=1,xticks=5,xminor=2,xtickname=replicate(' ',15), charsize=cs
      axis,xaxis=0,xticks=5,xminor=2,xrange=[xlow,xup],xstyle=1,$;font=-1,$;,xtickname=replicate(' ',6),$;,
      charsize=cs,xtitle=textoidl('\mu=cos(\theta)');,xthick=tickth,color=colors

      axis,yaxis=0,ytitle=textoidl('Polarization Degree'),yticks=9,yminor=2,$
             yrange=[ylow,yup],ystyle=1, charsize=cs
      axis,yaxis=1,ytickname=replicate(' ',12),yticks=9,yminor=2, charsize=cs
      dyy = 0.35
      y0 = 5
      for i=0,0 do begin
          oplot, vv(*), lv(*) , thick=4, color=black, linestyle=6;,psym=-4  

           oplot, cosq1,  abs(IQ_Q02 )/(IQ_I02 ) , thick=4, color=black, linestyle=2;,psym=-4 
           oplot, cosq1,  abs(IQ_Q05 )/(IQ_I05 ) , thick=4, color=black, linestyle=5;,psym=-4 
           oplot, cosq1,  abs(IQ_Q1 )/(IQ_I1 ) , thick=4, color=black, linestyle=3;,psym=-4 
           oplot, cosq1,  abs(IQ_Q2 )/(IQ_I2 ) , thick=4, color=black, linestyle=4;,psym=-4  
           oplot, cosq1,  abs(IQ_Q5 )/(IQ_I5 ) , thick=4, color=black, linestyle=1;,psym=-4  
           ;oplot, cosq1,  abs(IQ_Q20(0:ns-1))/(IQ_I20(0:ns-1)) , thick=6, color=black, linestyle=1;,psym=-4   
      endfor 
      
      cs = 1.4
      xyouts,0.82,0.158,textoidl('\tau=0.2'),charsize=cs 
      oplot, [0.67, 0.81],[0.16, 0.16], thick=4, color=black, linestyle=2
      xyouts,0.82,0.148,textoidl('\tau=0.5'),charsize=cs 
      oplot, [0.67, 0.81],[0.15, 0.15], thick=4, color=black, linestyle=5
      xyouts,0.82,0.138,textoidl('\tau=1.0'),charsize=cs
      oplot, [0.67, 0.81],[0.14, 0.14], thick=4, color=black, linestyle=3
      xyouts,0.82,0.128,textoidl('\tau=2.0'),charsize=cs
      oplot, [0.67, 0.81],[0.13, 0.13], thick=4, color=black, linestyle=4
      xyouts,0.82,0.118,textoidl('\tau=10.0'),charsize=cs
      oplot, [0.67, 0.81],[0.12, 0.12], thick=4, color=black, linestyle=1
      xyouts,0.82,0.108,textoidl('\tau=\infty'),charsize=cs
      oplot, [0.67, 0.81],[0.11, 0.11], thick=4, color=black, linestyle=6
      ;print, 3.e7 * (1e9)^0.25
      device,/close 
 
      set_plot,oldn

      end
 

;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PRO READFILE, infile = infile, IorQ, n 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    n = FILE_LINES(infile)
    IorQ = fltarr( n )
 
    print, 'file lines = ', n  
 
    OPENR, iunit, infile, /GET_LUN ;Open input file
    Point_lun, iunit, 0
  
    FOR j = 0, n-1 DO BEGIN 
        READF, iunit, x0 
        IorQ[j] = x0
    ENDFOR
    FREE_LUN, iunit
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
END
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




  
