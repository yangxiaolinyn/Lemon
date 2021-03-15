      pro pdi

      oldn=!D.name & set_plot,'ps'

      ns = 21
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Openr,lunAo10,'./polarizationDegree_ChanD.txt',/Get_Lun
      Point_lun,lunAo10,0
      lv=fltarr(ns)
      ReadF,lunAo10,lv
      free_lun,lunAo10
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Openr,lunAo10,'./ChandraI.txt',/Get_Lun
      Point_lun,lunAo10,0
      ChandI=fltarr(ns)
      ReadF,lunAo10,ChandI
      free_lun,lunAo10
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
 
      vv = fltarr(ns)
      for i=0,ns-1 do begin
          vv(i) = 0. + 1./20. * i
      endfor
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

      pi1=3.141592653589793 / 2.
      cosq1 = fltarr(ns)
      for i=0, ns-1 do begin
          cosq1(i) = ( i ) * 1. / (ns - 1);cos( pi1/(ns-1) * i )
      endfor
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ;ns = 101
      Openr,lunAo10,'./EIQ_I_tau=50.txt',/Get_Lun
      Point_lun, lunAo10, 0
      IQ_I=fltarr(ns)
      ReadF,lunAo10,IQ_I
      free_lun,lunAo10 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ;ns = 101
      Openr,lunAo10,'./EIQ_Q_tau=50.txt',/Get_Lun
      Point_lun, lunAo10, 0
      IQ_Q=fltarr(ns)
      ReadF,lunAo10,IQ_Q
      free_lun,lunAo10 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ;print, cosq1
 

      ratio=1.
      l=16 & xxss=l*(ratio) & yyss=l
      ;xoff=(LL-xxss)/2.,yoff=(3*LL/2.-yyss)/2.,
      !p.font = 0
      device,filename='./ChandraI_2.ps',xsize=xxss,ysize=yyss,bits_per_pixel=8,$
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
      ylow =0.
      yup =    1.5
      posup=[(1.-xlen)/2.,(1.-ylen)/2.,(1.+xlen)/2.,(1.+ylen)/2.] 
      plot,[xlow,xup],[ylow,yup],pos=[posup],/noerase,/nodata,/device,$
         xrange=[xlow,xup],yrange=[ylow,yup],/ynozero,/normal,xstyle=4+1,ystyle=4+1
 
      axis,xaxis=1,xticks=5,xminor=2,xtickname=replicate(' ',15)
      axis,xaxis=0,xticks=5,xminor=2,xrange=[xlow,xup],xstyle=1,$;font=-1,$;,xtickname=replicate(' ',6),$;,
      charsize=1,xtitle=textoidl('cos(\theta)');,xthick=tickth,color=colors

    axis,yaxis=0,ytitle=textoidl('Polarization Degree'),yticks=7,yminor=2,yrange=[ylow,yup],ystyle=1
      axis,yaxis=1,ytickname=replicate(' ',12),yticks=7,yminor=2
      dyy = 0.35
      y0 = 5
      for i=0,0 do begin
          ;oplot, vv(*), lv(*) , thick=2, color=blue, linestyle=6;,psym=-4 
          ;oplot, cosq1(*), q1(*) , thick=2, color=red, linestyle=6;,psym=-4 
          ;oplot, cosq1(*), qd1(*) , thick=2, color=green, linestyle=6;,psym=-4 
          ;oplot, cosq1(*), Qsp(*)/max(Isp)  , thick=2, color=red, linestyle=6;,psym=-4 
          ; oplot, cosq1(*), Usp(*)/max(Isp), thick=2, color=green, linestyle=6;,psym=-4
          ;oplot, cosq1(*), Isp(*)/max(Isp) * 1.26398 , thick=4, color=red, linestyle=6;,psym=-4  
          oplot, cosq1(*), IQ_I(*)/max(IQ_I) * 1.26398 , thick=4, color=green, linestyle=6;,psym=-4  
          oplot, vv, ChandI  , thick=2, color=blue, linestyle=6;,psym=-4  

          ;oplot, cosq1(*), Idsp(*)/max(Idsp)* 0.1 , thick=2, color=magenta, linestyle=6;,psym=-4  
      endfor 

      ;print, 3.e7 * (1e9)^0.25
      device,/close 

      free_lun,lunAo10
      set_plot,oldn

      end



   
