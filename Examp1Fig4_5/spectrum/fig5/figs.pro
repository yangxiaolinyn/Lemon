      pro figs

      oldn=!D.name & set_plot,'ps'

      ns = 500
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Openr,lunAo10,'./SARtobs=90.000pobs=0.000ne=1.0E+20.txt',/Get_Lun
      Point_lun,lunAo10,0
      Lv_Int=fltarr(ns)
      ReadF,lunAo10,Lv_Int
      free_lun,lunAo10

      Lv_Int_total = 0.
 
      vv = fltarr(ns)
      for i=0,ns-1 do begin
          vv(i) = 8.+(15.-8.)/ns*i
          Lv_Int_total = Lv_Int_total + Lv_Int(i)
      endfor
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      Openr,lunAo10,'./MCRtobs=90.000pobs=0.000ne=1.0E+20.txt',/Get_Lun
      Point_lun,lunAo10,0
      Lv_MC=fltarr(ns)
      ReadF,lunAo10,Lv_MC
      free_lun,lunAo10  

      Lv_MC_total = 0.
      for i=0,ns-1 do begin 
          Lv_MC_total = Lv_MC_total + Lv_MC(i)
      endfor
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ;total1=total(yy(*))

      ;print,max(yy),max(zz),total1/1.0^7
      ;print,alog10(lv)
      ;print,max(zz)
 

      ratio=1./0.67
      l=16 & xxss=l*(ratio) & yyss=l
      ;xoff=(LL-xxss)/2.,yoff=(3*LL/2.-yyss)/2.,
      !p.font = 0
      device,filename='vLv_01.ps',xsize=xxss,ysize=yyss,bits_per_pixel=8,$
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

      cs = 2.0

      xlen=0.7
      ylen=0.7
      xlow = 8.
      xup =  15. 
      ylow = -16
      yup = -1
      posup=[(1.-xlen)/2.,(1.-ylen)/2.,(1.+xlen)/2.,(1.+ylen)/2.] 
      !x.style = 0
      !y.style = 0
      plot,[xlow,xup],[ylow,yup],pos=[posup],/noerase,/nodata,/device,$
         xrange=[xlow,xup],yrange=[ylow,yup],/ynozero,/normal,xstyle=1+4,ystyle=1+4
 
      axis,xaxis=1,xticks=7,xminor=2,xtickname=replicate(' ',12), charsize=cs
      axis,xaxis=0,xticks=7,xminor=2,xrange=[xlow,xup],xstyle=1,$;font=-1,$;,xtickname=replicate(' ',6),$;,
      charsize=cs,xtitle=textoidl('Log(\nu)[Hz]');,xthick=tickth,color=colors

      ;axis,yaxis=0,ytitle=textoidl('Log(\nu L_\nu[erg s^{-1}])'),yticks=10,yminor=2,yrange=[ylow,yup],ystyle=1
      axis,yaxis=0,ytitle=textoidl('Log(\nu L_\nu[erg s^{-1}])'),yticks=5,$
          yminor=3,yrange=[ylow,yup],ystyle=1, charsize=cs
      axis,yaxis=1,ytickname=replicate(' ',20),yticks=5,yminor=3, charsize=cs
      dyy = 0.35
      y0 = 5
      for i=0,0 do begin
          oplot,vv(*),alog10(Lv_Int / Lv_Int_total),thick=2,color=black,linestyle=6;,psym=-4 
          oplot,vv(*),alog10(Lv_MC / Lv_MC_total)-0.02,thick=6,color=black,linestyle=1;,psym=-4 

          ;oplot,vv(*),alog10(Lv_Int2 / Lv_Int2_total),thick=2,color=black,linestyle=6;,psym=-4 
          ;oplot,vv(*),alog10(Lv_MC2 / Lv_MC2_total),thick=6,color=black,linestyle=1;,psym=-4 
          ;oplot,vv8(*),alog10(v_Lv2(*)/v_Lv2_total),thick=2,color=black,linestyle=6;,psym=-4
          ;oplot,vv8(*),alog10(vLv8(*)/Lv8_total),thick=6,color=black,linestyle=1;,psym=4
          ;oplot,vv(*)-3,vv(*)*1.5,thick=2,color=cyan,linestyle=6;,psym=-4
          ;oplot,vv(*)-1.6,vv(*)*1.,thick=2,color=red,linestyle=6;,psym=-4 
         ; oplot,vv(*),lv(*)/max(lv(*)),thick=2,color=blue,linestyle=6;,psym=-4
          ;oplot,tobs1(*),zz1(*),thick=2,color=blue,linestyle=6;,psym=-4   
          If i eq 3 then begin
              xyouts,-0.4,y0-i*dyy-0.02,textoidl('\vartheta_{*}=27^\circ'),charsize=1 
          endif  
      endfor

      cols = 243
      ;oplot,[-0.9,-0.1],[1.2,1.2],linestyle=0,color=cols,thick=2 
      
      ;oplot,[-0.,-0.],[1,6],linestyle=0,color=243,thick=2

      device,/close 

      free_lun,lunAo10
      set_plot,oldn

      end



  
