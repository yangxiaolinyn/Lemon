      pro figs

      oldn=!D.name & set_plot,'ps'

      ns = 500
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Openr,lunAo10,'./MCRTe=2.044tau=0.1000.dat',/Get_Lun
      Point_lun,lunAo10,0
      lv=fltarr(ns)
      ReadF,lunAo10,lv
      free_lun,lunAo10
 
      vv = fltarr(ns)
      a1 = 11.
      a2 = 22.
      for i=0,500-1 do begin
          vv(i) = (a2-a1)/499*i + a1
      endfor
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
 

      ;total1=total(yy(*))

      ;print,max(yy),max(zz),total1/1.0^7
      ;print,alog10(lv)
      ;print,max(zz)
 

      ratio=1.
      l=16 & xxss=l*(ratio) & yyss=l
      ;xoff=(LL-xxss)/2.,yoff=(3*LL/2.-yyss)/2.,
      !p.font = 0
      device,filename='vLv_1.ps',xsize=xxss,ysize=yyss,bits_per_pixel=8,$
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
      xlow = 11
      xup = 22
      ylow = -6
      yup = 1
      posup=[(1.-xlen)/2.,(1.-ylen)/2.,(1.+xlen)/2.,(1.+ylen)/2.] 
      plot,[xlow,xup],[ylow,yup],pos=[posup],/noerase,/nodata,/device,$
         xrange=[xlow,xup],yrange=[ylow,yup],/ynozero,/normal,xstyle=4+1,ystyle=4+1
 
      axis,xaxis=1,xticks=11,xminor=2,xtickname=replicate(' ',12)
      axis,xaxis=0,xticks=11,xminor=2,xrange=[xlow,xup],xstyle=1,$;font=-1,$;,xtickname=replicate(' ',6),$;,
      charsize=1,xtitle=textoidl('Log(v)[H_z]');,xthick=tickth,color=colors

      axis,yaxis=0,ytitle=textoidl('Log(\nu L_\nu[erg s^{-1}])'),yticks=7,yminor=5,$
                   yrange=[ylow,yup],ystyle=1
      axis,yaxis=1,ytickname=replicate(' ',12),yticks=7,yminor=5,ystyle=1
      dyy = 0.35
      y0 = 5
      for i=0,0 do begin
           oplot,vv(*),alog10(lv(*) /max(lv) )+0.8,thick=3,color=black,linestyle=6;,psym=-4
           ;oplot,vv(*),alog10(lv3(*) /max(lv3) )+0.5,thick=3,color=black,linestyle=6;,psym=-4
          ;oplot,vv(*),alog10(lv2(*) /max(lv2)),thick=3,color=black,linestyle=6;,psym=-4
          ;oplot,vv(*),alog10(lv4(*) /max(lv4)),thick=3,color=red,linestyle=3;,psym=-4
          ;print, alog10(lv2(0:200) /max(lv2)), vv(*)
          ;oplot,tobs1(*),zz1(*),thick=2,color=blue,linestyle=6;,psym=-4
          ;oplot,tobs21(*),dtdt1(*),thick=2,color=black,linestyle=6;,psym=-4
          ;oplot,tt03(*),rr03(*),thick=2,color=magenta,linestyle=6;,psym=-4
          ;oplot,tt04(*),rr04(*),thick=2,color=cyan,linestyle=6;,psym=-4
          ;oplot,tt05(*),rr05(*),thick=2,color=cyan,linestyle=6;,psym=-4 
          ;oplot,[-0.8,-0.5],[y0-i*dyy,y0-i*dyy],linestyle=0,color=254-i,thick=2
          ;temp = textoidl('\vartheta_{*}=');+string(i*9)+textoidl('^\circ')
          If i eq 0 then begin
             ; xyouts,-0.4,y0-i*dyy-0.02,textoidl('\vartheta_{*}=0^\circ'),charsize=1 
          endif
          If i eq 1 then begin
              xyouts,-0.4,y0-i*dyy-0.02,textoidl('\vartheta_{*}=9^\circ'),charsize=1 
          endif 
      endfor

      cols = 243
      ;oplot,[-0.9,-0.1],[1.2,1.2],linestyle=0,color=cols,thick=2
      ;oplot,[-0.1,-0.1],[1.2,3.4],linestyle=0,color=cols,thick=2
      ;oplot,[-0.9,-0.1],[3.4,3.4],linestyle=0,color=cols,thick=2
      ;oplot,[-0.9,-0.9],[1.2,3.4],linestyle=0,color=cols,thick=2
      
      ;oplot,[-0.,-0.],[1,6],linestyle=0,color=243,thick=2

      device,/close 

      free_lun,lunAo10
      set_plot,oldn

      end



 
