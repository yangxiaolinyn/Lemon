      pro fig

      oldn=!D.name & set_plot,'ps'
  
      infile = './dataFF_xx.dat'
      Openr, lunAo10, infile, /Get_Lun
      nn = FILE_LINES(infile) 
      Point_lun,lunAo10,0
      fi=fltarr(nn)
      xi=fltarr(nn)
      pi=fltarr(nn)
      FOR j = 0, nn-1 DO BEGIN
          ReadF,lunAo10, f, x, p
          fi[j] = f
          xi[j] = x 
          pi[j] = p
      ENDFOR
      free_lun,lunAo10
 
 
       ;print,max(pa)
  
      ;xx1 = xx1

      ;total1=total(yy(*))

      ;print,max(yy),max(zz),total1/1.0^7
      ;print,zz
      ;print,max(zz)
 
      l=10
      xxss = 15 
      yyss = l
      ratio=xxss / yyss 
      ;xoff=(LL-xxss)/2.,yoff=(3*LL/2.-yyss)/2.,
      !p.font = 0
      device,filename='samples_1.ps',xsize=xxss,ysize=yyss,bits_per_pixel=8,$
      /color,xoff=(2-xxss)/2.0,yoff=(2-yyss)/2.,$
      set_font='Times-Roman';, /tt_font

      loadct,30
      RRR=bytscl(findgen(256))
      GGG=bytscl(findgen(256))
      BBB=bytscl(findgen(256))
      RRR[243:255]=[0,255,0  ,0  ,0  ,255,255,200,0  ,0  ,200,200,0]
      GGG[243:255]=[0,0  ,255,0  ,255,0  ,255,0  ,200,0  ,200,0  ,200]
      BBB[243:255]=[0,0  ,0  ,255,255,255,0  ,0  ,0  ,200,0  ,200,200]
      black=243
      red=244
      blue=246
      green=245
      ;245=green,243=black,244=red,251=yellow,246=blue,253=cyan,254=magenta,255=white
      TVLCT,RRR,GGG,BBB

      xlen=0.7
      ylen=0.7
      xlow = 0.
      xup = !Pi * 2.0
      ylow = 0.0
      yup =  1.2;max(www/www(nn-2)) + 0.3
      posup=[(1.-xlen)/2.,(1.-ylen)/2.,(1.+xlen)/2.,(1.+ylen)/2.] 
      plot,[xlow,xup],[ylow,yup],pos=[posup],/noerase,/nodata,/device,$
         xrange=[xlow,xup],yrange=[ylow,yup],/ynozero,/normal,xstyle=4+1,ystyle=4+1
 
      axis,xaxis=1,xticks=6,xminor=2,xtickname=replicate(' ',12)
      ;axis,xaxis=0,xticks=6,xminor=2,xrange=[xlow,xup],xstyle=1,$;font=-1,$;,xtickname=replicate(' ',6),$;,
      ;charsize=1,xtitle=textoidl('x');,xthick=tickth,color=colors
              labels = ['0', textoidl('\pi/2'), textoidl('\pi'), textoidl('3\pi/4'), textoidl('2\pi  ')]
              axis,xaxis=0,xticks=4,xminor=4,xrange=[xlow,xup],xstyle=1,$
                        xtickname=labels,$; xtickname=replicate( textoidl('\pi'),5),$;,
                      charsize=1,xtitle=textoidl('x');,xthick=tickth,color=colors

      axis,yaxis=0,ytitle=textoidl('p(x)'),yticks=6,yminor=2,yrange=[ylow,yup],ystyle=1
      axis,yaxis=1,ytickname=replicate(' ',12),yticks=6,yminor=2
      dyy = 0.35
      y0 = 5

      xw = RANDOMU(seed, 100000)
      h = histogram(xw, min = 0.0, binsize = 0.01)
      ;print, x
      ;hist_plot, h, color=blue, /normal, xstyle=3
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   hist = xi/max(xi)
   max_value = max(fi)
   min_value = min(fi)
nhist = n_elements(xi)
  binsize = (max_value - min_value) / nn
bins = lindgen(nn) * binsize + min_value
x = fltarr(2 * nhist)
x[2 * lindgen(nhist)] = bins
x[2 * lindgen(nhist) + 1] = bins

y = fltarr(2 * nhist)
y[2 * lindgen(nhist)] = hist
y[2 * lindgen(nhist) + 1] = hist
y = shift(y, 1)

;- Plot the histogram 
fi[0]=fi[1]
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      for i=0,0 do begin
          myHist = HISTOGRAM(xi(*), min=0., BINSIZE=1. )
          ;print,  myHist[9]  
          ;oplot,xx1(*),myHist,thick=3,color=green,linestyle=0;,psym=-4 
          oplot, xi , pi /max(pi), thick=1, color=red,linestyle=6;, _extra = extra_keywords 
          oplot, xi , fi /max(fi), thick=4, color=black,linestyle=1;, _extra = extra_keywords 
 
      endfor

      cols = 243
      ;oplot,[-0.9,-0.1],[1.2,1.2],linestyle=0,color=cols,thick=2
      ;oplot,[-0.1,-0.1],[1.2,3.4],linestyle=0,color=cols,thick=2
      ;oplot,[-0.9,-0.1],[3.4,3.4],linestyle=0,color=cols,thick=2
      ;oplot,[-0.9,-0.9],[1.2,3.4],linestyle=0,color=cols,thick=2 

      device,/close 

      free_lun,lunAo10
      set_plot,oldn

      end




PRO HIST_PLOT, DATA, MIN=MIN_VALUE, MAX=MAX_VALUE, $
  BINSIZE=BINSlZE, NORMALIZE=NORMALIZE, FILL=FILL, $
  _EXTRA=EXTRA_KEYWORDS
;-Check arguments
if (n_params() ne 1) then message, 'Usage: Hist_Plot, Data'
if (n_elements(data) eq 0) then message, 'Data is undefined'
;- Check keywords
if (n_elements(min_value) eq 0) then min_value = min(data)
if (n_elements(max_value) eq 0) then max_value = max(data)
if (n_elements(binsize) eq 0) then $
  binsize = (max_value - min_value) * 0.01
binsize = binsize > ((max_value - min_value) * 1.0e-5)

;Compute h i s t o g r a m  

hist = histogram(float(data), binsize=binsize, $
   min=min_value, max=max_value)
hist = data;[hist, 0L]
nhist = n_elements(hist)

;- Normalize histogram if required
if keyword_set(normalize) then $
   hist = hist / float(n_elements(data))

;-Compute bin values
bins = lindgen(nhist) * binsize + min_value
 
;-Create plot arrays
x = fltarr(2 * nhist)
x[2 * lindgen(nhist)] = bins
x[2 * lindgen(nhist) + 1] = bins

y = fltarr(2 * nhist)
y[2 * lindgen(nhist)] = hist
y[2 * lindgen(nhist) + 1] = hist
y = shift(y, 1)

;- Plot the histogram
oplot, x, data, _extra = extra_keywords
;- F i l l the histogram i f r e q u i r e d
if keyword_set(fill) then $
    polyfill, [x, x[O]], [y, y[O]], _extra = extra_keywords
END

 
