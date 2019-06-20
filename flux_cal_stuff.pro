pro flux_cal_stuff

  .r /Users/james/idl/pro/remove.pro
  ;===== now compute the rough flux calibration for this target
  MOST_DIR = '~/Dropbox/research_projects/most-flares-kipping/'
  
  ; Bochanski M6 active template
  ;std = mrdfits('m6.active.ha.na.k.fits',0,hdr,/silent)
  std = mrdfits(MOST_DIR + 'm6.nactive.na.k.fits',0,hdr,/silent)
  std_flux = std[*,1]
  std_wave = sxpar(hdr, 'CRVAL1') + $
             findgen(n_elements(std_flux)) * sxpar(hdr, 'CD1_1')

  ; sorta...flux calibrated spectral of Prox Cen,
  ;    partial wavelength coverage of the MOST filter
  ; from DwarfArchives
  ;readcol, 'Gl_551.7328.txt', wave_cal, flux_cal


  ;--> need to get flux calibration spectrum instead
  ; from HST:
  ; http://archive.stsci.edu/proposal_search.php?id=6059&mission=hst
  cal2 = mrdfits(MOST_DIR + 'Y2WY0705T.fits',1,hdrcal2,/silent)
  cal1 = mrdfits(MOST_DIR + 'Y2WY0305T.fits',1,hdrcal1,/silent)

  wave_cal = [cal1.wavelength, $
              (cal2.wavelength)[where(cal2.wavelength gt max(cal1.wavelength))]]
  flux_cal = [cal1.flux, $
               (cal2.flux)[where(cal2.wavelength gt max(cal1.wavelength))]]

  remove,where(flux_cal lt 0), wave_cal, flux_cal

  ss = sort(wave_cal)
  wave_cal = wave_cal[ss]
  flux_cal = flux_cal[ss]


  ; now calibrate the template to this spectrum
  x1 = where(wave_cal ge 7300 and wave_cal le 7500)
  x2 = where(std_wave ge 7300 and std_wave le 7500)
  flux_convert = median(flux_cal[x1]) / median(std_flux[x2])


  ; the TESS filter curve
  readcol, 'tess-response-function-v1.0.csv',delim=',', w_f, trans_f,/silent
  wave_f = w_f * 10.


  ; - clip spectrum to filter wavelength range
  xok = where(std_wave ge min(wave_f) and $
              std_wave le max(wave_f))
  std_wave2 = std_wave[xok]
  std_flux2 = std_flux[xok] * flux_convert


  ;- resample filter on to spectrum wavelength grid
  trans_f2 = interpol(trans_f,wave_f,std_wave2)

  x50 = where(trans_f2 ge max(trans_f2)/2.)
  FWHM = max(std_wave2[x50]) - min(std_wave2[x50]) ; in Angstroms


  plx = 768.13 ; mas - from Lurie 2014
  dist = 1000. / plx * !pc

  print,'Distance (pc) = ',dist / !pc

  ;- convolve spectrum with filter
  ;; E_POINT = alog10(TSUM(std_wave2, std_flux2 * trans_f2 * std_wave2 * 1d-8) * $
  ;;                  (2. * !dpi * dist^2.) * FWHM)
  E_POINT = alog10(TSUM(std_wave2, std_flux2 * trans_f2) * $
                   (2. * !dpi * dist^2.))

  print,'E_POINT = ',e_point




return
end

