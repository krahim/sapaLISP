## test 1

(center&taper-time-series
 #(71 63 70 88 99 90 110) 
 :data-taper #'supplied-data-taper! 
 :data-taper-parameters (car (dpss-tapers-inverse-iteration 
                              7 
                              1 
                              :taper-parameter 1))
 :recenter-after-tapering-p t
 :restore-power-option-p t)

##test 2
(setf 20-pt-ts #(71.0  63.0  70.0  88.0  99.0  90.0 110.0 135.0 128.0 154.0
                   156.0 141.0 131.0 132.0 141.0 104.0 136.0 146.0 124.0 129.0))
(setf N-ts (length 20-pt-ts))
(setf sigma^2 (sample-variance 20-pt-ts))
(setf sampling-time 0.25) 
(setf result  (multiple-value-bind (dpss-NW-4-tapers eigenvalues)
                       (dpss-tapers-tri-diag
                        N-ts 7
                        :compute-true-eigenvalues-p t)
                (list dpss-NW-4-tapers eigenvalues)))
(setf  dpss-NW-4-tapers (car result))
(setf eigenvalues (nth 1 result))

(setf mt4Res (multiple-value-bind (mt-4 freqs N-f list-of-eigenspectra)
                         (multitaper-spectral-estimate
                          20-pt-ts
                          dpss-NW-4-tapers
                          :N-nonzero-freqs :Fourier
                          :sdf-transformation nil
                          :sampling-time sampling-time)
         (list mt-4 freqs N-f list-of-eigenspectra)))

(setf mt-4 (nth 0 mt4Res))
(setf freqs (nth 1 mt4Res))
(setf N-f (nth 2 mt4Res))
(setf list-of-eigenspectra  (nth 3 mt4Res))

(print (/ (+ (* 2 (sum mt-4 :end (1- N-f)))                   
             (svref mt-4 (1- N-f)))
          (* N-ts sampling-time sigma^2)))

(setf amtRes (multiple-value-bind (amt-4 dof max-iterations)
                 (eigenspectra->adaptive-multitaper-spectral-estimate 
                  list-of-eigenspectra eigenvalues sigma^2
                  :sdf-transformation nil)
               (list amt-4 dof max-iterations)))
(setf amt-4 (car amtRes))
(setf dof (nth 1 amtRes))
(setf max-iterations (nth 2 amtRes))


(print (/ (+ (* 2 (sum amt-4 :end (1- N-f)))
             (svref amt-4 (1- N-f)))
          (* N-ts sampling-time sigma^2)))
        
(transform-a-sequence! #'convert-to-dB amt-4)
(dotimes (i N-f)
  (format t "~&~6,4F: ~8,4F  ~5,1F"
          (svref freqs i)
          (svref amt-4 i)
          (svref dof i)))
#test 3
(setf ocean-noise-data (vector 0.36201143 2.4737856 -1.018699 1.236099
                                -0.38258758 3.590585 -0.9524991 -1.6762866
                                1.3824866 -0.34308758 -0.1589 0.1412
                                2.6574097 2.2194982 -2.7033095 2.0669982
                                -1.3256867 -2.020322 -0.0753 0.4515114
                                2.4969096 -0.26581144 -1.026599 0.31911144
                                0.5307876 1.2934105 -2.6631095 -3.0335972
                                -2.5025096 0.74761146 -0.26221144 -2.178198
                                -2.940797 -1.5878867 -4.16842 1.5515105
                                -1.136499 -1.879798 0.46968758 1.6587105
                                2.198098 -2.3408096 0.133 -1.4207866
                                -0.8249991 5.1213956 -0.64661145 -1.9440981 
                                -0.5948876 1.041699 1.6077867 -2.4599857
                                -2.989097 0.2540876 5.5098066 2.9239972 
                                0.7211114 -2.033898 -2.2643094 3.2546084
                                2.7131858 -0.896799 -0.8838991 2.9770973 
                                -0.67578757 -0.43901145 1.5423105 0.836399
                                0.0143 0.4059876 1.055099 -0.5244876 
                                -3.155997 0.83619905 0.2799876 0.0655
                                -2.3649857 -1.0694991 -1.3471105 1.3075105 
                                0.6634876 -2.0554981 2.4471095 -2.187498
                                0.6463876 1.2897105 2.4921856 -4.5647836 
                                -1.5162104 1.3197105 2.2606857 -2.0319982
                                -1.4144866 2.0942981 -1.5411105 1.4940104 
                                -1.2181991 -0.1618 0.828699 -0.79309905
                                2.6937857 -0.3346876 -1.5918106 -0.17650001 
                                1.068099 0.1367 -4.7230077 -2.5424857
                                -1.2936105 -4.000296 0.47421142 0.5830876 
                                -2.6602857 -0.46411142 -0.3632876 4.088496
                                -0.6841115 -2.8142972 -0.88159907 2.801497 
                                2.081398 -3.6499846 1.860398 -3.6201086
                                -0.6762114 2.5205097 -1.3579867 -1.4434105))
(Fisher-g-statistic ocean-noise-data)
;==> 0.1336096017468551      ;;; Fisher's g statistic
;    0.10876132671244287     ;;; g_F from Equation (491c)
;    :reject                 ;;; reject because g > g_F
(setf the-periodogram (periodogram ocean-noise-data
                                           :center-data T
                                           :N-nonzero-freqs :Fourier
                                           :return-est-for-0-freq-p nil
                                           :sdf-transformation nil
                                           :return-frequencies-p nil))


;;test direct
(setf time-series (x+b! (ranorms 16) 100))
;;#(99.90197046279583 100.1684510266307 97.92618795789346 100.50258297137479 100.19581476207976 100.77978398645772 99.44853098557564 100.45277812019387 98.68525209726057 100.51845380875595 100.17956956720474 101.6677900717272 99.75392378696179 100.06090819374168 99.29891411303511 100.56951871379512)
(setf time-series #(99.90197046279583 100.1684510266307 97.92618795789346 100.50258297137479 100.19581476207976 100.77978398645772 99.44853098557564 100.45277812019387 98.68525209726057 100.51845380875595 100.17956956720474 101.6677900717272 99.75392378696179 100.06090819374168 99.29891411303511 100.56951871379512))
(setf sampling-time 0.25)
(direct-spectral-estimate
                       time-series
                       :data-taper #'dpss-data-taper!
                       :data-taper-parameters 4.0
                       :N-nonzero-freqs :Fourier
                       :sdf-transformation nil
                       :sampling-time sampling-time 
                       :return-acvs-p T)