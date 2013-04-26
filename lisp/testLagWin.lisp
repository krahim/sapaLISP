(setf 20-pt-ts '#(71.0  63.0  70.0  88.0  99.0  90.0 110.0 135.0 128.0 154.0
                   156.0 141.0 131.0 132.0 141.0 104.0 136.0 146.0 124.0 129.0))
(setf sampling-time 0.25)
(setf the-acvs (acvs 20-pt-ts))
(setf Parzen-15-sdf-est (lag-window-spectral-estimate
                           the-acvs 
                           #'(lambda (tau)
                               (parzen-lag-window
                                tau 15))
                           :N-nonzero-freqs :Fourier
                           :sampling-time sampling-time))

(load (concatenate 'string path "nonparametric.lisp"))

 (lag-window-spectral-estimate
                        the-acvs 
                        #'(lambda (tau)
                            (parzen-lag-window
                             tau 5))
                        ;:max-lag 3
                        :N-nonzero-freqs :next-power-of-2
                        :sampling-time sampling-time)