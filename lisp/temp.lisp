(defun filter-time-series-fft
       (time-series
        the-filter
        &key
        (start 0)
        (end (length time-series))
        (fft-size (* 4 (next-power-of-2
                        (length the-filter))))
        (verbose-p nil)
        (result (make-array (- (- end start) (1- (length the-filter)))
                            :initial-element 0.0)
                result-supplied-p))
  "given
   [1] time-series (required)
       ==> a vector containing a time series
           x_0, x_1, ..., x_{N-1}
   [2] the-filter (required)
       ==> a vector containing the filter coefficients
           g_0, g_1, ..., x_{K-1}
   [3] start (keyword; 0)
       ==> start index of time-series to be used
   [4] end (keyword; length of time-series)
       ==> 1 + end index of time-series to be used
   [5] fft-size (keyword; 4 * power of 2 that is ceiling for filter length)
       ==> size of fft's to be used (must be a power of 2)
   [6] verbose-p (keyword; nil)
       ==> if t, prints line after each block of data is processed;
           if nil, prints nothing
   [7] result (keyword; vector of appropriate length)
       <== vector to contain filtered time series
                 K-1
           y_t = SUM g_k x_{t+K-1-k},  t = 0, ..., N-K+1
                 k=0
 returns
   [1] result, a vector containing the filtered time series
   [2] the number of values in the filtered time series
---
Note: result can be the same as time-series"
  (assert (power-of-2 fft-size))
  (let* ((N-time-series (- end start))
         (N-filter (length the-filter))
         (N-filter-1 (1- N-filter))
         (fft-of-filter (make-array fft-size :initial-element 0.0))
         fft-of-data
         (N-output (1+ (- N-time-series N-filter)))
         (N-block (1+ (- fft-size N-filter)))
         (i-result -1))
    (if result-supplied-p
      (fill result 0.0 :end N-output))
    (copy-vector the-filter fft-of-filter)
    (fft! fft-of-filter)
    (do* ((k start (+ k N-block))
          (still-to-go N-output (- still-to-go N-block))
          (m (+ still-to-go N-filter-1) (+ still-to-go N-filter-1)))
         ((not (plusp still-to-go)) . nil)
         (format t "inside do* k=~D,  fft-size=~D m=~D still-to-go=~D~%" k  fft-size m still-to-go)
      (setf fft-of-data
            (if (>= m fft-size)
              (subseq time-series k (+ k fft-size))
              (concatenate 'array
                           (subseq time-series k (+ k m))
                           (make-array (- fft-size m)
                                       :initial-element 0.0))))
      (if verbose-p (format t "~&k = ~D, still-to-go = ~D ~%" k still-to-go))
      (fft! fft-of-data)
      (dotimes (i fft-size)
        (setf (aref fft-of-data i)
              (/ (conjugate (* (aref fft-of-data i) (aref fft-of-filter i)))
                 fft-size)))
      (fft! fft-of-data)
      (dotimes (i (min still-to-go N-block))
        ;;(format t "dotimes i=~D still-to-go=~D N-block=~D i-result=~D ~%"
                ;;i still-to-go N-block i-result)
        (setf (aref result (incf i-result))
              (realpart (aref fft-of-data (+ N-filter-1 i))))
        (format t "i-result=~D aref result i-result= ~G + N-filter-1 i = ~D ~%" i-result (aref result i-result) (+ N-filter-1 i) )))
    (values result N-output)))
