(defun get-N-dft (N-nonzero-freqs sample-size)
  (case N-nonzero-freqs
    (:half-next-power-of-2
     (next-power-of-2 sample-size))
    (:next-power-of-2
     (* 2 (next-power-of-2 sample-size)))
    (:twice-next-power-of-2
     (* 4 (next-power-of-2 sample-size)))
    (:Fourier sample-size)
    (otherwise
     (if (and (integerp N-nonzero-freqs)
              (power-of-2 N-nonzero-freqs)
              (>= (* 2 N-nonzero-freqs) sample-size))
       (* 2 N-nonzero-freqs)
       (error "N-nonzero-freqs (~A) is should be either
[a] a power of 2 greater than or equal to ~F (half the sample size)
or
[b] one of the following keywords: half-next-power-of-2,
    next-power-of-2, twice-next-power-of-2 or Fourier"
              N-nonzero-freqs
              (float (/ sample-size 2))
              )))))

(defun get-N-freqs
       (N-nonzero-freqs
        sample-size
        return-est-for-0-freq-p)
  (let ((temp (truncate (get-N-dft N-nonzero-freqs sample-size) 2)))
    (if return-est-for-0-freq-p
      (1+ temp)
      temp)))