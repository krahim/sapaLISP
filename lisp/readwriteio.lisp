(defvar *ocean-wave-data*
  (let ((the-array (make-array 1024)))
    (with-open-file (in-stream "ocean-wave.dat"
                               :direction :input)
                    (dotimes (index (length the-array) the-array)
                      (setf (svref the-array index)
                            (float (read in-stream)))))))

(let ((output  *ocean-wave-data*))
  (with-open-file (out-stream "out22.dat" 
                              :direction :output
                              :if-exists :new-version
                              :if-does-not-exist :create)
                  (dotimes (index (length output) index)
                    (format out-stream "~D~%" (svref output index)))))


(defvar *of* (open "out22.dat" :direction :output))
(write 222 :stream *of*)
(close *of*)

(with-open-file (s "out223.dat"  :direction :output)
		(write 33333 s))