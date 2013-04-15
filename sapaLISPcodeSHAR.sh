#!/bin/sh
# This is a shell archive (produced by shar 3.50)
# To extract the files from this archive, save it to a file, remove
# everything above the "!/bin/sh" line above, and type "sh file_name".
#
# made 07/21/1993 16:54 UTC by dbp@vapor
# Source directory /home/vapor6/dbp/SAPA
#
# existing files will NOT be overwritten unless -c is specified
#
# This shar contains:
# length  mode       name
# ------ ---------- ------------------------------------------
#   6416 -rw-r--r-- acvs.lisp
#  38880 -rw-r--r-- basic-math.lisp
#  90436 -rw-r--r-- basic-statistics.lisp
#  18472 -rw-r--r-- dft-and-fft.lisp
#  73426 -rw-r--r-- examples.lisp
#  48477 -rw-r--r-- filtering.lisp
#   4605 -rw-r--r-- hacks.lisp
#   9544 -rw-r--r-- harmonic.lisp
#  10351 -rw-r--r-- index
#  30945 -rw-r--r-- matrix.lisp
#  90780 -rw-r--r-- multitaper.lisp
# 111652 -rw-r--r-- nonparametric.lisp
# 115199 -rw-r--r-- parametric.lisp
#  23412 -rw-r--r-- random.lisp
#   1822 -rw-r--r-- sapa-package.lisp
#  25572 -rw-r--r-- tapers.lisp
#  40336 -rw-r--r-- utilities.lisp
#
# ============= acvs.lisp ==============
if test -f 'acvs.lisp' -a X"$1" != X"-c"; then
	echo 'x - skipping acvs.lisp (File already exists)'
else
echo 'x - extracting acvs.lisp (Text)'
sed 's/^X//' << 'SHAR_EOF' > 'acvs.lisp' &&
;;;-*- Mode: LISP; Package: :SAPA; Syntax: COMMON-LISP -*-
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;
;  acvs.lisp
;
;  a collection of Lisp functions to compute the sample acvs and variogram ...
;  Note:  before compiling and loading acvs.lisp,
;         you should compile and load (in the order listed)
;            sapa-package.lisp, utilities.lisp,  basic-statistics.lisp and
;            dft-and-fft.lisp
;
;  6/3/93
;
;  SAPA, Version 1.0; Copyright 1993, Donald B. Percival, All Rights Reserved
;
;  Use and copying of this software and preparation of derivative works
;  based upon this software are permitted.  Any distribution of this
;  software or derivative works must comply with all applicable United
;  States export control laws.
; 
;  This software is made available AS IS, and no warranty -- about the
;  software, its performance, or its conformity to any
;  specification -- is given or implied. 
;
;  Comments about this software can be addressed to dbp@apl.washington.edu
;-------------------------------------------------------------------------------
;;; (compile-file "ccl:SAPA;acvs.lisp")
;;; (load "ccl:SAPA;acvs.fasl")
;-------------------------------------------------------------------------------
(if (not (find-package :SAPA))
X  (defpackage :SAPA (:USE :COMMON-LISP)))
X
(in-package :SAPA)
X
(export '(acvs
X          biased-acvs->unbiased-acvs
X          sample-variogram
X          ))
X
;-------------------------------------------------------------------------------
(defun acvs
X       (time-series
X        &key
X        (start 0)
X        (end (length time-series))
X        (center-data-p t)
X        (acs-p nil)
X        (result (make-array (- end start))))
X  "given
X   [1] time-series (required)
X       ==> a vector containing the time series
X   [2] start (keyword; 0)
X       ==> start index of time-series to be used
X   [3] end (keyword; length of time-series)
X       ==> 1 + end index of time-series to be used
X   [4] center-data-p (keyword; t)
X       ==> if t, subtract sample mean of time series
X           prior to computing the acvs;
X           if nil, do not subtract the sample mean
X   [5] acs-p (keyword; nil)
X       ==> return autocorrelation sequence
X           rather than autocovariance sequence
X   [6] result (keyword; vector of length n)
X       <== vector to hold acvs (or acs)
calculates autocovariance (or autocorrelation) sequence using fft's and
returns
X   [1] result, a vector with the sample acvs (or acs)
X   [2] sample variance of time series (with mean treated
X       as specified by center-data-p keyword)
---
Note: see Equations (191a) and (190b) of the SAPA book"
X  (let* ((N (- end start))
X         (ntot (next-power-of-2 (* 2 N)))
X         (N*ntot (float (* N ntot)))
X         (scratch (make-array ntot :initial-element 0.0))
X         sample-variance scale-factor)
X    #+mcl(declare (dynamic-extent scratch))
X    (copy-vector time-series scratch :start start :end end)
X    (if center-data-p
X      (let ((the-mean (sample-mean time-series :start start :end end)))
X        (dotimes (i N)
X          (decf (svref scratch i) the-mean))))
X    (fft! scratch)
X    (dotimes (i ntot)
X      (setf (svref scratch i) (expt (abs (svref scratch i)) 2)))
X    (fft! scratch)
X    (setf sample-variance (/ (realpart (svref scratch 0))
X                             N*ntot))
X    (setf scale-factor (if acs-p (* N*ntot sample-variance)
X                           N*ntot))
X    (dotimes (i N)
X      (setf (aref result i)
X            (/ (realpart (aref scratch i)) scale-factor)))
X    (values result sample-variance)))
X
;;; A test case for acvs is given along with the next function ...
X
;-------------------------------------------------------------------------------
(defun biased-acvs->unbiased-acvs
X       (biased-acvs)
X  "given a biased estimate of the acvs,
returns the corresponding unbiased estimate"
X  (let ((i -1)
X        (n (length biased-acvs)))
X    (map (type-of biased-acvs)
X         #'(lambda (x) (* x (/ n (- n (incf i)))))
X         biased-acvs)))
X        
#|
(let ((test (make-array 19 :initial-element 4.0)))
X  (multiple-value-bind (biased-acvs sample-var)
X                       (acvs test :center-data-p nil)
X    (print sample-var)
X    (let ((unbiased-acvs (biased-acvs->unbiased-acvs biased-acvs)))
X      (dotimes (i 19 (values))
X        (format t "~&~2D: ~10,7F  ~10,7F"
X                i
X                (aref biased-acvs i)
X                (aref unbiased-acvs i))))))
;==>
16.00000000000003 
X 0: 16.0000000  16.0000000
X 1: 15.1578947  16.0000000
X 2: 14.3157895  16.0000000
X 3: 13.4736842  16.0000000
X 4: 12.6315789  16.0000000
X 5: 11.7894737  16.0000000
X 6: 10.9473684  16.0000000
X 7: 10.1052632  16.0000000
X 8:  9.2631579  16.0000000
X 9:  8.4210526  16.0000000
10:  7.5789474  16.0000000
11:  6.7368421  16.0000000
12:  5.8947368  16.0000000
13:  5.0526316  16.0000000
14:  4.2105263  16.0000000
15:  3.3684211  16.0000000
16:  2.5263158  16.0000000
17:  1.6842105  16.0000000
18:  0.8421053  16.0000000
|#
X
;-------------------------------------------------------------------------------
(defun sample-variogram
X       (time-series
X        times)
X  "given:
X   [1] time-series (required)
X       ==> a vector containing the time series
X   [2] times (required)
X       ==> a vector containing the associated times
returns
X   [1] the sample variogram
X   [2] associated lags
---
Note: see Diggle's book"
X  (let* ((n (length time-series))
X         (m (/ (* n (1- n)) 2))
X         (lags (make-array m))
X         (variogram (make-array m))
X         (k 0)
X         (pairs '()))
X    (dotimes (i (1- n))
X      (dotimes (j (- n (1+ i)))
X        (setf (aref lags k)
X              (- (aref times (+ i j 1)) (aref times i)))
X        (setf (aref variogram k)
X              (* 0.5 (expt (- (aref time-series (+ i j 1)) (aref time-series i))
X                           2)))
X        (incf k)))
X    ;;; sort the variogram via size of lag ...
X    (dotimes (i m)
X      (push (list (aref lags i)
X                  (aref variogram i)) pairs))
X    (setf pairs (sort pairs #'< :key #'car))
X    (dotimes (i m (values variogram lags))
X      (setf (aref lags i) (car (elt pairs i))
X            (aref variogram i) (cadr (elt pairs i))))))
X
#|
(sample-variogram #(1 2 3 4 5) #(1 2 3 4 5))
;==> #(0.5 0.5 0.5 0.5 2.0 2.0 2.0 4.5 4.5 8.0)
;    #(1 1 1 1 2 2 2 3 3 4)
|#
SHAR_EOF
chmod 0644 acvs.lisp ||
echo 'restore of acvs.lisp failed'
Wc_c="`wc -c < 'acvs.lisp'`"
test 6416 -eq "$Wc_c" ||
	echo 'acvs.lisp: original size 6416, current size' "$Wc_c"
fi
# ============= basic-math.lisp ==============
if test -f 'basic-math.lisp' -a X"$1" != X"-c"; then
	echo 'x - skipping basic-math.lisp (File already exists)'
else
echo 'x - extracting basic-math.lisp (Text)'
sed 's/^X//' << 'SHAR_EOF' > 'basic-math.lisp' &&
;;;-*- Mode: LISP; Package: :SAPA; Syntax: COMMON-LISP -*-
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;
;  basic-math.lisp
;
;  a collection of Lisp functions for certain basic mathematical operations ...
;  Note:  before compiling and loading basic-math.lisp,
;         you should compile and load
;            sapa-package.lisp
;
;  6/3/93
;
;  SAPA, Version 1.0; Copyright 1993, Donald B. Percival, All Rights Reserved
;
;  Use and copying of this software and preparation of derivative works
;  based upon this software are permitted.  Any distribution of this
;  software or derivative works must comply with all applicable United
;  States export control laws.
; 
;  This software is made available AS IS, and no warranty -- about the
;  software, its performance, or its conformity to any
;  specification -- is given or implied. 
;
;  Comments about this software can be addressed to dbp@apl.washington.edu
;-------------------------------------------------------------------------------
;;; (compile-file "ccl:SAPA;basic-math.lisp")
;;; (load "ccl:SAPA;basic-math.fasl")
;-------------------------------------------------------------------------------
(if (not (find-package :SAPA))
X  (defpackage :SAPA (:USE :COMMON-LISP)))
X
(in-package :SAPA)
X
(export '(;;; functions related to the gamma function and its derivatives ...
X          log-of-gamma
X          factorial
X          digamma
X          trigamma
X
X          ;;; functions for dealing with polynomials ...
X          evaluate-polynomial
X          multiply-2-polynomials
X          multiply-polynomials
X          zeros-of-polynomial
X
X          ;;; functions for finding the root of a function ...
X          Newton-Raphson
X          bisection-with-Newton-Raphson
X          secant-method
X
X          ;;; functions for numerical integation ...
X          simple-numerical-integration
X          Gauss-Legendre-quadrature
X          ))
X
;-------------------------------------------------------------------------------
;;; used by the next set of functions ...
(defparameter +coefficents-for-log-of-gamma+
X  (list 76.18009172947146d0  -86.50532032941677d0      24.01409824083091d0
X        -1.231739572450155d0   0.1208650973866179d-2   -0.5395239384953d-5))
X
;-------------------------------------------------------------------------------
;;; used by the next set of functions ...
(defconstant +log-pi+ (log pi))
(defconstant +euler-constant+ 0.577215664901532860606512d0)
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  log-of-gamma
;;;                 factorial
;;;                 digamma
;;;                 trigamma
;;;  compute various quantities that are related to the gamma function
;;;  (log-of-gamma, factorial), its first derivative (digamma) and
;;;  its second derivative (trigamma).
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun log-of-gamma (xx)
X  "given xx (a real or complex valued number
whose real part is greater than 0),
returns the log (base e) of gamma(xx)
---
Note: based upon the discussion in Section 6.1,
Numerical Recipes, Second Edition"
X  (assert (plusp (realpart xx)))
X  (if (< (realpart xx) 1)
X    ;;; use reflection formula (6.1.4)
X    (let ((1-xx (- 1.0 xx)))
X      (- (+ +log-pi+ (log 1-xx))
X         (+ (log-of-gamma (1+ 1-xx)) (log (sin (* pi 1-xx))))))
X    ;;; since Re(xx) > 1, use approximation due to Lanczos
X    (let* ((x (1- xx))
X           (tmp-1 (+ x 5.5))
X           (tmp-2 (- (* (+ x 0.5) (log tmp-1)) tmp-1))
X           (ser 1.000000000190015d0))
X      (dolist (a-coefficient +coefficents-for-log-of-gamma+)
X        (incf x)
X        (incf ser (/ a-coefficient x)))
X      (+ tmp-2 (log (* 2.5066282746310005 ser))))))
X
;;; (log-of-gamma 6.0)  ;==> 4.787491742782046
;;; above should agree with
;;; (log (* 2 3 4 5))   ;==> 4.787491742782046
;;; and it does!
;;; (log-of-gamma 0.5)  ;==> 0.5723649429246563
;;; above should agree with
;;; (* 0.5 (log pi))    ;==> 0.5723649429247001
;;; and it does!
;;; (exp (log-of-gamma 1.755))  ;==> 0.9202092223790562
;;; above should agree with entry in row 2, column 2,
;;; page 270 of Abramowitz and Stegun, which has
;;;                                  0.9202092224
;;; and it does!
X
;-------------------------------------------------------------------------------
(defun factorial (k)
X  "given an integer k, returns k!"
X  (assert (and (integerp k) (not (minusp k))))
X  (cond
X   ((> k 15)  ; arbitrary choice ...
X    (exp (log-of-gamma (1+ k))))
X   ((or (= k 1) (= k 0))
X    1)
X   (t
X    (* k (factorial (1- k))))))
X
#|
(factorial 0)   ;==> 1
(factorial 1)   ;==> 1
(factorial 6)   ;==> 720
(factorial 25)  ;==> 1.5511210043610457E+25
|#
X
;-------------------------------------------------------------------------------
(defun digamma
X       (x
X        &key
X        (x-recursion 8.5)
X        (x-small 1.0e-5))
X  "given
X   [1] x (required)
X       ==> a positive number
X   [2] x-recursion (keyword; 8.5)
X       ==> for noninteger x,
X           if x-small < x < x-recursion,
X           recursive formula is used
X   [3] x-small (keyword; 1.0e-5)
X       ==> if x <= x-small,
X           small x approximation used (default )
returns
X   [1] value of digamma (psi, derivative of log(gamma)) function at x
---
Note: see Abramowitz and Stegun's equation 6.3.2;
expansion 6.3.18 plus recurrence 6.3.5;
the small x formula is Equation (5) from
Algorithm AS 103 -- Psi (Digamma) Function -- by J. M. Bernardo
in Applied Statistics, Vol. 25, No. 3, 1976, pp. 315--317;
note that better accuracy can be obtained by increasing
x-recursion --- 10.0 to 100.0 are probably useful values"
X  (assert (plusp x))
X  (cond
X   ((integerp x)
X    (let ((sum (- +euler-constant+)))
X      (dotimes (i (1- x) sum)
X        (incf sum (/ (1+ i))))))
X   ((<= x x-small)
X    (- (+ +euler-constant+ (/ x))))
X   ((< x x-recursion)
X    (- (digamma (1+ x) :x-recursion x-recursion) (/ x)))
X   (t
X    (let* ((x2 (* x x))
X           (x4 (* x2 x2)))
X      (- (log x)
X         (/ (* 2.0 x))
X         (/ (* 12.0 x2))
X         (/ (* -120.0 x4))
X         (/ (* 252.0 x2 x4)))))))
X
#|
(let ((x*1000 1500)
X      x)
X  (dotimes (i 15)
X    (setf x (float (/ x*1000 1000)))
X    (format t "~&~5,3F   ~12,10F" x (digamma x))
X    (incf x*1000 5)))
;==>
1.500   0.0364899738
1.505   0.0411536541
1.510   0.0457967894
1.515   0.0504195526
1.520   0.0550221144
1.525   0.0596046437
1.530   0.0641673072
1.535   0.0687102696
1.540   0.0732336935
1.545   0.0777377398
1.550   0.0822225674
1.555   0.0866883332
1.560   0.0911351924
1.565   0.0955632983
1.570   0.0999728023
;;; good agreement with 4th column of top of page 269, Abramowitz and Stegun
|#
X
;-------------------------------------------------------------------------------
(defun trigamma
X       (x
X        &key
X        (x-recursion 2.0))
X  "given
X   [1] x (required)
X       ==> a positive number
X   [2] x-recursion (keyword; 2.0)
X       ==> if  x < x-recursion,
X           recursive formula is used
returns
X   [1] value of trigamma function at x
---
Note: expansion 6.4.12 plus recurrence 6.4.6
of Abramowitz and Stegun;
better accuracy can be obtained by increasing
x-recursion ---10.0 to 30.0 is probably a useful upper limit"
X  (assert (plusp x))
X  (if (< x x-recursion)
X    (+ (trigamma (1+ x) :x-recursion x-recursion) (/ (* x x)))
X    (let* ((x2 (* x x))
X           (x3 (* x x2))
X           (x5 (* x2 x3))
X           (x7 (* x2 x5)))
X      (+ (/ x)
X         (/ (* 2.0 x2))
X         (/ (* 6.0 x3))
X         (/ (* -30.0 x5))
X         (/ (* 42.0 x7))
X         (/ (* -30.0 x2 x7))))))
#|
(let ((x*1000 1500)
X      x)
X  (dotimes (i 15)
X    (setf x (float (/ x*1000 1000)))
X    (format t "~&~5,3F   ~12,10F" x (trigamma x))
X    (incf x*1000 5)))
;==>
1.500   0.9348000492
1.505   0.9306736517
1.510   0.9265820503
1.515   0.9225248209
1.520   0.9185015462
1.525   0.9145118155
1.530   0.9105552242
1.535   0.9066313746
1.540   0.9027398746
1.545   0.8988803384
1.550   0.8950523864
1.555   0.8912556443
1.560   0.8874897441
1.565   0.8837543230
1.570   0.8800490239
;;; 5 place agreement with 5th column of top of page 269, Abramowitz and Stegun
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  evaluate-polynomial
;;;                 multiply-2-polynomials
;;;                 multiply-polynomials
;;;                 zeros-of-polynomial
;;;  evaluate, multiply and find the roots of polynomials whose coefficients
;;;  are represented by sequences of numbers.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun evaluate-polynomial
X       (polynomial
X        z
X        &key
X        (degree-of-polynomial (1- (length polynomial)))
X        (number-of-derivatives 0)
X        (complex-values (make-array (1+ number-of-derivatives)))
X        (bounds-to-be-computed nil)
X        (bounds (if bounds-to-be-computed
X                  (make-array (1+ number-of-derivatives)))))
X  "given
X   [1] polynomial (required)
X       ==> sequence of length n+1 with real or
X           complex-valued numbers giving the n+1
X           coefficients of a polynomial of nth order;
X           (elt polynomial i) = coefficient of z^(n-i)
X   [2] z (required)
X       ==> complex point at which the polynomial
X           to be evaluated
X   [3] degree-of-polynomial (keyword; (1- (length polynomial)))
X       ==> degree of polynomial
X   [4] number-of-derivatives (keyword; 0)
X       ==> number of derivatives to be evaluated
X   [5] complex-values (keyword; array of length (1+ number-of-derivatives))
X       <== sequence with values of polynomial
X       and derivatives
X   [6] bounds-to-be-computed (keyword; nil)
X       ==> if t, error bounds are computed
X       for values in complex-values
X   [7] bounds (keyword; array of length (1+ number-of-derivatives))
X       <== sequence with error bounds
returns
X   [1] the sequence complex-values containing the value
X       of the polynomial and its derivatives
X   [2] bounds, a sequence containing optional error bounds
---       
Note: this is a Lisp version of cmlib routine cpevl;
d1 was specified in cpevl as (expt 2 (- 1 (i1mach 11))),
where (i1mach 11) ==> THE NUMBER OF BASE-2 DIGITS (SINGLE PRECISION).
I have taken this to be equivalent to machine epsilon.
If this is in fact NOT the case, then the optional error bounds
might not be computed correctly."
X  (let ((d1 single-float-epsilon)
X        ci cim1)
X    (dotimes (j (1+ degree-of-polynomial))
X      (dotimes (i (min (1+ number-of-derivatives)
X                       (1+ (- degree-of-polynomial j))))
X        (setf ci (if (plusp j) (elt complex-values i) 0.0))
X        (setf cim1 (if (plusp i)
X                     (elt complex-values (1- i)) (elt polynomial j)))
X        (setf (elt complex-values i) (+ cim1 (* z ci)))
X        (if bounds-to-be-computed
X          (let* ((bi (if (plusp j) (elt bounds i) 0.0))
X                 (bim1 (if (plusp i) (elt bounds (1- i)) 0.0))
X                 (tf (+ bi (* (complex-of-absolutes ci)
X                             (+ (* 3.0 d1) (* 4.0 d1 d1)))))
X                 (r (realpart (* (complex-of-absolutes z)
X                                 (complex (realpart tf) (- (imagpart tf))))))
X                 (s (imagpart (* (complex-of-absolutes z) tf))))
X            (setf (elt bounds i)
X                  (if (plusp j)
X                    (* (1+ (* 8 d1))
X                       (+ bim1
X                          (* d1 (complex-of-absolutes cim1))
X                          (complex r s)))
X                    0.0)))))))
X  (values complex-values bounds))
X
#|
(evaluate-polynomial #(1 2 3) 1)
;==> #(6.0)
;    nil
(evaluate-polynomial #(1 2 3) 2)
;==> #(11.0)
;    nil
|#
X        
;-------------------------------------------------------------------------------
(defun multiply-2-polynomials
X       (polynomial-1
X        polynomial-2
X        &key               
X        (degree-of-1
X         (1- (length polynomial-1)))
X        (degree-of-2
X         (1- (length polynomial-2)))
X        (product-polynomial
X         (make-array (+ degree-of-1 degree-of-2 1))))
X  "given
X   [1] polynomial-1 (required)
X       ==> sequence of length n+1 with real or
X           complex-valued numbers giving the n+1
X           coefficients of a polynomial of nth order;
X           (elt polynomial-1 i) = coefficient of z^(n-i)
X   [2] polynomial-2 (required)
X       ==> another sequence representing another polynomial
X   [3] degree-of-1 (keyword; (1- (length polynomial-1)))
X       ==> degree of first polynomial
X   [4] degree-of-2 (keyword; (1- (length polynomial-2)))
X       ==> degree of second polynomial
X   [5] product-polynomial (keyword; array of length (+ degree-of-1 degree-of-2 1))
X       <== a sequence with product of two polynomials
returns
X   [1] product-polynomial, an sequence with coefficients
X       of polynomial given by the product of
X       polynomial-1 and polynomial-2
---       
Note: this routine works for a polynomial p represented either by
X        p(0) + p(1)*z + ... + p(degp)*z^degp
X                  or by
X        p(0)*z**degp + p(1)*z**(degp-1) + ... + p(degp)"
X  (let (k)
X    (dotimes (i (+ degree-of-1 degree-of-2 1) product-polynomial)
X      (setf (elt product-polynomial i) 0.0)
X      (dotimes (j (1+ degree-of-1))
X        (setf k (- i j))
X        (if (and (>= k 0) (<= k degree-of-2))
X          (incf (elt product-polynomial i)
X                (* (elt polynomial-1 j) (elt polynomial-2 k))))))))
X
#|
(multiply-2-polynomials #(1 2 3) #(4 3 2 1))
;==>#(4.0 11.0 20.0 14.0 8.0 3.0)
(multiply-2-polynomials #(4 3 2 1) #(1 2 3))
;==>#(4.0 11.0 20.0 14.0 8.0 3.0)
(multiply-2-polynomials #(4 3 2 1) #(2))
;==>#(8.0 6.0 4.0 2.0)
|#
X
;-------------------------------------------------------------------------------
(defun multiply-polynomials (&rest polys)
X  "given
X   [1] sequences representing any number of polynomials (required)
returns
X   [1] a sequence representing their product
---
Note: this function was written by Andrew G. Bruce"
X  (let ((n (length polys)))
X    (cond ((= n 0) polys)
X          ((= n 1) (first polys))
X          ((> n 1) (reduce #'multiply-2-polynomials polys)))))
X
#|
(multiply-2-polynomials #(1 2) #(4 3))
;==> #(4.0 11.0 6.0)
(multiply-2-polynomials #(4.0 11.0 6.0) #(1 2 3))
;==> #(4.0 19.0 40.0 45.0 18.0)
(multiply-polynomials #(1 2) #(4 3) #(1 2 3))
;==> #(4.0 19.0 40.0 45.0 18.0)
(multiply-polynomials #(1) #(4) #(2))
;==> #(8.0)
|#
X
;-------------------------------------------------------------------------------
(defun zeros-of-polynomial
X       (polynomial
X        &key
X        (degree-of-polynomial (1- (length polynomial)))
X        (the-roots (make-array degree-of-polynomial))
X        (maximum-number-of-iterations 25))
X  "given
X   [1] polynomial (required)
X       ==> sequence with coefficients of a polynomial;
X           (elt polynomial 0) must be nonzero
X   [2] degree-of-polynomial (keyword; (1- (length polynomial)))
X       ==> degree of polynomial
X   [3] the-roots (keyword; array of length degree-of-polynomial)
X       <== number of derivatives to be evaluated
X   [4] maximum-number-of-iterations (keyword; 25)
X       ==> maximum number of iterations
returns
X   [1] t or nil, where t indicates that all went well,
X       whereas nil indicates that convergence did not occur
X       at end of specificed number of iterations
X   [2] the-roots, a vector with the required roots
X       of the polynomial"
X  (cond
X   ;;; if this is a first degree, we're done
X   ((= degree-of-polynomial 1)
X    (setf (elt the-roots 0) (- (/ (elt polynomial 1) (elt polynomial 0))))
X    (values t the-roots))  ;t = all went well
X   ;;; If constant term is zero, we know one root and can get
X   ;;; the others by dealing with a polynomial of one less degree.
X   ((zerop (abs (elt polynomial degree-of-polynomial)))
X    (setf (elt the-roots (1- degree-of-polynomial)) 0.0)
X    (zeros-of-polynomial polynomial
X                         :degree-of-polynomial
X                         (1- degree-of-polynomial)
X                         :the-roots
X                         the-roots
X                         :maximum-number-of-iterations
X                         maximum-number-of-iterations))
X   ;;; Here we have to do some honest work: we first generate
X   ;;; some initial estimates for the roots and then the final estimates.
X   (t
X    (let ((initial-estimates (make-array degree-of-polynomial))
X          (scratch (make-array (* 2 (1+ degree-of-polynomial)))))
X      (initial-estimates-for-zeros-of-polynomial
X       polynomial
X       :degree-of-polynomial degree-of-polynomial
X       :initial-estimates initial-estimates
X       :scratch scratch)
X      ;;; If this function return nil, convergence was NOT found
X      (if (zeros-of-polynomial-with-initial-estimates
X           polynomial initial-estimates
X           :degree-of-polynomial degree-of-polynomial
X           :scratch scratch
X           :final-results the-roots
X           :maximum-number-of-iterations maximum-number-of-iterations)
X        (values t the-roots)    ;convergence occurred
X        (values nil the-roots)  ;did not converge
X        )))))
X
#|
(zeros-of-polynomial #(1 -2))
;==> t
;    #(2)
(zeros-of-polynomial #(1 -2 0))
;==> t
;    #(2 0.0)
(multiply-polynomials #(1 -2) #(1 -3) #(1 5))
;==> #(1.0 0.0 -19.0 30.0)
(zeros-of-polynomial #(1.0 0.0 -19.0 30.0))
;==>t
;   #(#c(2.9999999999999996 -2.2190775908547207E-22)
X      #c(-5.0                1.9700224365534028E-19)
X      #c(1.9999999999999998  2.219061434983382E-22))
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  Newton-Raphson
;;;                 bisection-with-Newton-Raphson
;;;                 secant-method
;;;  can be used to find the zero (root) of a function in a specfied interval.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun Newton-Raphson
X       (f
X        f-prime
X        x-left
X        x-right
X        &key
X        (accuracy (* 10.0 single-float-epsilon))
X        (maximum-number-of-iterations 20)
X        (prevent-bracket-jumping-p t))
X  "given
X   [1] f (required)
X       ==> a function with a single argument
X   [2] f-prime (required)
X       ==> another function with a single argument,
X           this one being the first derivative of f
X   [3] x-left (required)
X       ==> left-hand bracket for the desired root;
X           i.e., left-hand bracket <= desired root
X   [4] x-right (required)
X       ==> right-hand bracket for the desired root;
X           i.e., desired root <= right-hand bracket
X   [5] accuracy (keyword; (* 10.0 single-float-epsilon))
X       ==> desired relative accuracy for computed rood
X   [6] maximum-number-of-iterations (keyword; 20)
X       ==> maximum number of iterations
X   [7] prevent-bracket-jumping-p (keyword; t)
X       ==> if t, allows Newton-Raphson to continue
X           if it jumps out of the interval
X           [x-left, x-right];
X           if nil, jumping out of the interval
X           causes an error to be signaled
returns
X   [1] a root of f in [x-left, x-right];
X       i.e., a value x in that interval
X       such that f(x) = 0
X   [2] the number of iterations required
---
Note: this function is based loosely on rtnewt,
Section 9.4, Numerical Recipes, Second Edition"
X  (assert (< x-left x-right))
X  (let ((x (* 0.5 (+ x-left x-right)))
X        delta-x denom-for-accuracy-test)
X    (dotimes (j maximum-number-of-iterations
X                (if (not (cerror "returns solution so far"
X                                 "exceeding maximum number of iterations"))
X                  (values x maximum-number-of-iterations)))
X      (setf delta-x (/ (funcall f x)  (funcall f-prime x)))
X      (setf denom-for-accuracy-test (+ (abs x)
X                                       (abs (decf x delta-x))))
X      (cond
X       (prevent-bracket-jumping-p
X        (if (< x x-left) (setf x x-left))
X        (if (> x x-right) (setf x x-right))
X        (if (< (/ (abs delta-x) denom-for-accuracy-test) accuracy)
X          (return (values x (1+ j)))))
X       ((<= x-left x x-right)
X        (if (< (/ (abs delta-x) denom-for-accuracy-test) accuracy)
X          (return (values x (1+ j)))))
X       (t
X        (error "jumped out of brackets")
X        )))))
X  
#|
(defun exp-2 (x) (- (exp x) 2.0))
X
(Newton-Raphson
X #'exp-2
X #'exp        ;first derivative of exp-2
X 0.0
X 20.0)
X
;==> 0.6931471805599453   ;the root
;    15                   ;number of iterations needed
(exp-2 0.6931471805599453)
;==> 0.0
|#
X
;-------------------------------------------------------------------------------
(defun bisection-with-Newton-Raphson
X       (f
X        f-prime
X        x-left
X        x-right
X        &key
X        (accuracy (* 10.0 single-float-epsilon))
X        (maximum-number-of-iterations 100))
X  "given
X   [1] f (required)
X       ==> a function with a single argument
X   [2] f-prime (required)
X       ==> another function with a single argument,
X           this one being the first derivative of f
X   [3] x-left (required)
X       ==> left-hand bracket for the desired root;
X           i.e., left-hand bracket <= desired root
X   [4] x-right (required)
X       ==> right-hand bracket for the desired root;
X           i.e., desired root <= right-hand bracket
X   [5] accuracy (keyword; (* 10.0 single-float-epsilon))
X       ==> desired relative accuracy for computed rood
X   [6] maximum-number-of-iterations (keyword; 100)
X       ==> maximum number of iterations
returns
X   [1] a root of f in [x-left, x-right];
X       i.e., a value x in that interval
X       such that f(x) = 0
X   [2] the number of iterations required
---
Note: this function is based loosely on rtsafe,
Section 9.4, Numerical Recipes, Second Edition"
X  (let ((f-low (funcall f x-left))
X        (f-high (funcall f x-right))
X        df f-new
X        x-low x-high
X        rtsafe dxold dx temp)
X    (when (>= (* f-low f-high) 0.0)
X      (cond ((zerop f-low)
X             (values x-left 0))
X            ((zerop f-high)
X             (values x-right 0))
X            (t
X             (error "root not bracketed"))))
X    (cond ((< f-low 0.0)
X           (setf x-low x-left
X                 x-high x-right))
X          (t
X           (setf x-high x-left
X                 x-low x-right)))
X    (setf rtsafe (* 0.5 (+ x-left x-right))
X          dxold (abs (- x-right x-left))
X          dx dxold
X          f-new (funcall f rtsafe)
X          df (funcall f-prime rtsafe))
X    (dotimes (j maximum-number-of-iterations
X                (if (not (cerror "returns solution so far"
X                                 "exceeding maximum number of iterations"))
X                  (values rtsafe maximum-number-of-iterations)))
X      (cond ((or (>= (* (- (* (- rtsafe x-high) df) f-new)
X                        (- (* (- rtsafe x-low) df) f-new))
X                     0.0)
X                 (> (abs (* 2.0 f-new)) (abs (* dxold df))))
X             (setf dxold dx
X                   dx (* 0.5 (- x-high x-low))
X                   rtsafe (+ x-low dx))
X             (when (= x-low rtsafe) (return (values rtsafe (1+ j)))))
X            (t
X             (setf dxold dx
X                   dx (/ f-new df)
X                   temp rtsafe
X                   rtsafe (- rtsafe dx))
X             (when (= temp rtsafe) (return (values rtsafe (1+ j))))))
X      (when (< (abs dx) accuracy) (return (values rtsafe (1+ j))))
X      (setf f-new (funcall f rtsafe)
X            df (funcall f-prime rtsafe))
X      (if (< f-new 0.0)
X        (setf x-low rtsafe)
X        (setf x-high rtsafe)))))
X
#|
(defun exp-2 (x)
X  (- (exp x) 2.0))
X
(bisection-with-Newton-Raphson
X #'exp-2 #'exp 0.0 20.0
X :maximum-number-of-iterations 100)
;==> 0.6931471805599447
;    60
;;; good agreement with Newton-Raphson, but takes a lot longer!
|#
X
;-------------------------------------------------------------------------------
(defun secant-method
X       (f
X        x-left
X        x-right
X        &key
X        (accuracy (* 10.0 single-float-epsilon))
X        (maximum-number-of-iterations 50))
X  "given
X   [1] f (required)
X       ==> a function with a single argument
X   [2] x-left (required)
X       ==> left-hand bracket for the desired root;
X           i.e., left-hand bracket <= desired root
X   [3] x-right (required)
X       ==> right-hand bracket for the desired root;
X           i.e., desired root <= right-hand bracket
X   [4] accuracy (keyword; (* 10.0 single-float-epsilon))
X       ==> desired relative accuracy for computed rood
X   [5] maximum-number-of-iterations (keyword; 50)
X       ==> maximum number of iterations
returns
X   [1] a root of f in [x-left, x-right];
X       i.e., a value x in that interval
X       such that f(x) = 0
X   [2] the number of iterations required
---
Note: this function is based loosely on rtsec,
Section 9.2, Numerical Recipes, Second Edition"
X  (let ((f-left (funcall f x-left))
X        (f-right (funcall f x-right))
X        x-mid f-mid approx-f-mid approx-f-prime-mid delta-x x-new f-new
X        denom-for-accuracy-test)
X    (dotimes (j maximum-number-of-iterations
X                (if (not (cerror "returns solution so far"
X                                 "exceeding maximum number of iterations"))
X                  (values x-new maximum-number-of-iterations)))
X      (setf x-mid (* 0.5 (+ x-left x-right))
X            f-mid (funcall f x-mid)
X            approx-f-mid (* 0.5 (+ f-left f-right))
X            approx-f-prime-mid (/ (- f-right f-left)
X                                  (- x-right x-left))
X            delta-x (/ approx-f-mid  approx-f-prime-mid)
X            x-new (- x-mid delta-x)
X            f-new (funcall f x-new))
X      (setf denom-for-accuracy-test (+ (abs x-mid) (abs x-new)))
X      (if (or (zerop f-new)
X              (< (/ (abs delta-x) denom-for-accuracy-test) accuracy))
X        (return (values x-new (1+ j))))
X      (if (>= (* f-mid f-left) 0)
X        (setf x-left x-mid
X              f-left f-mid))
X      (if (>= (* f-mid f-right) 0)
X        (setf x-right x-mid
X              f-right f-mid))
X      (if (>= (* f-new f-left) 0)
X        (setf x-left x-new
X              f-left f-new)
X        (setf x-right x-new
X              f-right f-new))
X      )))
X  
#|
(defun exp-2 (x)
X  (- (exp x) 2.0))
X
(secant-method #'exp-2 0.0 20.0)
;==> 0.6931471805599453
;    13
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  simple-numerical-integration
;;;                 Gauss-Legendre-quadrature
;;;  are used for numerical integration.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun simple-numerical-integration
X       (f
X        a
X        b
X        &key
X        (accuracy 1.0e-6)
X        (maximum-number-of-iterations 20))
X  "given
X   [1] f (required)
X       ==> a function with a single argument
X   [2] a (required)
X       ==> left-hand limit for numerical integration
X   [3] b (required)
X       ==> right-hand limit for numerical integration
X           i.e., a < b
X   [4] accuracy (keyword; 1.0e-6)
X       ==> desired relative accuracy for computed rood
X   [5] maximum-number-of-iterations (keyword; 20)
X       ==> maximum number of iterations
returns
X   [1] the integral of f over the interval [a, b]
X   [2] the number of iterations required
---
Note: this function is based on qtrap,
Section 4.2, Numerical Recipes, Second Edition"
X  (let ((s-old nil)
X        s-new)
X    (dotimes (i maximum-number-of-iterations
X                (if (not (cerror "returns solution so far"
X                                 "exceeding maximum number of iterations"))
X                  (values s-old maximum-number-of-iterations)))
X      (setf s-new (trapezoidal-rule f a b s-old i))
X      (if (and s-old (< (abs (- s-new s-old))
X                        (* accuracy (abs s-old))))
X        (return (values s-new (1+ i))))
X      (setf s-old s-new))))
X
#|
;;; first test of trapezoidal-rule
(- (exp 1) 1) 
;==> 1.718281828459045
X
(let (s)
X  (dotimes (i 20)
X    (setf s (print (trapezoidal-rule #'exp 0.0 1.0 s i)))))
;==> 1.718281828459572  (great agreement, but takes a long time!)
X
;;; second test
(- (exp pi) (exp 1.0))
;==> 20.42241080432022
X
(let (s)
X  (dotimes (i 10)
X    (setf s (print (trapezoidal-rule #'exp 1.0 pi s i)))))
;==> 20.42244057984665 
X
;;; now test simple-numerical-integration on these two cases
(simple-numerical-integration #'exp 0.0 1.0)
;==> 1.7182823746860927  (good agreement)
;    10
X
(simple-numerical-integration #'exp 1.0 pi)
;==> 20.422412665291343
;    12
|#
X
;-------------------------------------------------------------------------------
(defun Gauss-Legendre-quadrature
X       (lower-limit
X        upper-limit
X        N)
X  "given
X   [1] lower-limit (required)
X       ==> lower limit of integration
X   [2] upper-limit (required)
X       ==> upper limit of integration
X   [3] N (required)
X       ==> number of points to be computed
X           in Gauss-Legendre quadrature
returns
X   [1] an N-dimensional vector of abscissas points
X   [2] an N-dimensional vector of weights
---
Note: this function is based on gauleg,
Section 4.5, Numerical Recipes, Second Edition"
X  (let ((abscissas (make-array N))
X        (weights (make-array N))
X        (eps 3.0d-14)
X        (m (truncate (+ N 1) 2))
X        (N+half (+ N 0.5))
X        (xm (* 0.5 (+ upper-limit lower-limit)))
X        (xl (* 0.5 (- upper-limit lower-limit)))
X        z pp)
X    ;;; Loop over desired roots.
X    (dotimes (i m (values abscissas weights))
X      (setf z (cos (* pi (/ (+ i 0.75) N+half))))
X      :label-1
X      (do ((p1 1.0 1.0)
X           (p2 0.0 0.0)
X           (z1 nil)
X            p3)
X          ;;; since z1 is initially nil,
X          ;;; we cannot exit before the first iteration
X          ((and z1 (<= (abs (- z z1)) eps)))
X        (dotimes (jm1 n)
X          (setf p3 p2
X                p2 p1
X                p1 (/ (- (* (1+ (* 2.0 jm1)) z p2) (* jm1 p3))
X                      (1+ jm1))))
X        (setf pp (* n (/ (- (* z p1) p2) (- (* z z) 1.0)))
X              z1 z
X              z (- z1 (/ p1 pp))))
X    (setf (aref abscissas i) (- xm (* xl z))
X          (aref abscissas (- n (1+ i))) (+ xm (* xl z))
X          (aref weights i) (* 2.0 (/ xl (* (- 1.0 (* z z)) pp pp)))
X          (aref weights (- n (1+ i))) (aref weights i)))))
X
#|
(multiple-value-bind (abscissas weights)
X                     (gauss-legendre-quadrature -1.0 1.0 32)
X  (dotimes (i (length abscissas))
X    (format t "~& ~2D  ~13,10F   ~13,10F"
X            (1+ i)
X            (aref abscissas i)
X            (aref weights i))))
X
X  1  -0.9972638618    0.0070186100
X  2  -0.9856115115    0.0162743947
X  3  -0.9647622556    0.0253920653
X  4  -0.9349060759    0.0342738629
X  5  -0.8963211558    0.0428358980
X  6  -0.8493676137    0.0509980593
X  7  -0.7944837960    0.0586840935
X  8  -0.7321821187    0.0658222228
X  9  -0.6630442669    0.0723457941
X 10  -0.5877157572    0.0781938958
X 11  -0.5068999089    0.0833119242
X 12  -0.4213512761    0.0876520930
X 13  -0.3318686023    0.0911738787
X 14  -0.2392873623    0.0938443991
X 15  -0.1444719616    0.0956387201
X 16  -0.0483076657    0.0965400885
X 17   0.0483076657    0.0965400885
X 18   0.1444719616    0.0956387201
X 19   0.2392873623    0.0938443991
X 20   0.3318686023    0.0911738787
X 21   0.4213512761    0.0876520930
X 22   0.5068999089    0.0833119242
X 23   0.5877157572    0.0781938958
X 24   0.6630442669    0.0723457941
X 25   0.7321821187    0.0658222228
X 26   0.7944837960    0.0586840935
X 27   0.8493676137    0.0509980593
X 28   0.8963211558    0.0428358980
X 29   0.9349060759    0.0342738629
X 30   0.9647622556    0.0253920653
X 31   0.9856115115    0.0162743947
X 32   0.9972638618    0.0070186100
X
;;; excellent agreement with top of page 917 of Abramowitz and Stegun
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  Everything below here consists of internal symbols in the SAPA package
;;;  and should be regarded as "dirty laundry" ...
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun zeros-of-polynomial-with-initial-estimates
X       (polynomial initial-estimates
X        &key
X        (degree-of-polynomial
X         (1- (length polynomial)))
X        (scratch (make-array (1+ degree-of-polynomial)))
X        (final-results initial-estimates)
X        (maximum-number-of-iterations 25))
X  (if (< (length initial-estimates) degree-of-polynomial)
X    (error (format nil "too few initial estimates (~D<~D)"
X                   (length initial-estimates) degree-of-polynomial)))
X  (cond
X   ((not (eq final-results initial-estimates))
X    (if (< (length final-results) degree-of-polynomial)
X      (error (format nil "not enough space for results (~D<~D)"
X                     (length final-results) degree-of-polynomial)))
X    (dotimes (i degree-of-polynomial)
X      (setf (elt final-results i) (elt initial-estimates i)))))
X  (let ((value-of-polynomial (make-array 1))
X        (error-bound-for-value (make-array 1))
X        (number-of-roots-found 0)
X        temp)
X    (dotimes (number-of-iterations (* maximum-number-of-iterations
X                                      degree-of-polynomial))
X      (dotimes (i degree-of-polynomial)
X        (cond
X         ((or (zerop number-of-iterations) (plusp (abs (elt scratch i))))
X          ;;; evaluate polynomial at point (elt final-results i)
X          ;;; and stuff result into (elt value-of-polynomial 0)
X          (evaluate-polynomial polynomial (elt final-results i)
X                               :degree-of-polynomial degree-of-polynomial
X                               :complex-values value-of-polynomial
X                               :bounds-to-be-computed t
X                               :bounds error-bound-for-value)
X          (cond
X           ((> (+ (abs (realpart (elt value-of-polynomial 0)))
X                  (abs (imagpart (elt value-of-polynomial 0))))
X               (+ (abs (realpart (elt error-bound-for-value 0)))
X                  (abs (imagpart (elt error-bound-for-value 0)))))
X            (setf temp (elt polynomial 0))
X            (dotimes (j degree-of-polynomial)
X              (if (not (= j i))
X                (setf temp (* temp (- (elt final-results i)
X                                      (elt final-results j))))))
X            (setf (elt scratch i) (/ (elt value-of-polynomial 0) temp)))
X           (t
X            (setf (elt scratch i) 0.0)
X            (incf number-of-roots-found))))))
X      (dotimes (i degree-of-polynomial)
X        (decf (elt final-results i) (elt scratch i)))
X      ;;; if exit is done here, the routine will return t;
X      ;;; otherwise, the routine will quit when the maximum number of
X      ;;; interations have been made and return nil.
X      (if (= number-of-roots-found degree-of-polynomial) (return t)))))
X
;-------------------------------------------------------------------------------
(defun initial-estimates-for-zeros-of-polynomial
X       (polynomial
X        &key
X        (degree-of-polynomial (1- (length polynomial)))
X        (initial-estimates (make-array degree-of-polynomial))
X        (scratch (make-array (* 2 (1+ degree-of-polynomial)))))
X  (let ((imax (1+ degree-of-polynomial))
X        (pn (make-array 1))
X        (scratch-offset
X         (make-array (+ degree-of-polynomial 2)
X                     :displaced-to scratch
X                     :displaced-index-offset degree-of-polynomial))
X        (temp (- (/ (elt polynomial 1)
X                    (* (elt polynomial 0) degree-of-polynomial))))
X        x u v)
X    (evaluate-polynomial polynomial temp
X                         :degree-of-polynomial degree-of-polynomial
X                         :number-of-derivatives degree-of-polynomial
X                         :complex-values scratch)
X    (setf (elt scratch degree-of-polynomial)
X          (abs (elt scratch degree-of-polynomial)))
X    (dotimes (i degree-of-polynomial)
X      (setf (elt scratch (+ degree-of-polynomial i 1))
X            (- (abs (elt scratch (- degree-of-polynomial i 1)))))
X      (if (< (realpart (elt scratch (+ degree-of-polynomial i 1)))
X             (realpart (elt scratch imax)))
X        (setf imax (+ degree-of-polynomial i 1))))
X    (setf x (expt (- (/ (realpart (elt scratch imax))
X                        (realpart (elt scratch degree-of-polynomial))))
X                  (/ (float (- imax degree-of-polynomial)))))
X    (loop
X      (setf x (* 2.0 x))
X      (evaluate-polynomial scratch-offset x
X                           :degree-of-polynomial degree-of-polynomial
X                           :complex-values pn)
X      (if (>= (realpart (elt pn 0)) 0.0) (return)))
X    (setf u (* 0.5 x))
X    (setf v x)
X    (loop
X      (setf x (* 0.5 (+ u v)))
X      (evaluate-polynomial scratch-offset x
X                           :degree-of-polynomial degree-of-polynomial
X                           :complex-values pn)
X      (if (plusp (realpart (elt pn 0)))
X        (setf v x)
X        (setf u x))
X      (if (<= (- v u) (* 0.001 (1+ v))) (return)))
X    (dotimes (i degree-of-polynomial)
X      (setf u (* (/ pi degree-of-polynomial) (+ 0.5 (* 2.0 i))))
X      (setf (elt initial-estimates i)
X            (+ temp (* (max x (* 0.001 (abs temp)))
X                       (complex (cos u) (sin u))))))))
X
;-------------------------------------------------------------------------------
(defun complex-of-absolutes (z)
X  (complex (abs (realpart z)) (abs (imagpart z))))
X
;-------------------------------------------------------------------------------
;;; trapezoidal-rule is based upon trapzd, Section 4.2, Numerical Recipes,
;;; Second Edition. It integrates f from a to b;
;;; on the first pass, should set n = 0, in which case s is ignored;
;;; for n = 1, 2, ..., calculates refined estimate of integral using s
;;; returned at previous stage; i.e.,
;;; at each stage n=m, function returns updated value of s,
;;; which should be used to call the function at stage n=m+1
(defun trapezoidal-rule (f a b s n)
X  (cond
X   ((zerop n)
X    (* 0.5 (- b a) (+ (funcall f a)
X                      (funcall f b))))
X   (t
X    (let* ((it (expt 2 (- n 1)))
X           (del (/ (- b a) it))
X           (x (+ a (* 0.5 del)))
X           (sum 0.0))
X      (dotimes (j it (* 0.5 (+ s (/ (* sum (- b a)) it))))
X        (incf sum (funcall f x))
X        (incf x del))))))
X
SHAR_EOF
chmod 0644 basic-math.lisp ||
echo 'restore of basic-math.lisp failed'
Wc_c="`wc -c < 'basic-math.lisp'`"
test 38880 -eq "$Wc_c" ||
	echo 'basic-math.lisp: original size 38880, current size' "$Wc_c"
fi
# ============= basic-statistics.lisp ==============
if test -f 'basic-statistics.lisp' -a X"$1" != X"-c"; then
	echo 'x - skipping basic-statistics.lisp (File already exists)'
else
echo 'x - extracting basic-statistics.lisp (Text)'
sed 's/^X//' << 'SHAR_EOF' > 'basic-statistics.lisp' &&
;;;-*- Mode: LISP; Package: :SAPA; Syntax: COMMON-LISP -*-
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;
;  basic-statistics.lisp
;
;  a collection of Lisp functions for certain basic statistical operations ...
;  Note:  before compiling and loading matrix.lisp,
;         you should compile and load (in the order listed)
;            sapa-package.lisp, utilities.lisp, basic-math.lisp and matrix.lisp
;
;  6/3/93
;
;  SAPA, Version 1.0; Copyright 1993, Donald B. Percival, All Rights Reserved
;
;  Use and copying of this software and preparation of derivative works
;  based upon this software are permitted.  Any distribution of this
;  software or derivative works must comply with all applicable United
;  States export control laws.
; 
;  This software is made available AS IS, and no warranty -- about the
;  software, its performance, or its conformity to any
;  specification -- is given or implied. 
;
;  Comments about this software can be addressed to dbp@apl.washington.edu
;-------------------------------------------------------------------------------
;;; (compile-file "ccl:SAPA;basic-statistics.lisp")
;;; (load "ccl:SAPA;basic-statistics.fasl")
;-------------------------------------------------------------------------------
(if (not (find-package :SAPA))
X  (defpackage :SAPA (:USE :COMMON-LISP)))
X
(in-package :SAPA)
X
(export '(;;; functions to compute various summary statistics ...
X          sample-mean
X          weighted-mean
X          sample-variance
X          sample-mean-and-variance
X          sample-skewness-and-kurtosis
X          sample-median
X          sample-correlation-coefficient
X          histogram
X          box-plot
X          symmetry-plot
X
X          ;;; quantiles ...
X          quantile-of-ordered-seq
X          quantile-plot
X          interquartile-range
X          quantile-of-normal-distribution
X          quantile-of-gamma-distribution
X          quantile-of-exponential-distribution
X          quantile-of-chi-square-2-distribution
X          quantile-of-chi-square-distribution
X          q-q-plot
X
X          ;;; the normal (Gaussian) distribution ...
X          standard-normal-pdf
X          tail-area-of-normal-distribution
X
X          ;;; robust statistics ...
X          median-absolute-deviation
X          Thomson-weight-function
X          Huber-weight-function
X          m-location-estimate
X
X          ;;; least squares calculations ...
X          ordinary-least-squares-Cholesky
X          ordinary-least-squares-Q-R
X          weighted-linear-least-squares
X          Durbin-Watson-test-statistic
X          predicted-y-at-x
X          var-predicted-y-at-x
X          var-predicted-mean-at-x
X
X          ;;; probability density function estimation ...
X          kernel-pdf-estimation
X          window-width-from-Silverman
X          ))
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  sample-mean
;;;                 weighted-mean
;;;                 sample-variance
;;;                 sample-mean-and-variance
;;;                 sample-skewness-and-kurtosis
;;;                 sample-median
;;;                 sample-correlation-coefficient
;;;                 histogram
;;;                 box-plot
;;;                 symmetry-plot
;;;  all take one or more sequences of numbers as input
;;;  and return the sample statistic(s) indicated by the function names.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun sample-mean
X       (the-seq
X        &key
X        (start 0)
X        (end (length the-seq)))
X  "given
X   [1] the-seq (required)
X       ==> a sequence of real or complex-valued numbers
X   [2] start (keyword; 0)
X       ==> start index of sequence to be used
X   [3] end (keyword; length of the-seq)
X       ==> 1 + end index of sequence to be used
returns
X   [1] sample mean of the specified numbers
X       in the-seq"
X  (/ (reduce #'+ the-seq :start start :end end) (- end start)))
X
#|
(sample-mean #(1 2 3 4 5 6 7 8 9 10))                    ;==> 11/2
(sample-mean #(1.0 2 3 4 5 6 7 8 9 10.0))                ;==> 5.5
(sample-mean #(1.0 2 3 #C(4 10) 5 6 7 8 9 10.0))         ;==> #c(5.5 1.0)
(sample-mean #(1 2 3 4 5 6 7 8 9 10) :start 5)           ;==> 8
(sample-mean #(1 2 3 4 5 6 7 8 9 10) :end 5)             ;==> 3
(sample-mean #(1 2 3 4 5 6 7 8 9 10) :start 0 :end 10)   ;==> 11/2
(sample-mean #(1 2 3 4 5 6 7 8 9 10) :start 1 :end 6)    ;==> 4
|#
X
;-------------------------------------------------------------------------------
(defun weighted-mean
X       (the-seq
X        &key
X        (weights nil)
X        (start 0)
X        (end (length the-seq)))
X  "given
X   [1] the-seq (required)
X       ==> a sequence of real or complex-valued numbers
X   [2] weights (keyword; nil)
X       ==> a sequence with (- end start) values
X           to be used as weights; if set to nil,
X           equal weights for all numbers
X           is assumed
X   [3] start (keyword; 0)
X       ==> start index of sequence to be used
X   [4] end (keyword; length of the-seq)
X       ==> 1 + end index of sequence to be used
returns
X   [1] the weighted average of the specified numbers
X       in the-seq"
X  (if weights
X    (/ (let ((weighted-sum 0)
X             (j (1- start)))
X         (dotimes (i (- end start) weighted-sum)
X           (incf weighted-sum (* (elt weights i)
X                                 (elt the-seq (incf j))))))
X       (sum weights :end (- end start)))
X    (sample-mean the-seq :start start :end end)))
X
#|
(weighted-mean '(1 2 3 4 5 6 7 8 9 10))                    ;==> 11/2
(weighted-mean '(1 2 3 4 5 6 7 8 9 10) :start 5)           ;==> 8
(weighted-mean '(1 2 3 4 5 6 7 8 9 10) :end 5)             ;==> 3
(weighted-mean '(1 2 3 4 5 6 7 8 9 10) :start 0 :end 10)   ;==> 11/2
(weighted-mean '(1 2 3 4 5 6 7 8 9 10) :start 1 :end 6)    ;==> 4
(weighted-mean '(1 2 3 4 5 6 7 8 9 10)
X               :weights '(0 1 0 1 0 1 0 1 0 1))            ;==> 6
(weighted-mean '(1 2 3 4 5 6 7 8 9 10)
X               :weights '(0 1 0 1 0 1 0 1 0 1)
X               :start 5)                                   ;==> 8
(weighted-mean '(1 2 3 4 5 6 7 8 9 10)
X               :weights '(0 1 0 1 0 1 0 1 0 1)
X               :end 5)                                     ;==> 3
(weighted-mean '(1 2 3 4 5 6 7 8 9 10)
X               :weights '(0 1 0 1 0 1 0 1 0 1)
X               :start 0 :end 10)                           ;==> 6
(weighted-mean '(1 2 3 4 5 6 7 8 9 10)
X               :weights '(0 1 0 1 0 1 0 1 0 1)
X               :start 1 :end 6)                            ;==> 4
|#
X
;-------------------------------------------------------------------------------
(defun sample-variance
X       (the-seq
X        &key
X        (start 0)
X        (end (length the-seq)))
X  "given
X   [1] the-seq (required)
X       ==> a sequence of real or complex-valued numbers
X   [2] start (keyword; 0)
X       ==> start index of sequence to be used
X   [3] end (keyword; length of the-seq)
X       ==> 1 + end index of sequence to be used
returns
X   [1] sample variance of the specified numbers
X       in the-seq"
X  (elt (multiple-value-list 
X        (sample-mean-and-variance the-seq :start start :end end))
X       1))
X
#|
(sample-variance #(1 2 3 4 5 6 7 8 9 10))                    ;==> 33/4
(sample-variance #(1.0 2 3 4 5 6 7 8 9 10.0))                ;==> 8.25
(sample-variance #(1.0 2 3 #C(4 10) 5 6 7 8 9 10.0))         ;==> 17.25
(sample-variance #(1 2 3 4 5 6 7 8 9 10) :start 5)           ;==> 2
(sample-variance #(1 2 3 4 5 6 7 8 9 10) :end 5)             ;==> 2
(sample-variance #(1 2 3 4 5 6 7 8 9 10) :start 0 :end 10)   ;==> 33/4
(sample-variance #(1 2 3 4 5 6 7 8 9 10) :start 1 :end 6)    ;==> 2
(sample-variance #(1 2 3 4 5 6 7 8 9 10) :start 2 :end 7)    ;==> 2
|#
X
;-------------------------------------------------------------------------------
(defun sample-mean-and-variance
X       (the-seq
X        &key
X        (start 0)
X        (end (length the-seq)))
X  "given
X   [1] the-seq (required)
X       ==> a sequence of real or complex-valued numbers
X   [2] start (keyword; 0)
X       ==> start index of sequence to be used
X   [3] end (keyword; length of the-seq)
X       ==> 1 + end index of sequence to be used
returns
X   [1] sample mean of the specified numbers
X       in the-seq
X   [2] sample variance of the specified numbers
X       in the-seq"
X  (let ((the-mean (sample-mean the-seq :start start :end end)))
X    (values the-mean (/ #-allegro
X                        (reduce #'+ the-seq
X                                :start start
X                                :end end
X                                :key #'(lambda (x)
X                                         (expt (abs (- x the-mean)) 2)))
X                        #+allegro
X                        (let ((SS 0.0)
X                              (j start))
X                          (dotimes (i (- end start) SS)
X                            (incf SS (expt (abs (- (elt the-seq j) the-mean)) 2))
X                            (incf j)))
X                        (- end start)))))
X
#|
(sample-mean-and-variance #(1 2 3 4 5 6 7 8 9 10))
;==> 11/2
;    33/4
(sample-mean-and-variance #(1.0 2 3 4 5 6 7 8 9 10.0))
;==> 5.5
;    8.25
(sample-mean-and-variance #(1.0 2 3 #C(4 10) 5 6 7 8 9 10.0))
;==> #c(5.5 1.0)
;    17.25
(sample-mean-and-variance #(1 2 3 4 5 6 7 8 9 10) :start 5)
;==> 8
;    2
(sample-mean-and-variance #(1 2 3 4 5 6 7 8 9 10) :end 5)
;==> 3
;    2
(sample-mean-and-variance #(1 2 3 4 5 6 7 8 9 10) :start 0 :end 10)
;==> 11/2
;    33/4
(sample-mean-and-variance #(1 2 3 4 5 6 7 8 9 10) :start 1 :end 6)
;==> 4
;    2
(sample-mean-and-variance #(1 2 3 4 5 6 7 8 9 10) :start 2 :end 7)
;==> 5
;    2
|#
X
;-------------------------------------------------------------------------------
(defun sample-skewness-and-kurtosis
X       (the-seq
X        &key
X        (start 0)
X        (end (length the-seq)))
X  "given
X   [1] the-seq (required)
X       ==> a sequence of real-valued numbers
X   [2] start (keyword; 0)
X       ==> start index of sequence to be used
X   [3] end (keyword; length of the-seq)
X       ==> 1 + end index of sequence to be used
returns
X   [1] sample skewness of the specified numbers
X       in the-seq
X   [2] sample kurtosis
X   [3] sample mean
X   [4] sample variance"
X  (multiple-value-bind (the-mean sample-variance)
X                       (sample-mean-and-variance
X                        the-seq :start start :end end)
X    (let ((n (- end start))
X          (the-sum-of-cubes 0)
X          (the-sum-of-quads 0)
X          (j (1- start)))
X      (dotimes (i n (values (/ (/ the-sum-of-cubes n)
X                               (expt sample-variance 3/2))
X                            (- (/ (/ the-sum-of-quads n)
X                                  (* sample-variance sample-variance)) 3)
X                            the-mean sample-variance))
X        (let ((centered-value (- (elt the-seq (incf j)) the-mean)))
X          (incf the-sum-of-cubes
X                (* centered-value centered-value centered-value))
X          (incf the-sum-of-quads
X                (* centered-value centered-value
X                   centered-value centered-value)))))))
X
#|
(sample-skewness-and-kurtosis #(1 2 3 4 5 6 7 8 9 10))
;==> 0.0
;    -202/165
;    11/2
;    33/4
(sample-skewness-and-kurtosis #(1 2 3 4 5 6 7 8 9 10) :start 1 :end 6)
;==> 0.0
;    -13/10
;    4
;    2
|#
X
;-------------------------------------------------------------------------------
(defun sample-median
X       (the-seq
X        &key
X        (start 0)
X        (end (length the-seq))
X        (the-seq-is-ordered-p nil))
X  "given
X   [1] the-seq (required)
X       ==> a sequence of real-valued numbers
X   [2] start (keyword; 0)
X       ==> start index of sequence to be used
X   [3] end (keyword; length of the-seq)
X       ==> 1 + end index of sequence to be used
X   [4] the-seq-is-ordered-p (keyword; nil)
X       ==> if t, the-seq is assumed to
X           be sorted from smallest to largest
X           value; if nil, the-seq will be sorted
returns
X   [1] sample median
X   [2] minimum value in sequence
X   [3] maximum value
X   [4] sorted sequence"
X  ;;; Note: subseq generates a new sequence ...
X  (let ((seq-to-be-sorted (subseq the-seq start end))
X        (n-effective (- end start)))
X    (if (not the-seq-is-ordered-p)
X      (setf seq-to-be-sorted (sort seq-to-be-sorted #'<)))
X    (values
X     (if (oddp n-effective)
X       (elt seq-to-be-sorted (/ (1- n-effective) 2))
X       (/ (+ (elt seq-to-be-sorted (/ n-effective 2))
X             (elt seq-to-be-sorted (/ (- n-effective 2) 2)))
X          2))
X     (elt seq-to-be-sorted 0)
X     (elt seq-to-be-sorted (1- n-effective))
X     seq-to-be-sorted)))
X
#|
(sample-median '(1 10 2 9 5))
;==> 5
;    1
;    10
;    (1 2 5 9 10)
(sample-median '(1 10 2 9 5 0 0 0 0 0) :end 5)
;==> 5
;    1
;    10
;    (1 2 5 9 10)
(sample-median #(1.0 2.0 3.0))
;==> 2.0
;    1.0
;    3.0
;    #(1.0 2.0 3.0)
(sample-median #(1.0 2.0 3.0 4.0))
;==> 2.5
;    1.0
;    4.0
;    #(1.0 2.0 3.0 4.0)
(sample-median #(1.0 2.0 3.0 4.0) :end 3)
;==> 2.0
;    1.0
;    3.0
;    #(1.0 2.0 3.0)
(sample-median (vector 1.0 2.0 3.0 4.0) :start 1 :end 3)
;==> 2.5
;    2.0
;    3.0
;    #(2.0 3.0)
|#
X
;-------------------------------------------------------------------------------
(defun sample-correlation-coefficient
X       (seq-1
X        seq-2
X        &key
X        (start 0)
X        (end (length seq-1)))
X  "given
X   [1] seq-1 (required)
X       ==> a sequence of real-valued numbers
X   [2] seq-2 (required)
X       ==> another sequence of real-valued numbers;
X           should have the same length as seq-1
X   [3] start (keyword; 0)
X       ==> start index of sequence to be used
X   [4] end (keyword; length of the-seq)
X       ==> 1 + end index of sequence to be used
returns
X   [1] sample correlation coefficient
X       between two sequences"
X  (multiple-value-bind (mean-1 variance-1)
X                       (sample-mean-and-variance
X                        seq-1 :start start :end end)
X    (multiple-value-bind (mean-2 variance-2)
X                         (sample-mean-and-variance
X                          seq-2 :start start :end end)
X      (let ((sum 0.0)
X            (n (- end start))
X            (j (1- start)))
X        (dotimes (i n (/ (/ sum n) (sqrt (* variance-1 variance-2))))
X          (incf sum (* (- (elt seq-1 (incf j)) mean-1)
X                       (- (elt seq-2 j) mean-2))))))))
X
#|
(sample-correlation-coefficient
X #(0.0 -2.0 -1.0 0.0  1.0  2.0 10.0)
X #(7.0  2.0  1.0 0.0 -1.0 -2.0 10.0))
;==> 0.6190943270844442
(sample-correlation-coefficient
X #(0.0 -2.0 -1.0 0.0  1.0  2.0 10.0)
X #(7.0  2.0  1.0 0.0 -1.0 -2.0 10.0)
X :start 1 :end 6)
;==> -1.0
|#
X
;-------------------------------------------------------------------------------
;;; The function histogram computes all the goodies needed to construct
;;; a histogram for a sequence of deviates.
(defun histogram
X       (sequence-of-deviates
X        &key
X        (number-of-bins 10)
X        (left-of-first-bin (min-of-seq sequence-of-deviates))
X        (right-of-last-bin (max-of-seq sequence-of-deviates))
X        (scale-as-density-p t)
X        (result (make-array number-of-bins)))
X  "given
X   [1] sequence-of-deviates (required)
X       ==> a sequence of real-valued numbers
X   [2] number-of-bins (keyword; 10)
X       ==> number of bins in histogram
X   [3] left-of-first-bin (keyword; min of sequence-of-deviates)
X       ==> leftmost value of first bin
X   [4] right-of-last-bin (keyword; max of sequence-of-deviates)
X       ==> rightmost value of last bin
X   [5] scale-as-density-p (keyword; t)
X       ==> if t, histogram is scaled as a density;
X           otherwise, cell counts are returned
X   [6] result (keyword; vector of size number-of-bins)
X       <== vector to hold the histogram
returns
X   [1] the vector result, which contains the values for the histogram
X       for sequence-of-deviates (one value for each bin); note that,
X       if a bin has boundaries a and b, all points x in sequence-of-deviates
X       for which a<=x<b are counted as being in that bin
X       EXCEPT for the last bin, where the rule is a<=x<=b.
X   [2] maximum value in the histogram (i.e., result)
X   [3] leftmost value of first bin
X   [4] bin-width of bins
X   [5] the number of unbinned elements in sequence-of-deviates
X       less than left-of-first-bin
X   [6] the number of unbinned elements in sequence-of-deviates
X       greater than right-of-last-bin
---
Note: the value in (svref result 0) is the histogram value
X      for the interval [leftmost value of first bin,
X                        leftmost value of first bin + bin-width of bins);
X      the value in (svref result 1) is the histogram value
X      for the interval [leftmost value of first bin + bin-width of bins,
X                        leftmost value of first bin + 2*bin-width of bins);
X      etc."
X  (assert (plusp (length sequence-of-deviates)))
X  (assert (plusp number-of-bins))
X  (assert (< left-of-first-bin right-of-last-bin))
X  (let* ((n (length sequence-of-deviates))
X         (low-elements-not-binned 0)
X         (high-elements-not-binned 0)
X         (spacing (/ (- right-of-last-bin left-of-first-bin) number-of-bins))
X         (density-factor (if scale-as-density-p
X                           (float (* n spacing))
X                           1))
X         (the-maximum 0.0)
X         j)
X    (fill result 0 :end number-of-bins)
X    (dotimes (i n)
X      (setf j (floor (- (elt sequence-of-deviates i) left-of-first-bin)
X                     spacing))
X      (cond
X       ((< -1 j number-of-bins)
X        (incf (elt result j)))
X       ((= (elt sequence-of-deviates i) right-of-last-bin)
X        (incf (elt result (1- number-of-bins))))
X       (t
X        (if (minusp j)
X          (incf low-elements-not-binned)
X          (incf high-elements-not-binned)))))
X    (dotimes (i number-of-bins (values result
X                                       the-maximum
X                                       left-of-first-bin
X                                       spacing
X                                       low-elements-not-binned
X                                       high-elements-not-binned))
X      (if (> (setf (elt result i) (/ (elt result i) density-factor))
X             the-maximum)
X        (setf the-maximum (elt result i))))))
X
#|
(histogram '(9 1 9 9 1 8 1 9 6 7 9 3 9)
X           :number-of-bins 10
X           :left-of-first-bin 0.5
X           :right-of-last-bin 10.5
X           :scale-as-density-p nil)
;==> #(3 0 1 0 0 1 1 1 6 0)     3 values in [0.5,1.5); 6 in [8.5,9.5); etc.
;    6                          maximum value in histogram vector
;    0.5                        left-hand side of first histogram bin 
;    1.0                        width of each histogram bin 
;    0                          number of values < 0.5
;    0                          number of values > 10.5
X
;;; Now we repeat the above, but now scale the histogram as a density:
(histogram '(9 1 9 9 1 8 1 9 6 7 9 3 9)
X           :number-of-bins 10
X           :left-of-first-bin 0.5
X           :right-of-last-bin 10.5)
;==> #(0.23076923076923078 0.0 0.07692307692307693 0.0 0.0 0.07692307692307693 0.07692307692307693 0.07692307692307693 0.46153846153846156 0.0)
;    0.46153846153846156
;    0.5
;    1.0
;    0
;    0
;;; Since the bin width is 1.0, the elements of the vector above
;;; sum to unity; had the bin width been 0.5, the elements would have
;;; summed to 2.0:
(sum (histogram '(9 1 9 9 1 8 1 9 6 7 9 3 9)
X                :number-of-bins 20
X                :left-of-first-bin 0.5
X                :right-of-last-bin 10.5))
;==> 2.0
|#
X
;-------------------------------------------------------------------------------
;;; The function box-plot computes all the goodies needed to construct
;;; a box plot for a sequence of deviates.
(defun box-plot
X       (a-seq
X        &key
X        (seq-sorted-p nil)
X        (get-bottom-line
X         #'(lambda (v)
X             (let* ((lower-quartile (quantile-of-ordered-seq v 0.25))
X                    (iqr (- (quantile-of-ordered-seq v 0.75)
X                            lower-quartile))
X                    (lower-limit (- lower-quartile (* 1.5 iqr)))
X                    (temp (binary-search lower-limit v))
X                    (lower-index (if temp (1+ temp) 0))
X                    (lower-adjacent-value (elt v lower-index)))
X               (values lower-adjacent-value lower-index))))
X        (get-bottom-of-box
X         #'(lambda (v)
X             (quantile-of-ordered-seq v 0.25)))
X        (get-line-through-box
X         #'(lambda (v)
X             (quantile-of-ordered-seq v 0.50)))
X        (get-top-of-box
X         #'(lambda (v)
X             (quantile-of-ordered-seq v 0.75)))
X        (get-top-line
X         #'(lambda (v)
X             (let* ((upper-quartile (quantile-of-ordered-seq v 0.75))
X                    (iqr (- upper-quartile
X                            (quantile-of-ordered-seq v 0.25)))
X                    (upper-limit (+ upper-quartile (* 1.5 iqr)))
X                    (temp (binary-search upper-limit v))
X                    (upper-index (if temp temp (1- (length v))))
X                    (upper-adjacent-value (elt v upper-index)))
X               (values upper-adjacent-value upper-index)))))
X  "given
X   [1] a-seq (required)
X       ==> a sequence of real-valued numbers
X   [2] seq-sorted-p (keyword; nil)
X       ==> if t, a-seq is already sorted;
X           if nil, a copy of a-seq is sorted
X           for use by the function
X   [3] get-bottom-line (keyword; computes lower adjacent value)
X       ==> a function for calculating
X           the line below the box (from sorted data)
X   [4] get-bottom-of-box (keyword; computes lower quartile)
X       ==> a function for calculating
X           bottom of the box
X   [5] get-line-through-box (keyword; computes median)
X       ==> a function for calculating
X           the line through the middle of the box
X   [6] get-top-of-box (keyword; computes upper quartile)
X       ==> a function for calculating
X           top of the box
X   [7] get-top-line (keyword; computes upper adjacent value)
X       ==> a function for calculating
X           the line above the box
returns
X   [1] bottom line below box
X   [2] bottom line of box
X   [3] line through middle of box
X   [4] top line of box
X   [5] top line above box
X   [6] a vector containing outside values
X       (note: this vector can be of length 0)
---
Note: the values that this function returns can
X      be used to contruct a box plot as defined in
X      Section 2.5 of ``Graphical Methods for Data Analysis''
X      by Chambers, Cleveland, Kleiner, Tukey"
X  (if (not seq-sorted-p)
X    (setf a-seq (sort (copy-seq a-seq) #'<)))
X  (multiple-value-bind (top-value top-index)
X                       (funcall get-top-line a-seq)
X    (multiple-value-bind (bottom-value bottom-index)
X                         (funcall get-bottom-line a-seq)
X      (let* ((n (length a-seq))
X             (n-upper (- n top-index 1))
X             (n-outside (+ bottom-index n-upper))
X             (lower-outside-values '())
X             (upper-outside-values '())
X             (outside-values (cond
X                              ((plusp n-outside)
X                               (dotimes (i bottom-index)
X                                 (push (elt a-seq i)
X                                       lower-outside-values))
X                               (dotimes (i (- n top-index 1))
X                                 (push (elt a-seq (- n i 1))
X                                       upper-outside-values))
X                               (make-array n-outside
X                                           :initial-contents
X                                           (append
X                                            (reverse lower-outside-values)
X                                            upper-outside-values)))
X                              (t
X                               (make-array 0)))))
X        (values
X         bottom-value                           ;;; bottom line below box
X         (funcall get-bottom-of-box a-seq)      ;;; bottom line of box
X         (funcall get-line-through-box a-seq)   ;;; line through middle of box
X         (funcall get-top-of-box a-seq)         ;;; top line of box
X         top-value                              ;;; top line above box
X         outside-values                         ;;; outside-values
X         )))))
X
#|
(box-plot '(9 5 4 5 6 5 4 5 6 -100 5))
;     4                    bottom line below box
;     4.25                 bottom line of box
;     5.0                  line through middle of box
;     5.75                 top line of box
;     6                    top line above box
;     #(-100 9)            outside values
|#
X
;-------------------------------------------------------------------------------
;;; The function symmetry-plot computes all the goodies needed to construct
;;; a symmetry plot for a sequence of deviates.
(defun symmetry-plot
X       (a-seq
X        &key
X        (a-seq-sorted-p nil))
X  "given
X   [1] a-seq (required)
X       ==> any sequence of real-valued numbers
X   [2] a-seq-sorted-p (keyword; nil)
X       ==> if t, a-seq is already sorted;
X           if nil, a copy of a-seq is sorted
X           for use by the function
returns
X   [1] y values for a symmetry plot for a-seq
X   [2] x values for a symmetry plot for a-seq"
X  (if (not a-seq-sorted-p)
X    (setf a-seq (sort (copy-seq a-seq) #'<)))
X  (let* ((n (length a-seq))
X         (n-1 (1- n))
X         (m (truncate n 2))
X         (m-1 (1- m))
X         (x-values (make-array m))
X         (y-values (make-array m))
X         (the-median (quantile-of-ordered-seq a-seq 0.5)))
X    (dotimes (i m (values y-values x-values))
X      (setf (elt x-values (- m-1 i)) (- the-median (elt a-seq i))
X            (elt y-values (- m-1 i)) (- (elt a-seq (- n-1 i)) the-median)))))
X
#|
(symmetry-plot '(9 5 4 5 6 5 4 5 6 -100 5))
;==>  #(0.0 0.0 1.0 1.0 4.0)
;     #(0.0 0.0 1.0 1.0 105.0)
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  quantile-of-ordered-seq
;;;                 quantile-plot
;;;                 interquartile-range
;;;                 quantile-of-normal-distribution
;;;                 quantile-of-gamma-distribution
;;;                 quantile-of-exponential-distribution
;;;                 quantile-of-chi-square-2-distribution
;;;                 quantile-of-chi-square-distribution
;;;                 q-q-plot
;;;  are all concerned with computing either sample or theoretical quantiles
;;;  for, respectively, either a sequence of numbers or a theoretical
;;;  distribution.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun quantile-of-ordered-seq
X       (the-seq
X        p)
X  "given
X   [1] the-seq (required)
X       ==> a sequence of real-valued numbers,
X           ordered from smallest to largest
X   [2] p (required)
X       ==> a percentile; i.e., 0 <= p <= 1
returns
X   [1] the sample quantile for p
---
Note: see Section 2.2 of ``Graphical Methods for Data Analysis''
X      by Chambers, Cleveland, Kleiner, and Tukey"
X  (assert (and (plusp (length the-seq)) (<= 0.0 p 1.0)))
X  (let* ((n (length the-seq))
X         (p-sub-0 (/ 0.5 n))
X         (p-sub-n-1 (- 1.0 (/ 0.5 n))))
X    (cond
X     ((<= p p-sub-0)
X      (elt the-seq 0))
X     ((>= p p-sub-n-1)
X      (elt the-seq (1- n)))
X     (t
X      (let* ((i (max 0 (truncate (- (* n p) 0.5))))
X             (p-sub-i (/ (+ i 0.5) n))
X             (i+1 (min (1+ i) (1- n)))
X             (p-sub-i+1 (/ (+ i+1 0.5) n))
X             (delta-p (- p-sub-i+1 p-sub-i))
X             (f (if (zerop delta-p)
X                  0.0
X                  (/ (- p p-sub-i) delta-p))))
X        (+ (* (- 1.0 f) (elt the-seq i))
X           (* f (elt the-seq i+1))))))))
X
#|
(quantile-of-ordered-seq
X #(1.0 1.1 2.0 2.1 3.0 3.1 4.0 4.1 5.0 5.1)
X 0.1)
;==> 1.05
(quantile-of-ordered-seq
X #(1.0 1.1 2.0 2.1 3.0 3.1 4.0 4.1 5.0 5.1)
X 0.9)
;==> 5.05
|#
X
;-------------------------------------------------------------------------------
;;; The function quantile-plot computes all the goodies needed to construct
;;; a quantile plot for a sequence of deviates.
(defun quantile-plot
X       (a-seq
X        &key
X        (a-seq-sorted-p nil))
X  "given
X   [1] a-seq (required)
X       ==> any sequence of real-valued numbers
X   [2] a-seq-sorted-p (keyword; nil)
X       ==> if t, a-seq is already sorted;
X           if nil, a copy of a-seq is sorted
X           for use by the function
returns
X   [1] sequence with sample quantiles for a-seq
X   [2] sequence with corresponding percentiles
---
Note: see Section 2.2 of ``Graphical Methods for Data Analysis''
X      by Chambers, Cleveland, Kleiner, and Tukey"
X  (if (not a-seq-sorted-p)
X    (setf a-seq (sort (copy-seq a-seq) #'<)))
X  (let* ((n (length a-seq))
X         (x-values (copy-seq a-seq)))
X    (dotimes (i n (values a-seq x-values))
X      (setf (elt x-values i) (/ (+ i 0.5) n)))))
X
#|
(quantile-plot '(9 5 4 5 6 5 4 6 -100 5))
;==> (-100 4 4 5 5 5 5 6 6 9)
;    (0.05 0.15 0.25 0.35 0.45 0.55 0.65 0.75 0.85 0.95)
|#
X
;-------------------------------------------------------------------------------
(defun interquartile-range
X       (the-seq
X        &key
X        (start 0)
X        (end (length the-seq))
X        (the-seq-is-ordered-p nil))
X  "given
X   [1] the-seq (required)
X       ==> a sequence of real-valued numbers
X   [2] start (keyword; 0)
X       ==> start index of sequence to be used
X   [3] end (keyword; length of the-seq)
X       ==> 1 + end index of sequence to be used
X   [4] the-seq-is-ordered-p (keyword; nil)
X       ==> if t, the-seq is assumed to
X           be sorted from smallest to largest
X           value; if nil, the-seq will be sorted
returns
X   [1] interquartile range; i.e,
X       0.75 quantile minus 0.25 quantile"
X  (let ((seq-to-be-sorted (subseq the-seq start end)))
X    (if (not the-seq-is-ordered-p)
X      (setf seq-to-be-sorted (sort seq-to-be-sorted #'<)))
X    (values (- (quantile-of-ordered-seq seq-to-be-sorted 0.75)
X               (quantile-of-ordered-seq seq-to-be-sorted 0.25)))))
X
#|
(quantile-of-ordered-seq
X #(1.0 1.1 2.0 2.1 3.0 3.1 4.0 4.1 5.0 5.1)
X 0.75)
;==> 4.1
(quantile-of-ordered-seq
X #(1.0 1.1 2.0 2.1 3.0 3.1 4.0 4.1 5.0 5.1)
X 0.25)
;==> 2.0
(interquartile-range #(1.0 1.1 2.0 2.1 3.0 3.1 4.0 4.1 5.0 5.1))
;==> 2.0999999999999996
|#
X
;-------------------------------------------------------------------------------
(defun quantile-of-normal-distribution
X       (p
X        &key
X        (accuracy-level :quick))
X  "given
X   [1] p (required)
X       ==> a number > 0 and < 1
X   [2] accuracy-level (keyword; :quick)
X       ==> if :quick, gives quick approximation;
X           if :better, gives slower, but better, approximation;
X           if :accurate, uses slowest, but more accurate,
X           approximation;
returns
X   [1] quantile of standard normal distribution;
---
Note: for :quick, see Table 6.5, page 227,
X      ``Graphical Methods for Data Analysis''
X      by Chambers, Cleveland, Kleiner, & Tukey;
X      for :better, see Sections 26.2.22 and 26.2.23
X      of Abramowitz and Stegun;
X      for :accurate, see AS 111 by Beasley and Springer,
X      Applied Statistics, 1977, vol. 26, p.118"
X  (case accuracy-level
X    (:quick
X     (let ((c (sqrt (* -2.0 (log (min p (- 1.0 p)))))))
X       (* (signum (- p 0.5))
X          (- c (/ (+ 2.30753 (* 0.27061 c))
X                  (+ 1.0 (* 0.99229 c) (* 0.04481 c c)))))))
X    (:better
X     (let ((c (sqrt (* -2.0 (log (min p (- 1.0 p)))))))
X       (* (signum (- p 0.5))
X          (- c (/ (+ 2.515517 (* 0.802853 c) (* 0.010328 c c))
X                  (+ 1.0 (* 1.432788 c) (* 0.189269 c c) (* 0.001308 c c c)))))))
X    (:accurate
X     (as-111 p))
X    (otherwise
X     (error "accuracy-level set to ~A; should be either quick, better or accurate"
X            accuracy-level))))
X
#|
(quantile-of-normal-distribution 0.975) 
;==> 1.960448273742351
(quantile-of-normal-distribution 0.975 :accuracy-level :better)
;==> 1.9603949169253396
(quantile-of-normal-distribution 0.975 :accuracy-level :accurate)
;==> 1.9599639822884356
|#
X
;-------------------------------------------------------------------------------
(defun quantile-of-gamma-distribution
X       (alpha
X        p)
X  "given
X   [1] alpha (required)
X       ==> shape parameter for gamma distribution
X   [2] p (required)
X       ==> a number > 0 and < 1
returns
X   [1] p-th quantile of standard gamma distribution
---
Note: Table 6.4, page 225, ``Graphical Methods for Data Analysis''
X      by Chambers, Cleveland, Kleiner, & Tukey"
X  (* alpha (expt (- 1.0
X                    (/ (* 9.0 alpha))
X                    (/ (quantile-of-normal-distribution p)
X                       ; Note: don't panic! --- -3.0 negates subtraction
X                       (* -3.0 (sqrt alpha))))
X                 3)))
X
#|
(quantile-of-gamma-distribution 1.0 0.975) 
;==> 3.6691637921673883
|#
X
;-------------------------------------------------------------------------------
(defun quantile-of-exponential-distribution
X       (p)
X  "given
X   [1] p (required)
X       ==> a number > 0 and < 1
returns
X   [1] p-th quantile of exponential distribution with mean 1
---
Note: Table 6.4, page 225, ``Graphical Methods for Data Analysis''
X      by Chambers, Cleveland, Kleiner, & Tukey"
X  (- (log (- 1.0 p))))
X
#|
(quantile-of-exponential-distribution 0.975)
;==> 3.6888794541139354
|#
X
;-------------------------------------------------------------------------------
(defun quantile-of-chi-square-2-distribution
X       (p)
X  "given
X   [1] p (required)
X       ==> a number > 0 and < 1
returns
X   [1] p-th quantile of chi-square distribution with 2 degrees of freedom
---
Note: Equation (31), P. 640 of Chave, Thomson, and Ander (1987)"
X  (* -2.0 (log (- 1.0 p))))
X
#|
(quantile-of-chi-square-2-distribution 0.975)
;==> 7.377758908227871
|#
X
;-------------------------------------------------------------------------------
(defun quantile-of-chi-square-distribution
X       (nu
X        p
X        &key
X        (accuracy-level :quick))
X  "given
X   [1] nu (required)
X       ==> degrees of freedom
X   [2] p (required)
X       ==> a number > 0 and < 1
X   [3] accuracy-level (keyword; :quick)
X       if :quick, gives quick approximation;
X       if :accurate, give accurate, but computationally
X       expensive, approximation
returns
X   [1] p-th quantile of chi-square distribution with nu degrees of freedom
---
Note: for :quick, see Table 6.4, page 225,
X      ``Graphical Methods for Data Analysis''
X      by Chambers, Cleveland, Kleiner, & Tukey;
X      if :accurate, see AS 91, Applied Statistics"
X  (case accuracy-level
X    (:quick
X     (* 2.0 (quantile-of-gamma-distribution (/ nu 2.0) p)))
X    (:accurate 
X     (as-91 p nu))
X    (otherwise
X     (error "accuracy-level set to ~A; should be either :quick or :accurate"
X            accuracy-level))))
X
#|
(quantile-of-chi-square-distribution 2 0.975)
;==> 7.338327584334777
(quantile-of-chi-square-distribution 2 0.975 :accuracy-level :accurate)
;==> 7.377758908227873
(quantile-of-chi-square-distribution 4 0.975)
;==> 11.130220797449477
(quantile-of-chi-square-distribution 4 0.975 :accuracy-level :accurate)
;==> 11.143286781877793
(quantile-of-chi-square-distribution 8 0.975)
;==> 17.533997167548357
(quantile-of-chi-square-distribution 8 0.975 :accuracy-level :accurate)
;==> 17.534546139484654
(quantile-of-chi-square-distribution 8.3 0.975)
;==> 17.984123910122445
(quantile-of-chi-square-distribution 8.3 0.975 :accuracy-level :accurate)
;==> 17.98423969597426
(quantile-of-chi-square-distribution 9 0.975)
;==> 19.023543945139263
(quantile-of-chi-square-distribution 9 0.975 :accuracy-level :accurate)
;==> 19.02276780221112
|#
X
;-------------------------------------------------------------------------------
;;; The function q-q-plot is organized so to give the x and y values
;;; needed to create a q-q plot for any of the the following combinations:
;;;   empirical (y) versus theoretical (x)
;;;   empirical     versus empirical
;;;   theoretical   versus theoretical
;;;
;;; For a theoretical component, you need to supply a function
;;; which, given a value from 0 to 1, calculates the corresponding quantile;
;;; for an empirical, you need to supply either a sequence of observations
;;; or a data object from which a vector is to be extracted (it is assumed
;;; to be the vector pointed to by y-values).
;;;
;;; OPERATIONAL NOTE:  if you are doing many different Q-Q plots of data sets
;;;                    with the same number of points versus the same
;;;                    theoretical distribution, it might be faster to
;;;                    calculate the quantiles just once pass them in as
;;;                    a sequence with [xy]-sorted-p set to t.
(defun q-q-plot
X       (y-vals
X        x-vals
X        &key
X        (n (if (typep y-vals 'sequence)
X             (length y-vals)
X             256))
X        (x-vals-sorted-p (functionp x-vals))
X        (y-vals-sorted-p (functionp y-vals)))
X  "given
X   [1] y-vals (required)
X       ==> a sequence or function
X   [2] x-vals (required)
X       ==> a sequence or function
X   [3] n (keyword; 256 if y-vals is a function, length of y-vals otherwise)
X       ==> positive integer
X   [4] x-vals-sorted-p (keyword; t if x-vals is a function, nil otherwise)
X       ==> a flag indicating whether x-vals is already sorted
X   [5] y-vals-sorted-p (keyword; t if y-vals is a function, nil otherwise)
X       ==> a flag indicating whether y-vals is already sorted
return
X   [1] y values needed to create a q-q-plot for y-vals
X   [2] x values
---
Note: For details, see Chapter 6 of ``Graphical Methods for Data Analysis''
X      by Chambers, Cleveland, Kleiner and Tukey"
X  (let* ((y-copy-already-made-p nil)
X         (x-copy-already-made-p nil)
X         (y-values (cond
X                    ((functionp y-vals)
X                     (setf y-copy-already-made-p
X                           (sample-from-a-function
X                            y-vals
X                            :x0 (/ 0.5 n)
X                            :delta-x (/ (float n))
X                            :n n)))
X                    ((arrayp y-vals)
X                     (if y-vals-sorted-p y-vals
X                         (setf y-copy-already-made-p
X                               (copy-seq y-vals))))
X                    (t
X                     (setf y-copy-already-made-p
X                           (make-array n
X                                       :initial-contents
X                                       y-vals)))))
X         (x-values (cond
X                    ((functionp x-vals)
X                     (setf x-copy-already-made-p
X                           (sample-from-a-function
X                            x-vals
X                            :x0 (/ 0.5 n)
X                            :delta-x (/ (float n))
X                            :n n)))
X                    ((arrayp x-vals)
X                     (if x-vals-sorted-p x-vals
X                         (setf x-copy-already-made-p
X                               (copy-seq x-vals))))
X                    (t
X                     (setf x-copy-already-made-p
X                           (make-array n
X                                       :initial-contents
X                                       x-vals))))))
X    (cond
X     ((not x-vals-sorted-p)
X      (if (not x-copy-already-made-p)
X        (setf x-copy-already-made-p (setf x-values (copy-seq x-values))))
X      (setf x-values (sort x-values #'<=))))
X    (cond
X     ((not y-vals-sorted-p)
X      (if (not y-copy-already-made-p)
X        (setf y-copy-already-made-p (setf y-values (copy-seq y-values))))
X      (setf y-values (sort y-values #'<=))))
X    (values y-values x-values)))
X
#|
(q-q-plot
X '(9 5 4 5 6 5 4 5 6 -100 5)
X '(8 4 3 4 5 4 3 4 5  -99 4))
;==> #(-100 4 4 5 5 5 5 5 6 6 9)
;     #(-99 3 3 4 4 4 4 4 5 5 8)
X
(q-q-plot
X #(9 5 4 5 6 5 4 5 6 -100 5)
X #(8 4 3 4 5 4 3 4 5  -99 4))
;==> #(-100 4 4 5 5 5 5 5 6 6 9)
;     #(-99 3 3 4 4 4 4 4 5 5 8)
X
(q-q-plot
X #(9 5 4 5 6 5 4 5 6 -100 5)
X  #'quantile-of-normal-distribution)
;==> #(-100 4 4 5 5 5 5 5 6 6 9)
;    #(-1.6903897928628724 -1.0948565688424345 -0.7451748604784388 -0.470072369736396 -0.22795082525859311 0.0 0.22795082525859311 0.470072369736396 0.7451748604784388 1.0948565688424352 1.6903897928628742)
X
(q-q-plot
X #'(lambda (p) (quantile-of-chi-square-distribution 100 p))
X #'quantile-of-normal-distribution
X :n 4)
;==> #(84.02560895964486 94.94772287048265 103.85500296347429 116.39899541971275)
;    #(-1.1485487296374308 -0.3163006861455602 0.3163006861455602 1.1485487296374308
|#
X
X
;-------------------------------------------------------------------------------
;;; This is needed for standard-normal-pdf below.
(defconstant +sqrt-of-two-pi+ (sqrt (* 2.0 pi)))
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  standard-normal-pdf
;;;                 tail-area-of-normal-distribution
;;;  compute the pdf and the tail area for the standard normal (Gaussian)
;;;  distribution.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun standard-normal-pdf
X       (x)
X  "given x, returns value of standard normal pdf at x"
X  (/ (exp (/ (* x x) -2.0)) +sqrt-of-two-pi+))
X
#|
(standard-normal-pdf  0.0)  ;==> 0.3989422804014327
(standard-normal-pdf  1.0)  ;==> 0.24197072451914337
(standard-normal-pdf -1.0)  ;==> 0.24197072451914337
|#
X
;-------------------------------------------------------------------------------
(defun tail-area-of-normal-distribution
X       (q
X        &key
X        (q-to-infinity-p nil))
X  "given
X   [1] q (required)
X       ==> a quantile
X   [2] q-to-infinity-p (keyword; nil)
X       ==> if t,   return integral from q to +infinity;
X           if nil, return integral from -infinity to q
return
X   tail area of the standard norm from
X   either q to +infinity (if q-to-infinity-p is true)
X   or     -infinity to q (if q-to-infinity-p is nil)
---
see  Algorithm AS-66, Applied Statistics,
1973, vol. 22, no. 3"
X  (if (minusp q)
X    (setf q-to-infinity-p (not q-to-infinity-p) 
X          q (- q)))
X  (cond
X   ((or (<= q 7.0d0)
X        (and q-to-infinity-p (<= q 18.66d0)))
X    (let* ((y (* 0.5 q q))
X           (alnorm (if (> q 1.28d0)
X                     (* 0.398942280385d0
X                        (/ (exp (- y))
X                           (+ q -3.8052d-8
X                              (/ 1.00000615302d0
X                                 (+ q 3.98064794d-4
X                                    (/ 1.98615381364d0
X                                       (+ q -0.151679116635d0
X                                          (/ 5.29330324926d0
X                                             (+ q 4.8385912808d0
X                                                (/ -15.1508972451d0
X                                                   (+ q 0.742380924027d0
X                                                      (/ 30.789933034d0 (+ q 3.99019417011d0)))))))))))))
X                     (- 0.5
X                        (* q
X                           (- 0.398942280444d0
X                              (* 0.39990348504d0
X                                 (/ y
X                                    (+ y 5.75885480458d0
X                                       (/ -29.8213557807d0
X                                          (+ y 2.62433121679d0 (/ 48.6959930692d0 (+ y 5.92885724438d0)))))))))))))
X      (if q-to-infinity-p alnorm (- 1.0 alnorm))))
X   ;;; too far out in one direction, so return 0.0 or 1.0 ...
X   (t
X    (if q-to-infinity-p 0.0 1.0))))
X
#|
;;; Here we compare tail-area-of-normal-distribution to values tabulated
;;; in Abramowitz and Stegun, Table 26.1, pages 966--72 -- the space
;;; in each number below indicates the point at which the computed value
;;; differs from the tabulated value.
(tail-area-of-normal-distribution 0.02)     ;==> 0.50797831371 76848
(tail-area-of-normal-distribution 1.14)     ;==> 0.8728568 399049156
(tail-area-of-normal-distribution 1.86)     ;==> 0.96855723701 83276
(tail-area-of-normal-distribution 2.54)     ;==> 0.99445737655 7139
(tail-area-of-normal-distribution 3.45)     ;==> 0.9997197067 231738
(tail-area-of-normal-distribution 4.85)     ;==> 0.9999993826 92628
X
;;; a few additional tests ...
(tail-area-of-normal-distribution
X (quantile-of-normal-distribution 0.025))
;==> 0.02497170911465468
(tail-area-of-normal-distribution
X (quantile-of-normal-distribution 0.975))
;==> 0.9750282908853453
(tail-area-of-normal-distribution
X (quantile-of-normal-distribution 0.975 :accuracy-level :better))
;==> 0.9750251752384552
(tail-area-of-normal-distribution
X (quantile-of-normal-distribution 0.975 :accuracy-level :accurate))
;==> 0.974999999867449
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  median-absolute-deviation
;;;                 Thomson-weight-function
;;;                 Huber-weight-function
;;;                 m-location-estimate
;;;  are all concerned with the calcuation of robust statistics.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun median-absolute-deviation
X       (the-seq
X        &key
X        (location-estimate (sample-median the-seq)))
X  "given
X   [1] the-seq (required)
X       ==> a sequence of real-valued numbers
X   [2] location-estimate (keyword; median of the-seq)
X       ==> robust estimate of location
returns
X   [1] the median absolute deviation about location-estimate
X   [2] location-estimate (median of the-seq by default)"
X  (let* ((n (length the-seq))
X         (scratch (make-array n)))
X    ;(declare (dynamic-extent scratch))
X    (dotimes (i n (values (sample-median scratch)
X                          location-estimate))
X      (setf (elt scratch i)
X            (abs (- (elt the-seq i) location-estimate))))))
X
#|
(median-absolute-deviation '(6.0 2.0 1.0 4.0 5.0 3.0))
;==> 1.5
;    3.5
X
(median-absolute-deviation ' (2.0 3.0 4.0 5.0 6.0))
;==> 1.0
;    4.0
|#
X
;-------------------------------------------------------------------------------
(defun Thomson-weight-function
X       (x
X        &key
X        (beta 1.0)
X        (symmetric-version t))
X  "given
X   [1] x (required)
X       ==> value at which to evaluate Thomson's
X           weight function
X   [2] beta (keyword; 1.0)
X       ==> tuning parameter
X   [3] symmetric-version (keyword; t)
X       ==> one-sided or two-sided down weightings
returns
X   [1] value of Thomson's weight function at x
---
Note: Equation (27), p. 637 of Chave, Thomson, and Ander (1987)"
X  (let* ((inner-argument (* beta (- (if symmetric-version (abs x) x)
X                                    beta)))
X         (inner-exp (if (< (abs inner-argument) 50.0)
X                      (exp inner-argument)
X                      50.0)))
X    (if (< inner-exp 50.0)
X      (exp (- inner-exp))
X      0.0)))
X
#|
(Thomson-weight-function 0.0)
;==> 0.6922006275553464
(Thomson-weight-function 0.0 :symmetric-version nil)
;==> 0.6922006275553464
|#
X
;-------------------------------------------------------------------------------
(defun Huber-weight-function
X       (x
X        &key
X        (parameter 1.5))
X  "given
X   [1] x (required)
X       ==> value at which to evaluate Huber's weight function
X   [2] parameter (keyword; 1.5)
X       ==> tuning parameter
returns
X   [1] value of Huber's weight function at x"
X  (if (> (abs x) parameter)
X      (/ parameter (abs x))
X      1.0))
X
#|
(Huber-weight-function -100.0)
;==> 0.015
|#
;;;
X
;-------------------------------------------------------------------------------
(defun m-location-estimate
X       (the-seq
X        &key
X        (weighting-functions
X         `(,#'(lambda (x)
X                (Huber-weight-function
X                 x :parameter 2.0))
X           ,#'(lambda (x)
X                (Thomson-weight-function
X                 x :beta (/ (- (length the-seq) 0.5)
X                            (length the-seq))
X                 :symmetric-version t))))
X        (iterations '(10 0)))
X  "given
X   [1] the-seq (required)
X       ==> a sequence of numbers
X   [2] weighting-functions (keyword; Huber and Thomson)
X       ==> a list of two weighting functions
X   [3] iterations (keyword; '(10 0))
X       ==> a list of nonegative integers
X           giving the number of iterations
X           for each weighting function
returns
X   [1] m-location estimate for sequence
X   [2] MAD scale estimate (based upon deviations
X       from m-location estimate)
---
Note: see article by Hogg in Launer and Wilkinson"
X  (let ((n (length the-seq))
X        (i -1)
X        scaled-deviation final-mad scale-estimate
X        weight sum-of-weights sum-of-weighted-scaled-deviations)
X    (multiple-value-bind (mad-scale-estimate the-answer)
X                         (median-absolute-deviation the-seq)
X      (setf scale-estimate (/ mad-scale-estimate 0.6745))
X      (dolist (a-weighting-function weighting-functions)
X        (dotimes (j (nth (incf i) iterations))
X          (setf sum-of-weights 0.0
X                sum-of-weighted-scaled-deviations 0.0)
X          (dotimes (k n)
X            (setf scaled-deviation (/ (- (elt the-seq k) the-answer)	
X                                      scale-estimate)
X                  weight (funcall a-weighting-function scaled-deviation))
X            (incf sum-of-weights weight)
X            (incf sum-of-weighted-scaled-deviations
X                  (* weight scaled-deviation))))
X        (incf the-answer (* (/ sum-of-weighted-scaled-deviations sum-of-weights)
X                            scale-estimate)))
X      (setf final-mad (/ (median-absolute-deviation
X                          the-seq :location-estimate the-answer)
X                         0.6745))
X      (values the-answer final-mad))))
X
#|
(m-location-estimate '(9 5 4 5 6 5 4 5 6 -100 5))
;==> 5.0
;    1.4825796886582654
X
(m-location-estimate '(9 5 4 5 6 5 4 5 6 -100 5)
X                     :iterations '(10 1))
;==> 5.00390370448126
;    1.4767921356838247
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  ordinary-least-squares-Cholesky
;;;                 ordinary-least-squares-Q-R
;;;                 weighted-linear-least-squares
;;;                 Durbin-Watson-test-statistic
;;;                 predicted-y-at-x
;;;                 var-predicted-y-at-x
;;;                 var-predicted-mean-at-x
;;;  handle linear regression analysis.  The first three of these functions
;;;  are the main computational engines.  Currently here is what distinguishes
;;;  these three routines:
;;;  [1] ordinary-least-squares-Cholesky uses a Cholesky decomposition;
;;;      it accepts an arbitrary list of independent variables, but all
;;;      variables in the list must be either functions or vectors; and
;;;      the function does not currently compute any error statistics
;;;      (other than the residuals themselves).
;;;  [2] ordinary-least-squares-Q-R uses a Q-R (i.e., modified Gram-Schmidt)
;;;      decomposition; it accepts an arbitrary list of independent variables,
;;;      any one of which can be either a function or a vector; and the
;;;      function optionally returns error statistics such as residuals,
;;;      variances of parameter estimates and the correlation matrix for
;;;      these estimates.  This function is somewhat slower than
;;;      ordinary-least-squares-Cholesky.
;;;  [3] weighted-linear-least-squares only works for the simple model
;;;          y_k = alpha + beta * x_k
;;;      and uses the method outlined in Chapter 2 of 
;;;      ``Applied Regression Analysis'' by Draper and Smith, 1966,
;;;      to estimate alpha and beta.  The independent variable can be
;;;      expressed as either a vector or a function.  This function optionally
;;;      returns error statistics.  It is faster than the other two functions.
;;;  The functions predicted-y-at-x, var-predicted-y-at-x and
;;;  var-predicted-mean-at-x are intended to be used in conjunction
;;;  with the values returned by weighted-linear-least-squares.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun ordinary-least-squares-Cholesky
X       (dependent-variable
X        list-of-independent-variables
X        &key
X        (compute-residuals-p nil)
X        (result (if compute-residuals-p
X                  (make-array
X                   (length dependent-variable)))))
X  "given
X   [1] dependent-variable (required)
X        ==> a vector
X   [2] list-of-independent-variables (required)
X       ==> a list of either vectors or functions
X           (must be ALL vectors or ALL functions)
X   [3] compute-residuals-p (keyword; nil)
X       ==> if t, residuals are computed
X   [4] result (keyword; new array if compute-residuals-p is t; otherwise nil)
X       <== storage space for residuals
X           (not used unless compute-residuals-p is true)
returns
X   [1] vector with estimated parameters
X   [2] residuals
---
Note: uses Lisp versions of Cholesky factorization routines
X      from linpack"
X  ;;; NOTE (5/7/93): the restriction that list-of-independent-variables be
X  ;;;                all vectors or all functions (and not a mixture thereof)
X  ;;;                was made for programming convenience.  To correct the code
X  ;;;                so that it works for a mixture (as is true for the Q-R
X  ;;;                version of this function), it is only necessary to modify
X  ;;;                the code in the cond clause below -- the code to compute
X  ;;;                the residuals can already handle a mixture.
X  (assert (plusp (length list-of-independent-variables)))
X  (let* ((N (length dependent-variable))
X         (p (length list-of-independent-variables))
X         (X-prime-X (make-array `(,p ,p)))
X         (X-prime-Y (make-array p))
X         (off-diag-to-do (1- p))
X         offset an-independent-variable another-independent-variable)
X    (cond
X     ((vectorp (car list-of-independent-variables))
X      ;;; independent variables are a list of vectors
X      (dotimes (i p)
X        (setf (aref X-prime-X i i) (sum-of-squares
X                                    (elt list-of-independent-variables i))
X              (aref X-prime-Y i) (dot-product
X                                  (elt list-of-independent-variables i)
X                                  dependent-variable)
X              offset (1+ i))
X        (dotimes (j off-diag-to-do)
X          (setf (aref X-prime-X i offset) (dot-product
X                                           (elt list-of-independent-variables i)
X                                           (elt list-of-independent-variables
X                                                (+ i j 1)))
X                (aref X-prime-X offset i) (aref X-prime-X i offset))
X          (incf offset))
X        (decf off-diag-to-do)))
X     (t
X      ;;; independent variables are a list of functions
X      (dotimes (i p)
X        (setf an-independent-variable (elt list-of-independent-variables i))
X        (setf (aref X-prime-X i i)
X              (let ((SS 0.0))
X                (dotimes (k N SS)
X                  (incf SS
X                        (expt (funcall an-independent-variable k) 2))))
X              (aref X-prime-Y i)
X              (let ((sum 0.0))
X                (dotimes (k N sum)
X                  (incf sum
X                        (* (funcall an-independent-variable k)
X                           (aref dependent-variable k)))))
X              offset (1+ i))
X        (dotimes (j off-diag-to-do)
X          (setf another-independent-variable
X                (elt list-of-independent-variables (+ i j 1)))
X          (setf (aref X-prime-X i offset)
X                (let ((sum 0.0))
X                  (dotimes (k N sum)
X                    (incf sum
X                          (* (funcall an-independent-variable k)
X                             (funcall another-independent-variable k)))))
X                (aref X-prime-X offset i) (aref X-prime-X i offset))
X          (incf offset))
X        (decf off-diag-to-do))))
X    (spofa! X-prime-X)
X    (sposl! X-prime-X X-prime-Y)
X    ;;; X-prime-Y now contains solution vector to normal equation
X    (when compute-residuals-p
X      (copy-vector dependent-variable result)
X      (dotimes (i p)
X        (setf an-independent-variable (elt list-of-independent-variables i))
X        (if (vectorp an-independent-variable)
X          (dotimes (j N)
X            (decf (aref result j)
X                  (* (aref X-prime-Y i)
X                     (aref an-independent-variable j))))
X          (dotimes (j N)
X            (decf (aref result j)
X                  (* (aref X-prime-Y i)
X                     (funcall an-independent-variable j)))))))
X    (values X-prime-Y result)))
X
#|
(ordinary-least-squares-Cholesky
X #(1.0 1.9 3.1 3.9 5.1 5.9 7.1 7.9 9.1 9.9)
X `(,(make-array 10 :initial-element 1.0)
X   ,#(0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5))
X :compute-residuals-p t)
;==> #(0.9927272727272731 1.9987878787878786)
;    #(0.0072727272727268755 -0.09212121212121249 0.1084848484848484 -0.09090909090909127 0.1096969696969694 -0.08969696969696894 0.1109090909090904 -0.08848484848484794 0.1121212121212114 -0.08727272727272783)
|#
X
;-------------------------------------------------------------------------------
(defun ordinary-least-squares-Q-R
X       (dependent-variable
X        list-of-independent-variables
X        &key
X        (compute-residuals-p nil)
X        (compute-standard-errors-p nil)
X        (compute-correlation-matrix-p nil)
X        (residuals (if (or compute-residuals-p
X                           compute-standard-errors-p
X                           compute-correlation-matrix-p)
X                     (make-array
X                      (length dependent-variable)))))
X  "given
X   [1] dependent-variable (required)
X       ==> a vector
X   [2] list-of-independent-variables (required)
X       ==> a list of either vectors, functions or mixture thereof
X   [3] compute-residuals-p (keyword; nil)
X       ==> if t, residuals are computed
X   [4] compute-standard-errors-p (keyword; nil)
X       ==> if t, standard errors for parameter estimates are computed,
X           along with residuals (even if compute-residuals-p is nil)
X   [5] compute-correlation-matrix-p (keyword; nil)
X       ==> if t, correlation matrix for parameter estimates is computed,
X           along with residuals and standard errors (even if either
X           of their associated keywords is nil) 
X   [6] residuals (keyword; nil or vector)
X       <== storage space for residuals
X           (nil unless either compute-residuals-p,
X           compute-standard-errors-p or
X           compute-correlation-matrix-p is true)
returns
X   [1] vector with parameter estimates
X   [2] vector with residuals
X   [2] standard error of residuals
X   [3] vector with standard errors of parameter estimates
X   [4] correlation matrix
---
Note: this routine is based on the discussion
X      in ``Nonlinear Regression Analysis and Its Applications''
X      by Bates and Watts, Wiley, 1988."
X  (assert (plusp (length list-of-independent-variables)))
X  ;;; To avoid some tricky code, we "stack up" various options:
X  (cond
X   (compute-correlation-matrix-p
X    (setf compute-residuals-p t
X          compute-standard-errors-p t))
X   (compute-standard-errors-p
X    (setf compute-residuals-p t)))
X  (let* ((N (length dependent-variable))
X         (p (length list-of-independent-variables))
X         (design-matrix (make-array `(,N ,p)))
X         (standard-errors (if compute-standard-errors-p
X                            (make-array p :initial-element 0.0)))
X         (s-estimate nil)
X         an-independent-variable)
X    ;;; (declare (dynamic-extent design-matrix))
X    ;;; create the design matrix ...
X    (dotimes (i p)
X      (setf an-independent-variable (elt list-of-independent-variables i))
X      (if (vectorp an-independent-variable)
X        (dotimes (j N)
X          (setf (aref design-matrix j i) (aref an-independent-variable j)))
X        (dotimes (j N)
X          (setf (aref design-matrix j i) (funcall an-independent-variable j)))))
X    ;;; Q-R factor the design matrix ...
X    (multiple-value-bind (Q R) (Q-R! design-matrix)
X      ;;; ... followed by computation of parameter estimates
X      (let ((parameter-estimates (upper-triangular-solve!
X                                  R
X                                  (multiply-matrix-and-vector
X                                   (transpose Q)
X                                   dependent-variable))))
X        (when compute-residuals-p
X          (copy-vector dependent-variable residuals)
X          (dotimes (i p)
X            (setf an-independent-variable (elt list-of-independent-variables i))
X            (if (vectorp an-independent-variable)
X              (dotimes (j N)
X                (decf (aref residuals j)
X                      (* (aref parameter-estimates i)
X                         (aref an-independent-variable j))))
X              (dotimes (j N)
X                (decf (aref residuals j)
X                      (* (aref parameter-estimates i)
X                         (funcall an-independent-variable j))))))
X          (setf s-estimate (sqrt (/ (sum-of-squares residuals)
X                                    (- N p)))))
X        (if compute-standard-errors-p
X          ;;; Here we produce standard errors for the parameter estimates
X          ;;; (along with a correlation matrix).
X          ;;; First, we need to produce the inverse of R:
X          (let ((R-inverse (make-array `(,p ,p)))
X                rhs-vector row-length)
X            (dotimes (i p)
X              (setf rhs-vector (make-array p :initial-element 0.0))
X              (setf (aref rhs-vector i) 1.0)
X              (upper-triangular-solve! R rhs-vector)
X              (dotimes (j p)
X                (setf (aref R-inverse j i) (aref rhs-vector j))))
X            ;;; We now have the inverse of R, from which we can produce
X            ;;; standard errors for the parameter estimates.
X            (dotimes (i p)
X              (dotimes (j p)
X                (incf (svref standard-errors i)
X                      (expt (aref R-inverse i j) 2)))
X              (setf row-length (sqrt (elt standard-errors i)))
X              (setf (svref standard-errors i) (* row-length s-estimate))
X              ;;; If correlation matrix is to be produced,
X              ;;; we need to normalize
X              (if compute-correlation-matrix-p
X                (let ((1-over-row-length (/ row-length)))
X                  (dotimes (j p)
X                    (multf (aref R-inverse i j) 1-over-row-length)))))
X            (values
X             parameter-estimates
X             residuals
X             s-estimate
X             standard-errors
X             (if compute-correlation-matrix-p
X               (multiply-two-matrices R-inverse (transpose R-inverse)))))
X          ;;; if compute-standard-errors-p is nil, here is what we return: 
X          (values
X             parameter-estimates
X             residuals
X             s-estimate))))))
X
#|
(ordinary-least-squares-Q-R
X #(1.0 1.9 3.1 3.9 5.1 5.9 7.1 7.9 9.1 9.9)
X `(,(make-array 10 :initial-element 1.0)
X   ,#(0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5))
X :compute-residuals-p t)
;==> #(0.9927272727272709 1.9987878787878792)
;    #(0.007272727272729096 -0.0921212121212106 0.10848484848484996 -0.09090909090908994 0.10969696969697074 -0.08969696969696805 0.11090909090909129 -0.08848484848484706 0.11212121212121229 -0.08727272727272606)
;    0.10545715775238804
X
(ordinary-least-squares-Q-R
X #(1.0 1.9 3.1 3.9 5.1 5.9 7.1 7.9 9.1 9.9)
X `(,(make-array 10 :initial-element 1.0)
X   ,#(0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5))
X :compute-correlation-matrix-p t)
;==> #(0.9927272727272709 1.9987878787878792)
;    #(0.007272727272729096 -0.0921212121212106 0.10848484848484996 -0.09090909090908994 0.10969696969697074 -0.08969696969696805 0.11090909090909129 -0.08848484848484706 0.11212121212121229 -0.08727272727272606)
;    0.10545715775238804
;    #(0.06198284664515573 0.023220901891718743)
;    #2a((0.9999999999999999 -0.8429272304235245) (-0.8429272304235245 1.0))
|#
X
;-------------------------------------------------------------------------------
(defun weighted-linear-least-squares
X       (y
X        &key
X        (versus #'(lambda (i) (float (- i (/ (1- (length y)) 2)))))
X        (weights nil)
X        (compute-residuals-p t)
X        (compute-covariance-matrix-p t)
X        (residuals (if compute-residuals-p (make-array (length y)))))
X  "given
X   [1] y (required)
X       ==> a sequence with values of y_k's, the dependent variable.
X   [2] versus (keyword; equally spaced & centered values)
X       ==> a sequence or function giving values of x_k's,
X           the independent variable
X   [3] weights (keyword; nil)
X       ==> if vector supplied, these are used as weights
X           on the observations (nil implies that all y_k's
X           are equally weighted)
X   [4] compute-residuals-p (keyword; nil)
X       ==> if t, residuals are computed
X   [5] compute-covariance-matrix-p (keyword; nil)
X       ==> if t, calculates and returns 2x2 var/covar matrix
X   [6] residuals (keyword; nil or vector)
X       <== a sequence to be stuffed with the residuals
X           (not used unless compute-residuals-p is true)
fits model y_k = alpha + beta * x_k, and
returns
X   [1] estimate of intercept alpha
X   [2] estimate of slope beta
X   [3] estimate of residual variance
X   [4] residuals (if compute-residuals-p is true)
X   [5] covariance matrix for parameter estimates
X       (if compute-covariance-matrix-p is true)
---
Note: see Chapter 2 of Draper and Smith, 1966"
X  (let ((sum-weights 0.0)
X        (sum-weights-x 0.0)
X        (sum-weights-x-sq 0.0)
X        (sum-weights-y 0.0)
X        (sum-weights-x-y 0.0)
X        (sum-weights-y-sq 0.0)
X        intercept slope det residual-variance loc-w loc-x loc-y
X        (n (length y))
X        (var-covar (if compute-covariance-matrix-p (make-array '(2 2))))
X        (get-x (if (typep versus 'sequence)
X                 #'(lambda (i)
X                     (elt versus i))
X                 versus)))
X    (cond
X     ((arrayp weights)
X      (dotimes (i n)
X        (setf loc-w (elt weights i)
X              loc-x (funcall get-x i)
X              loc-y (elt y i)
X              sum-weights (+ sum-weights loc-w)
X              sum-weights-x (+ sum-weights-x (* loc-w loc-x))
X              sum-weights-x-sq (+ sum-weights-x-sq (* loc-w loc-x loc-x))
X              sum-weights-y (+ sum-weights-y (* loc-w loc-y))
X              sum-weights-x-y (+ sum-weights-x-y (* loc-w loc-x loc-y))
X              sum-weights-y-sq (+ sum-weights-y-sq (* loc-w loc-y loc-y)))))
X     (t
X      (setf sum-weights (float n))
X      (dotimes (i n)
X        (setf loc-x (funcall get-x i)
X              loc-y (elt y i)
X              sum-weights-x (+ sum-weights-x loc-x)
X              sum-weights-x-sq (+ sum-weights-x-sq (* loc-x loc-x))
X              sum-weights-y (+ sum-weights-y loc-y)
X              sum-weights-x-y (+ sum-weights-x-y (* loc-x loc-y))
X              sum-weights-y-sq (+ sum-weights-y-sq (* loc-y loc-y))))))
X    (setf det (- (* sum-weights sum-weights-x-sq)
X                 (* sum-weights-x sum-weights-x)))
X    (cond
X     ((zerop det)
X      (cerror "slope, intercept, residual variance set to 0.0"
X              "determinant in linear least squares is 0")
X      (values 0.0 0.0 0.0))
X     (t
X      (setf intercept (/ (- (* sum-weights-y sum-weights-x-sq)
X                            (* sum-weights-x sum-weights-x-y))
X                         det))
X      (setf slope (/ (- (* sum-weights sum-weights-x-y)
X                        (* sum-weights-x sum-weights-y))
X                     det))
X      (if (arrayp residuals)
X        (dotimes (i n)
X          (setf (elt residuals i)
X                (- (elt y i)
X                   intercept (* slope (funcall get-x i))))))
X      (setf residual-variance (if (< n 3) 0.0
X                                  (/ (- sum-weights-y-sq
X                                        (* sum-weights-y intercept)
X                                        (* sum-weights-x-y slope))
X                                     (- n 2))))
X      (cond
X       ((arrayp var-covar)
X        (setf (aref var-covar 0 0)
X              (/ (* residual-variance sum-weights-x-sq) det))
X        (setf (aref var-covar 0 1)
X              (- (/ (* residual-variance sum-weights-x) det)))
X        (setf (aref var-covar 1 0) (aref var-covar 0 1))
X        (setf (aref var-covar 1 1)
X              (/ (* residual-variance sum-weights) det))))
X      (values intercept slope residual-variance residuals var-covar)))))
X
#|
(weighted-linear-least-squares
X #(1.0 1.9 3.1 3.9 5.1 5.9 7.1 7.9 9.1 9.9)
X :versus
X #(0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5))
;==> 0.9927272727272727       intercept
;    1.9987878787878788       slope
;    0.011121212121210533     residual variance
;    #(0.0072727272727273196 -0.09212121212121216 0.10848484848484863 -0.09090909090909083 0.10969696969696985 -0.08969696969696894 0.11090909090909129 -0.08848484848484794 0.11212121212121229 -0.08727272727272606)
;    #2a((0.003841873278236366 -0.0012132231404956945) (-0.0012132231404956945 5.392102846647531E-4))
|#
X        
;-------------------------------------------------------------------------------
(defun Durbin-Watson-test-statistic
X       (residuals)
X  "given a set of residuals in a seqeunce,
returns the Durbin-Watson test statistic"
X  (when residuals
X    (let ((top 0.0)
X          (i+1 0))
X      (dotimes (i (1- (length residuals))
X                  (/ top (sum-of-squares residuals)))
X        (incf top (expt (- (elt residuals (incf i+1))
X                           (elt residuals i))
X                        2))))))
X
#|
(Durbin-Watson-test-statistic
X #(0.0072727272727273196 -0.09212121212121216 0.10848484848484863 -0.09090909090909083 0.10969696969696985 -0.08969696969696894 0.11090909090909129 -0.08848484848484794 0.11212121212121229 -0.08727272727272606)
X )
;==> 3.7078028238791196
|#
X
;-------------------------------------------------------------------------------
(defun predicted-y-at-x
X       (x
X        intercept
X        slope)
X  "given
X   [1] x (required)
X       ==> a number
X   [2] intercept (required)
X       ==> intercept of a line
X   [3] slope (required)
X       ==> slope of a line
returns
X   [1] intercept + slope*x"
X  (+ intercept (* x slope)))
X
#|
(predicted-y-at-x 3.0 0.9927272727272727  1.9987878787878788)
;==> 6.989090909090908
(predicted-y-at-x 8.0 0.9927272727272727  1.9987878787878788)
;==> 16.983030303030304
|#
X
;-------------------------------------------------------------------------------
(defun var-predicted-y-at-x
X       (x
X        residual-variance
X        covariance-matrix)
X  "given
X   [1] x (required)
X       ==> a number
X   [2] residual-variance (required)
X       ==> the residual variance as returned
X           by weighted-linear-least-squares
X   [3] covariance-matrix (required)
X       ==> covariance matrix for estimated
X           intercept and slope as returned
X           by weighted-linear-least-squares
returns
X   [1] the variance of the predicted value of y
X       at x
---
Note: see Equation (1.4.8), p. 24, Draper and Smith, 1966"
X  (+ residual-variance (var-predicted-mean-at-x x covariance-matrix)))
X
#|
(var-predicted-y-at-x
X 3.0 0.011121212121210533
X #2a((0.003841873278236366 -0.0012132231404956945) (-0.0012132231404956945 5.392102846647531E-4))
X )
;==> 0.012536639118455508
(var-predicted-y-at-x
X 8.0 0.011121212121210533
X #2a((0.003841873278236366 -0.0012132231404956945) (-0.0012132231404956945 5.392102846647531E-4))
X )
;==> 0.030060973370059984
|#
X
;-------------------------------------------------------------------------------
(defun var-predicted-mean-at-x
X       (x
X        covariance-matrix)
X  "given
X   [1] x (required)
X       ==> a number
X   [2] covariance-matrix (required)
X       ==> covariance matrix for estimated
X           intercept and slope as returned
X           by weighted-linear-least-squares
returns
X   [1] the variance of the predicted mean of y
X       at x
---
Note: see page 56, Draper and Smith, 1966"
X  (+ (aref covariance-matrix 0 0)
X     (* 2 x (aref covariance-matrix 0 1))
X     (* x x (aref covariance-matrix 1 1))))
X
#|
(var-predicted-mean-at-x
X 3.0
X #2a((0.003841873278236366 -0.0012132231404956945) (-0.0012132231404956945 5.392102846647531E-4))
X )
;==> 0.0014154269972449754
(var-predicted-mean-at-x
X 8.0
X #2a((0.003841873278236366 -0.0012132231404956945) (-0.0012132231404956945 5.392102846647531E-4))
X )
;==> 0.018939761248849447
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  kernel-pdf-estimation
;;;                 window-width-from-Silverman
;;;  handle a simple case of probability density function estimation.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun kernel-pdf-estimation
X       (a-seq
X        start-x
X        increment-x
X        n-pdf
X        &key
X        (window-width (window-width-from-Silverman
X                       a-seq))
X        (result-pdf (make-array
X                     n-pdf
X                     :initial-element 0.0)
X                    result-pdf-supplied-p)
X        (result-x (make-array n-pdf)))
X  "given
X   [1] a-seq (required)
X       ==> a sequence of real-valued numbers
X   [2] start-x (required)
X       ==> first point at which pdf is
X           to be estimated
X   [3] increment-x (required)
X       ==> increment between points at which pdf is
X           to be estimated (i.e., second point is
X           at start-x + increment-x, third point is
X           at start-x + 2*increment-x, etc)
X   [4] n-pdf (required)
X       ==> number of points at which pdf is
X           to be estimated
X   [5] window-width (keyword; Silverman's formula)
X       ==> window width of kernel pdf estimator
X   [6] result-pdf (keyword; vector of length n-pdf)
X       <== vector to hold pdf estimate
X   [7] result-x (keyword; vector of length n-pdf)
X       <== vector to hold points at which
X           pdf was estimated
computes a pdf estimate using the normal (Gaussian) kernel and
returns
X   [1] the pdf estimate (in result-pdf)
X   [2] the points associated with the pdf estimates (in result-x)
X   [3] the window width"
X  (let* ((n (length a-seq))
X         (current-x start-x)
X         (normalization-factor (* n window-width)))
X    (if result-pdf-supplied-p
X      (fill result-pdf 0.0 :end n-pdf))
X    (dotimes (i n-pdf (values result-pdf result-x window-width))
X      (setf (elt result-x i) current-x)
X      (dotimes (j n)
X        (incf (elt result-pdf i)
X              (standard-normal-pdf
X               (/ (- (elt a-seq j) current-x)
X                  window-width))))
X      (setf (elt result-pdf i)
X            (/ (elt result-pdf i) normalization-factor))
X      (incf current-x increment-x))))
X
#|
;;; Now we repeat the above, but now scale the histogram as a density:
(kernel-pdf-estimation
X '(9 1 9 9 1 8 1 9 6 7 9 3 9)
X 0.5 1.0 11)
;==> #(0.040920458611761634 0.047495773125185846 0.053665548402843676 0.060012914127614225 0.06698186657888049 0.07432092485055371 0.08082908107969788 0.08464308352211156 0.08399125416137254 0.07803398106539552 0.06732806858217144)
;    #(0.5 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5 10.5)
;    3.308586641170241
X
(* (sum (kernel-pdf-estimation
X         '(9 1 9 9 1 8 1 9 6 7 9 3 9)
X         -10.0 0.1 300))
X   0.1)
;==> 0.9996684916122942
|#
X
;-------------------------------------------------------------------------------
(defun window-width-from-Silverman
X       (a-seq)
X  "given a-seq,
returns the window width for kernel pdf estimation
given on page 48, Equation (3.31), of Silverman, 1986"
X  (min (sqrt (sample-variance a-seq))
X       (/ (interquartile-range a-seq) 1.34)))
X
#|
(window-width-from-Silverman '(9 1 9 9 1 8 1 9 6 7 9 3 9))
;==> 3.308586641170241
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  Everything below here consists of internal symbols in the SAPA package
;;;  and should be regarded as "dirty laundry" ...
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defconstant +ppnd-a0+ 2.50662823884d0)
(defconstant +ppnd-a1+ -18.61500062529d0)
(defconstant +ppnd-a2+ 41.39119773534d0)
(defconstant +ppnd-a3+ -25.44106049637d0)
(defconstant +ppnd-b1+ -8.4735109309d0)
(defconstant +ppnd-b2+ 23.08336743743d0)
(defconstant +ppnd-b3+ -21.06224101826d0)
(defconstant +ppnd-b4+ 3.13082909833d0)
(defconstant +ppnd-c0+ -2.78718931138d0)
(defconstant +ppnd-c1+ -2.29796479134d0)
(defconstant +ppnd-c2+ 4.85014127135d0)
(defconstant +ppnd-c3+ 2.32121276858d0)
(defconstant +ppnd-d1+ 3.54388924762d0)
(defconstant +ppnd-d2+ 1.63706781897d0)
X
;-------------------------------------------------------------------------------
(defun as-111 (p)
X  "adapted from AS 111 --- FORTRAN routine ppnd"
X  (assert (and (plusp p) (< p 1.0)))
X  (let ((q (- p 0.5))
X         ppnd r)
X    (cond
X     ((<= (abs q) 0.42)
X      (setf r (* q q))
X      (* q (/ (+ (* (+ (* (+ (* +ppnd-a3+
X                                r)
X                             +ppnd-a2+)
X                          r)
X                       +ppnd-a1+)
X                    r)
X                 +ppnd-a0+)
X              (+ (* (+ (* (+ (* (+ (* +ppnd-b4+
X                                      r)
X                                   +ppnd-b3+)
X                                r)
X                             +ppnd-b2+)
X                          r)
X                       +ppnd-b1+)
X                    r)
X                 1.0))))
X                                     
X     (t
X      (setf r (sqrt (- (log (if (> q 0.0) (- 1.0 p) p))))
X            ppnd (/ (+ (* (+ (* (+ (* +ppnd-c3+
X                                      r)
X                                   +ppnd-c2+)
X                                r)
X                             +ppnd-c1+)
X                          r)
X                       +ppnd-c0+)
X                    (+ (* (+ (* +ppnd-d2+
X                                r)
X                             +ppnd-d1+)
X                          r)
X                       1.0)))
X      (if (< q 0.0) (- ppnd) ppnd)))))
X
;-------------------------------------------------------------------------------
(defconstant +ppchi2-c1+ 0.01)
(defconstant +ppchi2-c2+ 0.222222)
(defconstant +ppchi2-c3+ 0.32)
(defconstant +ppchi2-c4+ 0.4)
(defconstant +ppchi2-c5+ 1.24)
(defconstant +ppchi2-c6+ 2.2)
(defconstant +ppchi2-c7+ 4.67)
(defconstant +ppchi2-c8+ 6.66)
(defconstant +ppchi2-c9+ 6.73)
(defconstant +ppchi2-c10+ 13.32)
(defconstant +ppchi2-c11+ 60.0)
(defconstant +ppchi2-c12+ 70.0)
(defconstant +ppchi2-c13+ 84.0)
(defconstant +ppchi2-c14+ 105.0)
(defconstant +ppchi2-c15+ 120.0)
(defconstant +ppchi2-c16+ 127.0)
(defconstant +ppchi2-c17+ 140.0)
(defconstant +ppchi2-c18+ 175.0)
(defconstant +ppchi2-c19+ 210.0)
(defconstant +ppchi2-c20+ 252.0)
(defconstant +ppchi2-c21+ 264.0)
(defconstant +ppchi2-c22+ 294.0)
(defconstant +ppchi2-c23+ 346.0)
(defconstant +ppchi2-c24+ 420.0)
(defconstant +ppchi2-c25+ 462.0)
(defconstant +ppchi2-c26+ 606.0)
(defconstant +ppchi2-c27+ 672.0)
(defconstant +ppchi2-c28+ 707.0)
(defconstant +ppchi2-c29+ 735.0)
(defconstant +ppchi2-c30+ 889.0)
(defconstant +ppchi2-c31+ 932.0)
(defconstant +ppchi2-c32+ 966.0)
(defconstant +ppchi2-c33+ 1141.0)
(defconstant +ppchi2-c34+ 1182.0)
(defconstant +ppchi2-c35+ 1278.0)
(defconstant +ppchi2-c36+ 1740.0)
(defconstant +ppchi2-c37+ 2520.0)
(defconstant +ppchi2-c38+ 5040.0)
(defconstant +ppchi2-aa+ 0.6931471806d0)
(defconstant +ppchi2-e+ 0.0000005)
(defconstant +ppchi2-pmin+ 0.000002)
(defconstant +ppchi2-pmax+ 0.999998)
X
;-------------------------------------------------------------------------------
(defun as-91 (p degrees-of-freedom)
X  (let ((g (log-of-gamma (/ degrees-of-freedom 2.0)))
X        xx c ch)
X    ;;; p is outside range where approximation is valid ---
X    ;;; might want to change this to a cerror call with the continue option
X    ;;; of just going ahead and accepting a bad approximation ...
X    (assert (and (>= p +ppchi2-pmin+) (<= p +ppchi2-pmax+)))
X    ;;; the degrees of freedom should be positive ...
X    (assert (plusp degrees-of-freedom))
X    (setf xx (* 0.5 degrees-of-freedom)
X          c (- xx 1.0))
X    (cond
X     ((>= degrees-of-freedom (* (- +ppchi2-c5+) (log p)))
X      (setf ch
X            (if (> degrees-of-freedom +ppchi2-c3+)
X              (let* ((x (quantile-of-normal-distribution
X                         p
X                         :accuracy-level :accurate))
X                     (p1 (/ +ppchi2-c2+ degrees-of-freedom))
X                     (temp-ch (* degrees-of-freedom (expt (+ (* x (sqrt p1))
X                                                             (- 1.0 p1))
X                                                          3))))
X                (if (> temp-ch (+ (* +ppchi2-c6+ degrees-of-freedom) 6.0))
X                  (* (- 2.0) (+ (- (log (- 1.0 p))
X                                   (* c (log (* 0.5 temp-ch))))
X                                g))
X                  temp-ch))
X              (as-91-block-1 +ppchi2-c4+ (log (- 1.0 p)) g c)
X              ))
X      (as-91-block-4 ch xx c p g)
X      ;;; label 1 stuff goes here ...
X      )
X     (t
X      (setf ch (expt (* (* p xx)
X                        (exp (+ g (* xx +ppchi2-aa+))))
X                     (/ 1.0 xx)))
X      (if (< ch +ppchi2-e+) ch
X          (as-91-block-4 ch xx c p g))))))
X
;-------------------------------------------------------------------------------
(defun as-91-block-1 (ch a g c)
X  (let* ((q ch)
X         (p1 (+ 1.0 (* ch (+ +ppchi2-c7+ ch))))
X         (p2 (* ch (+ +ppchi2-c9+ (* ch (+ +ppchi2-c8+ ch)))))
X         (temp (+ (- 0.5) (- (/ (+ +ppchi2-c7+ (* 2.0 ch))
X                                  p1)
X                             (/ (+ +ppchi2-c9+
X                                   (* ch
X                                      (+ +ppchi2-aa+ 10
X                                         (* 3.0
X                                            ch))))
X                                p2)))))
X    (decf ch (/ (- 1.0 (* (exp (+ (+ (+ a g)
X                                     (* 0.5 ch))
X                                  (* c +ppchi2-aa+)))
X                          (/ p2 p1)))
X                temp))
X    (if (> (abs (- (/ q ch) 1.0)) +ppchi2-c1+)
X      (as-91-block-1 ch a g c))))
X
;-------------------------------------------------------------------------------
(defun as-91-block-4 (ch xx c p g)
X  (let* ((q ch)
X         (p1 (* 0.5 ch))
X         (p2 (- p (incomplete-gamma xx p1)))
X         (temp (* p2
X                  (exp
X                   (+ (+ (* xx +ppchi2-aa+) g)
X                      (- p1
X                         (* c (log ch)))))))
X         (b (/ temp ch))
X         (a (- (* 0.5 temp)
X               (* b c)))
X         (s1 (/ (+ +ppchi2-c19+
X                   (* a
X                      (+ +ppchi2-c17+
X                         (* a
X                            (+ +ppchi2-c14+
X                               (* a
X                                  (+ +ppchi2-c13+
X                                     (* a
X                                        (+ +ppchi2-c12+
X                                           (* +ppchi2-c11+
X                                              a))))))))))
X                +ppchi2-c24+))
X         (s2 (/ (+ +ppchi2-c24+
X                   (* a
X                      (+ +ppchi2-c29+
X                         (* a
X                            (+ +ppchi2-c32+
X                               (* a
X                                  (+ +ppchi2-c33+
X                                     (* +ppchi2-c35+
X                                        a))))))))
X                +ppchi2-c37+))
X         (s3 (/ (+ +ppchi2-c19+
X                   (* a
X                      (+ +ppchi2-c25+
X                         (* a
X                            (+ +ppchi2-c28+
X                               (* +ppchi2-c31+ a))))))
X                +ppchi2-c37+))
X         (s4 (/ (+ (+ +ppchi2-c20+
X                      (* a
X                         (+ +ppchi2-c27+ (* +ppchi2-c34+ a))))
X                   (* c
X                      (+ +ppchi2-c22+
X                         (* a
X                            (+ +ppchi2-c30+
X                               (* +ppchi2-c36+ a))))))
X                +ppchi2-c38+))
X         (s5 (/ (+ (+ +ppchi2-c13+ (* +ppchi2-c21+ a))
X                   (* c (+ +ppchi2-c18+ (* +ppchi2-c26+ a))))
X                +ppchi2-c37+))
X         (s6 (/ (+ +ppchi2-c15+
X                   (* c (+ +ppchi2-c23+ (* +ppchi2-c16+ c))))
X                +ppchi2-c38+)))
X    (setf ch (+ ch
X                (* temp
X                   (+ 1.0
X                      (- (* (* 0.5
X                               temp)
X                            s1)
X                         (* (* b c)
X                            (- s1
X                               (* b
X                                  (- s2
X                                     (* b
X                                        (- s3
X                                           (* b
X                                              (- s4
X                                                 (* b
X                                                    (- s5
X                                                       (* b
X                                                          s6))))))))))))))))
X    (if (> (abs (- (/ q ch) 1.0)) +ppchi2-e+)
X      (as-91-block-4 ch xx c p g)
X      ch)))
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;; incomplete-gamma is based upon the Numerical Recipes routine gammp
;;; with the following extension:
;;; (1) it returns three values instead of one
;;;     (see discussion on page 162, Numerical Recipes, First Edition)
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun incomplete-gamma (a x)
X  "arguments:
a  --- parameter value > 0
x  --- variable >= 0
returns three values:
value of incomplete gamma function at a and x;
the numerator (little gamma of a and x);
the denominator (gamma of a);
see Section 6.2 of NR"
X  (assert (and (plusp a) (or (plusp x) (zerop x))))
X  (if (< x (1+ a))
X    (series-representation-for-incomplete-gamma a x)
X    (continued-fraction-representation-for-incomplete-gamma a x)))
X
;-------------------------------------------------------------------------------
;;; incomplete-gamma-q is based upon the NR routine gammq with the following
;;; extension:
;;; (1) it returns three values instead of one (see discussion on page 162)
(defun incomplete-gamma-q (a x)
X  "arguments:
a  --- parameter value > 0
x  --- variable >= 0
returns three values:
value of incomplete gamma function at a and x;
the numerator (big gamma of a and x);
the denominator (gamma of a);
see Section 6.2 of NR"
X  (multiple-value-bind (ratio top bottom)
X                       (incomplete-gamma a x)
X    (values  (- 1.0 ratio)
X             (- bottom top)
X             bottom)))
X
;-------------------------------------------------------------------------------
;;; used internally by series-representation-for-incomplete-gamma and
;;;                    continued-fraction-representation-for-incomplete-gamma
(defvar +NR-gser-itmax+ 100)
(defvar +NR-gser-eps+ 3.E-7)
X
;-------------------------------------------------------------------------------
;;; series-representation-for-incomplete-gamma is based upon the NR routine
;;; gser and  is intended to be called from incomplete-gamma --- checking
;;; the validity of a and x is assumed to be already done.
(defun series-representation-for-incomplete-gamma (a x)
X  (let* ((gln (log-of-gamma a))
X         (gamma-a (exp gln))
X         ratio)
X    (if (zerop x)
X      (values 0.0 0.0 gamma-a)
X      (let* ((ap a)
X             (sum (/ a))
X             (del sum))
X        (dotimes (i +NR-gser-itmax+
X                    (error "a too large, +NR-gser-itmax+ too small"))
X          (incf ap)
X          (setf del (/ (* del x) ap))
X          (incf sum del)
X          (if (< (abs del) (* (abs sum) +NR-gser-eps+))
X            (return (values (setf ratio (* sum (exp (- (* a (log x))
X                                                       x
X                                                       gln))))
X                            (* ratio gamma-a)
X                            gamma-a))))))))
X
;-------------------------------------------------------------------------------
;;; continued-fraction-representation-for-incomplete-gamma is based upon
;;; the NR routine gcf (but it returns P(a,x) instead of Q(a,x) and is intended
;;; to be called from incomplete-gamma --- checking the validity of a and x
;;; is assumed to be already done.
(defun continued-fraction-representation-for-incomplete-gamma (a x)
X  (let* ((gln (log-of-gamma a))
X         (gamma-a (exp gln))
X         (gold 0.0)
X         (a0 1.0)
X         (a1 x)
X         (b0 0.0)
X         (b1 1.0)
X         (fac 1.0)
X         an ana anf g ratio)
X    (dotimes (i +NR-gser-itmax+
X                (error "a too large, +NR-gser-itmax+ too small"))
X      (setf an (float (1+ i))
X            ana (- an a)
X            a0 (* fac (+ a1 (* a0 ana)))
X            b0 (* fac (+ b1 (* b0 ana)))
X            anf (* fac an)
X            a1 (+ (* x a0) (* anf a1))
X            b1 (+ (* x b0) (* anf b1)))
X      (cond
X       ((not (zerop a1))
X        (setf fac (/ a1)
X              g (* b1 fac))
X        (if (< (abs (/ (- g gold) g)) +NR-gser-eps+)
X          (return (values (setf ratio (- 1.0
X                                         (* g (exp (- (* a (log x))
X                                                      x
X                                                      gln)))))
X                          (* ratio gamma-a)
X                          gamma-a)))
X        (setf gold g))))))
SHAR_EOF
chmod 0644 basic-statistics.lisp ||
echo 'restore of basic-statistics.lisp failed'
Wc_c="`wc -c < 'basic-statistics.lisp'`"
test 90436 -eq "$Wc_c" ||
	echo 'basic-statistics.lisp: original size 90436, current size' "$Wc_c"
fi
# ============= dft-and-fft.lisp ==============
if test -f 'dft-and-fft.lisp' -a X"$1" != X"-c"; then
	echo 'x - skipping dft-and-fft.lisp (File already exists)'
else
echo 'x - extracting dft-and-fft.lisp (Text)'
sed 's/^X//' << 'SHAR_EOF' > 'dft-and-fft.lisp' &&
;;;-*- Mode: LISP; Package: :SAPA; Syntax: COMMON-LISP -*-
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;
;  dft-and-fft.lisp
;
;  a collection of Lisp functions concerning the discrete Fourier transform ...
;  Note:  before compiling and loading dft-and-fft.lisp,
;         you should compile and load (in the order listed)
;            sapa-package.lisp and utilities.lisp
;
;  6/3/93
;
;  SAPA, Version 1.0; Copyright 1993, Donald B. Percival, All Rights Reserved
;
;  Use and copying of this software and preparation of derivative works
;  based upon this software are permitted.  Any distribution of this
;  software or derivative works must comply with all applicable United
;  States export control laws.
; 
;  This software is made available AS IS, and no warranty -- about the
;  software, its performance, or its conformity to any
;  specification -- is given or implied. 
;
;  Comments about this software can be addressed to dbp@apl.washington.edu
;-------------------------------------------------------------------------------
;;; (compile-file "ccl:SAPA;dft-and-fft.lisp")
;;; (load "ccl:SAPA;dft-and-fft.fasl")
;-------------------------------------------------------------------------------
(if (not (find-package :SAPA))
X  (defpackage :SAPA (:USE :COMMON-LISP)))
X
(in-package :SAPA)
X
(export '(dft!
X          inverse-dft!
X          fft!
X          inverse-fft!
X          dft-chirp!))
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;; The functions dft! and inverse-dft! implement the discrete Fourier
;;; transform and its inverse using, respectively, Equations (110a)
;;; and (111a) of the SAPA book.  They work on vectors of real
;;; or complex-valued numbers.  The vectors can be of any length.
;;; When the length of the input vector is a power of 2, these functions
;;; use fft!, a fast Fourier transform algorithm based upon one in Marple's
;;; 1987 book; when the length is not a power of 2, they use a chirp
;;; transform algorithm (as described in Section 3.10 of the SAPA book).
;;; Both functions return their results in the input vector and hence
;;; whip out whatever its contents were ---  hence we attach
;;; the stylistic "!" to the function names as a reminder
;;; that the function is mangling something given to it.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun dft!
X       (x
X        &key
X        (N (length x))
X        (sampling-time 1.0))
X  "given:
X   [1] x (required)
X       <=> a vector of real or complex-valued numbers
X   [2] N (keyword; length of x)
X       ==> number of points in dft
X   [3] sampling-time (keyword; 1.0)
returns:
X   [1] x, with contents replaced by
X      the discrete Fourier transform of x, namely,
X                            N-1
X      X(n) =  sampling-time SUM x(t) exp(-i 2 pi n t/N)
X                            t=0
---
Note: see Equation (110a) of the SAPA book"
X  (cond
X   ((power-of-2 N)
X    (fft! x :N N))
X   (t
X    (dft-chirp! x :N N)))
X  (when (not (= sampling-time 1.0))
X    (dotimes (i N)
X      (setf (aref x i) (* sampling-time (aref x i)))))
X  x)
X
;-------------------------------------------------------------------------------
(defun inverse-dft!
X       (X
X        &key
X        (N (length X))
X        (sampling-time 1.0))
X  "given:
X   [1] X (required)
X       <=> a vector of real or complex-valued numbers
X   [2] N (keyword; length of )
X       ==> number of points in inverse dft
X   [3] sampling-time (keyword; 1.0)
returns:
X   [1] X, with contents replaced by
X       the inverse discrete Fourier transform of X, namely,
X                                      N-1
X       x(t) =  (N sampling-time)^{-1} SUM X(n) exp(i 2 pi n t/N)
X                                      n=0
---
Note: see Equation (111a) of the SAPA book"
X  (let ((1-over-N*sampling-time (/ (* N sampling-time))))
X    (dotimes (i N)
X      (setf (aref x i) (conjugate (aref x i))))
X    (cond
X     ((power-of-2 N)
X      (fft! x :N N))
X     (t
X      (dft-chirp! x :N N)))
X    (dotimes (i N x)
X      (setf (aref x i) (* 1-over-N*sampling-time
X                          (conjugate (aref x i)))))))
X
#|
(inverse-dft! (dft! (vector 1.0 6.2 pi -7.0)))
;==>
#(#c(1.0 -0.0)
X  #c(6.199999999999999 3.996911308867813E-17)
X  #c(3.141592653589793 -0.0)
X  #c(-7.0 -3.996911308867813E-17))
X
(compare-seqs
X (inverse-dft! (dft! (vector 1.0 6.2 pi -7.0)))
X (vector 1.0 6.2 pi -7.0))
;==> 2.322616018446808E-16
;    8.89077294290045E-16
X
(inverse-dft! (dft! (vector 1.0 6.2 pi -7.0 (- pi))))
;==>
#(#c(1.0000000000000007 -3.552713678800501E-16)
X  #c(6.200000000000002 -2.486899575160351E-15)
X  #c(3.141592653589796 1.0658141036401504E-15)
X  #c(-7.000000000000002 2.131628207280301E-15)
X  #c(-3.1415926535897913 1.7763568394002506E-16))
X
(compare-seqs
X (inverse-dft! (dft! (vector 1.0 6.2 pi -7.0 (- pi))))
X (vector 1.0 6.2 pi -7.0 (- pi)))
;==> 2.2481756875131595E-15
;    3.0561598645713515E-15
X
;;; fft test from Singleton's IEEE article ...
(let ((c-v (vector #C(0.22925607 0.76687502)
X                   #C(0.68317685 0.50919111)
X                   #C(0.87455959 0.64464100)
X                   #C(0.84746840 0.35396343)
X                   #C(0.39889159 0.45709421)
X                   #C(0.23630936 0.13318189)
X                   #C(0.16605222 0.22602680)
X                   #C(0.66245903 0.25021174)
X                   #C(0.61769668 0.26246527)
X                   #C(0.51266762 0.93920734)
X                   #C(0.62402861 0.42238195)
X                   #C(0.93970599 0.28206823)
X                   #C(0.46921754 0.054879178)
X                   #C(0.51983086 0.39682690)
X                   #C(0.11315656 0.60751725)
X                   #C(0.70150672 0.88705479))))
X  (compare-seqs
X   (inverse-fft! (fft! (copy-seq c-v)))
X   c-v))
;==> 1.1741380797503165E-16
;    2.220446049250313E-16
|#
X
;-------------------------------------------------------------------------------
;;; used by the following routine fft!
(defparameter +use-pre-fft-with-cache-p+ t)
X
;-------------------------------------------------------------------------------
;;;  The following code is a Lisp version (with modifications) of
;;;  the FFT routines on pages 54--6 of ``Digital Spectral Analysis with
;;;  Applications'' by Marple, 1987.  The function fft! only works
;;;  for sample sizes which are powers of 2.
;;;  If the vector used as input to fft!  has elements x(0), ..., x(N-1),
;;;  the vector is modified on output to have elements X(0), ..., X(N-1),
;;;  where
;;;
;;;                 N-1
;;;         X(k) =  SUM x(n) exp(-i 2 pi n k/N)
;;;                 n=0
;;;
;;;  If +use-pre-fft-with-cache-p+ is true, the complex exponential table
;;;  for N=2^P is stored in an array of length N pointed to by the P-th element
;;;  of *sapa-cached-pre-fft-arrays*.  The intent here is to speed things up
;;;  and reduce garbage collection (at the expense of extra storage), but the
;;;  timing speed up in MCL is only about 10% using the cache.
(defun fft!
X       (complex-vector
X        &key
X        (N (length complex-vector)))
X  "given:
X   [1] complex-vector (required)
X       <=> a vector of real or complex-valued numbers
X   [2] N (keyword; length of complex-vector)
X       ==> number of points (must be a power of 2)
computes the discrete Fourier transform of complex-vector
using a fast Fourier transform algorithm and
returns:
X   [1] complex-vector, with contents replaced by
X       the discrete Fourier transform of the input, namely,
X                N-1
X       X(n) =   SUM x(t) exp(-i 2 pi n t/N)
X                t=0
---
Note: see Equation (110a) of the SAPA book
X      with the sampling time set to unity"
X  (let ((exponent (power-of-2 n))
X        (W (if +use-pre-fft-with-cache-p+
X             (pre-fft-with-cache n)
X             (pre-fft n)))
X        (MM 1)
X        (LL n)
X        NN JJ
X        c1 c2)
X    (dotimes (k exponent)
X      (setf NN (/ LL 2)
X            JJ MM)
X      (do* ((i 0 (+ i LL))
X            (kk NN (+ i NN)))
X           ((>= i N))
X        (setf c1 (+ (aref complex-vector i) (aref complex-vector kk))
X              (aref complex-vector kk)
X              (- (aref complex-vector i) (aref complex-vector kk))
X              (aref complex-vector i) c1))
X      (cond
X       ((> NN 1)
X        (do ((j 1 (1+ j)))
X            ((>= j NN))
X          (setf c2 (svref W JJ))
X          (do* ((i j (+ i LL))
X                (kk (+ j NN) (+ i NN)))
X               ((>= i N))
X            (setf c1 (+ (aref complex-vector i) (aref complex-vector kk))
X                  (aref complex-vector kk)
X                  (* (- (aref complex-vector i) (aref complex-vector kk))
X                     c2)
X                  (aref complex-vector i) c1))
X          (incf jj MM))
X        (setf LL NN)
X        (setf MM (* MM 2)))))
X    (let ((j 0)
X          (nv2 (/ n 2))
X          k)
X      (dotimes (i (1- N))
X        (if (< i j)
X          (setf c1 (aref complex-vector j)
X                (aref complex-vector j) (aref complex-vector i)
X                (aref complex-vector i) c1))
X        (setf k nv2)
X        (loop
X          (if (> k j) (return))
X          (decf j k)
X          (divf k 2))
X        (incf j k)))
X    complex-vector))
X
;-------------------------------------------------------------------------------
(defun inverse-fft!
X       (complex-vector
X        &key
X        (N (length complex-vector)))
X  "given:
X   [1] complex-vector (required)
X       <=> a vector of real or complex-valued numbers
X   [2] N (keyword; length of )
X       ==> number of points in inverse dft
X           (must be a power of 2)
computes the inverse discrete Fourier transform of complex-vector
using a fast Fourier transform algorithm and
returns:
X   [1] complex-vector, with contents replaced by
X       the inverse discrete Fourier transform of
X       the input X(n), namely,
X                   N-1
X       x(t) =  1/N SUM X(n) exp(i 2 pi n t/N)
X                   n=0
---
Note: see Equation (111a) of the SAPA book
X      with the sampling time set to unity"
X  (dotimes (i N)
X    (setf (aref complex-vector i)
X          (conjugate (aref complex-vector i))))
X  (fft! complex-vector :N N)
X  (dotimes (i N complex-vector)
X    (setf (aref complex-vector i)
X          (/ (conjugate (aref complex-vector i)) N))))
X
;-------------------------------------------------------------------------------
(defun dft-chirp!
X       (complex-vector
X        &key
X        (N (length complex-vector)))
X  "given:
X   [1] complex-vector (required)
X       <=> a vector of real or complex-valued numbers
X   [2] N (keyword; length of complex-vector)
X       ==> number of points (must be a power of 2)
computes the discrete Fourier transform of complex-vector
using a chirp transform algorithm and
returns:
X   [1] complex-vector, with contents replaced by
X       the discrete Fourier transform of the input, namely,
X                N-1
X       X(n) =   SUM x(t) exp(-i 2 pi n t/N)
X                t=0
---
Note: see Equation (110a) of the SAPA book
X      with the sampling time set to unity and
X      also Section 3.10 for a discussion on
X      the chirp transform algorithm"
X  (let* ((N-pot (next-power-of-2 (1- (* 2 N))))
X         (chirp-data (make-array N-pot :initial-element 0.0))
X         (chirp-factors (make-array N-pot :initial-element 0.0))
X         (N-1 (1- N))
X         (pi-over-N (/ pi N))
X         (j 1)
X         (N-pot-1-j (1- N-pot)))
X    ;;; copy data into chirp-data, where it will get multiplied
X    ;;; by the approriate chirp factors ...
X    (copy-vector complex-vector chirp-data :end N)
X    (setf (svref chirp-factors 0) 1.0
X          (aref complex-vector 0) 1.0)
X    ;;; note that, at the end of this dotimes form, complex-vector
X    ;;; contains the complex conjugate of the chirp factors
X    (dotimes (k N-1)
X      (multf (svref chirp-data j)
X             (setf (aref complex-vector j)
X                   (conjugate
X                    (setf (svref chirp-factors j)
X                          (setf (svref chirp-factors N-pot-1-j)
X                                (exp (complex 0.0 (* pi-over-N j j))))))))
X      (incf j)
X      (decf N-pot-1-j))
X    (fft! chirp-data)
X    (fft! chirp-factors)
X    (dotimes (i N-pot)
X      (multf (svref chirp-data i) (svref chirp-factors i)))
X    (inverse-fft! chirp-data)
X    (dotimes (i N complex-vector)
X      (multf (aref complex-vector i) (svref chirp-data i)))))
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  Everything below here consists of internal symbols in the SAPA package
;;;  and should be regarded as "dirty laundry" ...
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun pre-fft
X       (n
X        &key
X        (complex-exp-vector nil))
X  (if (power-of-2 n)
X    (let* ((the-vector (if complex-exp-vector complex-exp-vector
X                           (make-array n)))
X           (s (/ (* 2 pi) n))
X           (c1 (complex (cos s) (- (sin s))))
X           (c2 (complex 1.0 0.0)))
X      (dotimes (i n the-vector)
X        (setf (aref the-vector i) c2
X              c2 (* c2 c1))))))
X
;-------------------------------------------------------------------------------
;;; saves results up to size (expt 2 31) = 2147483648
(defvar *sapa-cached-pre-fft-arrays* (make-array 32 :initial-element nil))
X
;-------------------------------------------------------------------------------
(defun pre-fft-with-cache
X       (n)
X  (declare (special *sapa-cached-pre-fft-arrays*))
X  (let ((the-exponent (power-of-2 n)))
X    (if (numberp the-exponent)
X      (cond
X       ;;; if the following cond clause is true, the effect is
X       ;;; to return a pointer to the cached array
X       ((aref *sapa-cached-pre-fft-arrays* the-exponent))
X       (t
X        (setf (aref *sapa-cached-pre-fft-arrays* the-exponent)
X              (pre-fft n)))))))
X
#|
;-------------------------------------------------------------------------------
;;;  The following code is a Lisp version (with modifications) of
;;;  the FFT routines on pages 36 to 38 of ``Modern Spectrum Estimation:
;;;  Theory and Application'' by Kay, 1988.
;;;  It only works for sample sizes that are powers of 2.
;;;  If the vector used as input to fft! has elements x(0), ..., x(N-1),
;;;  the vector is modified on output to have elements X(0), ..., X(N-1),
;;;  where
;;;
;;;                 N-1
;;;         X(k) =  SUM x(n) exp(-i 2 pi n k/N)
;;;                 n=0
;;;
;;;  In MCL 2.0, it is slower than fft! (based upon Marple's
;;;  fft routine) by about a factor of 4.  This code can be speed up
;;;  by about 40% by creating a cache for the bit reversing operation, but
;;;  this still is not competitive with Marple's routine.
(defun fft-kay!
X       (complex-vector
X        &key
X        (N (length complex-vector)))
X  "given:
X   [1] complex-vector (required)
X       <=> a vector of real or complex-valued numbers
X   [2] N (keyword)
X       ==> length of vector (default = the length of complex-vector)
returns:
X   [1] same complex vector used as input,
X       but with contents altered to contain
X       its Fourier transform"
X  (let ((exponent (power-of-2 N))
X        (scratch (make-array N))
X        (N-1 (1- N))
X        num J L)
X    (declare (dynamic-extent scratch))
X    ;;; This assertion fails if N is not a power of two 
X    (assert exponent
X            ()
X            "vector input to fft! not a power of 2 (~A, ~D)"
X            complex-vector N)
X    (setf (svref scratch 0) (aref complex-vector 0)
X          (svref scratch N-1) (aref complex-vector N-1))
X    ;;; bit reverse the data (see Oppenheim and Schafer, 1989, p. 594)
X    (dotimes (i (1- N-1))
X      (setf num (1+ i)
X            J 0
X            L N)
X      (dotimes (k exponent)
X        (divf L 2)
X        (when (>= num L)
X          (incf j (expt 2 k))
X          (decf num L)))
X      (setf (svref scratch (1+ i)) (aref complex-vector j)))
X    ;;; main computational loop ...
X    (let ((ldft 1)
X          (ndft N)
X          (minus-two-pi (* -2 pi))
X          arg w NP Nq save)
X      (declare (dynamic-extent ldft ndft minus-two-pi
X                               arg w NP Nq save))
X      (dotimes (k exponent)
X        (multf ldft 2)
X        (divf ndft 2)
X        (dotimes (i ndft)
X          (dotimes (j (/ ldft 2))
X            (setf arg (/ (* minus-two-pi j) ldft)
X                  W (complex (cos arg) (sin arg))
X                  NP (+ j (* ldft i))
X                  NQ (+ NP (/ ldft 2))
X                  save (+ (svref scratch NP)
X                          (* W (svref scratch NQ)))
X                  (svref scratch NQ) (- (svref scratch NP)
X                                        (* W (svref scratch NQ)))
X                  (svref scratch NP) save)))))
X    ;;; overright original vector ...
X    (dotimes (i N complex-vector)
X      (setf (aref complex-vector i) (svref scratch i)))))
X
;;; Test case from Kay ...
(fft-kay! #(1.0 #C(0.0 1.0) -1.0 #C(0.0 -1.0)
X            1.0 #C(0.0 1.0) -1.0 #C(0.0 -1.0)))
;;; ... yields ...
#(#c(0.0 0.0) 
X  #c(0.0 0.0) 
X  #c(8.0 2.4492127076447545E-16) 
X  #c(0.0 0.0) 
X  #c(0.0 0.0) 
X  #c(0.0 0.0) 
X  #c(0.0 -2.4492127076447545E-16) 
X  #c(0.0 0.0))
;;; ... which agrees with Kay's result (to within stated precision)
X
;-------------------------------------------------------------------------------
(defun inverse-fft-kay!
X       (complex-vector
X        &key
X        (N (length complex-vector)))
X  "given:
X   [1] complex-vector (required)
X       <=> a vector of real or complex-valued numbers
X   [2] N (keyword)
X       ==> length of vector (default = the length of complex-vector)
returns:
X   [1] same complex vector used as input,
X       but with contents altered to contain
X       its inverse Fourier transform"
X  (dotimes (i N)
X    (setf (aref complex-vector i)
X          (conjugate (aref complex-vector i))))
X  (fft-kay! complex-vector)
X  (dotimes (i N complex-vector)
X    (setf (aref complex-vector i)
X          (/ (conjugate (aref complex-vector i)) N))))
X
;;; 
(let ((c-v #(1.0 #C(0.0 2.0) -3.0 #C(0.0 -4.0)
X             5.0 #C(0.0 6.0) -7.0 #C(0.0 -8.0))))
X  (compare-seqs
X   (inverse-fft-kay!
X    (fft-kay! (copy-seq c-v)))
X   c-v))
;==> 4.440892098500626E-16
;    1.0816271923258514E-16
|#
SHAR_EOF
chmod 0644 dft-and-fft.lisp ||
echo 'restore of dft-and-fft.lisp failed'
Wc_c="`wc -c < 'dft-and-fft.lisp'`"
test 18472 -eq "$Wc_c" ||
	echo 'dft-and-fft.lisp: original size 18472, current size' "$Wc_c"
fi
# ============= examples.lisp ==============
if test -f 'examples.lisp' -a X"$1" != X"-c"; then
	echo 'x - skipping examples.lisp (File already exists)'
else
echo 'x - extracting examples.lisp (Text)'
sed 's/^X//' << 'SHAR_EOF' > 'examples.lisp' &&
;;;-*- Mode: LISP; Package: :CL-USER; Syntax: COMMON-LISP -*-
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;
;  examples.lisp
;
;  examples showing how to repeat certain results in the SAPA book
;  using the functions in the SAPA package; to run these examples,
;  you should first evaluate the "in-package" and "use-package" forms
;  and then all the "defvar" forms, after which you can evaluate the
;  Lisp forms related to the various figures in the SAPA book
;  in any order you choose.
;
;  6/3/93
;
;  SAPA, Version 1.0; Copyright 1993, Donald B. Percival, All Rights Reserved
;
;  Use and copying of this software and preparation of derivative works
;  based upon this software are permitted.  Any distribution of this
;  software or derivative works must comply with all applicable United
;  States export control laws.
; 
;  This software is made available AS IS, and no warranty -- about the
;  software, its performance, or its conformity to any
;  specification -- is given or implied. 
;
;  Comments about this software can be addressed to dbp@apl.washington.edu
;-------------------------------------------------------------------------------
;;; (compile-file "ccl:SAPA;filtering.lisp")
;;; (load "ccl:SAPA;filtering.fasl")
;-------------------------------------------------------------------------------
(in-package :CL-USER)
X
(use-package :SAPA)
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;; We need to define the following special variables in order to run
;;; the test cases in this file.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defvar *rotation-of-earth-data*
X  (vector  71.0  63.0  70.0  88.0  99.0  90.0 110.0 135.0 128.0 154.0 156.0
X          141.0 131.0 132.0 141.0 104.0 136.0 146.0 124.0 129.0 110.0 119.0
X           97.0 115.0 119.0 113.0 115.0 127.0 122.0 147.0 123.0 130.0 149.0
X          183.0 186.0 185.0 189.0 203.0 191.0 217.0 220.0 245.0 213.0 227.0
X          247.0 266.0 221.0 253.0 231.0 236.0 243.0 246.0 238.0 261.0 248.0
X          264.0 270.0 282.0 255.0 286.0 277.0 260.0 261.0 260.0 289.0 297.0
X          319.0 311.0 314.0 309.0 309.0 310.0 315.0 310.0 296.0 264.0 284.0
X          260.0 286.0 271.0 271.0 259.0 279.0 268.0 296.0 280.0 310.0 265.0
X          277.0 257.0 296.0 307.0 306.0 266.0 285.0 279.0 271.0 257.0 270.0
X          232.0))
X
(defvar *centered-rotation-of-earth-data*
X  (x+b *rotation-of-earth-data* (- (sample-mean *rotation-of-earth-data*))))
X
(defvar *the-acvs*)
(defvar *C_h*)
X
(defvar *eigenvalues-NW-4*)
(defvar *eigenspectra-NW-4*)
(defvar *freqs*)
X
(defvar *ar-5-coeffs*)
(defvar *forward-pred-errors*)
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  Some of the examples below recreate analysis on time series as described
;;;  in the SAPA book.  These time series can be obtained from StatLib
;;;  by sending the message
;;;         send sapa from datasets
;;;  (for details, see page xxvii of the SAPA book).  Once these time
;;;  series have been obtained from StatLib, it is necessary to strip
;;;  out the pertinent numbers for each times series, place them into an
;;;  appropriate ASCII file, and then read each file to load the appropriate
;;;  time series values into Lisp.  The actual details of reading these files
;;;  depend on which implementation of Common Lisp is used.  Here we show
;;;  how it is done in Macintosh Common Lisp 2.0 (MCL 2.0) and in Symbolics
;;;  Genera 8.1.1 running on a Symbolics MacIvory model 3.  Evaluation of each
;;;  of the defvar forms below loads in a time series.  For example, in MCL 2.0,
;;;  "ccl:SAPA;ocean-wave" refers to a file named ocean-wave in a folder
;;;  called SAPA that lives in the folder containing the MCL 2.0 application.
;;;  The file "ccl:SAPA;ocean-wave" consists of the 1024 lines from StatLib
;;;  containing the ocean wave time series (one number each line).
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defvar *ocean-wave-data*
X  (let ((the-array (make-array 1024)))
X    (with-open-file (in-stream
X                     #+mcl "ccl:SAPA;ocean-wave"
X                     #+genera "L:>dbp>ocean-wave..newest"
X                     #+allegro "ocean-wave"
X                     :direction :input)
X      (dotimes (index (length the-array) the-array)
X        (setf (svref the-array index)
X              (float (read in-stream)))))))
X
;-------------------------------------------------------------------------------
(defvar *rough-ice-profile-data*
X  (let ((the-array (make-array 1121)))
X    (with-open-file (in-stream
X                     #+mcl "ccl:SAPA;rough-ice-profile"
X                     #+genera "L:>dbp>rough-ice-profile..newest"
X                     #+allegro "rough-ice-profile"
X                     :direction :input)
X      (dotimes (index (length the-array) the-array)
X        (setf (svref the-array index)
X              (float (read in-stream)))))))
X
;-------------------------------------------------------------------------------
(defvar *smooth-ice-profile-data*
X  (let ((the-array (make-array 288)))
X    (with-open-file (in-stream
X                     #+mcl "ccl:SAPA;smooth-ice-profile"
X                     #+genera "L:>dbp>smooth-ice-profile..newest"
X                     #+allegro "smooth-ice-profile"
X                     :direction :input)
X      (dotimes (index (length the-array) the-array)
X        (setf (svref the-array index)
X              (float (read in-stream)))))))
X
;-------------------------------------------------------------------------------
(defvar *Willamette-River-data*
X  (let ((the-array (make-array 395)))
X    (with-open-file (in-stream
X                     #+mcl "ccl:SAPA;Willamette-River"
X                     #+genera "L:>dbp>Willamette-River..newest"
X                     #+allegro "Willamette-River"
X                     :direction :input)
X      (dotimes (index (length the-array) the-array)
X        (setf (svref the-array index)
X              (float (read in-stream)))))))
X
;-------------------------------------------------------------------------------
(defvar *AR-2*
X  (let ((the-array (make-array 1024)))
X    (with-open-file (in-stream
X                     #+mcl "ccl:SAPA;AR-2"
X                     #+genera "L:>dbp>AR-2..newest"
X                     #+allegro "AR-2"
X                     :direction :input)
X      (dotimes (index (length the-array) the-array)
X        (setf (svref the-array index)
X              (float (read in-stream)))))))
X
#|
(svref *ocean-wave-data* 0)            ;==>  477.0
(svref *ocean-wave-data* 1023)         ;==> -113.0
X
(svref *rough-ice-profile-data* 0)     ;==>  -24.6
(svref *rough-ice-profile-data* 1120)  ;==>  -22.5
X
(svref *smooth-ice-profile-data* 0)    ;==> -9.1
(svref *smooth-ice-profile-data* 287)  ;==> -8.9
X
(svref *Willamette-River-data* 0)      ;==> 8.95361
(svref *Willamette-River-data* 394)    ;==> 9.06933
X
(svref *AR-2* 0)                       ;==> 1.619873110842541
(svref *AR-2* 1023)                    ;==> 1.081772237497269
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  Figure 172: filtering of rotation of earth data using simple
;;;              3 point filter with coefficients 1/4, 1/2 and 1/4.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(let* ((filtered-series (filter-time-series
X                         *centered-rotation-of-earth-data*
X                         #(1/4 1/2 1/4)))
X       (residuals (x-y (make-array
X                        (length filtered-series)
X                        :displaced-to *centered-rotation-of-earth-data*
X                        :displaced-index-offset 1)
X                       filtered-series)))
X  (format t "~&Figure 172 ...")
X  (dotimes (i 5)
X    (format t "~&~8,2F ~8,2F ~8,2F"
X            (svref *centered-rotation-of-earth-data* (1+ i))
X            (svref filtered-series i)
X            (svref residuals i)))
X  (format t "~&...")
X  (let ((N-fs (length filtered-series)))
X    (dotimes (i 5)
X      (format t "~&~8,2F ~8,2F ~8,2F"
X              (svref *centered-rotation-of-earth-data* (1+ (+ i (- N-fs 5))))
X              (svref filtered-series (+ i (- N-fs 5)))
X              (svref residuals (+ i (- N-fs 5))))))
X  (format t "~&... Figure 172")
X  (values))
X                        
#|
;;; Here is what is printed when the above form is evaluated:
Figure 172 ...
X -152.73  -148.98    -3.75
X -145.73  -142.98    -2.75
X -127.73  -129.48     1.75
X -116.73  -121.73     5.00
X -125.73  -118.48    -7.25
...
X   69.27    63.02     6.25
X   63.27    62.77     0.50
X   55.27    53.77     1.50
X   41.27    48.02    -6.75
X   54.27    41.52    12.75
... Figure 172
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  Figure 177: squared modulus of transfer function for
;;               the Kth order least squares approximation
;;;              to an ideal low-pass filter
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(let ((65-point-ls-filter
X       (create-least-squares-low-pass-filter 65 0.1)))
X  (format t "~&Figure 177 ...")
X  (dolist (i '(0 1))
X    (format t "~&~3D: ~F" (- i 32) (svref 65-point-ls-filter i)))
X  (format t "~&...")
X  (dolist (i '(31 32 33))
X    (format t "~&~3D: ~F" (- i 32) (svref 65-point-ls-filter i)))
X  (format t "~&...")
X  (dolist (i '(63 64))
X    (format t "~&~3D: ~F" (- i 32) (svref 65-point-ls-filter i)))
X  (format t "~& sum =  ~F" (sum 65-point-ls-filter))
X  (multiple-value-bind (mod-sq-trans-func freqs)
X                       (transfer-function-for-filter
X                        65-point-ls-filter
X                        :return-frequencies-p t)
X    (dolist (i '(0 1 2 3 4 5))
X      (format t "~&~6,4F: ~F"
X              (svref freqs i)
X              (svref mod-sq-trans-func i)))
X    (format t "~&...")
X    (dolist (i '(51 52 53))
X      (format t "~&~6,4F: ~F"
X              (svref freqs i)
X              (svref mod-sq-trans-func i)))
X    (format t "~&...")
X    (dolist (i '(254 255 256))
X      (format t "~&~6,4F: ~F"
X              (svref freqs i)
X              (svref mod-sq-trans-func i))))
X  (format t "~&... Figure 177")
X  (values))
X
#|
;;; Here is what is printed when the above form is evaluated:
Figure 177 ...
-32: 0.009474353340135194
-31: 0.00604435859161769
...
X -1: 0.18737511634014858
X  0: 0.2002963792180472
X  1: 0.18737511634014858
...
X 31: 0.00604435859161769
X 32: 0.009474353340135194
X sum =  1.0000000000000002
0.0000: 1.9286549331065735E-15
0.0020: -9.97589282841272E-4
0.0039: -0.003518706679958886
0.0059: -0.006268132747772469
0.0078: -0.0074626548387638734
0.0098: -0.005314313156846105
...
0.0996: -5.4415656590533095
0.1016: -7.797097385966063
0.1035: -10.914299436697938
...
0.4961: -42.85386826084243
0.4980: -40.44459647657306
0.5000: -39.73441540905796
... Figure 177
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  Figure 178: squared modulus of transfer function for
;;;              the Kth order least squares
;;;              + triangular convergence factors approximation
;;;              to an ideal low-pass filter
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(let ((65-point-ls-filter
X       (create-least-squares-low-pass-filter
X        65 0.1
X        :convergence-factors
X        #'(lambda (k) (triangular-convergence-factors k 65)))))
X  (format t "~&Figure 178 ...")
X  (dolist (i '(0 1))
X    (format t "~&~3D: ~F" (- i 32) (svref 65-point-ls-filter i)))
X  (format t "~&...")
X  (dolist (i '(31 32 33))
X    (format t "~&~3D: ~F" (- i 32) (svref 65-point-ls-filter i)))
X  (format t "~&...")
X  (dolist (i '(63 64))
X    (format t "~&~3D: ~F" (- i 32) (svref 65-point-ls-filter i)))
X  (format t "~& sum =  ~F" (sum 65-point-ls-filter))
X  (multiple-value-bind (mod-sq-trans-func freqs)
X                       (transfer-function-for-filter
X                        65-point-ls-filter
X                        :return-frequencies-p t)
X    (dolist (i '(0 1 2 3 4 5))
X      (format t "~&~6,4F: ~F"
X              (svref freqs i)
X              (svref mod-sq-trans-func i)))
X    (format t "~&...")
X    (dolist (i '(51 52 53))
X      (format t "~&~6,4F: ~F"
X              (svref freqs i)
X              (svref mod-sq-trans-func i)))
X    (format t "~&...")
X    (dolist (i '(254 255 256))
X      (format t "~&~6,4F: ~F"
X              (svref freqs i)
X              (svref mod-sq-trans-func i)))
X    )
X  (format t "~&... Figure 178")
X  (values))
X
#|
;;; Here is what is printed when the above form is evaluated:
Figure 178 ...
-32: 2.958988584176635E-4
-31: 3.7754952616136017E-4
...
X -1: 0.18726456497603497
X  0: 0.20643377319025294
X  1: 0.18726456497603497
...
X 31: 3.7754952616136017E-4
X 32: 2.958988584176635E-4
X sum =  0.9999999999999998
0.0000: -1.9286549331065743E-15
0.0020: 9.19973377882888E-4
0.0039: 0.0035297560721751367
0.0059: 0.007399570616709966
0.0078: 0.011879123556194419
0.0098: 0.016187305184551266
...
0.0996: -5.636216008278292
0.1016: -6.812966030344137
0.1035: -8.143286403474598
...
0.4961: -49.51447433638478
0.4980: -49.421175373585505
0.5000: -49.38851214695436
... Figure 178
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  Figure 181: squared modulus of transfer function for
;;;              dpss approximations to an ideal low-pass filter
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(let ((17-point-dpss
X       (dpss-data-taper!
X        (make-array 17 :initial-element 1.0)
X        :taper-parameter (* 0.1 17)))
X      (65-point-dpss
X       (dpss-data-taper!
X        (make-array 65 :initial-element 1.0)
X        :taper-parameter (* 0.1 65))))
X  (a*x! (/ (sum 17-point-dpss)) 17-point-dpss)
X  (a*x! (/ (sum 65-point-dpss)) 65-point-dpss)
X  (format t "~&Figure 181 ...")
X  (multiple-value-bind (mod-sq-trans-func-17 freqs)
X                       (transfer-function-for-filter
X                        17-point-dpss
X                        :return-frequencies-p t
X                        :N-fft 512)
X    (let ((mod-sq-trans-func-65
X           (transfer-function-for-filter
X            65-point-dpss
X            :return-frequencies-p t
X            :N-fft 512)))
X      (dolist (i '(0 1 2 3 4 5 6 7))
X        (format t "~&~6,4F: ~7,3F, ~7,3F"
X                (svref freqs i)
X                (svref mod-sq-trans-func-17 i)
X                (svref mod-sq-trans-func-65 i)))
X      (format t "~&...")
X      (dolist (i '(50 51 52 53 54))
X        (format t "~&~6,4F: ~7,3F, ~7,3F"
X                (svref freqs i)
X                (svref mod-sq-trans-func-17 i)
X                (svref mod-sq-trans-func-65 i)))
X      (format t "~&...")
X      (dolist (i '(254 255 256))
X        (format t "~&~6,4F: ~7,3F, ~7,3F"
X                (svref freqs i)
X                (svref mod-sq-trans-func-17 i)
X                (svref mod-sq-trans-func-65 i)))
X      ))
X  (format t "~&... Figure 181")
X  (values))
X
#|
;;; Here is what is printed when the above form is evaluated:
Figure 181 ...
0.0000:   0.000,  -0.0000
0.0020:  -0.008,  -0.033
0.0039:  -0.030,  -0.131
0.0059:  -0.068,  -0.294
0.0078:  -0.121,  -0.523
0.0098:  -0.189,  -0.818
0.0117:  -0.272,  -1.179
0.0137:  -0.371,  -1.607
...
0.0977: -26.612, -134.511
0.0996: -28.537, -176.340
0.1016: -30.737, -157.182
0.1035: -33.328, -172.494
0.1055: -36.530, -160.819
...
0.4961: -49.526, -184.556
0.4980: -49.369, -182.150
0.5000: -49.318, -181.440
... Figure 181
;;; Note: the above listing is from Macintosh Common Lisp.  When this form
;;;       is evaluated under Genera 8.1 or Allegro Common Lisp, the last
;;;       column differs somewhat at low dB levels:
0.0977: -26.612, -135.267
0.0996: -28.537, -157.133
0.1016: -30.737, -150.720
0.1035: -33.328, -157.869
0.1055: -36.530, -164.541
...
0.4961: -49.526, -163.929
0.4980: -49.369, -165.401
0.5000: -49.318, -165.964
;;;       This difference is evidently due to different default precisions
;;;       for numerical computations.
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  Figure 182: squared modulus of transfer function for
;;;              approximations to an ideal low-pass filter
;;;              using dpss as convergence factors to least squares
;;;              approximations
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(let ((delta-point-04
X       (create-dpss-low-pass-filter
X        65 0.04 0.1))
X      (delta-point-01
X       (create-dpss-low-pass-filter
X        65 0.01 0.1)))
X  (a*x! (/ (sum delta-point-04)) delta-point-04)
X  (a*x! (/ (sum delta-point-01)) delta-point-01)
X  (format t "~&Figure 182 ...")
X  (print (sum delta-point-04))
X  (print (sum delta-point-01))
X  (multiple-value-bind (mod-sq-trans-func-04 freqs)
X                       (transfer-function-for-filter
X                        delta-point-04
X                        :return-frequencies-p t
X                        :N-fft 512)
X    (let ((mod-sq-trans-func-01
X           (transfer-function-for-filter
X            delta-point-01
X            :return-frequencies-p t
X            :N-fft 512)))
X      (dolist (i '(0 1 2 3 4 5 6))
X        (format t "~&~6,4F: ~7,3F, ~7,3F"
X                (svref freqs i)
X                (svref mod-sq-trans-func-04 i)
X                (svref mod-sq-trans-func-01 i)))
X      (format t "~&...")
X      (dolist (i '(50 51 52 53 54))
X        (format t "~&~6,4F: ~7,3F, ~7,3F"
X                (svref freqs i)
X                (svref mod-sq-trans-func-04 i)
X                (svref mod-sq-trans-func-01 i)))
X      (format t "~&...")
X      (dolist (i '(254 255 256))
X        (format t "~&~6,4F: ~7,3F, ~7,3F"
X                (svref freqs i)
X                (svref mod-sq-trans-func-04 i)
X                (svref mod-sq-trans-func-01 i)))))
X  (format t "~&... Figure 182")
X  (values))
X
#|
;;; Here is what is printed when the above form is evaluated:
Figure 182 ...
0.9999999999999999 
1.0 
0.0000:   0.000,   0.000
0.0020:   0.0001,   0.0005
0.0039:   0.0002,   0.002
0.0059:   0.0004,   0.005
0.0078:   0.001,   0.009
0.0098:   0.001,   0.014
0.0117:   0.001,   0.021
...
0.0977:  -4.944,  -4.095
0.0996:  -5.830,  -5.601
0.1016:  -6.821,  -7.452
0.1035:  -7.923,  -9.733
0.1055:  -9.144, -12.575
...
0.4961: -96.614, -50.080
0.4980: -94.150, -47.668
0.5000: -93.426, -46.957
... Figure 182
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  Figure 200: Fejer's kernel for N = 4, 16 and 64
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(multiple-value-bind (the-spec-wind freqs N-f)
X                     (spectral-window-for-direct-spectral-estimate
X                      4)
X  (format t "~&Figure 200, top plot ...")
X  (dotimes (i N-f)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs i)
X            (svref the-spec-wind  i)))
X  (format t "~&... Figure 200, top plot")
X  (values))
X
#|
;;; Here is what is printed when the above form is evaluated:
Figure 200, top plot ...
0.0000:   6.0206
0.0625:   5.1644
0.1250:   2.3226
0.1875:  -3.9257
0.2500: -100.0000
0.3125:  -7.4278
0.3750:  -5.3329
0.4375:  -8.8624
0.5000: -100.0000
... Figure 200, top plot
|#
X
;;; Here we repeat the above calculation on a much finer grid of frequencies ...
(multiple-value-bind (the-spec-wind freqs N-f)
X                     (spectral-window-for-direct-spectral-estimate
X                      4
X                      :n-nonzero-freqs 512)
X  (format t "~&Figure 200, top plot ...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs i)
X            (svref the-spec-wind  i)))
X  (format t "~&...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs (+ i (- N-f 7)))
X            (svref the-spec-wind  (+ i (- N-f 7)))))
X  (format t "~&... Figure 200, top plot")
X  (values))
X
#|
;;; Here is what is printed when the above form is evaluated:
Figure 200, top plot ...
0.0000:   6.0206
0.0010:   6.0204
0.0020:   6.0198
0.0029:   6.0188
0.0039:   6.0173
0.0049:   6.0155
0.0059:   6.0132
...
0.4941: -28.6858
0.4951: -30.2674
0.4961: -32.2040
0.4971: -34.7016
0.4980: -38.2225
0.4990: -44.2426
0.5000: -100.0000
... Figure 200, top plot
|#
X
(multiple-value-bind (the-spec-wind freqs N-f)
X                     (spectral-window-for-direct-spectral-estimate
X                      16)
X  (format t "~&Figure 200, middle plot ...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs i)
X            (svref the-spec-wind  i)))
X  (format t "~&...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs (+ i (- N-f 7)))
X            (svref the-spec-wind  (+ i (- N-f 7)))))
X  (format t "~&... Figure 200, middle plot")
X  (values))
X
#|
;;; Here is what is printed when the above form is evaluated:
Figure 200, middle plot ...
0.0000:  12.0412
0.0156:  11.1326
0.0312:   8.1328
0.0469:   1.6181
0.0625: -100.0000
0.0781:  -2.7629
0.0937:  -1.2977
...
0.4062: -11.6589
0.4219: -14.7872
0.4375: -100.0000
0.4531: -14.9570
0.4687: -11.9993
0.4844: -15.0410
0.5000: -100.0000
... Figure 200, middle plot
|#
X
(multiple-value-bind (the-spec-wind freqs N-f)
X                     (spectral-window-for-direct-spectral-estimate
X                      64)
X  (format t "~&Figure 200, bottom plot ...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs i)
X            (svref the-spec-wind  i)))
X  (format t "~&...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs (+ i (- N-f 7)))
X            (svref the-spec-wind  (+ i (- N-f 7)))))
X  (format t "~&... Figure 200, bottom plot")
X  (values))
X
#|
;;; Here is what is printed when the above form is evaluated:
Figure 200, bottom plot ...
0.0000:  18.0618
0.0039:  17.1499
0.0078:  14.1403
0.0117:   7.6092
0.0156: -100.0000
0.0195:   3.1758
0.0234:   4.6048
...
0.4766: -18.0382
0.4805: -21.0557
0.4844: -100.0000
0.4883: -21.0662
0.4922: -18.0592
0.4961: -21.0714
0.5000: -100.0000
... Figure 200, bottom plot
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  Figure 211: dpss data tapers and spectral windows for N = 64,
;;;              NW = 1 and 2
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(let ((the-taper (dpss-data-taper!
X                  (make-array 64 :initial-element 1.0)
X                  :taper-parameter 1.0)))
X  (format t "~&sum of squares of taper elements =  ~8,4F"
X          (sum-of-squares the-taper))
X  (format t "~&Figure 211, top left-hand plot ...")
X  (dotimes (i 7)
X    (format t "~&~2D: ~8,4F"
X            (1+ i)
X            (svref the-taper i)))
X  (format t "~&...")
X  (dotimes (i 7)
X    (format t "~&~2D: ~8,4F"
X            (+ (1+ i) (- 64 7))
X            (svref the-taper (+ i (- 64 7)))))
X  (format t "~&... Figure 211, top left-hand plot")
X  (values))
X
#|
;;; Here is what is printed when the above form is evaluated:
sum of squares of taper elements =    1.0000
Figure 211, top left-hand plot ...
X 1:   0.0353
X 2:   0.0404
X 3:   0.0457
X 4:   0.0512
X 5:   0.0568
X 6:   0.0626
X 7:   0.0684
...
58:   0.0684
59:   0.0626
60:   0.0568
61:   0.0512
62:   0.0457
63:   0.0404
64:   0.0353
... Figure 211, top left-hand plot
|#
X
(let ((the-taper (dpss-data-taper!
X                  (make-array 64 :initial-element 1.0)
X                  :taper-parameter 2.0)))
X  (format t "~&sum of squares of taper elements =  ~8,4F"
X          (sum-of-squares the-taper))
X  (format t "~&Figure 211, top right-hand plot ...")
X  (dotimes (i 7)
X    (format t "~&~2D: ~8,4F"
X            (1+ i)
X            (svref the-taper i)))
X  (format t "~&...")
X  (dotimes (i 7)
X    (format t "~&~2D: ~8,4F"
X            (+ (1+ i) (- 64 7))
X            (svref the-taper (+ i (- 64 7)))))
X  (format t "~&... Figure 211, top right-hand plot")
X  (values))
X
#|
;;; Here is what is printed when the above form is evaluated:
sum of squares of taper elements =    1.0000
Figure 211, top right-hand plot ...
X 1:   0.0034
X 2:   0.0055
X 3:   0.0079
X 4:   0.0110
X 5:   0.0146
X 6:   0.0188
X 7:   0.0236
...
58:   0.0236
59:   0.0188
60:   0.0146
61:   0.0110
62:   0.0079
63:   0.0055
64:   0.0034
... Figure 211, top right-hand plot
|#
X
(multiple-value-bind (the-spec-wind freqs N-f)
X                     (spectral-window-for-direct-spectral-estimate
X                      64
X                      :data-taper #'dpss-data-taper!
X                      :data-taper-parameters 1.0)
X  (format t "~&Figure 211, bottom left-hand plot ...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs i)
X            (svref the-spec-wind  i)))
X  (format t "~&...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs (+ i (- N-f 7)))
X            (svref the-spec-wind  (+ i (- N-f 7)))))
X  (format t "~&... Figure 211, bottom left-hand plot")
X  (values))
X
#|
;;; Here is what is printed when the above form is evaluated:
Figure 211, bottom left-hand plot ...
0.0000:  17.4690
0.0039:  16.8738
0.0078:  15.0172
0.0117:  11.6357
0.0156:   6.0027
0.0195:  -4.8372
0.0234: -12.4591
...
0.4766: -29.6551
0.4805: -32.6327
0.4844: -78.3605
0.4883: -32.7070
0.4922: -29.6760
0.4961: -32.6802
0.5000: -100.0000
... Figure 211, bottom left-hand plot
|#
X
(multiple-value-bind (the-spec-wind freqs N-f)
X                     (spectral-window-for-direct-spectral-estimate
X                      64
X                      :data-taper #'dpss-data-taper!
X                      :data-taper-parameters 2.0)
X  (format t "~&Figure 211, bottom right-hand plot ...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs i)
X            (svref the-spec-wind  i)))
X  (format t "~&...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs (+ i (- N-f 7)))
X            (svref the-spec-wind  (+ i (- N-f 7)))))
X  (format t "~&... Figure 211, bottom right-hand plot")
X  (values))
X
#|
;;; Here is what is printed when the above form is evaluated:
Figure 211, bottom right-hand plot ...
0.0000:  16.3413
0.0039:  15.9770
0.0078:  14.8693
0.0117:  12.9709
0.0156:  10.1904
0.0195:   6.3636
0.0234:   1.1838
...
0.4766: -51.8112
0.4805: -54.6679
0.4844: -88.3687
0.4883: -54.9352
0.4922: -51.8305
0.4961: -54.8103
0.5000: -100.0000
... Figure 211, bottom right-hand plot
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  Figure 234: normalized cumulative periodogram for first 32 points
;;;              of AR(2) series
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;; For this example we use the first 64 points of the AR(2) time series
;;; shown in Figure 45 of the SAPA book.
(let* ((32-pt-ts (make-array 32 :displaced-to *AR-2*))
X       (cum-per-results (multiple-value-list
X                         (cumulative-periodogram 32-pt-ts)))
X       (cum-per (elt cum-per-results 0))
X       (cum-per-freqs (elt cum-per-results 1))
X       (KS-results (subseq cum-per-results 2))
X       (the-periodogram (periodogram 32-pt-ts
X                                     :N-nonzero-freqs :Fourier
X                                     :sdf-transformation nil)))
X  (format t "~&Figure 234, left-hand plot, right-hand plot ...")
X  (dotimes (i 5)
X    (format t "~&~6,4F:      ~8,4F         ~8,4F"
X            (svref cum-per-freqs i)
X            (svref the-periodogram i)
X            (svref cum-per i)))
X  (format t "~&...")
X  (let ((N-cum-per (length cum-per)))
X    (dotimes (i 5)
X      (format t "~&~6,4F:      ~8,4F         ~8,4F"
X              (svref cum-per-freqs (+ i (- N-cum-per 5)))
X              (svref the-periodogram (+ i (- N-cum-per 5)))
X              (svref cum-per (+ i (- N-cum-per 5))))))
X  (apply #'format t "~&... Figure 234
Kolmogorov test statistic             = ~A
index of maximum deviation            = ~A
frequency of maximum deviation        = ~A
quantile of Kolmogorov test statistic = ~A
reject/fail to reject null hypothesis = ~A
slope of upper and lower lines        = ~A
intercept of upper line               = ~A
intercept of lower line               = ~A"
X         KS-results)
X  (values))
X
#|
;;; Here is what is printed when the above form is evaluated:
Figure 234, left-hand plot, right-hand plot ...
0.0312:        0.2003           0.0055
0.0625:        2.9395           0.0866
0.0937:        7.8185           0.3022
0.1250:        7.9642           0.5219
0.1562:        5.6365           0.6773
...
0.3437:        2.0928           0.9507
0.3750:        0.6617           0.9689
0.4062:        0.1675           0.9735
0.4375:        0.7901           0.9953
0.4687:        0.1696           1.0000
... Figure 234
Kolmogorov test statistic             = 0.3916144875655616
index of maximum deviation            = 4
frequency of maximum deviation        = 0.15625
quantile of Kolmogorov test statistic = 0.349
reject/fail to reject null hypothesis = reject
slope of upper and lower lines        = 2.2857142857142856
intercept of upper line               = 0.2775714285714286
intercept of lower line               = -0.349
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  Figure 265: Parzen lag window, smoothing window and spectral windows
;;;              using default and dpss data tapers
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(progn
X  (format t "~&Figure 265, top left-hand plot ...")
X  (dotimes (i 7)
X    (format t "~&~2D: ~7,5F"
X            i
X            (parzen-lag-window i 37)))
X  (format t "~&...")
X  (dotimes (i 7)
X    (format t "~&~2D: ~7,5F"
X            (+ i 33)
X            (parzen-lag-window (+ i 33) 37)))
X  (format t "~&... Figure 265, top left-hand plot")
X  (values))
X
#|
;;; Here is what is printed when the above form is evaluated:
X 0: 1.00000
X 1: 0.99574
X 2: 0.98342
X 3: 0.96375
X 4: 0.93746
X 5: 0.90524
X 6: 0.86781
...
33: 0.00253
34: 0.00107
35: 0.00032
36: 0.00004
37: 0.00000
38: 0.00000
39: 0.00000
... Figure 265, top left-hand plot
|#
X
(multiple-value-bind (the-smooth-wind freqs N-f)
X                     (smoothing-window-for-lag-window-spectral-estimate
X                      64
X                      #'(lambda (tau)
X                          (parzen-lag-window tau 37)))
X  (format t "~&Figure 265, top right-hand plot ...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs i)
X            (svref the-smooth-wind  i)))
X  (format t "~&...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs (+ i (- N-f 7)))
X            (svref the-smooth-wind  (+ i (- N-f 7)))))
X  (format t "~&... Figure 265, top right-hand plot")
X  (values))
X
#|
;;; Here is what is printed when the above form is evaluated:
Figure 265, top right-hand plot ...
0.0000:  14.4326
0.0039:  14.2831
0.0078:  13.8316
0.0117:  13.0682
0.0156:  11.9756
0.0195:  10.5271
0.0234:   8.6825
...
0.4766: -46.8807
0.4805: -45.6696
0.4844: -44.6896
0.4883: -44.5244
0.4922: -45.2033
0.4961: -46.3784
0.5000: -47.0461
... Figure 265, top right-hand plot
|#
X
(multiple-value-bind (the-spec-wind freqs N-f)
X                     (spectral-window-for-lag-window-spectral-estimate
X                      64
X                      #'(lambda (tau)
X                          (parzen-lag-window tau 37)))
X  (format t "~&Figure 265, bottom left-hand plot ...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs i)
X            (svref the-spec-wind  i)))
X  (format t "~&...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs (+ i (- N-f 7)))
X            (svref the-spec-wind  (+ i (- N-f 7)))))
X  (format t "~&... Figure 265, bottom left-hand plot")
X  (values))
X
#|
;;; Here is what is printed when the above form is evaluated:
Figure 265, bottom left-hand plot ...
0.0000:  13.8038
0.0039:  13.6754
0.0078:  13.2890
0.0117:  12.6417
0.0156:  11.7287
0.0195:  10.5453
0.0234:   9.0893
...
0.4766: -21.0246
0.4805: -21.0307
0.4844: -21.0353
0.4883: -21.0396
0.4922: -21.0440
0.4961: -21.0475
0.5000: -21.0489
... Figure 265, bottom left-hand plot
|#
X
(multiple-value-bind (the-spec-wind freqs N-f)
X                     (spectral-window-for-lag-window-spectral-estimate
X                      64
X                      #'(lambda (tau)
X                          (parzen-lag-window tau 37))
X                      :data-taper #'dpss-data-taper!
X                      :data-taper-parameters 4.0)
X  (format t "~&Figure 265, bottom right-hand plot ...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs i)
X            (svref the-spec-wind  i)))
X  (format t "~&...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs (+ i (- N-f 7)))
X            (svref the-spec-wind  (+ i (- N-f 7)))))
X  (format t "~&... Figure 265, bottom right-hand plot")
X  (values))
X
#|
;;; Here is what is printed when the above form is evaluated:
Figure 265, bottom right-hand plot ...
0.0000:  13.2324
0.0039:  13.1422
0.0078:  12.8715
0.0117:  12.4186
0.0156:  11.7814
0.0195:  10.9566
0.0234:   9.9399
...
0.4766: -45.0465
0.4805: -45.2264
0.4844: -45.3652
0.4883: -45.4628
0.4922: -45.5266
0.4961: -45.5628
0.5000: -45.5747
... Figure 265, bottom right-hand plot
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  Figure 296: periodogram and other direct spectral estimates
;;;              for ocean wave data
;;;  Note: in what follows we assume that the 1024 values for
;;;        ocean wave data have been placed into the array pointed
;;;        to by *ocean-wave-data* -- at the beginning of this file
;;;        are examples of how to load numbers from an ASCII file
;;;        into a Lisp vector.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(sample-mean-and-variance *ocean-wave-data*)
;==> 209.103515625
;    143954.35647201538
X
(multiple-value-bind (the-periodogram freqs N-f)
X                     (periodogram
X                      *ocean-wave-data*
X                      :sampling-time 1/4)
X  (format t "~&Figure 296, top plot ...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs i)
X            (svref the-periodogram i)))
X  (format t "~&...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs (+ i (- N-f 7)))
X            (svref the-periodogram (+ i (- N-f 7)))))
X  (format t "~&... Figure 296, top plot")
X  (values))
X
#|
;;; Here is what is printed when the above form is evaluated:
Figure 296, top plot ...
0.0039:  47.8344
0.0078:  39.8411
0.0117:  43.0482
0.0156:  43.9995
0.0195:  50.3106
0.0234:  54.7701
0.0273:  54.0996
...
1.9766:  13.2352
1.9805:  11.5889
1.9844:  11.8137
1.9883:  12.3749
1.9922:  14.7094
1.9961:  12.9344
2.0000:  12.2424
... Figure 296, top plot
|#
X
(multiple-value-bind (the-direct-spectral-estimate freqs N-f)
X                     (direct-spectral-estimate
X                      *ocean-wave-data*
X                      :sampling-time 1/4
X                      :data-taper #'dpss-data-taper!
X                      :data-taper-parameters 1.0)
X  (format t "~&Figure 296, second plot ...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs i)
X            (svref the-direct-spectral-estimate i)))
X  (format t "~&...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs (+ i (- N-f 7)))
X            (svref the-direct-spectral-estimate (+ i (- N-f 7)))))
X  (format t "~&... Figure 296, second plot")
X  (values))
X
#|
;;; Here is what is printed when the above form is evaluated:
Figure 296, second plot ...
0.0039:  46.2041
0.0078:  23.2523
0.0117:  42.4188
0.0156:  47.8951
0.0195:  51.3775
0.0234:  54.3370
0.0273:  52.6524
...
1.9766:   2.3258
1.9805:  -3.7210
1.9844:  -0.6480
1.9883:  -1.8448
1.9922:   6.4357
1.9961:  -0.5324
2.0000:  -2.6338
... Figure 296, second plot
|#
X
(multiple-value-bind (the-direct-spectral-estimate freqs N-f)
X                     (direct-spectral-estimate
X                      *ocean-wave-data*
X                      :sampling-time 1/4
X                      :data-taper #'dpss-data-taper!
X                      :data-taper-parameters 2.0)
X  (format t "~&Figure 296, third plot ...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs i)
X            (svref the-direct-spectral-estimate i)))
X  (format t "~&...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs (+ i (- N-f 7)))
X            (svref the-direct-spectral-estimate (+ i (- N-f 7)))))
X  (format t "~&... Figure 296, third plot")
X  (values))
X
#|
;;; Here is what is printed when the above form is evaluated:
Figure 296, third plot ...
0.0039:  43.8021
0.0078:  36.5615
0.0117:  42.9830
0.0156:  49.5223
0.0195:  51.8054
0.0234:  53.3146
0.0273:  50.7869
...
1.9766: -11.0994
1.9805:  -8.6132
1.9844:  -5.4962
1.9883:  -6.0073
1.9922:  -0.1094
1.9961: -14.6655
2.0000: -16.8471
... Figure 296, third plot
|#
X
(multiple-value-bind (the-direct-spectral-estimate freqs N-f)
X                     (direct-spectral-estimate
X                      *ocean-wave-data*
X                      :sampling-time 1/4
X                      :data-taper #'dpss-data-taper!
X                      :data-taper-parameters 4.0)
X  (format t "~&Figure 296, bottom plot ...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs i)
X            (svref the-direct-spectral-estimate i)))
X  (format t "~&...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs (+ i (- N-f 7)))
X            (svref the-direct-spectral-estimate (+ i (- N-f 7)))))
X  (format t "~&... Figure 296, bottom plot")
X  (values))
X
#|
;;; Here is what is printed when the above form is evaluated:
Figure 296, bottom plot ...
0.0039:  40.9897
0.0078:  42.3745
0.0117:  46.2331
0.0156:  50.5033
0.0195:  51.8640
0.0234:  52.3855
0.0273:  50.8260
...
1.9766: -19.6861
1.9805:  -8.4228
1.9844:  -5.0683
1.9883:  -3.3882
1.9922:  -2.3835
1.9961:  -8.3321
2.0000: -22.7452
... Figure 296, bottom plot
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  Figure 298: periodogram for ocean wave data 
;;;              at a finer grid of frequencies ...
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(multiple-value-bind (the-periodogram freqs N-f)
X                     (periodogram
X                      *ocean-wave-data*
X                      :N-nonzero-freqs :next-power-of-2
X                      :sampling-time 1/4)
X  (format t "~&Figure 298 ...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs i)
X            (svref the-periodogram i)))
X  (format t "~&...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs (+ i (- N-f 7)))
X            (svref the-periodogram (+ i (- N-f 7)))))
X  (format t "~&... Figure 298")
X  (values))
X
#|
;;; Here is what is printed when the above form is evaluated:
Figure 298 ...
0.0020:  50.5514
0.0039:  47.8344
0.0059:  37.1331
0.0078:  39.8411
0.0098:  39.0109
0.0117:  43.0482
0.0137:  39.6363
...
1.9883:  12.3749
1.9902:   0.1867
1.9922:  14.7094
1.9941:  -5.7852
1.9961:  12.9344
1.9980: -10.8475
2.0000:  12.2424
... Figure 298
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  Figure 301: lag window spectral estimates for ocean wave data
;;;              based upon direct spectral estimate with NW=2 dpss data taper
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(multiple-value-bind (the-direct-spectral-estimate freqs N-f C_h the-acvs)
X                     (direct-spectral-estimate
X                      *ocean-wave-data*
X                      :sampling-time 1/4
X                      :return-acvs-p t
X                      :data-taper #'dpss-data-taper!
X                      :data-taper-parameters 2.0)
X  (declare (ignore the-direct-spectral-estimate freqs N-f))
X  (format t "~&C_h = ~5,3F" C_h)
X  (format t "~& lag      acvs")
X  (dotimes (i 7)
X    (format t "~&~4D: ~12,4F"
X            i
X            (svref the-acvs i)))
X  (format t "~&...")
X  (dotimes (i 7)
X    (format t "~&~4D: ~12,4F"
X            (+ i (- 1024 7))
X            (svref the-acvs (+ i (- 1024 7)))))
X  (setf *the-acvs* the-acvs
X        *C_h* C_h)
X  (values))
X
#|
;;; Here is what is printed when the above form is evaluated:
C_h = 2.005
X lag      acvs
X   0:  143954.3565
X   1:  138834.3641
X   2:  124334.1234
X   3:  102687.0485
X   4:   76732.5318
X   5:   49174.6362
X   6:   22178.5321
...
1017:      -0.3232
1018:      -0.2984
1019:      -0.2523
1020:      -0.1911
1021:      -0.1258
1022:      -0.0679
1023:      -0.0250
|#
X
(multiple-value-bind (Parzen-m-150 freqs N-f nu B_W)
X                     (lag-window-spectral-estimate
X                      *the-acvs*
X                      #'(lambda (lag)
X                          (parzen-lag-window lag 150))
X                      :sampling-time 1/4
X                      :C_h *C_h*)
X  (format t "~&Figure 301, top plot ...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs i)
X            (svref Parzen-m-150 i)))
X  (format t "~&...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs (+ i (- N-f 7)))
X            (svref Parzen-m-150 (+ i (- N-f 7)))))
X  (format t "~&... Figure 301, top plot")
X  (format t "~&equivalent degrees of freedom = ~5,1F" nu)
X  (format t "~&smoothing window bandwidth    = ~6,4F" B_W)
X  (values))
X
#|
;;; Here is what is printed when the above form is evaluated:
Figure 301, top plot ...
0.0039:  47.5847
0.0078:  48.2507
0.0117:  49.1616
0.0156:  50.1528
0.0195:  51.1108
0.0234:  51.9681
0.0273:  52.6864
...
1.9766:  -5.8887
1.9805:  -5.7822
1.9844:  -5.5922
1.9883:  -5.3665
1.9922:  -5.1577
1.9961:  -5.0117
2.0000:  -4.9595
... Figure 301, top plot
equivalent degrees of freedom =  12.6
smoothing window bandwidth    = 0.0494
|#
X
(multiple-value-bind (Parzen-m-55 freqs N-f nu B_W)
X                     (lag-window-spectral-estimate
X                      *the-acvs*
X                      #'(lambda (lag)
X                          (parzen-lag-window lag 55))
X                      :sampling-time 1/4
X                      :C_h *C_h*)
X  (format t "~&Figure 301, middle plot ...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs i)
X            (svref Parzen-m-55 i)))
X  (format t "~&...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs (+ i (- N-f 7)))
X            (svref Parzen-m-55 (+ i (- N-f 7)))))
X  (format t "~&... Figure 301, middle plot")
X  (format t "~&equivalent degrees of freedom = ~5,1F" nu)
X  (format t "~&smoothing window bandwidth    = ~6,4F" B_W)
X  (values))
X
#|
;;; Here is what is printed when the above form is evaluated:
Figure 301, middle plot ...
0.0039:  51.7292
0.0078:  51.7425
0.0117:  51.7636
0.0156:  51.7913
0.0195:  51.8242
0.0234:  51.8604
0.0273:  51.8985
...
1.9766:  -1.8860
1.9805:  -1.8645
1.9844:  -1.8520
1.9883:  -1.8457
1.9922:  -1.8433
1.9961:  -1.8428
2.0000:  -1.8428
... Figure 301, middle plot
equivalent degrees of freedom =  34.4
smoothing window bandwidth    = 0.1349
|#
X
(multiple-value-bind (Daniell-m-30 freqs N-f nu B_W)
X                     (lag-window-spectral-estimate
X                      *the-acvs*
X                      #'(lambda (lag)
X                          (daniell-lag-window lag 30))
X                      :sampling-time 1/4
X                      :C_h *C_h*)
X  (format t "~&Figure 301, bottom plot ...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs i)
X            (svref Daniell-m-30 i)))
X  (format t "~&...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs (+ i (- N-f 7)))
X            (svref Daniell-m-30 (+ i (- N-f 7)))))
X  (format t "~&... Figure 301, bottom plot")
X  (format t "~&equivalent degrees of freedom = ~5,1F" nu)
X  (format t "~&smoothing window bandwidth    = ~6,4F" B_W)
X  (values))
X
#|
;;; Here is what is printed when the above form is evaluated:
Figure 301, bottom plot ...
0.0039:  52.1040
0.0078:  52.2497
0.0117:  52.4965
0.0156:  52.5851
0.0195:  52.4657
0.0234:  52.0096
0.0273:  51.6623
...
1.9766:  -6.0432
1.9805:  -5.9070
1.9844:  -5.8789
1.9883:  -5.7424
1.9922:  -5.4674
1.9961:  -5.2132
2.0000:  -5.1368
... Figure 301, bottom plot
equivalent degrees of freedom =  34.1
smoothing window bandwidth    = 0.1337
;;; Note: page 302 of the SAPA book gives the degrees of freedom
;;;       for these three cases as 13, 35 and 35 -- the above
;;;       calculations differs slightly because we have obtained
;;;       C_h using Equation (251b) rather from Table 248.
|#
X
;;; Here we show how to compute and relocate the spectral window
;;; shown in the lower left-hand corner of the top plot ...
(multiple-value-bind (the-smooth-wind freqs N-f)
X                     (smoothing-window-for-lag-window-spectral-estimate
X                      1024
X                      #'(lambda (tau)
X                          (parzen-lag-window tau 150))
X                      :sampling-time 1/4
X                      :N-nonzero-freqs :half-next-power-of-2)
X  (let ((dB-offset (+ (- (svref the-smooth-wind 0)) -10.0)))
X    (setf the-smooth-wind (one-sided-sdf->two-sided-sdf the-smooth-wind)
X          freqs (one-sided-freq->two-sided-freq freqs)
X          N-f (1- (* 2 N-f)))
X    (x+b! the-smooth-wind dB-offset)
X    (x+b! freqs 0.25)
X  (dotimes (i N-f (values))
X    (if (<= 0.2 (svref freqs i) 0.3)
X      (format t "~&~6,4F: ~8,4F"
X              (svref freqs i)
X              (svref the-smooth-wind i))))))
X
#|
;;; Here is what is printed when the above form is evaluated:
0.2031: -44.8539
0.2070: -35.7970
0.2109: -29.5888
0.2148: -24.9154
0.2187: -21.2410
0.2227: -18.2942
0.2266: -15.9187
0.2305: -14.0171
0.2344: -12.5259
0.2383: -11.4022
0.2422: -10.6175
0.2461: -10.1536
0.2500: -10.0000
0.2539: -10.1536
0.2578: -10.6175
0.2617: -11.4022
0.2656: -12.5259
0.2695: -14.0171
0.2734: -15.9187
0.2773: -18.2942
0.2812: -21.2410
0.2852: -24.9154
0.2891: -29.5888
0.2930: -35.7970
0.2969: -44.8539
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  Figure 305: periodograms for rough and smooth ice profile data
;;;  Note: in what follows we assume that the 1121 values for
;;;        rough ice profile data have been placed into the array pointed
;;;        to by *rough-ice-profile-data* and that the 288 values for
;;;        smooth ice profile data have been placed into the array pointed
;;;        to by *smooth-ice-profile-data* -- at the beginning of this file
;;;        are examples of how to load numbers from an ASCII file
;;;        into a Lisp vector.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(sample-mean *rough-ice-profile-data*)   ;==> -19.916859946476382
X
(multiple-value-bind (the-periodogram freqs N-f)
X                     (periodogram
X                      *rough-ice-profile-data*
X                      :sampling-time 1.7712)
X  (format t "~&Figure 305, top plot ...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs i)
X            (svref the-periodogram i)))
X  (format t "~&...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs (+ i (- N-f 7)))
X            (svref the-periodogram (+ i (- N-f 7)))))
X  (format t "~&... Figure 305, top plot")
X  (values))
X
#|
;;; Here is what is printed when the above form is evaluated:
Figure 305, top plot ...
0.0003:  19.3206
0.0006:  29.8243
0.0008:  33.4726
0.0011:  33.9632
0.0014:  34.7072
0.0017:  35.0793
0.0019:  36.2898
...
0.2806:  -0.3305
0.2809:   0.2546
0.2812:  -1.4455
0.2815:   2.7054
0.2817:   4.2203
0.2820:  -5.3819
0.2823:  -0.1654
... Figure 305, top plot
|#
X
(sample-mean *smooth-ice-profile-data*)   ;==> -8.934027777777777
X
(multiple-value-bind (the-periodogram freqs N-f)
X                     (periodogram
X                      *smooth-ice-profile-data*
X                      :sampling-time 1.7712)
X  (format t "~&Figure 305, bottom plot ...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs i)
X            (svref the-periodogram i)))
X  (format t "~&...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs (+ i (- N-f 7)))
X            (svref the-periodogram (+ i (- N-f 7)))))
X  (format t "~&... Figure 305, bottom plot")
X  (values))
X
#|
;;; Here is what is printed when the above form is evaluated:
Figure 305, bottom plot ...
0.0011:   5.7414
0.0022:  -0.3703
0.0033:   0.9999
0.0044:   4.3387
0.0055:   7.5413
0.0066:   4.8875
0.0077:   1.8933
...
0.2757: -11.9119
0.2768: -17.4537
0.2779: -22.4416
0.2790: -15.1718
0.2801:  -9.8917
0.2812: -10.6674
0.2823: -16.0906
... Figure 305, bottom plot
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  Figure 307: lag window spectral estimates for rough and smooth
;;;              ice profile data
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(let ((rough-acvs (acvs *rough-ice-profile-data*)))
X  (multiple-value-bind (Parzen-m-88 freqs N-f nu B_W)
X                       (lag-window-spectral-estimate
X                        rough-acvs
X                        #'(lambda (lag)
X                            (parzen-lag-window lag 88))
X                        :sampling-time 1.7712)
X    (format t "~&Figure 307, top plot ...")
X    (dotimes (i 7)
X      (format t "~&~6,4F: ~8,4F"
X              (svref freqs i)
X              (svref Parzen-m-88 i)))
X    (format t "~&...")
X    (dotimes (i 7)
X      (format t "~&~6,4F: ~8,4F"
X              (svref freqs (+ i (- N-f 7)))
X              (svref Parzen-m-88 (+ i (- N-f 7)))))
X    (format t "~&... Figure 307, top plot")
X    (format t "~&equivalent degrees of freedom = ~5,1F" nu)
X    (format t "~&smoothing window bandwidth    = ~6,4F" B_W)
X    (format t "~&time series bandwidth         = ~6,4F"
X            (time-series-bandwidth rough-acvs :sampling-time 1.7712))))
#|
;;; Here is what is printed when the above form is evaluated:
Figure 307, top plot ...
0.0003:  32.2615
0.0006:  32.2518
0.0008:  32.2357
0.0011:  32.2130
0.0014:  32.1838
0.0017:  32.1480
0.0019:  32.1056
...
0.2806:   0.6920
0.2809:   0.6401
0.2812:   0.5964
0.2815:   0.5615
0.2817:   0.5360
0.2820:   0.5206
0.2823:   0.5154
... Figure 307, top plot
equivalent degrees of freedom =  47.2
smoothing window bandwidth    = 0.0119
time series bandwidth         = 0.0239
|#
X
(let ((smooth-acvs (acvs *smooth-ice-profile-data*)))
X  (multiple-value-bind (Parzen-m-20 freqs N-f nu B_W)
X                       (lag-window-spectral-estimate
X                        smooth-acvs
X                        #'(lambda (lag)
X                            (parzen-lag-window lag 20))
X                        :sampling-time 1.7712)
X    (format t "~&Figure 307, middle plot ...")
X    (dotimes (i 7)
X      (format t "~&~6,4F: ~8,4F"
X              (svref freqs i)
X              (svref Parzen-m-20 i)))
X    (format t "~&...")
X    (dotimes (i 7)
X      (format t "~&~6,4F: ~8,4F"
X              (svref freqs (+ i (- N-f 7)))
X              (svref Parzen-m-20 (+ i (- N-f 7)))))
X    (format t "~&... Figure 307, middle plot")
X    (format t "~&equivalent degrees of freedom = ~5,1F" nu)
X    (format t "~&smoothing window bandwidth    = ~6,4F" B_W)
X    (format t "~&time series bandwidth         = ~6,4F"
X            (time-series-bandwidth smooth-acvs :sampling-time 1.7712)))
X  (values))
X
#|
;;; Here is what is printed when the above form is evaluated:
Figure 307, middle plot ...
0.0011:   3.2732
0.0022:   3.2655
0.0033:   3.2527
0.0044:   3.2349
0.0055:   3.2120
0.0066:   3.1841
0.0077:   3.1513
...
0.2757: -13.2314
0.2768: -13.2507
0.2779: -13.2666
0.2790: -13.2789
0.2801: -13.2877
0.2812: -13.2930
0.2823: -13.2948
... Figure 307, middle plot
equivalent degrees of freedom =  53.4
smoothing window bandwidth    = 0.0523
time series bandwidth         = 0.1069
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  Figure 312: WOSA spectral estimates of AR(2) series
;;;  Note: in what follows we assume that the 1024 values for AR(20 series
;;;        have been placed into a vector pointed to by *AR-2* -- at the
;;;        beginning of this file are examples of how to load numbers
;;;        from an ASCII file into a Lisp vector.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;; block size = 8 points
(multiple-value-bind (wosa freqs dof)
X                     (wosa-spectral-estimate
X                      *AR-2*
X                      8
X                      :oversampling-factor 8)
X  (format t "~&Figure 307, top plot ...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs i)
X            (svref wosa i)))
X  (format t "~&...")
X  (let ((N-f (length freqs)))
X    (dotimes (i 7)
X      (format t "~&~6,4F: ~8,4F"
X              (svref freqs (+ i (- N-f 7)))
X              (svref wosa (+ i (- N-f 7))))))
X  (format t "~&... Figure 307, top plot")
X  (format t "~&equivalent degrees of freedom = ~5,1F" dof)
X  (values))
X
#|
;;; Here is what is printed when the above form is evaluated:
Figure 307, top plot ...
0.0000:   3.1104
0.0156:   3.1860
0.0312:   3.3992
0.0469:   3.7143
0.0625:   4.0842
0.0781:   4.4621
0.0937:   4.8075
...
0.4062:  -4.2055
0.4219:  -4.5725
0.4375:  -4.8838
0.4531:  -5.1387
0.4687:  -5.3309
0.4844:  -5.4517
0.5000:  -5.4930
... Figure 307, top plot
equivalent degrees of freedom = 453.2
|#
X
;;; block size = 16 points
(multiple-value-bind (wosa freqs dof)
X                     (wosa-spectral-estimate
X                      *AR-2*
X                      16
X                      :oversampling-factor 8)
X  (format t "~&Figure 307, middle plot ...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs i)
X            (svref wosa i)))
X  (format t "~&...")
X  (let ((N-f (length freqs)))
X    (dotimes (i 7)
X      (format t "~&~6,4F: ~8,4F"
X              (svref freqs (+ i (- N-f 7)))
X              (svref wosa (+ i (- N-f 7))))))
X  (format t "~&... Figure 307, middle plot")
X  (format t "~&equivalent degrees of freedom = ~5,1F" dof)
X  (values))
X
#|
;;; Here is what is printed when the above form is evaluated:
Figure 307, middle plot ...
0.0000:   2.7104
0.0078:   2.7175
0.0156:   2.7407
0.0234:   2.7860
0.0312:   2.8624
0.0391:   2.9807
0.0469:   3.1503
...
0.4531:  -5.7505
0.4609:  -5.9087
0.4687:  -6.0551
0.4766:  -6.1807
0.4844:  -6.2774
0.4922:  -6.3384
0.5000:  -6.3593
... Figure 307, middle plot
equivalent degrees of freedom = 233.8
|#
X
;;; block size = 64 points
(multiple-value-bind (wosa freqs dof)
X                     (wosa-spectral-estimate
X                      *AR-2*
X                      64
X                      :oversampling-factor 8)
X  (format t "~&Figure 307, bottom plot ...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs i)
X            (svref wosa i)))
X  (format t "~&...")
X  (let ((N-f (length freqs)))
X    (dotimes (i 7)
X      (format t "~&~6,4F: ~8,4F"
X              (svref freqs (+ i (- N-f 7)))
X              (svref wosa (+ i (- N-f 7))))))
X  (format t "~&... Figure 307, bottom plot")
X  (format t "~&equivalent degrees of freedom = ~5,1F" dof)
X  (values))
X
#|
;;; Here is what is printed when the above form is evaluated:
Figure 307, bottom plot ...
0.0000:   2.6662
0.0020:   2.6924
0.0039:   2.7645
0.0059:   2.8654
0.0078:   2.9718
0.0098:   3.0595
0.0117:   3.1077
...
0.4883:  -8.0203
0.4902:  -8.1908
0.4922:  -8.3616
0.4941:  -8.5197
0.4961:  -8.6494
0.4980:  -8.7351
0.5000:  -8.7650
... Figure 307, bottom plot
equivalent degrees of freedom =  58.5
;;; Note: the above values for the equivalent degrees of freedom
;;;       differ from what is stated on page 311 of the SAPA book
;;;       (the book has 483.3, 240.7 and 58.8 rather than 453.2, 233.8 and 58.5).
;;;       This discrepancy is due to the fact that wosa-spectral-estimate
;;;       uses Equation (292b), whereas the values quoted in the book
;;;       are based upon Equation (294) (an approximation to Equation (292b)).
;;;       Also, the spectra plotted in Figure 312 were actually computed
;;;       using a slightly different definition for the Hanning data taper
;;;       (essentially the one in Bloomfield's book).  The spectra computed
;;;       above thus are slightly different from what are plotted in the SAPA
;;;       book (the biggest changes are a 0.5 dB difference at f = 0 and
;;;       f = 0.5 for the 8 point block size, i.e., the top plot of Figure 312).
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  Figure 373: multitaper spectral estimates for ocean wave data
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(multiple-value-bind (dpss-NW-4 eigenvalues-NW-4)
X                     (dpss-tapers-tri-diag
X                      (length *ocean-wave-data*)
X                      7
X                      :taper-parameter 4.0
X                      :compute-true-eigenvalues-p t)
X  (setf *eigenvalues-NW-4* eigenvalues-NW-4)
X  (format t "~&order   eigenvalue")
X  (dotimes (k 7)
X    (format t "~&k = ~1D:  ~F" k (svref eigenvalues-NW-4 k)))
X  ;;; We compute 7 eigenspectra but use only 6 of them to form
X  ;;; the simple multitaper spectral estimate in plot (a).
X  ;;; Note: for Figure 373, the eigenspectra are computed by setting
X  ;;;       the keyword options recenter-after-tapering-p and
X  ;;;       restore-power-option-p to nil.  This produces a spectral
X  ;;;       estimate in agreement with the equations in Chapter 7,
X  ;;;       but arguably it might be better to just leave these keywords
X  ;;;       set to their default values of t.
X  (multiple-value-bind
X    (simple-multitaper-NW-4-0-to-6 freqs N-f eigenspectra-NW-4)
X    (multitaper-spectral-estimate
X     *ocean-wave-data*
X     dpss-NW-4
X     :sampling-time 1/4
X     :recenter-after-tapering-p nil
X     :restore-power-option-p nil)
X    (declare (ignore simple-multitaper-NW-4-0-to-6))
X    (setf *eigenspectra-NW-4* eigenspectra-NW-4
X          *freqs* freqs)
X    (let ((simple-multitaper-NW-4-0-to-5
X           (eigenspectra->multitaper-spectral-estimate
X            eigenspectra-NW-4 :N-eigenspectra 6)))
X      (format t "~&Figure 373, plot (a) ...")
X      (dotimes (i 7)
X        (format t "~&~6,4F: ~8,4F"
X                (svref freqs i)
X                (svref simple-multitaper-NW-4-0-to-5 i)))
X      (format t "~&...")
X      (dotimes (i 7)
X        (format t "~&~6,4F: ~8,4F"
X                (svref freqs (+ i (- N-f 7)))
X                (svref simple-multitaper-NW-4-0-to-5 (+ i (- N-f 7)))))
X      (format t "~&... Figure 373, plot (a)")
X      (values))))
X
#|
;;; Here is what is printed when the above form is evaluated:
order   eigenvalue
k = 0:  0.9999999997056523
k = 1:  0.9999999723287881
k = 2:  0.9999987902598706
k = 3:  0.9999675626065103
k = 4:  0.9994101803916527
k = 5:  0.9925053051988856
k = 6:  0.9366554082073846
Figure 373, plot (a) ...
0.0039:  45.0265
0.0078:  46.7926
0.0117:  48.5376
0.0156:  50.5297
0.0195:  51.1797
0.0234:  53.3410
0.0273:  53.9492
...
1.9766:  -5.0610
1.9805:  -4.9186
1.9844:  -4.7007
1.9883:  -3.5230
1.9922:  -3.0773
1.9961:  -3.3919
2.0000:  -3.3314
... Figure 373, plot (a)
|#
X
(multiple-value-bind (adaptive-multitaper-NW-4 dof)
X                     (eigenspectra->adaptive-multitaper-spectral-estimate
X                      *eigenspectra-NW-4*
X                      *eigenvalues-NW-4*
X                      (sample-variance *ocean-wave-data*)
X                      :sampling-time 1/4)
X  (multiple-value-bind (upper lower)
X                       (create-ci-for-amt-sdf-estimate
X                        adaptive-multitaper-NW-4
X                        dof)
X    (format t "~&Figure 373, plot (c) ...               plot (d) ...")
X    (dotimes (i 7)
X      (format t "~&~6,4F: ~8,4F ~8,4F ~8,4F    ~8,4F"
X              (svref *freqs* i)
X              (svref upper i)
X              (svref adaptive-multitaper-NW-4 i)
X              (svref lower i)
X              (svref dof i)))
X    (format t "~&...")
X    (let ((N-f (length *freqs*)))
X      (dotimes (i 7)
X        (format t "~&~6,4F: ~8,4F ~8,4F ~8,4F    ~8,4F"
X                (svref *freqs* (+ i (- N-f 7)))
X                (svref upper (+ i (- N-f 7)))
X                (svref adaptive-multitaper-NW-4 (+ i (- N-f 7)))
X                (svref lower (+ i (- N-f 7)))
X                (svref dof (+ i (- N-f 7))))))
X    (format t "~&... Figure 373, plots (c) and (d)")
X    (values)))
X
#|
;;; Here is what is printed when the above form is evaluated:
Figure 373, plot (c) ...               plot (d) ...
0.0039:  48.5041  44.5444  41.8350     13.9860
0.0078:  50.9711  47.0137  44.3053     13.9985
0.0117:  53.7760  49.8188  47.1104     13.9996
0.0156:  54.1360  50.1787  47.4704     13.9994
0.0195:  55.6164  51.6589  48.9505     13.9982
0.0234:  56.6806  52.7230  50.0145     13.9973
0.0273:  57.6636  53.7058  50.9973     13.9966
...
1.9766:   1.0163  -5.8131  -9.6216      6.0346
1.9805:  -0.5008  -7.4550 -11.3023      5.8786
1.9844:   2.9011  -3.7747  -7.5345      6.2375
1.9883:   2.5575  -4.1425  -7.9100      6.2047
1.9922:   3.3604  -3.2706  -7.0160      6.2991
1.9961:   2.2524  -4.4749  -8.2511      6.1681
2.0000:   4.0624  -2.5117  -6.2389      6.3789
... Figure 373, plots (c) and (d)
|#
X
(multiple-value-bind (simple-multitaper-NW-6-0-to-9 freqs N-f)
X                     (multitaper-spectral-estimate
X                      *ocean-wave-data*
X                      (dpss-tapers-tri-diag
X                       (length *ocean-wave-data*)
X                       10
X                       :taper-parameter 6.0)
X                      :sampling-time 1/4
X                      :recenter-after-tapering-p nil
X                      :restore-power-option-p nil)
X  (format t "~&Figure 373, plot (b) ...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs i)
X            (svref simple-multitaper-NW-6-0-to-9 i)))
X  (format t "~&...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs (+ i (- N-f 7)))
X            (svref simple-multitaper-NW-6-0-to-9 (+ i (- N-f 7)))))
X  (format t "~&... Figure 373, plot (b)")
X  (values))
X
#|
Figure 373, plot (b) ...
0.0039:  48.0430
0.0078:  49.0666
0.0117:  49.7126
0.0156:  51.5444
0.0195:  51.9033
0.0234:  53.5943
0.0273:  53.7669
...
1.9766:  -4.2660
1.9805:  -4.8461
1.9844:  -4.4768
1.9883:  -4.8493
1.9922:  -3.9439
1.9961:  -2.6031
2.0000:  -1.7041
... Figure 373, plot (b)
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  Figures 440-1: designing a prewhitening filter for ocean wave data
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;; estimation of AR(5) model using the Yule-Walker method,
;;; from which we compute a parametric sdf estimate ...
(multiple-value-bind (ar-5-coeffs innovations-variance)
X                     (yule-walker-algorithm-given-data
X                      *ocean-wave-data*
X                      5)
X  (multiple-value-bind (Y-W-5-sdf freqs N-f)
X                       (ar-coeffs->sdf ar-5-coeffs
X                                       innovations-variance
X                                       :sampling-time 1/4
X                                       :return-frequencies-p t)
X    (format t "~&Figure 440, thin curve, top plot ...")
X    (dotimes (i 7)
X      (format t "~&~6,4F: ~8,4F"
X              (svref freqs i)
X              (svref Y-W-5-sdf i)))
X    (format t "~&...")
X    (dotimes (i 7)
X      (format t "~&~6,4F: ~8,4F"
X              (svref freqs (+ i (- N-f 7)))
X              (svref Y-W-5-sdf (+ i (- N-f 7)))))
X    (format t "~&... Figure 440, thin curve, top plot")
X    (values)))
X
#|
;;; Here is what is printed when the above form is evaluated:
Figure 440, thin curve, top plot ...
0.0000:  51.2480
0.0078:  51.2659
0.0156:  51.3197
0.0234:  51.4097
0.0312:  51.5365
0.0391:  51.7008
0.0469:  51.9034
...
1.9531:  10.5844
1.9609:  10.5978
1.9687:  10.6089
1.9766:  10.6175
1.9844:  10.6237
1.9922:  10.6275
2.0000:  10.6287
... Figure 440, thin curve, top plot
|#
X
;;; estimation of AR(27) model using Burg's algorithm,
;;; from which we compute a parametric sdf estimate ...
(multiple-value-bind (ar-27-coeffs innovations-variance)
X                     (Burg-algorithm
X                      *ocean-wave-data*
X                      27)
X  (multiple-value-bind (Burg-27-sdf freqs N-f)
X                       (ar-coeffs->sdf ar-27-coeffs
X                                       innovations-variance
X                                       :sampling-time 1/4
X                                       :return-frequencies-p t)
X    (format t "~&Figure 440, thin curve, bottom plot ...")
X    (dotimes (i 7)
X      (format t "~&~6,4F: ~8,4F"
X              (svref freqs i)
X              (svref Burg-27-sdf i)))
X    (format t "~&...")
X    (dotimes (i 7)
X      (format t "~&~6,4F: ~8,4F"
X              (svref freqs (+ i (- N-f 7)))
X              (svref Burg-27-sdf (+ i (- N-f 7)))))
X    (format t "~&... Figure 440, thin curve, bottom plot")
X    (values)))
X
#|
;;; Here is what is printed when the above form is evaluated:
Figure 440, thin curve, bottom plot ...
0.0000:  51.8188
0.0078:  51.8771
0.0156:  52.0414
0.0234:  52.2796
0.0312:  52.5356
0.0391:  52.7342
0.0469:  52.7986
...
1.9531:  -6.2274
1.9609:  -5.8871
1.9687:  -5.5495
1.9766:  -5.2402
1.9844:  -4.9885
1.9922:  -4.8231
2.0000:  -4.7653
... Figure 440, thin curve, bottom plot
|#
X
;;; estimation of AR(5) model using Burg's algorithm,
;;; from which we [a] compute a parametric sdf estimate
;;;               [b] compute the periodogram for the forward prediction errors
;;;               [c] smooth the periodogram using a Parzen lag window and
;;;                   then postcolor this lag window estimate using the
;;;                   AR prewhitening filter
(multiple-value-bind (ar-5-coeffs innovations-variance
X                                  junk-1 junk-2
X                                  forward-pred-errors)
X                     (Burg-algorithm *ocean-wave-data* 5)
X  (declare (ignore junk-1 junk-2))
X  (setf *ar-5-coeffs* ar-5-coeffs
X        *forward-pred-errors* forward-pred-errors)
X  (multiple-value-bind (Burg-5-sdf freqs N-f)
X                       (ar-coeffs->sdf ar-5-coeffs
X                                       innovations-variance
X                                       :sampling-time 1/4
X                                       :return-frequencies-p t)
X    (format t "~&Figure 440, thick curve, both plots ...")
X    (dotimes (i 7)
X      (format t "~&~6,4F: ~8,4F"
X              (svref freqs i)
X              (svref Burg-5-sdf i)))
X    (format t "~&...")
X    (dotimes (i 7)
X      (format t "~&~6,4F: ~8,4F"
X              (svref freqs (+ i (- N-f 7)))
X              (svref Burg-5-sdf (+ i (- N-f 7)))))
X    (format t "~&... Figure 440, thick curve, both plots")
X    (values)))
X
#|
;;; Here is what is printed when the above form is evaluated:
Figure 440, thick curve, both plots ...
0.0000:  51.2758
0.0078:  51.2919
0.0156:  51.3403
0.0234:  51.4215
0.0312:  51.5361
0.0391:  51.6849
0.0469:  51.8690
...
1.9531:  -7.9316
1.9609:  -7.9411
1.9687:  -7.9489
1.9766:  -7.9549
1.9844:  -7.9592
1.9922:  -7.9618
2.0000:  -7.9627
... Figure 440, thick curve, both plots
|#
X
(progn
X  (format t "~&Figure 441, plot (a) ...")
X  (dotimes (i 7)
X    (format t "~&~4D: ~8,4F"
X            (+ i 6)
X            (svref *forward-pred-errors* i)))
X  (format t "~&...")
X  (let ((N-fpe (length *forward-pred-errors*)))
X    (dotimes (i 7)
X      (format t "~&~4D: ~8,4F"
X              (+  i 6 (- N-fpe 7))
X              (svref *forward-pred-errors* (+ i (- N-fpe 7))))))
X  (format t "~&... Figure 441, plot (a)")
X  (values))
X
#|
;;; Here is what is printed when the above form is evaluated:
Figure 441, plot (a) ...
X   6:   6.6713
X   7:  14.7491
X   8:  -5.1810
X   9:   2.9906
X  10:  -9.3390
X  11:  21.8462
X  12:  -0.5174
...
1018: -12.8457
1019:   0.9300
1020:  -6.7758
1021:   3.0929
1022:   3.7940
1023:   3.0392
1024:   0.2890
... Figure 441, plot (a)
|#
X
(multiple-value-bind (the-periodogram freqs N-f)
X                     (periodogram
X                      *forward-pred-errors*
X                      :sampling-time 1/4)
X  (format t "~&Figure 441, plot (b) ...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs i)
X            (svref the-periodogram i)))
X  (format t "~&...")
X  (dotimes (i 7)
X    (format t "~&~6,4F: ~8,4F"
X            (svref freqs (+ i (- N-f 7)))
X            (svref the-periodogram (+ i (- N-f 7)))))
X  (format t "~&... Figure 441, plot (b)")
X  (values))
X
#|
;;; Here is what is printed when the above form is evaluated:
Figure 441, plot (b) ...
0.0039:  12.3693
0.0078:   4.4997
0.0117:   7.5193
0.0156:   8.1898
0.0195:  14.8114
0.0234:  19.1417
0.0273:  18.2599
...
1.9766:  17.3194
1.9805:  17.9742
1.9844:  18.6913
1.9883:  10.7566
1.9922:  25.4055
1.9961:  10.0394
2.0000:   6.1584
... Figure 441, plot (b)
|#
X
(multiple-value-bind (Parzen-m-55 freqs N-f nu B_W)
X                     (lag-window-spectral-estimate
X                      (acvs *forward-pred-errors*)
X                      #'(lambda (lag)
X                          (parzen-lag-window lag 55))
X                      :sampling-time 1/4
X                      :sdf-transformation nil)
X  (let ((postcolored-Parzen-m-55 (postcolor-spectral-estimate
X                                  Parzen-m-55
X                                  (ar-coeffs->prewhitening-filter *ar-5-coeffs*)
X                                  (length *forward-pred-errors*))))
X    (format t "~&Figure 441, thick curve, plot (c) ...")
X    (dotimes (i 7)
X      (format t "~&~6,4F: ~8,4F"
X              (svref freqs i)
X              (svref postcolored-Parzen-m-55 i)))
X    (format t "~&...")
X    (dotimes (i 7)
X      (format t "~&~6,4F: ~8,4F"
X              (svref freqs (+ i (- N-f 7)))
X              (svref postcolored-Parzen-m-55 (+ i (- N-f 7)))))
X    (format t "~&... Figure 441, thick curve, plot (c)")
X    (format t "~&equivalent degrees of freedom = ~5,1F" nu)
X    (format t "~&smoothing window bandwidth    = ~6,4F" B_W))
X  (values))
X
#|
;;; Here is what is printed when the above form is evaluated:
Figure 441, thick curve, plot (c) ...
0.0039:  51.6943
0.0078:  51.7057
0.0117:  51.7240
0.0156:  51.7487
0.0195:  51.7789
0.0234:  51.8134
0.0273:  51.8512
...
1.9766:  -5.6470
1.9805:  -5.5772
1.9844:  -5.5188
1.9883:  -5.4726
1.9922:  -5.4392
1.9961:  -5.4190
2.0000:  -5.4122
... Figure 441, thick curve, plot (c)
equivalent degrees of freedom =  68.7
smoothing window bandwidth    = 0.1349
|#
SHAR_EOF
chmod 0644 examples.lisp ||
echo 'restore of examples.lisp failed'
Wc_c="`wc -c < 'examples.lisp'`"
test 73426 -eq "$Wc_c" ||
	echo 'examples.lisp: original size 73426, current size' "$Wc_c"
fi
# ============= filtering.lisp ==============
if test -f 'filtering.lisp' -a X"$1" != X"-c"; then
	echo 'x - skipping filtering.lisp (File already exists)'
else
echo 'x - extracting filtering.lisp (Text)'
sed 's/^X//' << 'SHAR_EOF' > 'filtering.lisp' &&
;;;-*- Mode: LISP; Package: :SAPA; Syntax: COMMON-LISP -*-
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;
;  filtering.lisp
;
;  a collection of Lisp functions to filter a time series ...
;  Note:  before compiling and loading filtering.lisp,
;         you should compile and load (in the order listed)
;            sapa-package.lisp, utilities.lisp, basic-statistics.lisp,
;            dft-and-fft.lisp and tapers.lisp
;
;  6/3/93
;
;  SAPA, Version 1.0; Copyright 1993, Donald B. Percival, All Rights Reserved
;
;  Use and copying of this software and preparation of derivative works
;  based upon this software are permitted.  Any distribution of this
;  software or derivative works must comply with all applicable United
;  States export control laws.
; 
;  This software is made available AS IS, and no warranty -- about the
;  software, its performance, or its conformity to any
;  specification -- is given or implied. 
;
;  Comments about this software can be addressed to dbp@apl.washington.edu
;-------------------------------------------------------------------------------
;;; (compile-file "ccl:SAPA;filtering.lisp")
;;; (load "ccl:SAPA;filtering.fasl")
;-------------------------------------------------------------------------------
(if (not (find-package :SAPA))
X  (defpackage :SAPA (:USE :COMMON-LISP)))
X
(in-package :SAPA)
X
(export '(;;; general functions for filtering a time series ...
X          filter-time-series-direct
X          filter-time-series-fft
X          filter-time-series
X
X          ;;; creation of filters ...
X          ideal-low-pass-filter-irs
X          ideal-high-pass-filter-irs
X          ideal-band-pass-filter-irs
X          create-least-squares-low-pass-filter
X          triangular-convergence-factors
X          create-dpss-low-pass-filter
X          compose-symmetric-filters
X
X          ;;; calculation of the transfer function for a filter ...
X          transfer-function-for-filter
X
X          ;;; specialized functions for filtering a time series ...
X          three-point-smoother
X          n-applications-of-three-point-smoother
X          exponential-smoothing
X          running-median
X
X          ;;; specialized functions for filtering a time series ...
X          center&prewhiten&taper-time-series
X          ))
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  filter-time-series-direct
;;;                 filter-time-series-fft
;;;                 filter-time-series
;;;  all take a time series and a filter and return a filtered
;;;  time series.  The function filter-time-series-direct uses
;;;  the direct ``time domain'' definition of filtering to compute
;;;  the filtered series; filter-time-series-fft uses fft operations;
;;;  and filter-time-series call one of these two functions based upon
;;;  a crude assessment of which of the two methods is the fastest.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;; Note: filter-time-series-direct is primarily intended to be used
;;;       with short filters for which use of fft's would not be effecient.
(defun filter-time-series-direct
X       (time-series
X        the-filter
X        &key
X        (start 0)
X        (end (length time-series))
X        (result (make-array (- (- end start) (1- (length the-filter))))))
X  "given
X   [1] time-series (required)
X       ==> a vector containing a time series
X           x_0, x_1, ..., x_{N-1}
X   [2] the-filter (required)
X       ==> a vector containing the filter coefficients
X           g_0, g_1, ..., x_{K-1}
X   [3] start (keyword; 0)
X       ==> start index of time-series to be used
X   [4] end (keyword; length of time-series)
X       ==> 1 + end index of time-series to be used
X   [5] result (keyword; vector of appropriate length)
X       <== vector to contain filtered time series
X                 K-1
X           y_t = SUM g_k x_{t+K-1-k},  t = 0, ..., N-K+1
X                 k=0
X returns
X   [1] result, a vector containing the filtered time series
X   [2] the number of values in the filtered time series
---
Note: result can be the same as time-series"
X  (let* ((N-filter (length the-filter))
X         (N-filter-1 (1- N-filter))
X         (N-output (- (- end start) N-filter-1)))
X    ;;; reverse needed to conform to convention for
X    ;;; filtering (a convolution)
X    (setf the-filter (reverse the-filter))
X    (dotimes (i N-output (values result N-output))
X      (setf (aref result i)
X            (* (aref time-series (+ i start))
X               (aref the-filter 0)))
X      (dotimes (j N-filter-1)
X        (incf (aref result i)
X              (* (aref time-series (+ i start j 1))
X                 (aref the-filter (1+ j))))))))
X
#|
(filter-time-series-direct #(1 2 3 4 5 -5 -7 -9) #(1 -1))
;==> #(1 1 1 1 -10 -2 -2)
;    7
X
(filter-time-series-direct #(1 2 3 4 5 -5 -7 -9) #(1 -1) :start 2 :end 6)
;==> #(1 1 -10)
;    3
X
(let ((test #(1.1 2.2 3.3 4.4 5.4 -5 -7 -9)))
X  (filter-time-series-direct test #(1 -1) :result test))
;==> #(1.1 1.0999999999999996 1.1000000000000005 1.0 -10.4 -2 -2 -9)
;    7
|#
X
;-------------------------------------------------------------------------------
(defun filter-time-series-fft
X       (time-series
X        the-filter
X        &key
X        (start 0)
X        (end (length time-series))
X        (fft-size (* 4 (next-power-of-2
X                        (length the-filter))))
X        (verbose-p nil)
X        (result (make-array (- (- end start) (1- (length the-filter)))
X                            :initial-element 0.0)
X                result-supplied-p))
X  "given
X   [1] time-series (required)
X       ==> a vector containing a time series
X           x_0, x_1, ..., x_{N-1}
X   [2] the-filter (required)
X       ==> a vector containing the filter coefficients
X           g_0, g_1, ..., x_{K-1}
X   [3] start (keyword; 0)
X       ==> start index of time-series to be used
X   [4] end (keyword; length of time-series)
X       ==> 1 + end index of time-series to be used
X   [5] fft-size (keyword; 4 * power of 2 that is ceiling for filter length)
X       ==> size of fft's to be used (must be a power of 2)
X   [6] verbose-p (keyword; nil)
X       ==> if t, prints line after each block of data is processed;
X           if nil, prints nothing
X   [7] result (keyword; vector of appropriate length)
X       <== vector to contain filtered time series
X                 K-1
X           y_t = SUM g_k x_{t+K-1-k},  t = 0, ..., N-K+1
X                 k=0
X returns
X   [1] result, a vector containing the filtered time series
X   [2] the number of values in the filtered time series
---
Note: result can be the same as time-series"
X  (assert (power-of-2 fft-size))
X  (let* ((N-time-series (- end start))
X         (N-filter (length the-filter))
X         (N-filter-1 (1- N-filter))
X         (fft-of-filter (make-array fft-size :initial-element 0.0))
X         fft-of-data
X         (N-output (1+ (- N-time-series N-filter)))
X         (N-block (1+ (- fft-size N-filter)))
X         (i-result -1))
X    (if result-supplied-p
X      (fill result 0.0 :end N-output))
X    (copy-vector the-filter fft-of-filter)
X    (fft! fft-of-filter)
X    (do* ((k start (+ k N-block))
X          (still-to-go N-output (- still-to-go N-block))
X          (m (+ still-to-go N-filter-1) (+ still-to-go N-filter-1)))
X         ((not (plusp still-to-go)) . nil)
X      (setf fft-of-data
X            (if (>= m fft-size)
X              (subseq time-series k (+ k fft-size))
X              (concatenate 'array
X                           (subseq time-series k (+ k m))
X                           (make-array (- fft-size m)
X                                       :initial-element 0.0))))
X      (if verbose-p (format t "~&k = ~D, still-to-go = ~D" k still-to-go))
X      (fft! fft-of-data)
X      (dotimes (i fft-size)
X        (setf (aref fft-of-data i)
X              (/ (conjugate (* (aref fft-of-data i) (aref fft-of-filter i)))
X                 fft-size)))
X      (fft! fft-of-data)
X      (dotimes (i (min still-to-go N-block))
X        (setf (aref result (incf i-result))
X              (realpart (aref fft-of-data (+ N-filter-1 i))))))
X    (values result N-output)))
X
#|
(filter-time-series-fft #(1 2 3 4 5 -5 -7 -9) #(1 -1))
;==> #(0.9999999999999991 1.0000000000000018 0.9999999999999993 1.0
X       -10.0 -2.0000000000000018 -1.9999999999999993)
;    7
X
(filter-time-series-fft #(1 2 3 4 5 -5 -7 -9) #(1 -1) :start 2 :end 6)
;==> #(1.0000000000000004 0.9999999999999996 -10.0)
;    3
X
(let* ((a-cosine (sample-from-a-function 
X                  #'cos :delta-x (/ pi 100) :N 788))
X       (direct-filtered
X        (filter-time-series-direct a-cosine #(1 -1)))
X       (fft-filtered
X        (filter-time-series-fft a-cosine #(1 -1))))
X  (compare-seqs direct-filtered fft-filtered))
;==> 4.894807734342944E-17
;    1.97758476261356E-16
|#
X
;-------------------------------------------------------------------------------
(defun filter-time-series
X       (time-series
X        the-filter
X        &key
X        (start 0)
X        (end (length time-series))
X        (technique :fastest)
X        (result (make-array (- (- end start) (1- (length the-filter))))))
X  "given
X   [1] time-series (required)
X       ==> a vector containing a time series
X           x_0, x_1, ..., x_{N-1}
X   [2] the-filter (required)
X       ==> a vector containing the filter coefficients
X           g_0, g_1, ..., x_{K-1}
X   [3] start (keyword; 0)
X       ==> start index of time-series to be used
X   [4] end (keyword; length of time-series)
X       ==> 1 + end index of time-series to be used
X   [5] technique (keyword; :fastest)
X       ==> if :fastest, tries to pick fastest method;
X           if :fft, uses fft-based method;
X           if :direct, uses direct method
X   [6] result (keyword; vector of appropriate length)
X       <== vector to contain filtered time series
X                 K-1
X           y_t = SUM g_k x_{t+K-1-k},  t = 0, ..., N-K+1
X                 k=0
X returns
X   [1] result, a vector containing the filtered time series
X   [2] the number of values in the filtered time series
---
Note: result can be the same as time-series"
X  ;;; If fastest method is chosen, we do a ROUGH count of the number
X  ;;; of floating point operations for the direct and fft methods.
X  (if (eq technique :fastest)
X    (setf technique (which-is-faster? (length the-filter)
X                                      (- end start))))
X  (case technique
X    (:direct (filter-time-series-direct
X              time-series
X              the-filter
X              :start start
X              :end end
X              :result result))
X    (:fft (filter-time-series-fft
X           time-series
X           the-filter
X           :start start
X           :end end
X           :result result))
X    (otherwise
X     (error "technique (~A) should be either fft, direct or fastest"
X            technique))))
X
#|
(let* ((a-cosine (sample-from-a-function 
X                  #'cos :delta-x (/ pi 100) :N 788))
X       (direct-filtered
X        (filter-time-series-direct a-cosine #(1 -1)))
X       (fft-filtered
X        (filter-time-series-fft a-cosine #(1 -1)))
X       (filtered
X        (filter-time-series a-cosine #(1 -1))))
X  (print (compare-seqs direct-filtered fft-filtered))
X  (print (compare-seqs direct-filtered filtered))
X  (print (compare-seqs filtered fft-filtered))
X  (values))
;==> 4.894807734342944E-17 
;    4.894807734342944E-17 
;    0.0 
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  ideal-low-pass-filter-irs
;;;                 ideal-high-pass-filter-irs
;;;                 ideal-band-pass-filter-irs
;;;                 create-least-squares-low-pass-filter
;;;                 triangular-convergence-factors
;;;                 create-dpss-low-pass-filter
;;;                 compose-symmetric-filters
;;;  can be used to create a filter (by which we mean a vector containing
;;;  the filter coefficients).  The first three of these functions
;;;  return a single member of the impulse response sequence for an ideal
;;;  low-pass, high-pass or band-pass filter.  The next three functions
;;;  can be used to create one of the approximations to an ideal low-pass
;;;  filter discussed in Sections 5.8 and 5.9 of the SAPA book.  The
;;;  final function takes any number of symmetric filters of odd length
;;;  and returns the equivalent composite filter.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun ideal-low-pass-filter-irs (k W)
X  "given
X   [1] k (required)
X       ==> index of member of impulse response sequence (irs)
X           to be calculated (must be an integer)
X   [2] W (required)
X       ==> the cutoff frequency, standardized such that
X            0 < W < 0.5 = Nyquist frequency
returns
X   [1] kth member of the impulse response sequence
X       for an ideal low-pass filter with cutoff frequency W
---
Note: see Section 5.8 of the SAPA book"
X  ;(assert (and (integerp k) (plusp W) (< W 0.5)))
X  (if (zerop k)
X    (* 2 W)
X    (/ (sin (* 2 pi W k)) (* pi k))))
X
#|
(ideal-low-pass-filter-irs 0 0.1)   ;==> 0.2
(ideal-low-pass-filter-irs 1 0.1)   ;==> 0.1870978567577278
(ideal-low-pass-filter-irs -1 0.1)  ;==> 0.1870978567577278
|#
X
;-------------------------------------------------------------------------------
(defun ideal-high-pass-filter-irs (k W)
X  "given
X   [1] k (required)
X       ==> index of member of impulse response sequence (irs)
X           to be calculated (must be an integer)
X   [2] W (required)
X       ==> the cutoff frequency, standardized such that
X            0 < W < 0.5 = Nyquist frequency
returns
X   [1] kth member of the impulse response sequence
X       for an ideal low-pass filter with cutoff frequency W"
X  (assert (and (integerp k) (plusp W) (< W 0.5)))
X  (if (zerop k)
X    (- 1.0 (* 2 W))
X    (- (/ (sin (* 2 pi W k)) (* pi k)))))
X
#|
(ideal-high-pass-filter-irs 0 0.1)   ;==> 0.8
(ideal-high-pass-filter-irs 1 0.1)   ;==> -0.1870978567577278
(ideal-high-pass-filter-irs -1 0.1)  ;==> -0.1870978567577278
|#
X
;-------------------------------------------------------------------------------
;;; Note: by setting W-high = 0.5, this routine produces the irs for
;;;       a high-pass filter;
;;;       by setting W-low = 0.5, this routine produces the irs for
;;;       a low-pass filter
(defun ideal-band-pass-filter-irs (k W-low W-high)
X  "given
X   [1] k (required)
X       ==> index of member of impulse response sequence (irs)
X           to be calculated (must be an integer)
X   [2] W-low (required)
X       ==> the low frequency cutoff (in standardized units)
X   [3] W-high (required)
X       ==> the high frequency cutoff (in standardized units
X           so that 0 <= W-low < W-high <= 0.5, the assumed
X           Nyquist frequency).
returns
X   [1] k-th member of the impulse response sequence
X       for an ideal band-pass filter"
X  (assert (and (integerp k) (/= W-low W-high) (<= 0.0 W-low W-high 0.5)))
X  (let ((width (- W-high W-low)))
X    (if (zerop k)
X      (* 2 width)
X      (/ (* 2 (cos (* pi (+ W-high W-low) k)) (sin (* pi width k)))
X         (* pi k)))))
X
#|
;;; low-pass
(ideal-band-pass-filter-irs  0 0.0 0.1)   ;==>  0.2
(ideal-band-pass-filter-irs  1 0.0 0.1)   ;==>  0.1870978567577278
(ideal-band-pass-filter-irs -1 0.0 0.1)   ;==>  0.1870978567577278
;;; high-pass
(ideal-band-pass-filter-irs  0 0.3 0.5)   ;==>  0.4
(ideal-band-pass-filter-irs  1 0.3 0.5)   ;==> -0.3027306914562628
(ideal-band-pass-filter-irs -1 0.3 0.5)   ;==> -0.3027306914562628
|#
X
;-------------------------------------------------------------------------------
(defun create-least-squares-low-pass-filter
X       (filter-length
X        W
X        &key
X        (convergence-factors nil)
X        (Nyquist-frequency 0.5)
X        (result (make-array filter-length)))
X  "given
X   [1] filter-length (required)
X       ==> an odd positive integer = 2L +1
X   [2] W (required)
X       ==> a cutoff frequency greater than 0
X           and less than the Nyquist frequency
X   [3] convergence-factors (keyword; nil)
X       ==> a one-argument function that maps
X           an integer to a convergence factor;
X           nil is also acceptable, in which case
X           no convergence factors are used
X   [4] Nyquist-frequency (keyword; 0.5)
X       ==> the Nyquist frequency
X   [5] result (keyword; vector of length filter-length)
X       <== vector of length filter-length
X           into which filter coefficients
X           are placed (returned by the function)
uses a least squares approximation to a low-pass filter and
returns
X   [1] a symmetric low-pass filter of length 2L+1
X       and a gain of unity at zero frequency;
X       i.e., result(0) = result(2L)
X             result(1) = result(2L-1)
X             etc., and
X             (sum result) = 1.0
---
Note: see Section 5.8 of the SAPA book;
X      (aref result 0)       corresponds to element -L of the filter;
X      (aref result L)       corresponds to element  0 of the filter;
X      (aref result (* 2 L)) corresponds to element  L of the filter"
X  (assert (and (plusp filter-length)
X               (oddp filter-length)
X               (plusp W)
X               (plusp Nyquist-frequency)
X               (< 0.0 W Nyquist-frequency)))
X  ;;; convert W from user units to standardized units ...
X  (if (not (= Nyquist-frequency 0.5))
X    (divf W (* 2 Nyquist-frequency)))
X  (let* ((L (/ (1- filter-length) 2))
X         (minus-L-to-L (iota (- L) L)))
X    (transform-a-sequence!
X     #'(lambda (k) (ideal-low-pass-filter-irs k W))
X     minus-L-to-L
X     :result result)
X    (if convergence-factors (x*y! result (map-into
X                                          minus-L-to-L
X                                          convergence-factors
X                                          minus-L-to-L)))
X    (a*x! (/ (sum result)) result)))
X
#|
(create-least-squares-low-pass-filter 5 0.1)
;==> #(0.1726089497020141 0.21335639535653148 0.22806930988290877 0.21335639535653148 0.1726089497020141)
X
(sum (create-least-squares-low-pass-filter 5 0.1))
;==> 1.0
|#
X
;-------------------------------------------------------------------------------
(defun triangular-convergence-factors (k filter-length)
X  (let ((abs-k (abs k)))
X    (if (> abs-k (truncate filter-length 2))
X      0.0
X      (- 1.0 (/ (* 2 abs-k) (1+ filter-length))))))
X
#|
(triangular-convergence-factors -2 3)  ;==> 0.0
(triangular-convergence-factors -1 3)  ;==> 0.5
(triangular-convergence-factors  0 3)  ;==> 1.0
(triangular-convergence-factors  1 3)  ;==> 0.5
(triangular-convergence-factors  2 3)  ;==> 0.0
X
(create-least-squares-low-pass-filter
X 5 0.1
X :convergence-factors
X #'(lambda (k) (triangular-convergence-factors k 5)))
;==> #(0.09167422811028572 0.22663115545827048 0.3633892328628876 0.22663115545827048 0.09167422811028572)
X
(sum (create-least-squares-low-pass-filter
X      5 0.1
X      :convergence-factors
X      #'(lambda (k) (triangular-convergence-factors k 5))))
;==> 1.0
|#
X
;-------------------------------------------------------------------------------
(defun create-dpss-low-pass-filter
X       (filter-length
X        delta
X        W
X        &key
X        (Nyquist-frequency 0.5)
X        (result (make-array filter-length)))
X  "given
X   [1] filter-length (required)
X       ==> an odd positive integer = 2L +1
X   [2] delta (required)
X       ==> ``W'' parameter for dpss in user
X           units (see page 182 of the SAPA book)
X   [2] W (required)
X       ==> a cutoff frequency greater than 0
X           and less than the Nyquist frequency
X   [3] Nyquist-frequency (keyword; 0.5)
X       ==> the Nyquist frequency
X   [4] result (keyword; vector of length filter-length)
X       <== vector of length filter-length
X           into which filter coefficients
X           are placed (returned by the function)
uses a dpss as convergence factors in least squares approximation
to a low-pass filter and
returns
X   [1] a symmetric low-pass filter of length 2L+1
X       and a gain of unity at zero frequency;
X       i.e., result(0) = result(2L)
X             result(1) = result(2L-1)
X             etc., and
X             (sum result) = 1.0
---
Note: see Section 5.9 of the SAPA book;
X      (aref result 0)       corresponds to element -L of the filter;
X      (aref result L)       corresponds to element  0 of the filter;
X      (aref result (* 2 L)) corresponds to element  L of the filter"
X  (assert (and (plusp filter-length)
X               (oddp filter-length)
X               (plusp delta)               
X               (plusp W)
X               (plusp Nyquist-frequency)
X               (< 0.0 W Nyquist-frequency)))
X  ;;; convert delta and W from user units to standardized units ...
X  (when (not (= Nyquist-frequency 0.5))
X    (divf delta (* 2 Nyquist-frequency))
X    (divf W (* 2 Nyquist-frequency)))
X  (let* ((L (/ (1- filter-length) 2))
X         (minus-L-to-L (iota (- L) L)))
X    (transform-a-sequence!
X     #'(lambda (k) (ideal-low-pass-filter-irs k W))
X     minus-L-to-L
X     :result result)
X    (x*y! result (dpss-data-taper! (make-array filter-length
X                                               :initial-element 1.0)
X                                   :taper-parameter
X                                   (* filter-length delta)))
X    (a*x! (/ (sum result)) result)))
X
#|
(create-dpss-low-pass-filter 5 0.04 0.1)
;==> #(0.16887701282533799 0.21504897804256012 0.23214801826420375 0.21504897804256012 0.16887701282533799)
X
(sum (create-dpss-low-pass-filter 5 0.04 0.1))
;==> 0.9999999999999999
|#
X
;-------------------------------------------------------------------------------
(defun compose-symmetric-filters
X       (&rest filters)
X  "given any number of symmetric filters
(each of odd length and each represented
by a vector of coefficients),
returns the composite filter
that will produce an output
identical to that obtained
by cascading the individual filters"
X  (let* ((first-filter (car filters))
X         (second-filter (if (= (length filters) 2)
X                          (nth 1 filters)
X                          (apply 'compose-symmetric-filters (cdr filters))))
X         (half-length-first-filter
X          (/ (1- (array-total-size first-filter)) 2))
X         (half-length-second-filter
X          (/ (1- (array-total-size second-filter)) 2)))
X    (assert (integerp half-length-first-filter))
X    (assert (integerp half-length-second-filter))
X    (let* ((K (min half-length-first-filter half-length-second-filter))
X           (L (max half-length-first-filter half-length-second-filter))
X           (K+L (+ K L))
X           (2-times-K+L (* 2 K+L))
X           (short-length (1+ (* 2 K)))
X           (index-of-last-short (1- short-length))
X           (composite-filter (make-array (1+ 2-times-K+L)
X                                         :initial-element 0.0))
X           (short-filter (if (<= half-length-first-filter
X                                 half-length-second-filter)
X                           first-filter second-filter))
X           (long-filter (if (eq short-filter first-filter)
X                          second-filter first-filter)))
X      (dotimes (j (1+ K+L) composite-filter)
X        (setf (aref composite-filter (- 2-times-K+L j))
X              (dotimes (m (min (1+ j) short-length) (aref composite-filter j))
X                (incf (aref composite-filter j)
X                      (* (aref short-filter (- index-of-last-short m))
X                         (aref long-filter (- j m))))))))))
X
#|
(compose-symmetric-filters #(1/4 1/2 1/4) #(1/4 1/2 1/4))
;==> #(0.0625 0.25 0.375 0.25 0.0625)
X
(compose-symmetric-filters #(1/8 1/4 1/4 1/4 1/8) #(1/4 1/2 1/4))
;==> #(0.03125 0.125 0.21875 0.25 0.21875 0.125 0.03125)
X
(compose-symmetric-filters #(1/4 1/2 1/4) #(1/4 1/2 1/4) #(1/4 1/2 1/4))
;==> #(0.015625 0.09375 0.234375 0.3125 0.234375 0.09375 0.015625)
X
(compose-symmetric-filters #(0.0625 0.25 0.375 0.25 0.0625) #(1/4 1/2 1/4))
;==> #(0.015625 0.09375 0.234375 0.3125 0.234375 0.09375 0.015625)
X
(compose-symmetric-filters #(1/4 1/2 1/4) #(0.0625 0.25 0.375 0.25 0.0625))
;==> #(0.015625 0.09375 0.234375 0.3125 0.234375 0.09375 0.015625)
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The function   transfer-function-for-filter
;;;  can be used to compute the transfer function for a filter
;;;  represented by a vector of coefficents.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun transfer-function-for-filter
X       (the-filter
X        &key
X        (tf-transformation #'(lambda (x)
X                               (careful-convert-to-dB
X                                (expt (abs x) 2)
X                                -100.0)))
X        (actual-index-of-1st-filter-coeff 0)
X        (N-fft (* 4 (next-power-of-2 (length the-filter))))
X        (return-frequencies-p nil)
X        (Nyquist-frequency 0.5)
X        (result-tf (make-array (1+ (/ N-fft 2))))
X        (result-freq (if return-frequencies-p (make-array (1+ (/ N-fft 2))))))
X  "given
X   [1] the-filter (required)
X       ==> a vector of filter coefficients
X   [2] tf-transformation (keyword; mod squared in dB with 0 mapped to -100 dB)
X       ==> a function of one argument or nil;
X           if bound to a function, the function is used
X           to transform all of the elements of the transfer function
X   [3] actual-index-of-1st-filter-coeff (keyword; 0)
X       ==> the filter coefficient in (aref the-filter 0) is assumed
X           to have an index equal to whatever is given here
X           (this is only needed if the phase of the transfer function
X           is needed)
X   [4] N-fft (keyword; 4 * (next-power-of-2 (length the-filter)))
X       ==> the length of vector to be used in the fast Fourier transform 
X           -- the larger this is, the finer the grid of frequencies
X           over which the transfer function is computed; this number
X           MUST be a power of 2.
X   [5] return-frequencies-p (keyword; nil)
X       ==> if t, the frequencies associated with the transfer function
X           are computed and returned in result-freq
X   [6] Nyquist-frequency (keyword; 0.5)
X       ==> the Nyquist frequency
X   [7] result-tf (keyword; vector of length (1+ (/ N-fft 2)))
X       <== vector of length (1+ (/ N-fft 2))
X           into which the properly transformed transfer function
X           is placed (returned by the function)
X   [8] result-freq (keyword; vector of length (1+ (/ N-fft 2)) if return-frequencies-p t)
X       <== vector of length (1+ (/ N-fft 2))
X           into which the frequencies associated with the values
X           in result-tf are placed if return-frequencies-p is true
X           (returned by the function)
returns
X   [1] result-tf, a vector holding
X       the properly transformed transfer function
X   [2] nil (if return-frequencies-p is nil) or
X       result-freq (if return-frequencies-p is t),
X       where result-freq is a vector holding
X       the frequencies associated with values
X       in  result-tf
---
Note: see Section 5.3 of the SAPA book"
X  (let* ((fft-scratch (make-array N-fft :initial-element 0.0))
X         (N-freq (1+ (/ N-fft 2))))
X    ;(declare (dynamic-extent fft-scratch))
X    (copy-vector the-filter fft-scratch)
X    (fft! fft-scratch)
X    (if (not (zerop actual-index-of-1st-filter-coeff))
X      (let ((mult-factor
X             (exp (complex
X                   0.0
X                   (/ (* -2 pi actual-index-of-1st-filter-coeff) N-fft))))
X            (freq-factor 1.0))
X        (dotimes (i N-freq)
X          (multf (svref fft-scratch i) freq-factor)
X          (multf freq-factor mult-factor))))
X    (copy-vector fft-scratch result-tf :end N-freq)
X    (if tf-transformation
X      (map-into result-tf tf-transformation result-tf))
X    (cond
X     (return-frequencies-p
X      (let ((fund-freq (float (/ (* N-fft (Nyquist-frequency->sampling-time
X                                           Nyquist-frequency))))))
X        (dotimes (i N-freq)
X          (setf (aref result-freq i) (* i fund-freq))))
X      (values result-tf result-freq))
X     (t
X      (values result-tf)))))
X
#|
;;; This examples uses the 3 point filter discussed in Section 5.7
;;; of the SAPA book -- a plot of the squared modulus of the transfer function
;;; is given in Figure 172.
;;; First example: squared modulus of transfer function in dB
;;;                (associated frequences are j/16, j = 0, 1, ..., 8)
(transfer-function-for-filter
X #(1/4 1/2 1/4))
;==> #( 0.0
X       -0.33704242622473096
X       -1.3753861631621742
X       -3.2061447321250065
X       -6.020599913279622
X      -10.210441257596415
X      -16.686413576676696
X      -28.390570977011226
X       -100.0)
X
;;; Second example: squared modulus of transfer function
(transfer-function-for-filter
X #(1/4 1/2 1/4)
X :tf-transformation #'(lambda (x) (expt (abs x) 2)))
;==> #(1.0
X       0.9253281139039617
X       0.7285533905932737
X       0.4779533685342263
X       0.25 0.09526993616913669
X       0.02144660940672623
X       0.0014485813926750574
X       0.0)
X
;;; Third example: transfer function itself (untransformed)
(transfer-function-for-filter
X #(1/4 1/2 1/4)
X :actual-index-of-1st-filter-coeff -1
X :tf-transformation nil)
;==> #(  1.0
X      #c(0.9619397662556433 0.0)
X      #c(0.8535533905932737 5.551115123125783E-17)
X      #c(0.6913417161825448 2.7755575615628914E-17)
X      #c(0.5 0.0)
X      #c(0.30865828381745514 0.0)
X      #c(0.14644660940672624 -8.326672684688674E-17)
X      #c(0.03806023374435655 -3.2959746043559335E-17)
X      #c(-0.0 -0.0))
X
;;; Page 171 of the SAPA book says that this transfer function
;;; is equal to the following:
(sample-from-a-function #'(lambda (f) (expt (cos (* pi f)) 2))
X                        :delta-x 1/16
X                        :n 9)
;==> #(1.0
X       0.9619397662556434
X       0.8535533905932737
X       0.6913417161825449
X       0.5000000000000001
X       0.3086582838174552
X       0.1464466094067263
X       0.038060233744356645
X       3.7491518045553436E-33)
;    #(0.0 0.0625 0.125 0.1875 0.25 0.3125 0.375 0.4375 0.5)
;;; ... quite good agreement
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  three-point-smoother
;;;                 n-applications-of-three-point-smoother
;;;                 exponential-smoothing
;;;                 running-median
;;;  implement some specialized filters, all of a (more or less)
;;;  low-pass form (the running median filter is nonlinear).
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun three-point-smoother (a-seq)
X  "given a sequence of length N,
returns a sequence of length N-2 formed by filtering
the input sequence with a three-point filter
with coefficients 1/4, 1/2 and 1/4
---
Note: see Section 5.7 of the SAPA book"
X  (let* ((n (length a-seq))
X         (n-2 (- n 2))
X         (smoothed-seq (subseq a-seq 1 (1- n))))
X    (dotimes (i n-2 smoothed-seq)
X      (multf (elt smoothed-seq i) 0.5)
X      (incf (elt smoothed-seq i)
X            (* 0.25 (+ (elt a-seq i) (elt a-seq (+ i 2))))))))
X
#|
(three-point-smoother #(71.0 63.0 70.0 88.0 99.0 90.0 110.0))
;==> #(66.75 72.75 86.25 94.0 97.25)
(three-point-smoother #(66.75 72.75 86.25 94.0 97.25))
;==> #(74.625 84.8125 92.875)
|#
X
;-------------------------------------------------------------------------------
(defun n-applications-of-three-point-smoother (a-seq n)
X  "given
X   [1] a-seq (required)
X       ==> a sequence of length 2*n + 1 or greater
X   [2] n (required)
X       ==> positive integer indicating the number
X           of times the three-point smoother is
X           to be applied
returns
X   [1] the result of repetitively smoothing a-seq
X       using the function three-point-smoother"
X  (dotimes (i n a-seq)
X    (setf a-seq (three-point-smoother a-seq))))
X
#|
(n-applications-of-three-point-smoother
X  #(71.0 63.0 70.0 88.0 99.0 90.0 110.0)
X  2)
;==> #(74.625 84.8125 92.875)
|#
X
;-------------------------------------------------------------------------------
(defun exponential-smoothing
X       (time-series
X        alpha
X        &key
X        (initial-prediction (elt time-series 0))
X        (return-smoothed-values-p t)
X        (n-predictions 0)
X        (result-smoothed
X         (if return-smoothed-values-p
X           (make-array (length time-series))))
X        (result-predictions
X         (if (plusp n-predictions)
X           (make-array n-predictions))))
X  "given
X   [1] time-series (required)
X       ==> a sequence of numbers
X   [2] alpha (required)
X       ==> parameter controlling the degree of exponential smoothing
X           (usually  0 < alpha < 1)
X   [3] initial-prediction (keyword; first element of time-series)
X       ==> to get the smoother going, we need a ``prediction'' for
X           the first element in the sequence time-series -- two
X           common hacks are the first element itself (the default)
X           or 0 (if the time series has a zero sample mean)
X   [4] return-smoothed-values-p (keyword; t)
X       ==> if t, returns smoothed (i.e., filtered) time series
X   [5] n-predictions (keyword; 0)
X       ==> number of predictions desired
X   [6] result-smoothed (keyword; vector of length of time-series or nil)
X       <== a sequence to hold smoothed (filtered) time series -- used
X           only if return-smoothed-values-p is true
X   [7] result-predictions (keyword; vector of length n-predictions or nil)
X       <== a sequence to hold predicted values of the time series --
X           used only if n-predictions is positive
returns
X   [1] sum-of-squares of prediction errors
X   [2] either one step ahead prediction (if n-predictions is 0)
X       or vector of length n-predictions if n-predictions > 0)
X   [3] vector with smoothed values
X       (nil if return-smoothed-values-p is true)
---
Note: see Section 7.3 of ``Time Series: A Biostatistical Introduction''
by Diggle, 1990"
X  (let ((n-ts (length time-series))
X        (sum-of-squares-of-prediction-errors 0.0)
X        (current-prediciton initial-prediction)
X        current-prediction-error)
X    (dotimes (i n-ts)
X      (if return-smoothed-values-p
X        (setf (elt result-smoothed i) current-prediciton))
X      (setf current-prediction-error (- (elt time-series i) current-prediciton)
X            current-prediciton (+ current-prediciton
X                                  (* alpha current-prediction-error)))
X      (incf sum-of-squares-of-prediction-errors
X            (* current-prediction-error current-prediction-error)))
X    (if (plusp n-predictions)
X      (let ((last-value-in-time-series (elt time-series (1- n-ts)))
X            (1-alpha (- 1.0 alpha)))
X        (dotimes (i n-predictions (values
X                                   sum-of-squares-of-prediction-errors
X                                   result-predictions
X                                   result-smoothed))
X          (setf (elt result-predictions i) current-prediciton
X                current-prediciton (+ (* alpha last-value-in-time-series)
X                                      (* 1-alpha current-prediciton)))))
X      (values
X       sum-of-squares-of-prediction-errors
X       current-prediciton
X       result-smoothed))))
X
#|
(exponential-smoothing 
X #(71.0 63.0 70.0 88.0 99.0 90.0 110.0)
X 0.2)
;==> 2064.469974937598
;    86.558592
;    #(71.0 71.0 69.4 69.52000000000001 73.21600000000001 78.37280000000001 80.69824000000001)
X
(exponential-smoothing 
X #(71.0 63.0 70.0 88.0 99.0 90.0 110.0)
X 0.9)
;==> 1043.8179815444005
;    108.077138
;    #(71.0 71.0 63.8 69.38 86.138 97.7138 90.77138)
X
(exponential-smoothing 
X #(71.0 63.0 70.0 88.0 99.0 90.0 110.0)
X 0.9
X :n-predictions 10)
;==> 1043.8179815444005
;    #(108.077138 109.8077138 109.98077138 109.998077138 109.9998077138 109.99998077138 109.999998077138 109.99999980771379 109.99999998077138 109.99999999807713)
;    #(71.0 71.0 63.8 69.38 86.138 97.7138 90.77138)
|#
X
;-------------------------------------------------------------------------------
(defun running-median
X       (time-series
X        K
X        &key
X        (result (make-array (- (length time-series) (1- K)))))
X  "given
X   [1] time-series (required)
X       ==> a sequence of numbers
X   [2] K (required)
X       ==> positive integer giving the number of points
X           in the running median
X   [3] result (keyword; vector of appropriate length)
X       <== a sequence to hold running medians
X           of time series
return
X   [1] result, the vector of running medians
X        of K consecutive points"
X  (dotimes (i (- (length time-series) (1- K)) result)
X    (setf (aref result i)
X          (sample-median (subseq time-series i (+ i K))))))
X
#|
(running-median #(71.0 63.0 70.0 88.0 99.0 90.0 110.0) 3)
;==> #(70.0 70.0 88.0 90.0 99.0)
X
(running-median #(71.0 63.0 70.0 88.0 99.0 90.0 110.0) 4)
;==> #(70.5 79.0 89.0 94.5)
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The function   center&prewhiten&taper-time-series
;;;  handles centering, prewhitening and tapering of a time series
;;;  (if just centering and tapering are needed, it is probably
;;;  best to use center&taper-time-series in tapering.lisp).
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun center&prewhiten&taper-time-series
X       (time-series
X        &key
X        (center-data t)  ;t, nil or value to be subtracted off ...
X        (start 0)
X        (end (length time-series))
X        (prewhitening-filter nil)
X        (recenter-after-prewhitening-p t)
X        (data-taper nil)
X        (data-taper-parameters)
X        (recenter-after-tapering-p t)
X        (restore-power-option-p t)
X        (result (make-array
X                 (if prewhitening-filter
X                   (- end start (1- (length prewhitening-filter)))
X                   (- end start)))))
X  "given
X   [1] time-series (required)
X       ==> a sequence of time series values
X   [2] center-data (keyword; t)
X       ==> if t, subtract sample mean from time series;
X           if a number, subtracts that number from time series;
X           if nil, time-series is not centered
X   [3] start (keyword; 0)
X       ==> start index of sequence to be used
X   [4] end (keyword; length of time-series)
X       ==> 1 + end index of sequence to be used
X   [5] prewhitening-filter (keyword; nil)
X       ==> vector with coeffients of prewhitening filter
X           or nil (if no prewhitening is to be done)
X   [6] recenter-after-prewhitening-p (keyword; t)
X       ==> if t, prewhitened series is centered using sample mean
X           (not used if  prewhitening-filter is nil)
X   [7] data-taper (keyword; nil)
X       ==> nil or a tapering function
X   [8] data-taper-parameters (keyword)
X       ==> parameters for tapering function (not used
X           if data-taper is nil)
X   [9] recenter-after-tapering-p (keyword; t)
X       ==> if t and data-taper is a function,
X           centers tapered series by subtracting
X           off its sample mean
X  [10] restore-power-option-p (keyword; t)
X       ==> if t and data-taper is a function,
X           normalizes tapered series to have same
X           sum of squares as before tapering
X  [11] result (keyword; vector of size (- end start))
X       <== sequence to hold centered, prewhitened
X           and tapered time series
returns
X   [1] result, the sequence containing the centered, prewhitened
X       and tapered time series (time-series are unaltered unless
X       it is bound to result)
X   [2] the number used to center the time series
X       (this is nil if center-data is nil)
X   [3] C_h, the variance inflation factor due to the data taper
---
Note: see also center&taper-time-series in tapers.lisp"
X  (cond
X   ((not prewhitening-filter)
X    (center&taper-time-series
X     time-series
X     :center-data center-data
X     :start start
X     :end end
X     :data-taper data-taper
X     :data-taper-parameters data-taper-parameters
X     :recenter-after-tapering-p recenter-after-tapering-p
X     :restore-power-option-p restore-power-option-p
X     :result result))
X   (t
X    ;;; Because result might not be long enough to hold the
X    ;;; the centered time series, we must create an intermediate
X    ;;; array (since we don't want to trash time-series. The call
X    ;;; to subseq thus does two things at once for us:
X    ;;;       [1] subsets the time series and
X    ;;;       [2] creates a NEW sequence so that
X    ;;;           time-series is NOT trashed.
X    (setf time-series (subseq time-series start end))
X    (let ((center-factor nil)
X          (C_h 1.0)
X          (N (- end start (1- (length prewhitening-filter)))))
X      ;;; subtract the sample mean from time series if required
X      (when center-data
X        (setf center-factor (if (numberp center-data)
X                              center-data
X                              (sample-mean time-series)))
X        (x+b! time-series (- center-factor)))
X      ;;; prewhiten time series ...
X      (filter-time-series time-series prewhitening-filter :result result)
X      ;;; Note: prewhitening can reintroduce a nonzero sample mean,
X      ;;;       so here we remove it if so desired.
X      (if recenter-after-prewhitening-p
X        (let ((recenter-factor (sample-mean result :end N)))
X          (dotimes (i N)
X            (decf (elt result i) recenter-factor))))
X      ;;; taper the time series if required
X      (when data-taper
X        (multiple-value-bind (junk more-junk local-C_h)
X                             (center&taper-time-series
X                              result
X                              :center-data nil
X                              :end N
X                              :data-taper data-taper
X                              :data-taper-parameters data-taper-parameters
X                              :recenter-after-tapering-p recenter-after-tapering-p
X                              :restore-power-option-p restore-power-option-p
X                              :result result)
X          (declare (ignore junk more-junk))
X          (setf C_h local-C_h)))
X      (values result center-factor C_h)))))
X
#|
(defvar *20-point-time-series*
X  #(71.0 63.0 70.0 88.0 99.0 90.0 110.0 135.0 128.0 154.0
X    156.0 141.0 131.0 132.0 141.0 104.0 136.0 146.0 124.0 129.0))
(sample-mean-and-variance *20-point-time-series*)
;==> 117.4
;    782.6399999999998
X
(center&prewhiten&taper-time-series
X *20-point-time-series*)
;==> #(-46.400000000000006 -54.400000000000006 -47.400000000000006
X       -29.400000000000006 -18.400000000000006 -27.400000000000006
X        -7.400000000000006 17.599999999999994 10.599999999999994
X        36.599999999999994 38.599999999999994 23.599999999999994
X        13.599999999999994 14.599999999999994 23.599999999999994
X       -13.400000000000006 18.599999999999994 28.599999999999994
X         6.599999999999994 11.599999999999994)
;    117.4
;    1.0
(sum (center&prewhiten&taper-time-series *20-point-time-series*))
;==> -1.1368683772161603E-13
X
(center&prewhiten&taper-time-series
X *20-point-time-series*
X :start 5 :end 15)
;==> #(-41.80000000000001 -21.80000000000001 3.1999999999999886
X       -3.8000000000000114 22.19999999999999 24.19999999999999
X       9.199999999999989 -0.8000000000000114 0.19999999999998863
X       9.199999999999989)
;     131.8
;     1.0
(sum (center&prewhiten&taper-time-series
X      *20-point-time-series*
X      :start 5 :end 15))
;==> -1.1368683772161603E-13
X
(center&prewhiten&taper-time-series
X *20-point-time-series*
X :start 5 :end 15 :prewhitening-filter #(-1 1)
X :recenter-after-prewhitening-p nil)
;==> #(-20.0 -25.0 7.0 -26.0 -2.0 15.0 10.0 -1.0 -9.0)
;    131.8
;    1.0
(sample-mean (center&prewhiten&taper-time-series
X              *20-point-time-series*
X              :start 5 :end 15 :prewhitening-filter #(-1 1)
X              :recenter-after-prewhitening-p nil))
;==> -5.666666666666667
X
(center&prewhiten&taper-time-series
X *20-point-time-series*
X :start 5 :end 15 :prewhitening-filter #(-1 1))
;==> #(-14.333333333333332 -19.333333333333332 12.666666666666668
X       -20.333333333333332 3.666666666666667 20.666666666666668
X        15.666666666666668 4.666666666666667 -3.333333333333333)
;    131.8
;    1.0
(sample-mean-and-variance (center&prewhiten&taper-time-series
X                           *20-point-time-series*
X                           :start 5 :end 15 :prewhitening-filter #(-1 1)))
;==> 5.921189464667501E-16
;    208.0
X
(center&prewhiten&taper-time-series
X *20-point-time-series*
X :start 5 :end 15 :prewhitening-filter #(-1 1)
X :data-taper #'dpss-data-taper!
X :data-taper-parameters 1.0)
;==> #(-8.029535056948914 -15.939710264707552 11.111941111527308
X       -25.99193293330027 2.981287578952983 23.07865115777016
X       14.13595435534983 1.7916496140314773 -3.13830556267502)
;    131.8
;    1.327454606904972
(sample-mean-and-variance (center&prewhiten&taper-time-series
X                           *20-point-time-series*
X                           :start 5 :end 15 :prewhitening-filter #(-1 1)
X                           :data-taper #'dpss-data-taper!
X                           :data-taper-parameters 1.0))
;==> 0.0
;    207.99999999999997
X
(center&prewhiten&taper-time-series
X *20-point-time-series*
X :start 5 :end 15 :prewhitening-filter #(-1 1)
X :data-taper #'dpss-data-taper!
X :data-taper-parameters 1.0
X :restore-power-option-p nil)
;==> #(-8.525293352136908 -16.92385735797956 11.798012838297177
X       -27.596722783305598 3.1653577694589425 24.503569620927365
X       15.00873423393055 1.902269363733894 -3.3320703329258587)
;    131.8
;    1.327454606904972
(sample-mean-and-variance (center&prewhiten&taper-time-series
X                           *20-point-time-series*
X                           :start 5 :end 15 :prewhitening-filter #(-1 1)
X                           :data-taper #'dpss-data-taper!
X                           :data-taper-parameters 1.0
X                           :restore-power-option-p nil))
;==> 5.427757009278544E-16
;    234.4775142934783
X
(center&prewhiten&taper-time-series
X *20-point-time-series*
X :start 5 :end 15 :prewhitening-filter #(-1 1)
X :data-taper #'dpss-data-taper!
X :data-taper-parameters 1.0
X :recenter-after-tapering-p nil
X :restore-power-option-p nil)
;==> #(-6.766926964426519 -15.165490970269172 13.556379226007566
X       -25.83835639559521 4.923724157169332 26.261936008637754
X        16.76710062164094 3.6606357514442833 -1.5737039452154695)
;    131.8
;    1.327454606904972
(sample-mean-and-variance (center&prewhiten&taper-time-series
X                           *20-point-time-series*
X                           :start 5 :end 15 :prewhitening-filter #(-1 1)
X                           :data-taper #'dpss-data-taper!
X                           :data-taper-parameters 1.0
X                           :recenter-after-tapering-p nil
X                           :restore-power-option-p nil))
;==> 1.7583663877103892
;    234.4775142934783
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  Everything below here consists of internal symbols in the SAPA package
;;;  and should be regarded as "dirty laundry" ...
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;; Here we attempt to ROUGHLY count the number of operations to do
;;; filtering using the direct and fft methods for different combinations
;;; of filter length and sample size.  This function is used
;;; by filter-time-series to choose between the fft and direct way
;;; of computing a filtered time series.
(defun which-is-faster? (N-filter N-ts)
X  (let* ((N-output (-  N-ts (1- N-filter)))
X         (fft-size (* 4 (next-power-of-2 N-filter)))
X         (fft-cost (* 1.5 fft-size (log fft-size 2))))
X    (if (< (* N-output N-filter)
X           (+ fft-cost
X              (* 2 fft-cost
X                 (/ (1+ (- fft-size N-filter)) N-output))))
X      :direct
X      :fft)))
SHAR_EOF
chmod 0644 filtering.lisp ||
echo 'restore of filtering.lisp failed'
Wc_c="`wc -c < 'filtering.lisp'`"
test 48477 -eq "$Wc_c" ||
	echo 'filtering.lisp: original size 48477, current size' "$Wc_c"
fi
# ============= hacks.lisp ==============
if test -f 'hacks.lisp' -a X"$1" != X"-c"; then
	echo 'x - skipping hacks.lisp (File already exists)'
else
echo 'x - extracting hacks.lisp (Text)'
sed 's/^X//' << 'SHAR_EOF' > 'hacks.lisp' &&
;;;-*- Mode: LISP; Package: :SAPA; Syntax: COMMON-LISP -*-
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;
;  hacks.lisp
;
;  a hack to correct a deficiency in map-into as implemented in Genera 8.1
;  and in Allegro Common Lisp ...
;  Note:  before compiling and loading hacks.lisp,
;         you should compile and load
;            sapa-package.lisp
;
;  6/3/93
;
;  SAPA, Version 1.0; Copyright 1993, Donald B. Percival, All Rights Reserved
;
;  Use and copying of this software and preparation of derivative works
;  based upon this software are permitted.  Any distribution of this
;  software or derivative works must comply with all applicable United
;  States export control laws.
; 
;  This software is made available AS IS, and no warranty -- about the
;  software, its performance, or its conformity to any
;  specification -- is given or implied. 
;
;  Comments about this software can be addressed to dbp@apl.washington.edu
;-------------------------------------------------------------------------------
;;; (compile-file "ccl:SAPA;hacks.lisp")
;;; (load "ccl:SAPA;hacks.fasl")
;-------------------------------------------------------------------------------
(if (not (find-package :SAPA))
X  (defpackage :SAPA (:USE :COMMON-LISP)))
X
(in-package :SAPA)
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;; For some unknown reason, the function map-into is NOT defined
;;; in the package Common-Lisp in the version of Genera under which
;;; the SAPA package was tested (Symbolics Genera 8.1.1 on a Symbolics
;;; MacIvory model 3).  In that Lisp environment, there is a function
;;; called map-into in the packages Symbolics-Common-Lisp
;;; and Future-Common-Lisp, but these do not act in all respects like
;;; map-into as defined in Steele2 (and as implemented in Macintosh Common
;;; Lisp) -- an example of its deficiency is ilustrated in the commented-out
;;; material following the function definition.  Here is a hack that attempts
;;; to provide a map-into that is close enough to the definition in Steele2
;;; so that the functions in SAPA package will run properly.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
#+genera
(defun map-into (result function &rest rest)
X  (if (= (length result) (length (car rest)))
X    (apply #'scl:map-into result function rest)
X    (let ((temp (make-array (length (car rest)))))
X      (apply #'scl:map-into temp function rest)
X      (dotimes (i (min (length result) (length temp)) result)
X        (setf (elt result i) (svref temp i))))))
#|
;;; Evaluation of the following form in Genera causes an error
;;; because this version of map-into does not like the sequences
;;; to have different lengths -- Steele2 explicitly says that
;;; this is allowed.
(scl:map-into #(1 2 3) #'log #(1 10 100 1000))
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;; For some unknown reason, the function map-into is incorrectly defined
;;; in the package Common-Lisp in the version of Allegro Common Lisp under which
;;; the SAPA package was tested -- an example of its deficiency is ilustrated
;;; in the commented-out material.  Unfortunately, Allegro Common Lisp
;;; does not allow a redefinition of map-into, so I have defined the following
;;; hack to be used in situations where the incorrect definition makes a
;;; difference (currently just in spectral-window-for-direct-spectral-estimate).
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
#+allegro
(defun sapa-map-into (result function &rest rest)
X  (if (= (length result) (length (car rest)))
X    (apply #'map-into result function rest)
X    (let ((temp (make-array (length (car rest)))))
X      (apply #'map-into temp function rest)
X      (dotimes (i (min (length result) (length temp)) result)
X        (setf (elt result i) (svref temp i))))))
X
#|
;;; Evaluation of the following form in Allegro Common Lisp causes an error
;;; because this version of map-into does not allow the first sequence to be
;;; shorter than the second sequence -- Steele2 explicitly says that
;;; this is ok.
(map-into #(1 2 3) #'log #(1 10 100 1000))
|#
SHAR_EOF
chmod 0644 hacks.lisp ||
echo 'restore of hacks.lisp failed'
Wc_c="`wc -c < 'hacks.lisp'`"
test 4605 -eq "$Wc_c" ||
	echo 'hacks.lisp: original size 4605, current size' "$Wc_c"
fi
# ============= harmonic.lisp ==============
if test -f 'harmonic.lisp' -a X"$1" != X"-c"; then
	echo 'x - skipping harmonic.lisp (File already exists)'
else
echo 'x - extracting harmonic.lisp (Text)'
sed 's/^X//' << 'SHAR_EOF' > 'harmonic.lisp' &&
;;;-*- Mode: LISP; Package: :SAPA; Syntax: COMMON-LISP -*-
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;
;  harmonic.lisp
;
;  a collection of Lisp functions for harmonic analysis ...
;  Note:  before compiling and loading harmonic.lisp,
;         you should compile and load (in the order listed)
;            sapa-package.lisp, utilities.lisp, basic-math.lisp,
;            basic-statistics.lisp, dft-and-fft.lisp, tapers.lisp
;            and nonparametric.lisp
;
;  6/3/93
;
;  SAPA, Version 1.0; Copyright 1993, Donald B. Percival, All Rights Reserved
;
;  Use and copying of this software and preparation of derivative works
;  based upon this software are permitted.  Any distribution of this
;  software or derivative works must comply with all applicable United
;  States export control laws.
; 
;  This software is made available AS IS, and no warranty -- about the
;  software, its performance, or its conformity to any
;  specification -- is given or implied. 
;
;  Comments about this software can be addressed to dbp@apl.washington.edu
;-------------------------------------------------------------------------------
;;; (compile-file "ccl:SAPA;harmonic.lisp")
;;; (load "ccl:SAPA;harmonic.fasl")
;-------------------------------------------------------------------------------
(if (not (find-package :SAPA))
X  (defpackage :SAPA (:USE :COMMON-LISP)))
X
(in-package :SAPA)
X
(export '(;;; function to compute cosine and sine components of periodogram ...
X          periodogram-at-one-frequency
X
X          ;;; function to compute Fisher's g statistic
X          Fisher-g-statistic
X          ))
X
;-------------------------------------------------------------------------------
(defun periodogram-at-one-frequency
X       (time-series
X        freq
X        &key
X        (center-data t)
X        (start 0)
X        (end (length time-series))
X        (sampling-time 1.0))
X  "given
X   [1] time-series (required)
X       ==> a vector of time series values
X   [2] center-data (keyword; t)
X       ==> if t, subtract sample mean from time series;
X           if a number, subtracts that number from time series;
X           if nil, time-series is not centered
X   [3] start (keyword; 0)
X       ==> start index of vector to be used
X   [4] end (keyword; length of time-series)
X       ==> 1 + end index of vector to be used
X   [5] sampling-time (keyword; 1.0)
X       ==> sampling time (called delta t in the SAPA book)
returns
X   [1] value of periodogram at freq
X   [2] approximate conditional least squares estimate for A,
X       the amplitude of cosine term in Equation (461a)
X       of the SAPA book
X   [3] approximate conditional least squares estimate for B
X       the amplitude of sine term
---
Note: see Section 10.2 of the SAPA book"
X  (let ((N (- end start))
X        (j start)
X        (center-factor (cond
X                        ((numberp center-data) center-data)
X                        (center-data (sample-mean time-series
X                                                  :start start
X                                                  :end end))
X                        (t 0.0)))
X        (cos-sum 0.0)
X        (sin-sum 0.0)
X        (2-pi-f-deltt (* 2.0 pi freq sampling-time))
X        (t-index 1))
X    (dotimes (i N (values (* (/ sampling-time N)
X                             (+ (expt cos-sum 2)
X                                (expt sin-sum 2)))
X                          (* (/ 2.0 N) cos-sum)
X                          (* (/ 2.0 N) sin-sum)))
X      (incf cos-sum (* (- (aref time-series j) center-factor)
X                       (cos (* 2-pi-f-deltt t-index))))
X      (incf sin-sum (* (- (aref time-series j) center-factor)
X                       (sin (* 2-pi-f-deltt t-index))))
X      (incf j)
X      (incf t-index))))
X
#|
;;; For this example we use the 20 point time series of Section 6.16
;;; of the SAPA book:
(let ((20-pt-ts #(71.0  63.0  70.0  88.0  99.0  90.0 110.0 135.0 128.0 154.0
X                  156.0 141.0 131.0 132.0 141.0 104.0 136.0 146.0 124.0 129.0)))
X  (dolist (a-freq '(0.2 0.4 0.6 1.6 1.8 1.9 2.0) (values))
X    (multiple-value-bind (pgrm A B)
X                         (periodogram-at-one-frequency
X                          20-pt-ts
X                          a-freq
X                          :sampling-time 0.25)
X      (format t "~&~6,4F: ~8,4F ~8,4F  ~8,4F"
X              a-freq
X              (convert-to-dB pgrm)
X              A
X              B))))
;==>
0.2000:  30.8586 -19.3691  -24.4889
0.4000:  24.9217  10.7329  -11.5441
0.6000:  21.0531   0.0949  -10.0967
1.6000:  19.9132   8.8552    0.0457
1.8000:  11.6759   1.6742   -2.9941
1.9000:   9.1024   2.4743    0.6196
2.0000:   5.0515   1.6000    0.00000
|#
X
;-------------------------------------------------------------------------------
(defun Fisher-g-statistic
X       (time-series
X        &key
X        (center-data t)
X        (start 0)
X        (end (length time-series))
X        (alpha 0.05))
X    "given
X   [1] time-series (required)
X       ==> a vector of time series values
X   [2] center-data (keyword; t)
X       ==> if t, subtract sample mean from time series;
X           if a number, subtracts that number from time series;
X           if nil, time-series is not centered
X   [3] start (keyword; 0)
X       ==> start index of vector to be used
X   [4] end (keyword; length of time-series)
X       ==> 1 + end index of vector to be used
X   [5] alpha (keyword; 0.05)
X       ==> critical level at which Fisher's g test
X           is performed
returns
X   [1] Fisher's g statistic
X   [2] approximation to critical level of Fisher's g test
X       from Equation (491c)
X   [3] either :reject or :fail-to-reject, depending on whether
X       or not we reject or fail to reject the null hypothesis
X       of white noise
---
Note: see Section 10.9 of the SAPA book"
X  (let* ((N (- end start))
X         (m (if (oddp N)
X              (/ (1- N) 2)
X              (/ (- N 2) 2))))
X    (cond
X     ((<= m 1) (values 1.0 1.0 :fail-to-reject))
X     (t
X      (let* ((g_F (- 1.0 (expt (/ alpha m) (/ (1- m)))))
X             (the-periodogram (periodogram time-series
X                                           :center-data center-data
X                                           :start start
X                                           :end end
X                                           :N-nonzero-freqs :Fourier
X                                           :return-est-for-0-freq-p nil
X                                           :sdf-transformation nil
X                                           :return-frequencies-p nil))
X             (max-Shp (max-of-seq the-periodogram :end m))
X             (sum-Shp (sum the-periodogram :end m))
X             (g-statistic (/ max-Shp sum-Shp)))
X        (values g-statistic g_F (if (> g-statistic g_F)
X                                  :reject
X                                  :fail-to-reject)))))))
X
#|
;;; example discussed in Section 10.10 of the SAPA book ...
(let ((ocean-noise-data (vector 0.36201143 2.4737856 -1.018699 1.236099
X                                -0.38258758 3.590585 -0.9524991 -1.6762866
X                                1.3824866 -0.34308758 -0.1589 0.1412
X                                2.6574097 2.2194982 -2.7033095 2.0669982
X                                -1.3256867 -2.020322 -0.0753 0.4515114
X                                2.4969096 -0.26581144 -1.026599 0.31911144
X                                0.5307876 1.2934105 -2.6631095 -3.0335972
X                                -2.5025096 0.74761146 -0.26221144 -2.178198
X                                -2.940797 -1.5878867 -4.16842 1.5515105
X                                -1.136499 -1.879798 0.46968758 1.6587105
X                                2.198098 -2.3408096 0.133 -1.4207866
X                                -0.8249991 5.1213956 -0.64661145 -1.9440981 
X                                -0.5948876 1.041699 1.6077867 -2.4599857
X                                -2.989097 0.2540876 5.5098066 2.9239972 
X                                0.7211114 -2.033898 -2.2643094 3.2546084
X                                2.7131858 -0.896799 -0.8838991 2.9770973 
X                                -0.67578757 -0.43901145 1.5423105 0.836399
X                                0.0143 0.4059876 1.055099 -0.5244876 
X                                -3.155997 0.83619905 0.2799876 0.0655
X                                -2.3649857 -1.0694991 -1.3471105 1.3075105 
X                                0.6634876 -2.0554981 2.4471095 -2.187498
X                                0.6463876 1.2897105 2.4921856 -4.5647836 
X                                -1.5162104 1.3197105 2.2606857 -2.0319982
X                                -1.4144866 2.0942981 -1.5411105 1.4940104 
X                                -1.2181991 -0.1618 0.828699 -0.79309905
X                                2.6937857 -0.3346876 -1.5918106 -0.17650001 
X                                1.068099 0.1367 -4.7230077 -2.5424857
X                                -1.2936105 -4.000296 0.47421142 0.5830876 
X                                -2.6602857 -0.46411142 -0.3632876 4.088496
X                                -0.6841115 -2.8142972 -0.88159907 2.801497 
X                                2.081398 -3.6499846 1.860398 -3.6201086
X                                -0.6762114 2.5205097 -1.3579867 -1.4434105)))
X  (Fisher-g-statistic ocean-noise-data))
;==> 0.1336096017468551      ;;; Fisher's g statistic
;    0.10876132671244287     ;;; g_F from Equation (491c)
;    :reject                 ;;; reject because g > g_F
|#
SHAR_EOF
chmod 0644 harmonic.lisp ||
echo 'restore of harmonic.lisp failed'
Wc_c="`wc -c < 'harmonic.lisp'`"
test 9544 -eq "$Wc_c" ||
	echo 'harmonic.lisp: original size 9544, current size' "$Wc_c"
fi
# ============= index ==============
if test -f 'index' -a X"$1" != X"-c"; then
	echo 'x - skipping index (File already exists)'
else
echo 'x - extracting index (Text)'
sed 's/^X//' << 'SHAR_EOF' > 'index' &&
;;;-*- Mode: LISP; Package: :CL-USER; Syntax: COMMON-LISP -*-
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;
;  index for the sapaclisp contribution to StatLib, Version 1.0, 3 June 1993
;
;  sapaclisp is a collection of Common Lisp functions that can be used
;  to carry out many of the computations described in the book
;     "Spectral Analysis for Physical Applications: Multitaper and
;      Conventional Univariate Techniques",
;     by Donald B. Percival and Andrew T. Walden,
;     Cambridge University Press, Cambridge, England, 1993
;  (we refer to this book later on as "the SAPA book").
;  Questions concerning these Lisp functions can be addressed to
;  Don Percival via electronic mail at the Internet address
;             dbp@apl.washington.edu
;  or via traditional mail at the address
;             Donald B. Percival
;             Applied Physics Laboratory
;             HN-10
;             University of Washington
;             Seattle, WA  98195
;  The SAPA book uses a number of time series as examples
;  of various spectral analysis techniques.  The most important of
;  these series are also available from StatLib by sending the command
;             send sapa from datasets
;  to the Internet address
;             statlib@lib.stat.cmu.edu
; 
;  Version 1.0 of the sapaclisp contribution consists of this index
;  and the following 16 files:
;   [ 1]  sapa-package
;   [ 2]  hacks
;   [ 3]  utilities
;   [ 4]  basic-math
;   [ 5]  matrix
;   [ 6]  basic-statistics
;   [ 7]  dft-and-fft
;   [ 8]  tapers
;   [ 9]  filtering
;   [10]  random
;   [11]  acvs
;   [12]  parametric
;   [13]  nonparametric
;   [14]  multitaper
;   [15]  harmonic
;   [16]  examples
;  All 16 files can be obtaining by sending the command
;             send everything from sapaclisp
;  to the Internet address
;             statlib@lib.stat.cmu.edu
;  Individual files can be obtained by sending a command such as
;             send dft-and-fft from sapaclisp
;  to the same address.
;
;  Here is a brief description of the contents of each file:
;   [ 1]  sapa-package contains Lisp forms needed to create
;         the package SAPA (packages are mysterious and tricky
;         creatures in Common Lisp that serve a noble purpose,
;         namely, to make sure that different collections of
;         Lisp functions can gracefully live together).  All
;         of the Lisp functions in the sapaclisp contribution
;         live in the SAPA package.  You should compile and load
;         this file before trying to do ANYTHING from within Lisp
;         with the other files (including viewing them with a
;         Lisp-oriented text editor!).
;   [ 2]  hacks at present patches up faulty definitions for map-into
;         in Allegro Common Lisp and in the Genera 8.1 version of Common Lisp
;         (if you are using Macintosh Common Lisp, you can ignore this file).
;   [ 3]  utilities contains a collection of utility functions
;         that are used extensively throughout SAPA.  These are mostly
;         functions for manipulating various sequences of data, but
;         there are also functions for converting back and from
;         decibels, a Lisp implementation of the Fortran sign function, etc.
;   [ 4]  basic-math contains functions for supporting certain
;         basic mathematical operations such as computing the log
;         of the gamma function, manipulating polynomials, root finding
;         and simple numerical integration.
;   [ 5]  matrix contains functions for manipulating matrices, including
;         the Cholesky and modified Gram-Schmidt (i.e., Q-R) decompositions.
;   [ 6]  basic-statistics contains functions to carry out basic
;         statistical operations, such as sample means and variances,
;         sample medians, computation of quantiles from various distributions,
;         linear least squares, etc.
;   [ 7]  dft-and-fft contains functions for computing the discrete Fourier
;         transform of a vector of numbers via a fast Fourier transform
;         algorithm or a chirp transform algorithm.
;   [ 8]  tapers contains functions for computing the cosine and dpss
;         data tapers and applying them to a time series.
;   [ 9]  filtering contains functions that approximate ideal low-pass,
;         high-pass and band-pass filters using techniques outlined
;         in Chapter 5 of the SAPA book.  It also has functions for
;         some simple smoothers.
;   [10]  random contains functions for generating realizations of
;         various stationary processes.
;   [11]  acvs contains functions for computing the sample autocovariance
;         sequence and sample variogram for a time series.
;   [12]  parametric contains functions for computing autoregressive
;         spectral estimates using a variety of techniques, including
;         the Yule-Walker method, Burg's algorithm, forward least squares
;         and forward/backward least squares (these are discussed in
;         Chapter 9 of the SAPA book).
;   [13]  nonparametric contains functions for computing nonparametric
;         spectral estimates, including the periodogram, direct spectral
;         estimates, lag window spectral estimates and WOSA spectral
;         estimates.  There are also functions for computing the sample
;         cepstrum, the time series bandwidth and the cumulative periodogram
;         test statistic for white noise (all of these are discussed
;         in Chapter 6 of the SAPA book).
;   [14]  multitaper contains functions for computing the orthonormal
;         dpss data tapers and both the simple and adaptive multitaper
;         spectral estimates (discussed in Chapters 7 and 8 of the SAPA book).
;   [15]  harmonic contains functions for computing the cosine and
X          sine component of the periodogram and Fisher's g statistic
X          (discussed in Chapter 10 of the SAPA book).
;   [16]  examples contains just that: some examples of how to use
;         the functions in SAPA to reproduce some of the results in
;         the SAPA book.
;  All exported functions in SAPA have a documentation line that tells
;  how to use them.  In addition, almost all functions have a brief
;  example of their use immediately following the definition of the
;  the function.  You can repeat these examples by evaluating the
;  appropriate Lisp forms.  Note: some of these Lisp forms return
;  a short vector, the contents of which can be examined easily
;  by setting the Common Lisp variable *print-array* to t prior to
;  evaluating the Lisp form (see page 565 of the second edition of
;  "Common Lisp: The Language" by Guy Steele, Digital Press, 1990 --
;  we refer to this book later on as "Steele2").
;
;  Here are the steps to take to prepare the SAPA files for use.
;   [ 1]  Retrieve all 16 files from StatLib, and place them on your
;         computer in 16 individual files.  As an example, let's
;         assume that the files are named
;              sapa-package.lisp
;              hacks.lisp
;              utilities.lisp
;                  ...
;              harmonic.lisp
;              examples.lisp
;         The first line in each file should have some gunk on it
;         concerning Mode, Package and Syntax.
;   [ 2]  Bring up Lisp on your computer.
;   [ 3]  Compile sapa-package.lisp, and then load it into your
;         Lisp world.  Here are examples of two Lisp forms that
;         will do this if you evaluate them:
;              (compile-file "sapa-package.lisp")
;              (load "sapa-package.fasl")
;         Here the function compile-file compiles the file sapa-package.lisp,
;         and the results of the compilation are placed in sapa-package.fasl.
;         The compiled version is then loaded into the Lisp world using
;         the load function.  Note that the ".lisp" and ".fasl" extensions
;         in this example may need to be changed to match the naming
;         conventions for source files and compiled binary files
;         in Lisp systems other than Macintosh Common Lisp 2.0.
;   [ 4]  Compile and load all the remaining files EXCEPT examples.lisp
;         in the order in which they are listed above (you can actually
;         compile and load these remaining files in any order if you
;         are willing to put up with snippy comments from certain compilers
;         about undefined functions).
;   [ 5]  To use the functions in SAPA, you must evaluate the following
;         Lisp form
;              (use-package :SAPA)
;         from within whatever package you plan to use the functions in the
;         SAPA package (usually this will be the package CL-USER).
;         This Lisp form is the second form in the file examples.lisp,
;         so you can just evaluate it there (that file does things in the
;         package CL-USER).  You can evaluate various forms given
;         in examples.lisp to learn how to use some of the SAPA functions.
;
;  For the record, the functions in the SAPA package have been tested
;  successfully in the following versions of Common Lisp:
;   [a]  Macintosh Common Lisp 2.0p2
;   [b]  Symbolics Genera 8.1.1 on a Symbolics MacIvory model 3
;   [c]  Allegro Common Lisp running under SunOS Release 4.1.2
;
;
;  Finally, here is some legal nonsense, shamelessly stolen from someone else
;  and placed here in an attempt to fend off lawyers (I don't know what it
;  means; I don't care what it means; I wasn't in my right mind when I did
;  all of this; I've NEVER been in my right mind; I'm not responsible
;  for anything I've ever done; the Devil made me do it; ... ):
;
;  SAPA, Version 1.0; Copyright 1993, Donald B. Percival, All Rights Reserved
;
;  Use and copying of this software and preparation of derivative works
;  based upon this software are permitted.  Any distribution of this
;  software or derivative works must comply with all applicable United
;  States export control laws.
; 
;  This software is made available AS IS, and no warranty -- about the
;  software, its performance, or its conformity to any
;  specification -- is given or implied. 
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
SHAR_EOF
chmod 0644 index ||
echo 'restore of index failed'
Wc_c="`wc -c < 'index'`"
test 10351 -eq "$Wc_c" ||
	echo 'index: original size 10351, current size' "$Wc_c"
fi
# ============= matrix.lisp ==============
if test -f 'matrix.lisp' -a X"$1" != X"-c"; then
	echo 'x - skipping matrix.lisp (File already exists)'
else
echo 'x - extracting matrix.lisp (Text)'
sed 's/^X//' << 'SHAR_EOF' > 'matrix.lisp' &&
;;;-*- Mode: LISP; Package: :SAPA; Syntax: COMMON-LISP -*-
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;
;  matrix.lisp
;
;  a collection of Lisp functions for some basic matrix operations ...
;  Note:  before compiling and loading matrix.lisp,
;         you should compile and load
;         sapa-package.lisp
;
;  6/3/93
;
;  SAPA, Version 1.0; Copyright 1993, Donald B. Percival, All Rights Reserved
;
;  Use and copying of this software and preparation of derivative works
;  based upon this software are permitted.  Any distribution of this
;  software or derivative works must comply with all applicable United
;  States export control laws.
; 
;  This software is made available AS IS, and no warranty -- about the
;  software, its performance, or its conformity to any
;  specification -- is given or implied. 
;
;  Comments about this software can be addressed to dbp@apl.washington.edu
;-------------------------------------------------------------------------------
;;; (compile-file "ccl:SAPA;matrix.lisp")
;;; (load "ccl:SAPA;matrix.fasl")
;-------------------------------------------------------------------------------
(if (not (find-package :SAPA))
X  (defpackage :SAPA (:USE :COMMON-LISP)))
X
(in-package :SAPA)
X
(export '(;;; functions for printing arrays ...
X          print-1-or-2-d-array
X          print-rows-of-nxm-matrix
X
X          ;;; functions for converting to/from linearized storage ...
X          linearized-upper-trangular->2d-matrix
X          2d-matrix->linearized-upper-trangular
X
X          ;;; functions to do elementary matrix operations ...
X          transpose
X          Hermitian-transpose
X          zero-strict-lower-diagonal!
X          multiply-two-matrices
X          multiply-matrix-and-vector
X          multiply-matrix-and-scalar
X          subtract-two-matrices
X          trace-matrix
X          2d-matrix-move!
X          
X          ;;; functions to do Cholesky and Gram-Schmidt (Q-R) decompositions ...
X          cholesky!
X          spofa!
X          sposl!
X          modified-Gram-Schmidt!
X          Q-R!
X
X          ;;; functions to solve triangular system of equations ...
X          upper-triangular-solve!
X          lower-triangular-solve!
X          ))
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  print-1-or-2-d-array
;;;                 print-rows-of-nxm-matrix
;;;  can be used to print out the elements of one and two dimensional array.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun print-1-or-2-d-array
X       (an-array
X        &key
X        (tag "array")
X        (format-control-for-array-element "~F"))
X  "given
X   [1] an-array (required)
X       ==> a one-dimensional or two-dimensional array
X   [2] tag (keyword; `array')
X       ==> an optional tag to be printed
X           along with array elements
X   [3] format-control-for-array-element (keyword; `~F')
X       ==> format control for a single array element
prints the elements of an-array and
returns
X   [1] an-array"
X  (let ((dims (array-dimensions an-array))
X        (for-format (concatenate 'string
X                                 "~&~A, element ~A: "
X                                 format-control-for-array-element)))
X    (cond
X     ((= (length dims) 1)
X      ;;; a vector ...
X      (dotimes (i (car dims) (values an-array))
X        (format t for-format tag i (aref an-array i))))
X     ((= (length dims) 2)
X      (dotimes (i (nth 0 dims) (values an-array))
X        (dotimes (j (nth 1 dims))
X          (format t for-format tag (list i j) (aref an-array i j)))))
X     (t
X      (error "can't print ~A" an-array)))))
X
#|
(print-1-or-2-d-array #(1 2 3 4))
;==>
array, element 0: 1.0
array, element 1: 2.0
array, element 2: 3.0
array, element 3: 4.0
#(1 2 3 4)
X
(print-1-or-2-d-array (make-array '(3 2)
X                                  :initial-contents
X                                  '((1 2) (3 4) (5 6))))
;==>
array, element (0 0): 1.0
array, element (0 1): 2.0
array, element (1 0): 3.0
array, element (1 1): 4.0
array, element (2 0): 5.0
array, element (2 1): 6.0
#2a((1 2) (3 4) (5 6))
|#
X
;-------------------------------------------------------------------------------
(defun print-rows-of-nxm-matrix
X       (an-nxm-matrix
X        &key
X        (format-control-for-array-element "~F "))
X  "given
X   [1] an-nxm-matrix (required)
X       ==> a two-dimensional array
X   [3] format-control-for-array-element (keyword; `~F')
X       ==> format control for a single array element
prints the elements of the 2d array and
returns
X   [1] an-nxm-matrix"
X  (let ((dimensions (array-dimensions an-nxm-matrix)))
X    (assert (= (length dimensions) 2))
X    (let ((m-columns (elt dimensions 1)))
X      (dotimes (i-row (elt dimensions 0))
X        (format t "~&")
X        (dotimes (j-col m-columns)
X          (format t format-control-for-array-element
X                  (aref an-nxm-matrix i-row j-col))))
X      (format t "~&")
X      (values an-nxm-matrix))))
X
#|
(print-rows-of-nxm-matrix (make-array '(3 2)
X                                      :initial-contents
X                                      '((1 2) (3 4) (5 6))))
;==>
1.0 2.0 
3.0 4.0 
5.0 6.0
#2a((1 2) (3 4) (5 6))
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  linearized-upper-trangular->2d-matrix
;;;                 2d-matrix->linearized-upper-trangular
;;;  convert back and forth between a straightforward representation for 
;;;  a 2d square matrix and a scheme for storing the upper triangular
;;;  portion of a two dimensional matrix in a vector.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun linearized-upper-trangular->2d-matrix
X       (A
X        &key
X        (result (let ((N (round (/ (1- (sqrt (1+ (* 8 (length A)))))
X                                   2))))
X                  (make-array `(,N ,N)))))
X  "given
X   [1] A (required)
X       ==> a vector with an upper triangular matrix
X           stored in linearized form
X   [2] result (keyword; new 2d array of appropriate size)
X       <== a 2d hermitian matrix filled according
X           to contents of A
converts linearized form A into 2d hermitian array and
returns
X   [1] result, a 2d hermitian matrix"
X  (let ((L -1))
X    (dotimes (i (car (array-dimensions result)) result)
X      (dotimes (j (1+ i))
X        (setf (aref result i j)
X              (conjugate (setf (aref result j i)
X                               (aref A (incf L)))))))))
X
#|
(linearized-upper-trangular->2d-matrix
X #(1 #C(1 1) 1))
;==> #2a((1 #c(1 1)) (#c(1 -1) 1))
|#
X
;-------------------------------------------------------------------------------
(defun 2d-matrix->linearized-upper-trangular
X       (A
X        &key
X        (result (let ((N (car (array-dimensions A))))
X                  (make-array (/ (* N (1+ N)) 2)))))
X  "given
X   [1] A (required)
X       ==> a 2d hermitian matrix
X   [2] result (keyword; new vector of appropriate size)
X       <== a vector filled linearly with elements
X           of A
converts 2d hermitian array into linearized form and
returns
X   [1] result, a vector with 2d hermitian array
stored in linearized form"
X  (let ((L -1))
X    (dotimes (i (car (array-dimensions A)) result)
X      (dotimes (j (1+ i))
X        (setf (aref result (incf L))
X              (aref A j i))))))
X
#|
(2d-matrix->linearized-upper-trangular
X (linearized-upper-trangular->2d-matrix
X   #(1 #C(1 1) 1)))
;==> #(1 #c(1 1) 1)
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  transpose
;;;                 Hermitian-transpose
;;;                 zero-strict-lower-diagonal!
;;;                 multiply-two-matrices
;;;                 multiply-matrix-and-vector
;;;                 multiply-matrix-and-scalar
;;;                 subtract-two-matrices
;;;                 trace-matrix
;;;                 2d-matrix-move!
;;;  perform fairly simple operations on matrices and vectors.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun transpose
X       (a-matrix
X        &key
X        (result
X         (make-array
X          (reverse (array-dimensions a-matrix)))))
X  "given
X   [1] A (required)
X       ==> a 2d matrix
X   [2] result (keyword; new 2d array of appropriate size)
X       <== a 2d matrix to contain transpose of a-matrix
returns
X   [1] transpose of a-matrix (placed in result)"
X  (let ((list-of-two-integers (array-dimensions a-matrix)))
X    (dotimes (i (nth 0 list-of-two-integers) result)
X      (dotimes (j (nth 1 list-of-two-integers))
X        (setf (aref result j i)
X              (aref a-matrix i j))))))
X
#|
(transpose #2a((1 2) (3 4) (5 6)))
;==> #2a((1 3 5) (2 4 6))
|#
X
;-------------------------------------------------------------------------------
(defun Hermitian-transpose
X       (a-matrix
X        &key
X        (result
X         (make-array
X          (reverse (array-dimensions a-matrix)))))
X  "given
X   [1] A (required)
X       ==> a 2d matrix
X   [2] result (keyword; new 2d array of appropriate size)
X       <== a 2d matrix to contain Hermitian transpose of a-matrix
returns
X   [1] Hermitian transpose of a-matrix (placed in result)"
X  (let ((list-of-two-integers (array-dimensions a-matrix)))
X    (dotimes (i (nth 0 list-of-two-integers) result)
X      (dotimes (j (nth 1 list-of-two-integers))
X        (setf (aref result j i)
X              (conjugate (aref a-matrix i j)))))))
X
#|
(Hermitian-transpose #2a((1 2) (3 4) (5 6)))
;==> #2a((1 3 5) (2 4 6))
(Hermitian-transpose #2a((1 #c(2 1)) (#c(3 -1) 1)))
;==>                 #2a((1 #c(3 1)) (#c(2 -1) 1))
|#
X
;-------------------------------------------------------------------------------
(defun zero-strict-lower-diagonal!
X       (a-matrix)
X  "given a square matrix,
zeros its lower diagonal and returns
the modified square matrix"
X  (let* ((n (nth 0 (array-dimensions a-matrix)))
X         (n-in-column-to-zap (1- n))
X         (row-start 1))
X    (dotimes (j-column (1- n) a-matrix)
X      (dotimes (i n-in-column-to-zap)
X        (setf (aref a-matrix (+ i row-start) j-column ) 0.0))
X      (incf row-start)
X      (decf n-in-column-to-zap))))
X
#|
(zero-strict-lower-diagonal! #2a((1 #c(2 1)) (#c(3 -1) 1)))
;==>                         #2a((1 #c(2 1)) (0.0 1))
(zero-strict-lower-diagonal! #2a((1 2 3)   (4 5 6)     (7 8 9)))
;==>                         #2a((1 2 3) (0.0 5 6) (0.0 0.0 9))
|#
X
X
;-------------------------------------------------------------------------------
(defun multiply-two-matrices
X       (a-matrix
X        b-matrix
X        &key
X        (result
X         (make-array
X          (list (nth 0 (array-dimensions a-matrix))
X                (nth 1 (array-dimensions b-matrix))))))
X  "given
X   [1] a-matrix (required)
X       ==> a 2d matrix
X   [2] b-matrix (required)
X       ==> another 2d matrix, with dimensions such that
X           the product of a-matrix and b-matrix is defined
X   [3] result (keyword; new 2d array of appropriate size)
X       <== a 2d matrix to contain product of two matrices
returns
X   [1] product of two matrices (placed in result)"
X  (let ((m (nth 0 (array-dimensions a-matrix)))
X        (n (nth 1 (array-dimensions b-matrix)))
X        (common (nth 0 (array-dimensions b-matrix))))
X    (dotimes (i m result)
X      (dotimes (j n)
X        (setf (aref result i j) 0.0)
X        (dotimes (k common)
X          (incf (aref result i j)
X                (* (aref a-matrix i k) (aref b-matrix k j))))))))
X
#|
(multiply-two-matrices #2a((0 0 1) (0 1 0) (1 0 0))
X                       #2a((10 9) (8 7) (6 5)))
;==> #2a((6.0 5.0) (8.0 7.0) (10.0 9.0))
|#
X
;-------------------------------------------------------------------------------
(defun multiply-matrix-and-vector
X       (a-matrix
X        b-vector
X        &key
X        (result
X         (make-array
X          (nth 0 (array-dimensions a-matrix)))))
X  "given
X   [1] a-matrix (required)
X       ==> a 2d matrix
X   [2] b-vector (required)
X       ==> a vector, with dimensions such that
X           the product of a-matrix and b-vector is defined
X   [3] result (keyword; new vector of appropriate size)
X       <== a vector to contain product of a-matrix and b-vector
returns
X   [1] product of a-matrix and b-vector (placed in result)"
X  (let ((m (nth 0 (array-dimensions a-matrix)))
X        (n (length b-vector)))
X    (dotimes (i m result)
X      (setf (aref result i) 0.0)
X      (dotimes (j n)
X        (incf (aref result i)
X              (* (aref a-matrix i j) (aref b-vector j)))))))
X
#|
(multiply-matrix-and-vector #2a((0 0 1) (0 1 0) (1 0 0))
X                            #(10 9 8))
;==> #(8.0 9.0 10.0)
|#
X
;-------------------------------------------------------------------------------
(defun multiply-matrix-and-scalar
X       (a-matrix
X        scalar
X        &key
X        (result
X         (make-array
X          (array-dimensions a-matrix))))
X  "given
X   [1] a-matrix (required)
X       ==> a 2d matrix
X   [2] scalar (required)
X       ==> an arbitrary number
X   [3] result (keyword; new matrix of same size as a-matrix)
X       <== a matrix to contain product of a-matrix and scalar
returns
X   [1] product of a-matrix and scalar (placed in result)"
X  (let ((m (nth 0 (array-dimensions a-matrix)))
X        (n (nth 1 (array-dimensions a-matrix))))
X    (dotimes (i m result)
X      (dotimes (j n)
X        (setf (aref result i j)
X              (* scalar (aref a-matrix i j)))))))
X
#|
(multiply-matrix-and-scalar #2a((0 0 1) (0 1 0) (1 0 0))
X                            1/3)
;==> #2a((0 0 1/3) (0 1/3 0) (1/3 0 0))
|#
X
;-------------------------------------------------------------------------------
(defun subtract-two-matrices
X       (a-matrix
X        b-matrix
X        &key
X        (result
X         (make-array (array-dimensions a-matrix))))
X  "given
X   [1] a-matrix (required)
X       ==> a 2d matrix
X   [2] b-matrix (required)
X       ==> a 2d matrix, with dimensions the same
X           as a-matrix
X   [3] result (keyword; new vector of appropriate size)
X       <== a matrix to contain result of subtracting
X           b-matrix from a-matrix
returns
X   [1] a-matrix minus b-matrix (placed in result)"
X  (let ((m (nth 0 (array-dimensions a-matrix)))
X        (n (nth 1 (array-dimensions a-matrix))))
X    (dotimes (i m result)
X      (dotimes (j n)
X        (setf (aref result i j)
X              (- (aref a-matrix i j) (aref b-matrix i j)))))))
X
#|
(subtract-two-matrices #2a((1 0 0) (0 1 0) (0 0 1))
X                       #2a((0 0 1) (0 1 0) (1 0 0)))
;==> #2a((1 0 -1) (0 0 0) (-1 0 1))
|#
X
;-------------------------------------------------------------------------------
(defun trace-matrix (a-matrix)
X  "given a square matrix,
returns its trace"
X  (let ((m (nth 0 (array-dimensions a-matrix)))
X        (n (nth 1 (array-dimensions a-matrix))))
X    (assert (= m n))
X    (let ((sum 0.0))
X      (dotimes (i m sum)
X        (incf sum (aref a-matrix i i))))))
X
;;; (trace-matrix #2a((1 2 3) (4 5 6) (7 8 9)))  ;==> 15.0
X
;-------------------------------------------------------------------------------
(defun 2d-matrix-move!
X       (from-this
X        to-this)
X  "given
X   [1] from-this (required)
X       ==> a 2d matrix
X   [2] to-this (required)
X       <== another 2d matrix
transfer contents of from-this to corresponding
locations in to-this, and
returns
X   [1] to-this"
X  (let* ((temp (array-dimensions from-this))
X         (n-rows (nth 0 temp))
X         (n-columns (nth 1 temp)))
X    (dotimes (i n-rows to-this)
X      (dotimes (j n-columns)
X        (setf (aref to-this i j) (aref from-this i j))))))
X
#|
(2d-matrix-move! #2a((1 2 3) (4 5 6) (7 8 9)) (make-array '(3 3)))
;==> #2a((1 2 3) (4 5 6) (7 8 9))
(2d-matrix-move! #2a((1 2 3) (4 5 6) (7 8 9)) (make-array '(4 4)))
;==> #2a((1 2 3 nil) (4 5 6 nil) (7 8 9 nil) (nil nil nil nil))
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  cholesky!
;;;                 spofa!
;;;                 sposl!
;;;                 modified-Gram-Schmidt!
;;;                 Q-R!
;;;  carry out Cholesky and modified-Gram-Schmidt factorizations of a two
;;;  dimensional matrix.  These functions destroy whatever matrices
;;;  are given to them, and hence all have names that end with a ``!''.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun cholesky! (A b &key (eps single-float-epsilon))
X  "given
X   [1] A (required)
X       ==> vector of length N*(N+1)/2 representing
X           a two dimensional N by N positive definite matrix
X           with elements stored columnwise; i.e, elements with indices
X           0     1     2     3     4     5    ...
X           correspond to matrix positions
X           1,1   1,2   2,2   1,3   2,3   3,3  ...
X   [2] b (required)
X       <=> right-hand vector on input;
X           solution X on output
X   [3] eps (keyword; single-float-epsilon)
X       number used to test for a computationally singular matrix
solves A X = b using the Cholesky decomposition method, and
returns
X   [1] X, the solution vector (stored in b)"
X  (let* ((N (length b))
X         (N-1 (1- N))
X         (XL (make-array `(,N ,N)))
X         (Y (make-array N))
X         (D (make-array N))
X         (L 0))
X    ;(declare (dynamic-extent XL Y D))
X    (assert (= (length A) (/ (* N (1+ N)) 2)))
X    ;;; Factor into triangular and diagonal form ...
X    (setf (aref D 0) (realpart (aref A 0)))
X    (dotimes (i N-1)
X      (dotimes (j (1+ i))
X        (incf L)
X        (setf (aref XL (1+ i) j) (/ (conjugate (aref A L)) (aref D j)))
X        (when (not (zerop j))
X          (dotimes (k j)
X            (decf (aref XL (1+ i) j)
X                  (/ (* (aref XL (1+ i) k)
X                        (conjugate (aref XL j k))
X                        (aref D k))
X                     (aref D j))))))
X      (incf L)
X      (setf (aref D (1+ I)) (realpart (aref A L)))
X      (dotimes (k (1+ I))
X        (decf (aref D (1+ I))
X              (* (aref D k) (expt (abs (aref XL (1+ i) k)) 2))))
X      ;;; Test for nonpositive value (i.e., matrix is too close to
X      ;;; being singular)
X      (if (< (aref D (1+ I)) eps)
X        (error "(aref D ~D) = ~F < ~F = esp~&         matrix is computationally singular"
X               (1+ I) (aref D (1+ I)) eps)))
X    ;;; Solve for intermediate column vector solution ...
X    (setf (aref Y 0) (aref b 0))
X    (dotimes (k N-1)
X      (setf (aref Y (1+ k)) (aref b (1+ k)))
X      (dotimes (j (1+ k))
X        (decf (aref Y (1+ k))
X              (* (aref XL (1+ k) j) (aref Y j)))))
X    ;;; Solve for final column vector solution ...
X    (setf (aref b N-1) (/ (aref Y N-1) (aref D N-1)))
X    (let ((k-index N-1))
X      (dotimes (k-count N-1 (values b))
X        (decf k-index)
X        (setf (aref b k-index) (/ (aref Y k-index) (aref D k-index)))
X        (dotimes (j (- N k-index 1))
X          (decf (aref b k-index)
X                (* (conjugate (aref XL (+ k-index j 1) k-index))
X                   (aref b (+ k-index j 1)))))))))
X
#|
;;; Example inspired by Yule-Walker equations for AR(1) model
;;; with phi = 0.9
(2d-matrix->linearized-upper-trangular
X #2a((1.0 0.9 0.81) (0.9 1.0 0.9) (0.81 0.9 1.0)))
;==> #(1.0 0.9 1.0 0.81 0.9 1.0)
(cholesky! #(1.0 0.9 1.0 0.81 0.9 1.0)
X           #(0.9 0.81 0.729))
;==> #(0.9 5.258951169277058E-16 -5.84327907697451E-16)
|#
X
;-------------------------------------------------------------------------------
(defun spofa!
X       (A
X        &key
X        (n (car (array-dimensions a))))
X  "given
X   [1] A (required)
X       <=> real symmetric positive definite matrix
X           (this gets trashed)
X   [2] n (keyword; (array-dimensions a))
X       ==> size of matrix to be factored
calculates the upper triangular matrix U
of the Cholesky decomposition A = LU,
where L = U^T, and
returns
X   [1] U, stored in upper triangular part of A"
X  (let (info temp sum)
X    (dotimes (j n (values A 0))
X      (setf info j
X            sum 0.0)
X      (when (plusp j)
X        (dotimes (k j)
X          (setf temp (/ (- (aref A k j)
X                           (let ((dot-product 0.0))
X                             (dotimes (m k dot-product)
X                               (incf dot-product
X                                     (* (aref A m j) (aref A m k))))))
X                        (aref A k k))
X                (aref A k j) temp)
X          (incf sum (* temp temp))))
X      (setf sum (- (aref A j j) sum))
X      (if (not (plusp sum)) (return (values A (1+ info))))
X      (setf (aref A j j) (sqrt sum)))))
X
#|
(setf *U* (spofa! #2a((1.0 0.9 0.81) (0.9 1.0 0.9) (0.81 0.9 1.0))))
;==> #2a((1.0 0.9 0.81) (0.9 0.4358898943540673 0.39230090491866054) (0.81 0.9 0.4358898943540673))
(zero-strict-lower-diagonal! *U*)
;==> #2a((1.0 0.9 0.81) (0.0 0.4358898943540673 0.39230090491866054) (0.0 0.0 0.4358898943540673))
(multiply-two-matrices (transpose *U*) *U*)
;==> #2a((1.0 0.9 0.81) (0.9 1.0 0.9) (0.81 0.9 1.0))
|#
X
;-------------------------------------------------------------------------------
(defun sposl!
X       (A
X        b
X        &key
X        (N (length b)))
X  "given
X   [1] A (required)
X       ==> real symmetric positive definite matrix
X           AFTER it has been crunched by spofa!
X   [2] b (required)
X       <=> the right-hand side vector on input;
X            on output, this gets replaced by the solution X
X   [3] N (keyword; length of b)
X       ==> order of matrix A (default is length of b)
returns
X   [1] X, the solution to A X = b; note that X is the same as
X       the vector which contained b"
X  (let ((kb (1- N))
X        temp)
X    (dotimes (k-col N)
X      (setf (aref b k-col) (/ (- (aref b k-col)
X                                 (if (zerop k-col) 0.0
X                                     (let ((dot-product 0.0))
X                                       (dotimes (j-row k-col dot-product)
X                                         (incf dot-product
X                                               (* (aref A j-row k-col)
X                                                  (aref b j-row)))))))
X                              (aref a k-col k-col))))
X    (dotimes (k N b)
X      (setf (aref b kb) (/ (aref b kb) (aref a kb kb)))
X      (setf temp (- (aref b kb)))
X      (dotimes (j kb)
X        (setf (aref b j) (+ (* temp (aref a j kb))
X                            (aref b j))))
X      (decf kb))))
X
#|
(setf *U* (spofa! #2a((1.0 0.9 0.81) (0.9 1.0 0.9) (0.81 0.9 1.0))))
;==> #2a((1.0 0.9 0.81) (0.9 0.4358898943540673 0.39230090491866054) (0.81 0.9 0.4358898943540673))
(sposl! *U* #(0.9 0.81 0.729))
;==> #(0.9 5.25895116927706E-16 -5.8432790769745105E-16)
|#
X
;-------------------------------------------------------------------------------
(defun modified-Gram-Schmidt! (A-matrix)
X  "given
X   [1] A (required)
X       <=> a matrix of real or complex-valued numbers with
X           size m by n with rank N
computes the factorization A = QR, where
X  Q is of size m by n and has orthonormal columns and
X  R is of size n by n and is upper triangular, and
returns
X   [1] Q
X   [2] R
---
Note: see Algorithm 5.2.5, p. 219, of Golub and Van Loan,
Second Edition, with simple changes to handle matrices
with real or complex-valued numbers"
X  (let* ((dim-of-A-matrix (array-dimensions A-matrix))
X         (m (car dim-of-A-matrix))
X         (n (cadr dim-of-A-matrix))
X         (Q (make-array dim-of-A-matrix))
X         (R (make-array `(,n ,n) :initial-element 0.0))
X         j)
X    (dotimes (k n (values Q R))
X      (setf (aref R k k)
X            (let ((sum 0.0))
X              (dotimes (i m (sqrt sum))
X                (incf sum (realpart
X                           (* (conjugate (aref A-matrix i k))
X                              (aref A-matrix i k)))))))
X      (dotimes (i m)
X        (setf (aref Q i k) (/ (aref A-matrix i k)
X                              (aref R k k))))
X      (setf j k)
X      (dotimes (j-shifted (1- (- n k)))
X        (incf j)
X        (setf (aref R k j)
X              (let ((sum 0.0))
X                (dotimes (i m sum)
X                  (incf sum (* (conjugate (aref Q i k))
X                               (aref A-matrix i j))))))
X        (dotimes (i m)
X          (setf (aref A-matrix i j) (- (aref A-matrix i j)
X                                       (* (aref Q i k)
X                                          (aref R k j)))))))))
X
#|
;;; example from pages 219--20 of Golub and Van Loan:
(multiple-value-bind (Q R)
X                     (modified-Gram-Schmidt!
X                      (make-array '(3 2)
X                                  :initial-contents
X                                  '((1.0   1.0)
X                                    (0.001 0.0)
X                                    (0.0   0.001))))
X  (print-rows-of-nxm-matrix Q
X                            :format-control-for-array-element "~A ")
X  (terpri)
X  (print-rows-of-nxm-matrix (multiply-two-matrices Q R)
X                            :format-control-for-array-element "~A ")
X  (terpri)
X  (print-rows-of-nxm-matrix (multiply-two-matrices
X                             (transpose Q) Q)
X                            :format-control-for-array-element "~A ")
X  (values))
X
;==> 0.999999500000375    7.0710625089215E-4 
X     9.99999500000375E-4 -0.7071062508568814 
X     0.0                  0.7071069579631323 
X
X    1.0   1.0            ;recovers original matrix perfectly!
X    0.001 0.0 
X    0.0   0.001
X
X    1.0                     3.5268662990084465E-14 ;Q^T Q is close to identity
X    3.5268662990084465E-14  0.9999999999999998 
X
;;; results look quite reasonable and agree fairly well with Golub and Van Loan
X
;;; complex-valued test case ...
(multiple-value-bind (Q R)
X                     (modified-Gram-Schmidt!
X                      (make-array '(3 2)
X                                  :initial-contents
X                                  `(,(list #C( 0.40586 -1.24797)
X                                           #C(-0.21314 -1.13773))
X                                    ,(list #C( 0.28417 -0.82029)
X                                           #C( 0.72524 -0.73988))
X                                    ,(list #C(-0.27198 -1.11286)
X                                           #C( 1.48759 -1.290507)))))
X  (print-rows-of-nxm-matrix Q
X                            :format-control-for-array-element "~A ")
X  (terpri)
X  (print-rows-of-nxm-matrix (multiply-two-matrices Q R)
X                            :format-control-for-array-element "~A ")
X  (terpri)
X  (print-rows-of-nxm-matrix (multiply-two-matrices
X                             (Hermitian-transpose Q) Q)
X                            :format-control-for-array-element "~A ")
X  (values))
;==> #c(0.2085255208117019 -0.6411905440481438) #c(-0.6474034073977328 -0.15810480125045623) 
X     #c(0.14600280207229419 -0.421454194714017) #c(0.0744424854739143 -0.10374014288729083) 
X     #c(-0.13973974067502756 -0.5717728061166673) #c(0.7211687432221876 -0.1395838200540448) 
X
X     #c(0.40586 -1.24797)  #c(-0.2131400000000001 -1.13773) 
X     #c(0.28417 -0.82029)  #c(0.72524 -0.73988) 
X     #c(-0.27198 -1.11286) #c(1.48759 -1.290507) 
X
X     #c(1.0 0.0) #c(-1.249000902703301E-16 -2.220446049250313E-16) 
X     #c(-1.249000902703301E-16 2.220446049250313E-16)  #c(1.0 0.0)
|#
X
X
;-------------------------------------------------------------------------------
(defun Q-R! (A-matrix)
X  "same as modified-Gram-Schmidt!"
X  (modified-Gram-Schmidt! A-matrix))
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  upper-triangular-solve!
;;;                 lower-triangular-solve!
;;;  solve upper and lower triangular systems of equations.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun upper-triangular-solve!
X       (A
X        b
X        &key
X        (n (length b)))
X  "given
X   [1] A (required)
X       ==> upper n by n triangular matrix
X           (note: lower subdiagonal part of this matrix is ignored)
X   [2] b (required)
X       <=> on input, the right-hand side vector;
X           on output, the solution X to A X = b
X   [3] n (keyword; length of b)
X       ==> order of matrix A
returns
X   [4] X, the solution to A X = b, where X is the same as vector
which contained b (note that A is left unchanged)"
X  (let ((current-row (1- n))
X        (k-total 0))
X    (dotimes (i n b)
X      (setf (aref b current-row)
X            (/ (- (aref b current-row)
X                  (let ((sum 0.0)
X                        (k (1+ current-row)))
X                    ;(declare (dynamic-extent sum))
X                    (dotimes (j k-total sum)
X                      (incf sum (* (aref b k) (aref A current-row k)))
X                      (incf k))))
X               (aref A current-row current-row)))
X      (decf current-row)
X      (incf k-total))))
X
#|
(upper-triangular-solve!
X #2A((1 2 3) (0 4 5) (0 0 6))
X #(8 9 10))
;==> #(2.666666666666667 0.16666666666666652 1.6666666666666667)
X
(multiply-matrix-and-vector
X #2A((1 2 3) (0 4 5) (0 0 6))
X #(2.666666666666667 0.16666666666666652 1.6666666666666667))
;==> #(8.0 9.0 10.0)
|#
X
;-------------------------------------------------------------------------------
(defun lower-triangular-solve!
X       (A
X        b
X        &key
X        (n (length b)))
X  "given
X   [1] A (required)
X       ==> lower n by n triangular matrix
X           (note: upper superdiagonal part of this matrix is ignored)
X   [2] b (required)
X       <=> on input, the right-hand side vector;
X           on output, the solution X to A X = b
X   [3] n (keyword; length of b)
X       ==> order of matrix A
returns
X   [4] X, the solution to A X = b, where X is the same as vector
which contained b (note that A is left unchanged)"
X  (let ((current-row 0)
X        (k-total 0))
X    (dotimes (i n b)
X      (setf (aref b current-row)
X            (/ (- (aref b current-row)
X                  (let ((sum 0.0))
X                    ;(declare (dynamic-extent sum))
X                    (dotimes (k k-total sum)
X                      (incf sum (* (aref b k) (aref A current-row k))))))
X               (aref A current-row current-row)))
X      (incf current-row)
X      (incf k-total))))
X
#|
(lower-triangular-solve!
X (transpose #2A((1 2 3) (0 4 5) (0 0 6)))
X #(8 9 10))
;==> #(8.0 -1.75 -0.875)
X
(multiply-matrix-and-vector
X (transpose #2A((1 2 3) (0 4 5) (0 0 6)))
X #(8.0 -1.75 -0.875))
;==> #(8.0 9.0 10.0)
|#
SHAR_EOF
chmod 0644 matrix.lisp ||
echo 'restore of matrix.lisp failed'
Wc_c="`wc -c < 'matrix.lisp'`"
test 30945 -eq "$Wc_c" ||
	echo 'matrix.lisp: original size 30945, current size' "$Wc_c"
fi
# ============= multitaper.lisp ==============
if test -f 'multitaper.lisp' -a X"$1" != X"-c"; then
	echo 'x - skipping multitaper.lisp (File already exists)'
else
echo 'x - extracting multitaper.lisp (Text)'
sed 's/^X//' << 'SHAR_EOF' > 'multitaper.lisp' &&
;;;-*- Mode: LISP; Package: :SAPA; Syntax: COMMON-LISP -*-
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;
;  multitaper.lisp
;
;  a collection of Lisp functions for multitaper spectral estimation ...
;  Note:  before compiling and loading multitaper.lisp,
;         you should compile and load (in the order listed)
;            sapa-package.lisp, utilities.lisp, basic-math.lisp,
;            basic-statistics.lisp, dft-and-fft.lisp, tapers.lisp,
;            random.lisp, acvs.lisp and nonparametric.lisp
;
;  6/3/93
;
;  SAPA, Version 1.0; Copyright 1993, Donald B. Percival, All Rights Reserved
;
;  Use and copying of this software and preparation of derivative works
;  based upon this software are permitted.  Any distribution of this
;  software or derivative works must comply with all applicable United
;  States export control laws.
; 
;  This software is made available AS IS, and no warranty -- about the
;  software, its performance, or its conformity to any
;  specification -- is given or implied. 
;
;  Comments about this software can be addressed to dbp@apl.washington.edu
;-------------------------------------------------------------------------------
;;; (compile-file "ccl:SAPA;multitaper.lisp")
;;; (load "ccl:SAPA;multitaper.fasl")
;-------------------------------------------------------------------------------
(if (not (find-package :SAPA))
X  (defpackage :SAPA (:USE :COMMON-LISP)))
X
(in-package :SAPA)
X
(export '(;;; utility function to check orthonormality of tapers ...
X          check-orthonormality
X          
X          ;;; functions to compute a set of dpss data tapers ...
X          dpss-tapers-tri-diag
X          dpss-tapers-Thomson-approx
X          dpss-tapers-inverse-iteration
X          
X          ;;; function to compute an approximation to dpss data tapers ...
X          trig-prolate-tapers
X          
X          ;;; functions to compute multitaper spectral estimates ...
X          multitaper-spectral-estimate
X          eigenspectra->multitaper-spectral-estimate
X          eigenspectra->adaptive-multitaper-spectral-estimate
X
X          ;;; function to create confidence intervals for adaptive estimate ...
X          create-ci-for-amt-sdf-estimate
X          ))
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The function  check-orthonormality
;;;  checks the orthonormality of a list of vectors.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun check-orthonormality (list-of-vectors)
X  "given a list of vectors, computes and prints
their sum of squares and pairwise dot products, and 
returns the maximum absolute deviation from 0
of the pairwise dot products"
X  (let ((the-rest-of-them (cdr list-of-vectors))
X        (max-dev 0.0)
X        temp
X        kth-order-vector)
X    (dotimes (k (length list-of-vectors))
X      (format t "~&k = ~D: sum of squares = ~F"
X              k (sum-of-squares (nth k list-of-vectors))))
X    (dotimes (k (1- (length list-of-vectors)) max-dev)
X      (setf kth-order-vector (nth k list-of-vectors))
X      (format t "~&k = ~D:" k)
X      (dolist (another-vector the-rest-of-them)
X        (format t "~&      dot product = ~F" 
X                (setf temp (dot-product kth-order-vector another-vector))))
X      (if (> (abs temp) max-dev)
X        (setf max-dev (abs temp)))
X      (setf the-rest-of-them (cdr the-rest-of-them)))))
X
#|
single-float-epsilon  ;==> 1.1107651257113995E-16
(check-orthonormality `(#(1 0 0) #(0 1 0) #(,single-float-epsilon 0 1)))
;==>
k = 0: sum of squares = 1.0
k = 1: sum of squares = 1.0
k = 2: sum of squares = 1.0
k = 0:
X      dot product = 0.0
X      dot product = 1.1107651257113995E-16
k = 1:
X      dot product = 0.0
1.1107651257113995E-16
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  dpss-tapers-tri-diag
;;;                 dpss-tapers-Thomson-approx
;;;                 dpss-tapers-inverse-iteration
;;;  each compute a set of orthonormal dpss data tapers for a specified
;;;  sample size N and taper parameter NW (the product of the duration N
;;;  and half bandwidth W -- note that W is taken to be in standardized units
;;;  so that 0 < W < 1/2).  In general, dpss-tapers-tri-diag should be used
;;;  because it is the fastest method and is generally quite accurate.
;;;  The function dpss-tapers-Thomson-approx implements the scheme
;;;  in Thomson's 1982 paper (with some modifications), but,
;;;  although it is about as fast as dpss-tapers-tri-diag, it is less
;;;  accurate.  The function dpss-tapers-inverse-iteration is potentially
;;;  the most accurate of the three methods, but it is much slower.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun dpss-tapers-tri-diag
X       (N
X        number-of-tapers
X        &key
X        (taper-parameter 4.0)  ;NW
X        (print-progress-p nil)
X        (compute-true-eigenvalues-p nil))
X  "given
X   [1] N (required)
X       the sample size
X   [2] number-of-tapers (required)
X       number of orthonormal dpss data tapers
X       to be computed
X   [3] taper-parameter (keyword; 4.0)
X       NW, the duration-half-bandwidth product
X       (must be such that 0 < NW/N < 1/2)
X   [4] print-progress-p (keyword; nil)
X       if t, prints a dot after each eigenvalue
X       and eigenvector has been computed
X   [5] compute-true-eigenvalues-p (keyword; nil)
X       if t, returns eigenvalues for eigenproblem
X       of Equation (378) of the SAPA book;
X       if nil, returns eigenvalues for tridiagonal
X       formulation
returns
X   [1] a list of length number-of-tapers of N dimensional vectors
X       with orthonormal dpss's of orders 0, 1, ..., number-of-tapers - 1;
X   [2] a vector of length number-of-tapers with
X       eigenvalues as specified by compute-true-eigenvalues-p
---
Note: computes the dpss tapers using the tridiagonal
X      formulation (see Section 8.3 of the SAPA book)"
X  (assert (< 0 (/ taper-parameter N) 0.5))
X  (let* ((Nm1 (1- N))
X         (Nm1o2 (float (/ Nm1 2.0)))
X         (ctpiW (cos (* 2 pi (/ taper-parameter N))))
X         (diag (make-array N))
X         (off-diag (make-array Nm1))
X         (results '()))
X    ;;; generate diagonal elements of symmetric tridiagonal system ...
X    (dotimes (i N)
X      (setf (aref diag i) (* ctpiW (expt (- Nm1o2 i) 2))))
X    ;;; generate off-diagonal elements ...
X    (dotimes (i Nm1)
X      (setf (aref off-diag i) (* 0.5 (1+ i) (- N (1+ i)))))
X    (if print-progress-p
X      (format t "~&finding eigenvalues  "))
X    ;;; get eigenvalues ...
X    (let ((eigenvalues (largest-eigenvalues-of-tridiagonal-matrix
X                        diag off-diag number-of-tapers
X                        :print-progress-p print-progress-p)))
X      (if print-progress-p
X        (format t "~&finding eigenvectors "))
X      ;;; get eigenvectors (the dpss's) ...
X      (dotimes (k number-of-tapers (values (reverse results) eigenvalues))
X        (push (fast-tridiag-eigenvalue->dpss!
X               (aref eigenvalues k) k diag off-diag)
X              results)
X        (if compute-true-eigenvalues-p
X          (setf (aref eigenvalues k)
X                (dpss->eigenvalue (car results) taper-parameter)))
X        (if print-progress-p (format t "."))))))
X
#|
(multiple-value-bind (list-of-tapers eigenvalues)
X                     (dpss-tapers-tri-diag
X                      16 4
X                      :taper-parameter 4
X                      :compute-true-eigenvalues-p t)
X  (dotimes (i 4)
X    (format t "~&~2D: ~16,13F" i (svref eigenvalues i)))
X  (format t "~&---")
X  (dotimes (i 16)
X    (format t "~&~2D: ~11,8F ~11,8F ~11,8F ~11,8F"
X            i
X            (svref (elt list-of-tapers 0) i)
X            (svref (elt list-of-tapers 1) i)
X            (svref (elt list-of-tapers 2) i)
X            (svref (elt list-of-tapers 3) i)))
X  (check-orthonormality list-of-tapers))
;==>
X 0:  0.9999999999819
X 1:  0.9999999970607
X 2:  0.9999997876512
X 3:  0.9999910301796
---
X 0:  0.00083944  0.00442231  0.01607017  0.04634116
X 1:  0.00655153  0.02829405  0.08186426  0.18053656
X 2:  0.02694253  0.09460907  0.21480029  0.35196200
X 3:  0.07617175  0.21249814  0.36208333  0.39776040
X 4:  0.16388770  0.34799492  0.40188579  0.19828108
X 5:  0.28236294  0.42176848  0.24234799 -0.13646373
X 6:  0.40070313  0.35562544 -0.05975464 -0.31466700
X 7:  0.47568617  0.14005351 -0.30328392 -0.16191173
X 8:  0.47568617 -0.14005351 -0.30328392  0.16191173
X 9:  0.40070313 -0.35562544 -0.05975464  0.31466700
10:  0.28236294 -0.42176848  0.24234799  0.13646373
11:  0.16388770 -0.34799492  0.40188579 -0.19828108
12:  0.07617175 -0.21249814  0.36208333 -0.39776040
13:  0.02694253 -0.09460907  0.21480029 -0.35196200
14:  0.00655153 -0.02829405  0.08186426 -0.18053656
15:  0.00083944 -0.00442231  0.01607017 -0.04634116
k = 0: sum of squares = 1.0
k = 1: sum of squares = 1.0
k = 2: sum of squares = 1.0
k = 3: sum of squares = 1.0
k = 0:
X      dot product = 9.17082571992231E-18
X      dot product = -2.802662615875029E-16
X      dot product = 5.319366908757006E-18
k = 1:
X      dot product = -2.847385955490056E-17
X      dot product = 1.2495430037895439E-17
k = 2:
X      dot product = -3.859759734048396E-17
3.859759734048396E-17
|#
X
;-------------------------------------------------------------------------------
(defun dpss-tapers-Thomson-approx
X       (N
X        number-of-tapers
X        &key
X        (taper-parameter 4.0)  ;NW
X        (print-progress-p nil)
X        (compute-true-eigenvalues-p nil)
X        (abscissas (let ()
X                     (declare (special *abscissas-32-point*))
X                     *abscissas-32-point*))
X        (weights (let ()
X                   (declare (special *weights-32-point*))
X                   *weights-32-point*)))
X  "given
X   [1] N (required)
X       the sample size
X   [2] number-of-tapers (required)
X       number of orthonormal dpss data tapers
X       to be computed
X   [3] taper-parameter (keyword; 4.0)
X       NW, the duration-half-bandwidth product
X       (must be such that 0 < NW/N < 1/2)
X   [4] print-progress-p (keyword; nil)
X       if t, prints a dot after each eigenvalue
X       and eigenvector has been computed
X   [5] compute-true-eigenvalues-p (keyword; nil)
X       if t, returns eigenvalues for eigenproblem
X       of Equation (378) of the SAPA book;
X       if nil, returns nil
X   [6] abscissas (keyword; *abscissas-32-point*)
X       a vector of abscissas points used
X       in Gauss-Legendre quadrature
X   [7] weights (keyword; *abscissas-32-point*)
X       a vector of weights used
X       in Gauss-Legendre quadrature
returns
X   [1] a list of length number-of-tapers of N dimensional vectors
X       with orthonormal dpss's of orders 0, 1, ..., number-of-tapers - 1;
X   [2] a vector of length number-of-tapers with
X       eigenvalues if compute-true-eigenvalues-p is t;
X       nil if if compute-true-eigenvalues-p is nil
---
Note: computes the dpss tapers using Thomson's numerical
X      integration scheme (see Section 8.2 of the SAPA book)"
X  (declare (special *abscissas-32-point* *weights-32-point*))
X  (assert (< 0 (/ taper-parameter N) 0.5))
X  (let* ((c (* pi taper-parameter))
X         (Cap-J (length abscissas))
X         (sqrt-weights (transform-a-sequence #'sqrt weights))
X         (Cap-Psi (make-array `(,Cap-J ,Cap-J))))
X    ;;; Compute elements of Cap-Psi matrix
X    (dotimes (i Cap-J)
X      (dotimes (j (1+ i))
X        (setf (aref Cap-Psi j i)
X              (setf (aref Cap-Psi i j)
X                    (if (= i j)
X                      (/ (* c (aref weights i)) pi)
X                      (/ (* (aref sqrt-weights i)
X                            (aref sqrt-weights j)
X                            (sin (* c (- (aref abscissas i)
X                                         (aref abscissas j)))))
X                         (* pi (- (aref abscissas i)
X                                  (aref abscissas j)))))))))
X    (multiple-value-bind (diagonal-elements off-diagonal-elements)
X                         (Householder-reduction-of-real-sym-matrix! Cap-Psi)
X      (multiple-value-bind (eigenvalues array-with-eigenvectors)
X                           (eigenvalues-and-vectors-of-sym-tridiag-matrix
X                            diagonal-elements off-diagonal-elements
X                            :z Cap-Psi)
X        (let ((sorted-eigenvalues (sort (copy-seq eigenvalues) #'>=))
X              (true-eigenvalues (if compute-true-eigenvalues-p
X                                  (make-array number-of-tapers)))
X              (list-of-dpss '())
X              (N-minus-1 (1- N))
X              (N-over-2 (/ N 2))
X              index-unsorted a-dpss x i-down factor)
X          (if print-progress-p
X            (format t "~&finding eigenvectors "))
X          (dotimes (k number-of-tapers (values (reverse list-of-dpss)
X                                               true-eigenvalues))
X            (setf index-unsorted (position (aref sorted-eigenvalues k)
X                                           eigenvalues))
X            (setf a-dpss (make-array N :initial-element 0.0))
X            (push a-dpss list-of-dpss)
X            (dotimes (i N)
X              (setf x (float (1- (/ (1+ (* 2 i)) N))))
X              (dotimes (j Cap-J)
X                (incf (aref a-dpss i)
X                      (* (aref sqrt-weights j)
X                         (if (= x (aref abscissas j))
X                           (/ c pi)
X                           (/ (sin (* c (- x (aref abscissas j))))
X                              (* pi (- x (aref abscissas j)))))
X                         (aref array-with-eigenvectors j index-unsorted)))))
X            ;;; force symmetry ...
X            (setf i-down N-minus-1)
X            (if (evenp k)
X              (dotimes (i-up N-over-2)
X                (setf (aref a-dpss i-up)
X                      (setf (aref a-dpss i-down)
X                            (* 0.5 (+ (aref a-dpss i-up)
X                                      (aref a-dpss i-down)))))
X                (decf i-down))
X              (dotimes (i-up N-over-2)
X                (setf (aref a-dpss i-up)
X                      (- (setf (aref a-dpss i-down)
X                               (* 0.5 (- (aref a-dpss i-up)
X                                         (aref a-dpss i-down))))))
X                (decf i-down)))
X            ;;; normalize properly ...
X            (setf factor (/ (sqrt (sum-of-squares a-dpss))))
X            (if (or (and (evenp k) (minusp (sum a-dpss)))
X                    (and (oddp k) (minusp (let ((the-sum 0.0))
X                                            (dotimes (j N the-sum)
X                                              (incf the-sum
X                                                    (* (- N-minus-1 (* 2 j))
X                                                       (aref a-dpss j))))))))
X              (setf factor (- factor)))
X            (a*x! factor a-dpss)
X            (if compute-true-eigenvalues-p
X              (setf (svref true-eigenvalues k)
X                    (dpss->eigenvalue a-dpss taper-parameter)))
X            (if print-progress-p (format t "."))))))))
X
#|
(multiple-value-bind (list-of-tapers eigenvalues)
X                     (dpss-tapers-Thomson-approx
X                      16 4
X                      :taper-parameter 4
X                      :compute-true-eigenvalues-p t)
X  (dotimes (i 4)
X    (format t "~&~2D: ~16,13F" i (svref eigenvalues i)))
X  (format t "~&---")
X  (dotimes (i 16)
X    (format t "~&~2D: ~11,8F ~11,8F ~11,8F ~11,8F"
X            i
X            (svref (elt list-of-tapers 0) i)
X            (svref (elt list-of-tapers 1) i)
X            (svref (elt list-of-tapers 2) i)
X            (svref (elt list-of-tapers 3) i)))
X  (check-orthonormality list-of-tapers))
;==>
X 0:  0.9999999996987
X 1:  0.9999999691381
X 2:  0.9999988186339
X 3:  0.9999757483772
---
X 0:  0.00029156  0.00199509  0.00917919  0.03251189
X 1:  0.00391148  0.01960816  0.06457599  0.15858593
X 2:  0.02019310  0.07831855  0.19367378  0.34072886
X 3:  0.06496358  0.19446610  0.35286192  0.41188210
X 4:  0.15154808  0.33902060  0.41377799  0.22832948
X 5:  0.27490901  0.42755451  0.26579779 -0.11427479
X 6:  0.40237312  0.36919277 -0.04345084 -0.31005778
X 7:  0.48467585  0.14703590 -0.29995992 -0.16346730
X 8:  0.48467585 -0.14703590 -0.29995992  0.16346730
X 9:  0.40237312 -0.36919277 -0.04345084  0.31005778
10:  0.27490901 -0.42755451  0.26579779  0.11427479
11:  0.15154808 -0.33902060  0.41377799 -0.22832948
12:  0.06496358 -0.19446610  0.35286192 -0.41188210
13:  0.02019310 -0.07831855  0.19367378 -0.34072886
14:  0.00391148 -0.01960816  0.06457599 -0.15858593
15:  0.00029156 -0.00199509  0.00917919 -0.03251189
k = 0: sum of squares = 0.9999999999999999
k = 1: sum of squares = 1.0
k = 2: sum of squares = 1.0000000000000002
k = 3: sum of squares = 0.9999999999999997
k = 0:
X      dot product = 1.889836384442751E-18
X      dot product = -3.781650095408522E-9
X      dot product = -4.313091767418897E-18
k = 1:
X      dot product = -3.18890963982299E-17
X      dot product = -3.870056405406281E-8
k = 2:
X      dot product = -2.9761349634727097E-17
3.870056405406281E-8
|#
X
;-------------------------------------------------------------------------------
(defun dpss-tapers-inverse-iteration
X       (N 
X        number-of-tapers
X        &key
X        (taper-parameter 4.0)
X        (print-progress-p nil)
X        (eps 0.5e-4))
X  "given
X   [1] N (required)
X       the sample size
X   [2] number-of-tapers (required)
X       number of orthonormal dpss data tapers
X       to be computed
X   [3] taper-parameter (keyword; 4.0)
X       NW, the duration-half-bandwidth product
X       (must be such that 0 < NW/N < 1/2)
X   [4] print-progress-p (keyword; nil)
X       if t, prints a dot after each eigenvalue
X       and eigenvector has been computed
X   [5] eps (keyword;  0.5e-4)
X       controls accuracy
returns
X   [1] a list of length number-of-tapers of N dimensional vectors
X       with orthonormal dpss's of orders 0, 1, ..., number-of-tapers - 1;
X   [2] a vector of length number-of-tapers with the
X       corresponding eigenvalues
---
Note: computes the dpss tapers using inverse
X      iteration (see Section 8.1 of the SAPA book)"
X  (assert (< 0 (/ taper-parameter N) 0.5))
X  (let* ((list-of-vectors nil)
X         (eigenvalues (make-array number-of-tapers))
X         (sines (make-array (1- (* 2 N))))
X         (W (/ taper-parameter N))
X         (two-pi-W (* 2.0 pi W))
X         (root-N (sqrt N))
X         (displacement N)
X         (v-old (make-array N))
X         1-over-norm sum diff)
X    #+mcl(declare (dynamic-extent v-old))
X    ;;; set up the Toeplitz vector of sines
X    (dotimes (m (1- N))
X      (decf displacement)
X      (setf (aref sines m)
X            (setf (aref sines (- (* 2 n) m 2))
X                  (/ (sin (* two-pi-W displacement))
X                     (* pi displacement)))))
X    ;;; loop once for each sequence
X    (if print-progress-p
X      (format t "~&finding eigenvalues and eigenvectors "))
X    (dotimes (k number-of-tapers)
X      (setf (aref sines (1- N))
X            (if (zerop k)
X              (1- (* 2 W))
X              (- (* 2 W) (1+ (aref eigenvalues (1- k))))))
X      (let ((v (generate-initial-guess-at-dpss N k)))
X        ;;; The `dotimes' form returns nil if the convergence
X        ;;; criterion was satisfied --- a value of t indicates
X        ;;; non-convergence.
X        (if (dotimes (it (ceiling (* root-N (+ k 3))) t)
X              (replace v-old v)
X              (toeplitz sines v-old v)
X              (dolist (v-lower-order list-of-vectors)
X                (a*x+y! (- (dot-product v-lower-order v)) v-lower-order v))
X              (setf 1-over-norm (/ (euclidean-norm v)))
X              (a*x! 1-over-norm v)
X              (a*x+y! -1.0 v v-old)
X              (setf diff (euclidean-norm v-old))
X              (a*x+y! 2.0 v v-old)
X              (setf sum (euclidean-norm v-old))
X              (setf (aref eigenvalues k)
X                    (+ (if (plusp k) (aref eigenvalues (1- k)) 0.0)
X                       (if (< sum diff)
X                         (- 1-over-norm) 1-over-norm)))
X              (if (<= (min diff sum) eps) (return nil)))
X          (format t "~&dpss: order ~D did not converge, err = ~F"
X                  k (min diff sum)))
X        ;;; normalize dpss properly ...
X        (if (minusp (sum v :end (max 2 (truncate N (1+ k)))))
X          (a*x! -1 v))
X        (push v list-of-vectors)
X        (if print-progress-p (format t "."))))
X    ;;; convert from sigma to actual eigenvalue of interest
X    (x+b! eigenvalues 1)
X    (values (reverse list-of-vectors) eigenvalues)))
X
#|
(multiple-value-bind (list-of-tapers eigenvalues)
X                     (dpss-tapers-inverse-iteration
X                      16 4
X                      :taper-parameter 4)
X  (dotimes (i 4)
X    (format t "~&~2D: ~16,13F" i (svref eigenvalues i)))
X  (format t "~&---")
X  (dotimes (i 16)
X    (format t "~&~2D: ~11,8F ~11,8F ~11,8F ~11,8F"
X            i
X            (svref (elt list-of-tapers 0) i)
X            (svref (elt list-of-tapers 1) i)
X            (svref (elt list-of-tapers 2) i)
X            (svref (elt list-of-tapers 3) i)))
X  (check-orthonormality list-of-tapers))
;==>
X 0:  0.9999999999819
X 1:  0.9999999970607
X 2:  0.9999997876512
X 3:  0.9999910301796
---
X 0:  0.00083944  0.00442231  0.01607017  0.04634116
X 1:  0.00655153  0.02829405  0.08186426  0.18053656
X 2:  0.02694253  0.09460907  0.21480029  0.35196201
X 3:  0.07617175  0.21249814  0.36208333  0.39776040
X 4:  0.16388770  0.34799492  0.40188579  0.19828107
X 5:  0.28236294  0.42176848  0.24234798 -0.13646374
X 6:  0.40070313  0.35562544 -0.05975464 -0.31466699
X 7:  0.47568617  0.14005351 -0.30328392 -0.16191172
X 8:  0.47568617 -0.14005351 -0.30328392  0.16191173
X 9:  0.40070313 -0.35562544 -0.05975464  0.31466699
10:  0.28236294 -0.42176848  0.24234798  0.13646373
11:  0.16388770 -0.34799492  0.40188579 -0.19828108
12:  0.07617175 -0.21249814  0.36208333 -0.39776040
13:  0.02694253 -0.09460907  0.21480029 -0.35196200
14:  0.00655153 -0.02829405  0.08186426 -0.18053656
15:  0.00083944 -0.00442231  0.01607017 -0.04634116
k = 0: sum of squares = 0.9999999999999999
k = 1: sum of squares = 1.0000000000000002
k = 2: sum of squares = 1.0000000000000002
k = 3: sum of squares = 0.9999999999999999
k = 0:
X      dot product = 3.260695682102792E-17
X      dot product = 2.659683454378503E-19
X      dot product = 1.9203930980149497E-17
k = 1:
X      dot product = 2.2386064356394453E-16
X      dot product = -6.830473686658678E-17
k = 2:
X      dot product = -1.283261691353843E-15
1.283261691353843E-15
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The function  trig-prolate-tapers
;;;  computes the trig-prolate data tapers, an orthonormal set of tapers
;;;  that are a useful approximation to the dpss data tapers.  Their chief
;;;  advantage is the speed with which they can be computed.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun trig-prolate-tapers
X       (N
X        number-of-tapers
X        &key
X        (taper-parameter 4)  ;NW
X        (print-progress-p nil)
X        (compute-true-eigenvalues-p nil))
X  "given
X   [1] N (required)
X       the sample size
X   [2] number-of-tapers (required)
X       number of orthonormal trig prolate tapers
X       to be computed; currently restricted to
X       one of the following maximum values:
X       2 if taper-parameter is 2;
X       4 if taper-parameter is 3;
X       5 if taper-parameter is 4;
X       7 if taper-parameter is 5
X   [3] taper-parameter (keyword; 4.0)
X       NW, the duration-half-bandwidth product;
X       currently restricted to one of the
X       following integers: 2, 3, 4, 5
X       (must also be such that 0 < NW/N < 1/2)
X   [4] print-progress-p (keyword; nil)
X       if t, prints a dot after each taper
X       has been computed
X   [5] compute-true-eigenvalues-p (keyword; nil)
X       if t, returns eigenvalues for eigenproblem
X       of Equation (378) of the SAPA book;
X       if nil, returns nil
returns
X   [1] a list of length number-of-tapers of N dimensional vectors
X       with orthonormal trig prolate data tapers of orders
X       0, 1, ..., number-of-tapers - 1;
X   [2] a vector of length number-of-tapers with
X       eigenvalues if compute-true-eigenvalues-p is t;
X       nil if if compute-true-eigenvalues-p is nil
---
Note: computes the trig prolate approximation to
X      the dpss tapers (see Section 8.4 of the SAPA book)"
X  (assert (and (integerp taper-parameter) (<= 2 taper-parameter 5)))
X  (case taper-parameter
X    (2 (assert (<= number-of-tapers 2)))
X    (3 (assert (<= number-of-tapers 4)))
X    (4 (assert (<= number-of-tapers 5)))
X    (5 (assert (<= number-of-tapers 7))))
X  (if print-progress-p
X    (format t "~&computing trig prolates "))
X  (let ((results '())
X        (true-eigenvalues (if compute-true-eigenvalues-p
X                            (make-array number-of-tapers))))
X    (dotimes (k number-of-tapers (values (reverse results)
X                                         true-eigenvalues))
X      (push (generate-tri-prolate-taper
X             N
X             (produce-upsilon-func
X              (trig-prolate-sign-cosine-coefficients
X               taper-parameter k)))
X            results)
X      (if compute-true-eigenvalues-p
X        (setf (svref true-eigenvalues k)
X              (dpss->eigenvalue (car results) taper-parameter)))
X      (if print-progress-p (format t ".")))))
X
#|
(multiple-value-bind (list-of-tapers eigenvalues)
X                     (trig-prolate-tapers
X                      16 4
X                      :taper-parameter 4
X                      :compute-true-eigenvalues-p t)
X  (dotimes (i 4)
X    (format t "~&~2D: ~16,13F" i (svref eigenvalues i)))
X  (format t "~&---")
X  (dotimes (i 16)
X    (format t "~&~2D: ~11,8F ~11,8F ~11,8F ~11,8F"
X            i
X            (svref (elt list-of-tapers 0) i)
X            (svref (elt list-of-tapers 1) i)
X            (svref (elt list-of-tapers 2) i)
X            (svref (elt list-of-tapers 3) i)))
X  (check-orthonormality list-of-tapers))
;==>
X 0:  0.9999999577008
X 1:  1.0000000927627
X 2:  0.9999979984130
X 3:  0.9999612860661
---
X 0:  0.00025830  0.00192627  0.00826665  0.03228915
X 1:  0.00369660  0.01882494  0.06316808  0.15240000
X 2:  0.01950040  0.07672601  0.19103732  0.33537529
X 3:  0.06367965  0.19247292  0.35125873  0.41454498
X 4:  0.15001960  0.33788630  0.41520495  0.23422290
X 5:  0.27391468  0.42811787  0.26873512 -0.11173293
X 6:  0.40252728  0.37068759 -0.04152554 -0.31135663
X 7:  0.48578549  0.14781044 -0.29954212 -0.16468134
X 8:  0.48578549 -0.14781044 -0.29954212  0.16468134
X 9:  0.40252728 -0.37068759 -0.04152554  0.31135663
10:  0.27391468 -0.42811787  0.26873512  0.11173293
11:  0.15001960 -0.33788630  0.41520495 -0.23422290
12:  0.06367965 -0.19247292  0.35125873 -0.41454498
13:  0.01950040 -0.07672601  0.19103732 -0.33537529
14:  0.00369660 -0.01882494  0.06316808 -0.15240000
15:  0.00025830 -0.00192627  0.00826665 -0.03228915
k = 0: sum of squares = 0.99999995816752
k = 1: sum of squares = 1.00000013603332
k = 2: sum of squares = 1.0000001262701201
k = 3: sum of squares = 0.99999995576222
k = 0:
X      dot product = -3.1223116743451206E-17
X      dot product = -3.757875999225208E-8
X      dot product = -2.6359665318553827E-18
k = 1:
X      dot product = -1.1868625656927256E-17
X      dot product = -2.1731819964144203E-8
k = 2:
X      dot product = -3.6483403104137224E-17
2.1731819964144203E-8
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  multitaper-spectral-estimate
;;;                 eigenspectra->multitaper-spectral-estimate
;;;                 eigenspectra->adaptive-multitaper-spectral-estimate
;;;  implement multitaper spectral estimation.  The simpliest formulation
;;;  is implemented by multitaper-spectral-estimate, which computes
;;;  the simple average of direct spectral estimates given by Equation (333)
;;;  of the SAPA book.  This function also optionally returns a list
;;;  of the individual eigenspectra, all or part of which can then be presented
;;;  to either eigenspectra->multitaper-spectral-estimate
;;;  or        eigenspectra->adaptive-multitaper-spectral-estimate
;;;  to compute, respectively, a simple average with a few number of terms
;;;  or the adaptive weighted multitaper spectral estimate given by
;;;  Equation (370a) of the SAPA book.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun multitaper-spectral-estimate
X       (time-series
X        list-of-data-tapers
X        &key
X        (center-data t)
X        (start 0)
X        (end (length time-series))
X        (N-tapers (length list-of-data-tapers))
X        (N-nonzero-freqs :half-next-power-of-2)
X        (return-est-for-0-freq-p nil)
X        (sampling-time 1.0)
X        (scratch-dft (make-array (get-N-dft N-nonzero-freqs (- end start))))
X        (recenter-after-tapering-p t)
X        (restore-power-option-p t)
X        (sdf-transformation #'convert-to-dB)
X        (result-sdf (make-array (get-N-freqs N-nonzero-freqs (- end start)
X                                             return-est-for-0-freq-p)))
X        (return-eigenpspectra-p t)
X        (return-frequencies-p t)
X        (freq-transformation nil)
X        (result-freq (if return-frequencies-p
X                       (make-array (get-N-freqs N-nonzero-freqs (- end start)
X                                                return-est-for-0-freq-p)))))
X  "given
X   [1] time-series (required)
X       ==> a sequence of time series values
X   [2] list-of-data-tapers (required)
X       ==> a list of orthonormal data tapers,
X           each of length (- end start)
X   [3] center-data (keyword; t)
X       ==> if t, subtract sample mean from time series;
X           if a number, subtracts that number from time series;
X           if nil, time-series is not centered
X   [4] start (keyword; 0)
X       ==> start index of sequence to be used
X   [5] end (keyword; length of time-series)
X       ==> 1 + end index of sequence to be used
X   [6] N-tapers (keyword; length of list-of-data-tapers)
X       ==> number of data tapers to be used in list-of-data-tapers 
X   [7] N-nonzero-freqs (keyword; :half-next-power-of-2)
X       ==> specifies at how many nonzero frequencies
X           multitaper spectral estimate is to be computed -- choices are:
X           :half-next-power-of-2
X            ==> 1/2 * power of two >= sample size;
X           :next-power-of-2
X            ==> power of two >= sample size;
X           :twice-next-power-of-2
X            ==> 2 * power of two >= sample size;
X           :Fourier
X            ==> just at Fourier frequencies
X            -- or --
X           any power of 2 >= 1/2 * [power of two >= sample size]
X   [8] return-est-for-0-freq-p (keyword; nil)
X       ==> if t, sdf is computed at zero frequency;
X           otherwise, it is not computed.
X   [9] sampling-time (keyword; 1.0)
X       ==> sampling time (called delta t in the SAPA book)
X  [10] scratch-dft (keyword; vector of correct length)
X       ==> vector in which the in-place dft is done
X  [11] recenter-after-tapering-p (keyword; t)
X       ==> if t and data-taper is a function,
X           centers tapered series by subtracting
X           off its sample mean
X  [12] restore-power-option-p (keyword; t)
X       ==> if t and data-taper is a function,
X           normalizes tapered series to have same
X           sum of squares as before tapering
X  [13] sdf-transformation (keyword; #'convert-to-dB)
X       ==> a function of one argument or nil;
X           if bound to a function, the function is used
X           to transform all elements of result-sdf
X  [14] result-sdf (keyword; vector of correct length)
X       <== vector into which multitaper spectral estimate is placed;
X           it must be exactly of the length dictated
X           by N-nonzero-freqs and return-est-for-0-freq-p
X  [15] return-eigenpspectra-p (keyword; t)
X       ==> if t, individual eigenspectra are returned in a list
X           (each eigenspectrum in the list is associated with
X           the corresponding taper in list-of-data-tapers)
X  [16] return-frequencies-p (keyword; t)
X       ==> if t, the frequencies associated with the spectral estimate
X           are computed and returned in result-freq
X  [17] freq-transformation (keyword; nil)
X       ==> a function of one argument or nil;
X           if bound to a function, the function is used
X           to transform all elements of result-freq
X           (ignored unless return-frequencies-p is true)
X  [18] result-freq (keyword; nil or vector of correct length)
X       <== not used if return-frequencies-p nil; otherwise,
X           vector of length N-nonzero-freqs (if return-est-for-0-freq-p is nil)
X           or N-nonzero-freqs +1 (if return-est-for-0-freq-p is t)
X           into which the frequencies associated with the values
X           in result-sdf are placed
returns
X   [1] result-sdf, a vector holding
X       the properly transformed multitaper spectral estimate
X   [2] result-freq (if return-frequencies-p is t),
X       where result-freq is a vector holding
X       the properly transformed frequencies
X       associated with values in  result-sdf,
X       or
X       nil (if return-frequencies-p is nil)
X   [3] the length of the vector result-sdf
X   [4] a list of untransformed eigenspectra (return-eigenpspectra-p is t),
X       or
X       nil (if return-eigenpspectra-p is nil)
---
Note: see Section 7.1 of the SAPA book"
X  (let* ((sample-size (- end start))
X         (N-freqs (get-N-freqs N-nonzero-freqs sample-size
X                               return-est-for-0-freq-p))
X         (N-dft (get-N-dft N-nonzero-freqs sample-size))
X         (offset (if return-est-for-0-freq-p 0 1))
X         (fiddle-factor-sdf (/ sampling-time sample-size))
X         (fiddle-factor-freq (/ (* N-dft sampling-time)))
X         (list-of-eigenspectra '()))
X    (fill result-sdf 0.0)
X    ;;; loop once for each data taper ...
X    (dotimes (k N-tapers)
X      (center&taper-time-series time-series
X                                :center-data center-data
X                                :start start
X                                :end end
X                                :data-taper #'supplied-data-taper!
X                                :data-taper-parameters
X                                (elt list-of-data-tapers k)
X                                :recenter-after-tapering-p
X                                recenter-after-tapering-p
X                                :restore-power-option-p
X                                restore-power-option-p
X                                :result scratch-dft)
X        ;;; zero the rest of scratch-dft ...
X        (fill scratch-dft 0.0 :start sample-size :end N-dft)
X        (dft! scratch-dft :N N-dft)
X        (if return-eigenpspectra-p
X          (let ((an-eigenpspectrum (make-array N-freqs)))
X            (dotimes (i N-freqs)
X              (setf (aref an-eigenpspectrum i)
X                    (* fiddle-factor-sdf
X                       (expt (abs (aref scratch-dft (+ i offset))) 2))))
X            (push an-eigenpspectrum list-of-eigenspectra)
X            (x+y! result-sdf an-eigenpspectrum))
X          (dotimes (i N-freqs)
X            (incf (aref result-sdf i)
X                  (* fiddle-factor-sdf
X                     (expt (abs (aref scratch-dft (+ i offset))) 2))))))
X    (a*x! (/ (float N-tapers)) result-sdf)
X    (if sdf-transformation
X      (transform-a-sequence! sdf-transformation result-sdf))
X    (when return-frequencies-p
X      (dotimes (i N-freqs)
X        (setf (aref result-freq i)
X              (* fiddle-factor-freq (float (+ i offset)))))
X      (if freq-transformation
X        (transform-a-sequence! freq-transformation result-freq)))
X    (values result-sdf result-freq N-freqs (reverse list-of-eigenspectra))))
X
#|
;;; Here we test the multitaper spectral estimate using the Parseval result
;;; analogous to the one we used to check the periodogram ... 
(dolist (N '(22 32 63 64 65) (values))
X  (let* ((time-series (x+b! (ranorms N) 100))
X         (sigma^2 (sample-variance time-series))
X         (sampling-time 0.25)
X         (the-sdf-est (multitaper-spectral-estimate
X                       time-series
X                       (dpss-tapers-tri-diag N 4)
X                       :N-nonzero-freqs :Fourier
X                       :sdf-transformation nil
X                       :sampling-time sampling-time))
X         (N-f (length the-sdf-est))
X         (Parseval-sum (if (evenp N)
X                         (/ (+ (* 2 (sum the-sdf-est  :end (1- N-f)))
X                               (svref the-sdf-est  (1- N-f)))
X                            (* N sampling-time))
X                         (/ (* 2 (sum the-sdf-est ))
X                            (* N sampling-time)))))
X    (format t "~& N = ~2D, N-f = ~2D: ~11,8F ~11,8F ~F"
X            N N-f sigma^2 Parseval-sum (/ Parseval-sum sigma^2))))
;==>
X N = 22, N-f = 11:  0.55906682  0.55906682 1.0000000000000064
X N = 32, N-f = 16:  0.63422723  0.63422723 0.9999999999999991
X N = 63, N-f = 31:  1.20061865  1.20061865 1.0000000000000042
X N = 64, N-f = 32:  0.77272584  0.77272584 1.0000000000000024
X N = 65, N-f = 32:  1.01794147  1.01794147 1.000000000000031
;;; Note: the important thing is that the last column of numbers
;;;       be close to unity -- since we are using random numbers,
;;;       these three columns of numbers will change each time
;;;       the above Lisp form is evaluated.
|#
X
;-------------------------------------------------------------------------------
(defun eigenspectra->multitaper-spectral-estimate
X       (list-of-eigenspectra
X        &key
X        (N-eigenspectra (length list-of-eigenspectra))
X        (sdf-transformation #'convert-to-dB)
X        (result-sdf (make-array (length (car list-of-eigenspectra)))))
X  "given
X   [1] list-of-eigenspectra (required)
X       ==> list of eigenspectra (such as optionally returned
X           by multitaper-spectral-estimate); each eigenspectrum
X           is assumed to be untransformed (e.g., not expressed
X           in dB) and represented by a vector
X   [2] N-eigenspectra (keyword; length of list-of-eigenspectra)
X       ==> number of eigenspectra to be used (must be less than
X           or equal to the length of list-of-eigenspectra)
X   [3] sdf-transformation (keyword; #'convert-to-dB)
X       ==> a function of one argument or nil;
X           if bound to a function, the function is used
X           to transform all elements of result-sdf
X   [4] result-sdf (keyword; vector of correct length)
X       <== vector into which multitaper spectral estimate is placed;
X           it must be exactly the same length as each of the vectors
X           in list-of-eigenspectra
returns
X   [1] result-sdf, a vector holding
X       the properly transformed multitaper spectral estimate
---
Note: see Section 7.1 of the SAPA book"
X  (replace result-sdf (car list-of-eigenspectra))
X  (let ((the-rest (cdr list-of-eigenspectra)))
X    (dotimes (k (1- N-eigenspectra))
X      (x+y! result-sdf (elt the-rest k))))
X  (a*x! (/ (float N-eigenspectra)) result-sdf)
X  (if sdf-transformation
X    (transform-a-sequence! sdf-transformation result-sdf))
X  (values result-sdf))
X
#|
;;; Here we create two multitaper spectral estimates, one using 4 tapers
;;; and the other using 3 tapers.  We check the Parseval result in
;;; both cases.
(let* ((time-series (x+b! (ranorms 32) 100))
X       (sigma^2 (sample-variance time-series))
X       (sampling-time 0.25))
X  (multiple-value-bind (mt-4 junk N-f list-of-eigenspectra)
X                       (multitaper-spectral-estimate
X                        time-series
X                        (dpss-tapers-tri-diag 32 4)
X                        :N-nonzero-freqs :Fourier
X                        :sdf-transformation nil
X                        :sampling-time sampling-time)
X    (declare (ignore junk))
X    (print (/ (+ (* 2 (sum mt-4 :end (1- N-f)))
X                 (svref mt-4 (1- N-f)))
X              (* 32 sampling-time sigma^2)))
X    (let ((mt-3 (eigenspectra->multitaper-spectral-estimate
X                 list-of-eigenspectra
X                 :N-eigenspectra 3
X                 :sdf-transformation nil)))
X      (print (/ (+ (* 2 (sum mt-3 :end (1- N-f)))
X                   (svref mt-3 (1- N-f)))
X                (* 32 sampling-time sigma^2)))
X      (values))))
;==>
0.9999999999999991 
0.9999999999999993
|#
X
;-------------------------------------------------------------------------------
(defun eigenspectra->adaptive-multitaper-spectral-estimate
X       (list-of-eigenspectra
X        eigenvalues
X        variance-of-time-series
X        &key
X        (sampling-time 1.0)
X        (N-eigenspectra (length list-of-eigenspectra))
X        (maximum-number-of-iterations 100)
X        (result-dof (make-array (length (car list-of-eigenspectra))))
X        (sdf-transformation #'convert-to-dB)
X        (result-sdf (make-array (length (car list-of-eigenspectra)))))
X  "given
X   [1] list-of-eigenspectra (required)
X       ==> list of eigenspectra (such as optionally returned
X           by multitaper-spectral-estimate); each eigenspectrum
X           is assumed to be untransformed (e.g., not expressed
X           in dB)
X   [2] eigenvalues (required)
X       ==> vector of eigenvalues corresponding to the dpss's
X           used to create the eigenspectra (the length of eigenvalues 
X           should be at least as large as the length specified
X           by N-eigenspectra)
X   [3] variance-of-time-series (required)
X       ==> variance of time series
X           from which eigenspectra were computed
X   [4] sampling-time (keyword; 1.0)
X       ==> sampling time (called delta t in the SAPA book)
X   [5] N-eigenspectra (keyword; length of list-of-eigenspectra)
X       ==> number of eigenspectra to be used (must be less than
X           or equal to the length of list-of-eigenspectra)
X   [6] maximum-number-of-iterations (keyword; 100)
X       ==> maximum number of iterations
X   [7] result-dof (keyword; vector of correct length)
X       <== vector into which degrees of freedom for each value
X           in the adaptive multitaper spectral estimate is placed;
X           it must be exactly the same length as each of the vectors
X           in list-of-eigenspectra
X   [8] sdf-transformation (keyword; #'convert-to-dB)
X       ==> a function of one argument or nil;
X           if bound to a function, the function is used
X           to transform all elements of result-sdf
X   [9] result-sdf (keyword; vector of correct length)
X       <== vector into which adaptive multitaper spectral estimate is placed;
X           it must be exactly the same length as each of the vectors
X           in list-of-eigenspectra
returns
X   [1] result-sdf, a vector holding
X       the properly transformed adaptive multitaper spectral estimate
X   [2] result-dof, a vector holding corresponding degrees of freedom
X   [3] the maximum number of iterations required to reach convergence
X       for all of the values in result-sdf
---
Note: see Section 7.4 of the SAPA book"
X  (assert (> (length list-of-eigenspectra) 1))
X  (let* ((N-f (length (car list-of-eigenspectra)))
X         (sig2*Delta-t (* variance-of-time-series sampling-time))
X         (lambda-0 (elt eigenvalues 0))
X         (lambda-1 (elt eigenvalues 1))
X         (lambda-0+lambda-1 (+ lambda-0 lambda-1))
X         (weights (make-array N-eigenspectra))
X         previous-est new-est-top new-est-bot
X         (max-iterations 0)
X         local-max-iterations
X         (first-eigenspectrum  (elt list-of-eigenspectra 0))
X         (second-eigenspectrum (elt list-of-eigenspectra 1)))
X    ;;; determine multitaper estimator, freq by freq ...
X    (dotimes (i N-f)
X      (setf previous-est (/ (+ (* lambda-0 (aref first-eigenspectrum i))
X                               (* lambda-1 (aref second-eigenspectrum i)))
X                            lambda-0+lambda-1))
X      (setf local-max-iterations 1)
X      (dotimes (j maximum-number-of-iterations
X                  (format t "~&nonconvergence at index ~D" i))
X        ;;; compute current set of weights and sums needed to new est ...
X        (setf new-est-top 0.0
X              new-est-bot 0.0)
X        (dotimes (k N-eigenspectra)
X          (setf (aref weights k)
X                (/ previous-est
X                   (+ (* (elt eigenvalues k) previous-est)
X                      (* (- 1.0 (elt eigenvalues k)) sig2*Delta-t))))
X          (incf new-est-top (* (aref weights k)
X                               (aref weights k)
X                               (elt eigenvalues k)
X                               (aref (elt list-of-eigenspectra k) i)))
X          (incf new-est-bot (* (aref weights k)
X                               (aref weights k)
X                               (elt eigenvalues k))))
X        (setf (aref result-sdf i) (/ new-est-top new-est-bot))
X        (if (< (/ (abs (- (aref result-sdf i) previous-est))
X                  previous-est) 0.05) (return))
X        (setf previous-est (aref result-sdf i))
X        (incf local-max-iterations))
X      (if (> local-max-iterations max-iterations)
X        (setf max-iterations local-max-iterations))
X      (setf (aref result-dof i) (/ (* 2.0 new-est-bot new-est-bot)
X                            (let ((the-sum 0.0))
X                              (dotimes (k n-eigenspectra the-sum)
X                                (incf the-sum
X                                      (* (expt (aref weights k) 4)
X                                         (elt eigenvalues k)
X                                         (elt eigenvalues k))))))))
X    (if sdf-transformation
X      (transform-a-sequence! sdf-transformation result-sdf))
X    (values result-sdf result-dof max-iterations)))
X
#|
;;; For this example we use the 20 point time series of Section 6.16
;;; of the SAPA book:
(let* ((20-pt-ts #(71.0  63.0  70.0  88.0  99.0  90.0 110.0 135.0 128.0 154.0
X                   156.0 141.0 131.0 132.0 141.0 104.0 136.0 146.0 124.0 129.0))
X       (N-ts (length 20-pt-ts))
X       (sigma^2 (sample-variance 20-pt-ts))
X       (sampling-time 0.25))
X  (multiple-value-bind (dpss-NW-4-tapers eigenvalues)
X                       (dpss-tapers-tri-diag
X                        N-ts 7
X                        :compute-true-eigenvalues-p t)
X    (multiple-value-bind (mt-4 freqs N-f list-of-eigenspectra)
X                         (multitaper-spectral-estimate
X                          20-pt-ts
X                          dpss-NW-4-tapers
X                          :N-nonzero-freqs :Fourier
X                          :sdf-transformation nil
X                          :sampling-time sampling-time)
X      (print (/ (+ (* 2 (sum mt-4 :end (1- N-f)))
X                   (svref mt-4 (1- N-f)))
X                (* N-ts sampling-time sigma^2)))
X      (multiple-value-bind (amt-4 dof)
X                           (eigenspectra->adaptive-multitaper-spectral-estimate 
X                            list-of-eigenspectra eigenvalues sigma^2
X                            :sdf-transformation nil)
X        (print (/ (+ (* 2 (sum amt-4 :end (1- N-f)))
X                     (svref amt-4 (1- N-f)))
X                  (* N-ts sampling-time sigma^2)))
X        (transform-a-sequence! #'convert-to-dB amt-4)
X        (dotimes (i N-f)
X          (format t "~&~6,4F: ~8,4F  ~5,1F"
X                  (svref freqs i)
X                  (svref amt-4 i)
X                  (svref dof i)))
X        (values)))))
;==>
1.0000000000000062 
0.9989962635770837 
0.2000:  27.3962   14.0
0.4000:  27.0953   14.0
0.6000:  25.3427   13.9
0.8000:  23.6051   13.9
1.0000:  19.4616   13.4
1.2000:  17.1247   13.0
1.4000:  18.2214   13.2
1.6000:  16.8903   13.0
1.8000:  16.9013   13.0
2.0000:  15.4482   12.7
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The function  create-ci-for-amt-sdf-estimate
;;;  takes the spectral estimate and vector of degrees of freedom returned by
;;;  eigenspectra->adaptive-multitaper-spectral-estimate and creates
;;;  a confidence interval for the true sdf at each frequency.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun create-ci-for-amt-sdf-estimate
X       (sdf-dB
X        dofs
X        &key
X        (confidence-level 0.95))
X    "given
X   [1] sdf-dB (required)
X       ==> vector containing an adaptive multitaper spectral
X           estimate expressed in decibels
X   [2] dofs (required)
X       ==> vector containing the degrees of freedom associated
X           with the values in sdf-dB
X   [3] confidence-level (keyword; 0.95)
X       ==> the level of the confidence intervals to be created
returns
X   [1] a vector containing the upper confidence interval
X   [2] a vector containing the lower confidence interval
---
Note: see Section 7.4 of the SAPA book"
X  (let* ((p (/ (- 1.0 confidence-level) 2.0))
X         (one-minus-p (- 1.0 p))
X         (n-f (length sdf-dB))
X         (upper-ci (copy-seq sdf-dB))
X         (lower-ci (copy-seq sdf-dB)))
X    (dotimes (i n-f (values upper-ci lower-ci))
X      (incf (aref upper-ci i)
X            (convert-to-dB (/ (aref dofs i)
X                              (quantile-of-chi-square-distribution
X                               (aref dofs i)
X                               p
X                               :accuracy-level :accurate))))
X      (incf (aref lower-ci i)
X            (convert-to-dB (/ (aref dofs i)
X                              (quantile-of-chi-square-distribution
X                               (aref dofs i)
X                               one-minus-p
X                               :accuracy-level :accurate)))))))
X
#|
;;; For this example we use the 20 point time series of Section 6.16
;;; of the SAPA book:
(let* ((20-pt-ts #(71.0  63.0  70.0  88.0  99.0  90.0 110.0 135.0 128.0 154.0
X                   156.0 141.0 131.0 132.0 141.0 104.0 136.0 146.0 124.0 129.0))
X       (N-ts (length 20-pt-ts))
X       (sigma^2 (sample-variance 20-pt-ts))
X       (sampling-time 0.25))
X  (multiple-value-bind (dpss-NW-4-tapers eigenvalues)
X                       (dpss-tapers-tri-diag
X                        N-ts 7
X                        :compute-true-eigenvalues-p t)
X    (multiple-value-bind (mt-4-dB freqs N-f list-of-eigenspectra)
X                         (multitaper-spectral-estimate
X                          20-pt-ts
X                          dpss-NW-4-tapers
X                          :N-nonzero-freqs :Fourier
X                          :sampling-time sampling-time)
X      (multiple-value-bind (amt-4-dB dof)
X                           (eigenspectra->adaptive-multitaper-spectral-estimate 
X                            list-of-eigenspectra eigenvalues sigma^2)
X        (multiple-value-bind (upper-dB lower-dB)
X                             (create-ci-for-amt-sdf-estimate
X                              amt-4-dB dof)
X          (dotimes (i N-f)
X            (format t "~&~6,4F: ~8,4F    ~8,4F ~8,4F ~8,4F"
X                    (svref freqs i)
X                    (svref mt-4-dB i)
X                    (svref lower-dB i)
X                    (svref amt-4-dB i)
X                    (svref upper-dB i)))
X        (values))))))
;==>
0.2000:  27.3532     24.6867  27.3962  31.3562
0.4000:  27.1112     24.3855  27.0953  31.0559
0.6000:  25.3879     22.6297  25.3427  29.3099
0.8000:  23.6973     20.8867  23.6051  27.5840
1.0000:  19.8763     16.7062  19.4616  23.5201
1.2000:  16.7150     14.3342  17.1247  21.2598
1.4000:  17.8006     15.4454  18.2214  22.3247
1.6000:  16.5783     14.0954  16.8903  21.0347
1.8000:  16.7305     14.1055  16.9013  21.0478
2.0000:  16.3822     12.6271  15.4482  19.6505
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  Everything below here consists of internal symbols in the SAPA package
;;;  and should be regarded as "dirty laundry" ...
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun dpss->eigenvalue (dpss NW)
X  "given dpss (a vector of length N) and NW,
computes the corresponding eigenvalue using
the method of Exercise [8.1], page 390,
of the SAPA book"
X  (let* ((eigenvalue 0.0)
X         (dpss-acvs
X          (acvs dpss :center-data-p nil))
X         (N (length dpss))
X         (W (/ NW N))
X         (vector-of-ratios (let ((ratios (make-array N))
X                                 (j 0))
X                             (setf (aref ratios 0) (* 2 W))
X                             (dotimes (i (1- N) ratios)
X                               (incf j)
X                               (setf (aref ratios j)
X                                     (/ (sin (* 2 pi W j))
X                                        (* pi j))))))
X         
X         (j (1- N)))
X    #+mcl(declare (dynamic-extent dpss-acvs vector-of-ratios))
X    (dotimes (i N)
X      (setf (aref dpss-acvs i)
X            (* N (aref dpss-acvs i))))
X    ;;; Note: both vector of ratios and dpss-acvs
X    ;;;       roughy decrease in amplitude with increasing index,
X    ;;;       so we sum things up in reverse order.
X    (dotimes (i (1- N) (+ (* 2 eigenvalue)
X                          (* 2 W (aref dpss-acvs 0))))
X      (incf eigenvalue
X            (* (aref dpss-acvs j)
X               (aref vector-of-ratios j)))
X      (decf j))))
X
#|
(multiple-value-bind (list-of-dpss eigenvalues)
X                     (dpss-tapers-inverse-iteration
X                      63 6 :taper-parameter 3.0)
X  (dotimes (k 6)
X    (let ((temp (dpss->eigenvalue (elt list-of-dpss k) 3.0)))
X      (format t "~&k = ~D: ~14,11F ~14,11F ~14,11F"
X              k
X              (svref eigenvalues k)
X              temp
X              (- (svref eigenvalues k) temp)))))
;==>
k = 0:  0.99999987329  0.99999987329 -0.000000000000
k = 1:  0.99999116408  0.99999116408 -0.000000000000
k = 2:  0.99972363504  0.99972363504 -0.000000000000
k = 3:  0.99500772991  0.99500779999 -0.00000007008
k = 4:  0.94627164450  0.94659849352 -0.00032684902
k = 5:  0.70835607249  0.70835638986 -0.00000031737
|#
X
;-------------------------------------------------------------------------------
(defun largest-eigenvalues-of-tridiagonal-matrix
X       (diag
X        off-diag
X        number-of-tapers
X        &key
X        (squared-off-diag (map 'vector #'(lambda (x) (* x x))  off-diag))
X        (macheps single-float-epsilon)
X        (print-progress-p nil))
X  (let* ((n (length diag))
X         (n-1 (1- n))
X         (n-2 (1- n-1)))
X    ;;; Set to zero all elements of squared-off-diag that correspond to
X    ;;; small elements of off-diag (see do-loop 40 of tridib) ...
X    (let ((previous-abs (aref diag 0)))
X      (dotimes (i n-1)
X        (if (<= (aref off-diag i) (* macheps (+ previous-abs
X                                                (setf previous-abs
X                                                      (aref diag (1+ i))))))
X          (setf (aref squared-off-diag i) 0.0))))
X    ;;; Use Equation (6) of Barth, Martin and Wilkinson to find
X    ;;; upper and lower bounds for all eigenvalues ...
X    (let* ((abs-off-diag-behind (abs (aref off-diag 0)))
X           (abs-off-diag-ahead  (abs (aref off-diag n-2)))
X           sum-of-abs
X           (lower-bound-all-eigenvalues
X            (min (- (aref diag 0) abs-off-diag-behind)
X                 (- (aref diag n-1) abs-off-diag-ahead)))
X           (upper-bound-all-eigenvalues
X            (max (+ (aref diag 0) abs-off-diag-behind)
X                 (+ (aref diag n-1) abs-off-diag-ahead)))
X           (i 0))
X      (dotimes (j n-2)
X        (setf abs-off-diag-ahead (abs (aref off-diag (incf i)))
X              sum-of-abs (+ abs-off-diag-behind abs-off-diag-ahead))
X        (if (> lower-bound-all-eigenvalues (- (aref diag i) sum-of-abs))
X          (setf lower-bound-all-eigenvalues (- (aref diag i) sum-of-abs)))
X        (if (< upper-bound-all-eigenvalues (+ (aref diag i) sum-of-abs))
X          (setf upper-bound-all-eigenvalues (+ (aref diag i) sum-of-abs)))
X        (setf abs-off-diag-behind abs-off-diag-ahead))
X      ;;; Expand upper and lower bounds a little (evidently to avoid
X      ;;; numerical problems -- see code following do-loop 40 of tridib) ...
X      (let ((eigenvalues (make-array number-of-tapers))
X            (upper-bounds (make-array number-of-tapers
X                                      :initial-element
X                                      upper-bound-all-eigenvalues))
X            (lower-bounds (make-array number-of-tapers
X                                      :initial-element
X                                      lower-bound-all-eigenvalues)))
X        (dotimes (k number-of-tapers (values eigenvalues))
X          (setf (aref eigenvalues k)
X                ;;; use bisection to isolate eigenvalues, keeping track
X                ;;; of bisections at each k to speed up subsequent k's
X                (do* ((update-other-bounds (< k (1- number-of-tapers)))
X                      (upper-target-count (- n k))
X                      (lower-target-count (1- upper-target-count))
X                      current-count
X                      (U (svref upper-bounds k))
X                      (L (svref lower-bounds k))
X                      (mid-point (/ (+ U L) 2.0) (/ (+ U L) 2.0)))
X                     ((<= (abs (- U L)) (* macheps (+ (abs U) (abs L))))
X                      mid-point)
X                  (setf current-count (sturm-sequence-count
X                                       mid-point
X                                       diag
X                                       off-diag
X                                       squared-off-diag
X                                       macheps))
X                  (when update-other-bounds
X                    (let* ((N-current-count (- n current-count))
X                           (j-to-inc N-current-count))
X                      (dotimes (j (min N-current-count number-of-tapers))
X                        (if (> mid-point (svref lower-bounds j))
X                          (setf (svref lower-bounds j) mid-point)))
X                      (dotimes (j (- number-of-tapers N-current-count))
X                        (if (< mid-point (svref upper-bounds j-to-inc))
X                          (setf (svref upper-bounds j-to-inc) mid-point))
X                        (incf j-to-inc)))
X                    (if (>= current-count lower-target-count)
X                      (setf update-other-bounds nil)))
X                  (if (<= current-count lower-target-count)
X                    (setf L mid-point)
X                    (setf U mid-point))))
X          (if print-progress-p (format t ".")))))))
X
;-------------------------------------------------------------------------------
(defun sturm-sequence-count
X       (test-lambda
X        diag
X        off-diag
X        squared-off-diag
X        machep
X        &key
X        (start 0)
X        (end (length diag)))
X  (let ((count start)
X        (bottom 1.0))
X    (dotimes (i (- end start) count)
X      (let ((ratio (cond
X                    ((zerop start)
X                     0.0)
X                    ((zerop bottom)
X                     (if (zerop (aref squared-off-diag (1- start)))
X                       0.0
X                       (/ (abs (aref off-diag (1- start))) machep)))
X                    (t
X                     (/ (aref squared-off-diag (1- start))
X                        bottom)))))
X        (setf bottom
X              (- (aref diag start) test-lambda ratio))
X        (incf start)
X        (if (minusp bottom) (incf count))))))
X
;-------------------------------------------------------------------------------
(defun symmetric-tridiagonal-solve! (diag off-diag b)
X  "given
X   [1] diag (required)
X       <=> diagonal part of symmetric tridiagonal matrix A;
X           trash on output
X   [1] off-diag (required)
X       <=> off-diagonal part of symmetric tridiagonal matrix;
X           trash on output
X   [2] b (required)
X       <=> on input, the right-hand side vector;
X           on output, the solution X to A X = b
returns
X   [4] X, the solution to A X = b, where X is the vector
X       that contained b on input
---
Note: this is an implementation of Algorithm 4.3.6,
p. 156, Golub and Van Loan, 1989, with modifications
to avoid divides by zero"
X  (let* (temp
X         (n (length diag))
X         (n-1 (1- n))
X         (zero-replacement (* single-float-epsilon
X                              (symmetric-tridiagonal-infinity-norm
X                               diag off-diag))))
X    (if (< (abs (aref diag 0)) single-float-epsilon)
X      (setf (aref diag 0) zero-replacement))
X    (dotimes (k n-1)
X      (setf temp (aref off-diag k)
X            (aref off-diag k) (/ temp (aref diag k)))
X      (decf (aref diag (1+ k)) (* temp (aref off-diag k)))
X      (if (< (abs (aref diag (1+ k))) single-float-epsilon)
X        (setf (aref diag (1+ k)) zero-replacement)))
X    (dotimes (k n-1)
X      (decf (aref b (1+ k)) (* (aref b k) (aref off-diag k))))
X    (setf (aref b n-1) (/ (aref b n-1) (aref diag n-1)))
X    (let ((k n-1))
X      (dotimes (j n-1 (values b))
X        (decf k)
X        (setf (aref b k) (- (/ (aref b k) (aref diag k))
X                            (* (aref off-diag k) (aref b (1+ k)))))))))
X
;-------------------------------------------------------------------------------
(defun symmetric-tridiagonal-infinity-norm (diag off-diag)
X  (let* ((n (length diag))
X         (n-2 (- n 2))
X         (infinity-norm
X          (max (+ (abs (aref diag 0)) (abs (aref off-diag 0)))
X               (+ (abs (aref diag (1+ n-2))) (abs (aref off-diag n-2)))))
X         (j 0))
X    (dotimes (i n-2 infinity-norm)
X      (if (> (+ (abs (aref off-diag j))
X                (abs (aref off-diag (incf j)))
X                (abs (aref diag j)))
X             infinity-norm)
X        (setf infinity-norm
X              (+ (abs (aref off-diag (1- j)))
X                 (abs (aref off-diag j))
X                 (abs (aref diag j))))))))
X           
;-------------------------------------------------------------------------------
(defun fast-tridiag-eigenvalue->dpss!
X       (eigenvalue
X        order-of-taper
X        diag
X        off-diag
X        &key
X        (eps (* 10 single-float-epsilon))
X        (maximum-number-of-iterations 25)
X        (b (generate-initial-guess-at-dpss (length diag) order-of-taper
X                                           :half-size-p t)))
X  (let* ((N (length diag))
X         (M (/ (if (evenp N) N (1- N)) 2))
X         (shorter-diag (make-array M))
X         (shorter-off-diag (make-array M))
X         (scratch-diag (make-array M))
X         (scratch-off-diag (make-array M)))
X    #+mcl(declare (dynamic-extent shorter-diag shorter-off-diag
X                                  scratch-diag scratch-off-diag))
X    (replace shorter-diag diag)
X    (replace shorter-off-diag off-diag)
X    (if (evenp N)
X      ;;; even N
X      (if (evenp order-of-taper)
X        (incf (svref shorter-diag (1- M))
X              (aref off-diag (1- M)))
X        (decf (svref shorter-diag (1- M))
X              (aref off-diag (1- M))))
X      ;;; odd N
X      (if (evenp order-of-taper)
X        (incf (svref shorter-diag (1- M))
X              (/ (* 2.0 (expt (aref off-diag (1- M)) 2)) eigenvalue))))
X    (x+b! shorter-diag (- eigenvalue))
X    (let ((b-old (copy-seq b))
X          (end-sum (max 2 (truncate N (* 2 (1+ order-of-taper)))))
X          iter-total)
X      (setf iter-total
X            (dotimes (i maximum-number-of-iterations)
X              (symmetric-tridiagonal-solve!
X               (replace scratch-diag shorter-diag)
X               (replace scratch-off-diag shorter-off-diag)
X               b)
X              (a*x! (* (/ (sqrt (sum-of-squares b)))
X                       (if (plusp (sum b :end end-sum))
X                         1.0
X                         -1.0))
X                    b)
X              (if (< (compare-seqs b b-old) eps) (return (values (1+ i)))
X                  (replace b-old b))))
X      (let ((result-dpss (make-array N)))
X        (replace result-dpss b)
X        (cond
X         ((evenp N)
X          (let ((i-rev N))
X            (if (evenp order-of-taper)
X              (dotimes (i M)
X                (setf (svref result-dpss (decf i-rev))
X                      (svref result-dpss i)))
X              (dotimes (i M)
X                (setf (svref result-dpss (decf i-rev))
X                      (* -1 (svref result-dpss i)))))))
X         ;;; N is odd - branch according to order of taper
X         ((evenp order-of-taper)
X          (let ((i-rev N))
X            (dotimes (i M)
X              (setf (svref result-dpss (decf i-rev))
X                    (svref result-dpss i)))
X            (setf (svref result-dpss M)
X                  (/ (* 2
X                        (svref result-dpss (1- M))
X                        (aref off-diag (1- M)))
X                     eigenvalue))))
X         ;;; N is odd, and order of taper is odd
X         (t 
X          (let ((i-rev N))
X            (dotimes (i M)
X              (setf (svref result-dpss (decf i-rev))
X                    (* -1 (svref result-dpss i))))
X            (setf (svref result-dpss M) 0.0))))
X        (values (a*x! (/ (sqrt (sum-of-squares result-dpss))) result-dpss)
X                iter-total)))))
X
;-------------------------------------------------------------------------------
(defun generate-initial-guess-at-dpss
X       (N
X        order-of-taper
X        &key
X        (half-size-p nil)
X        (result (make-array (if half-size-p
X                              (/ (if (evenp N) N (1- N)) 2)
X                              N))))
X  (iota 1 (length result) :result result)
X  (a*x! (float (/ (1+ N))) result)
X  (map-into result #'(lambda (x)
X                       (* x x (expt (1- x) 2)
X                          (let ((prod 1.0)
X                                (order+1 (1+ order-of-taper)))
X                            (dotimes (j order-of-taper prod)
X                              (multf prod (- (/ (1+ j) order+1) x))))))
X            result)
X  (let ((factor (/ (sqrt (sum-of-squares result)))))
X    (a*x! factor result)))
X
;-------------------------------------------------------------------------------
;;; not currently used, but useful for pulling out a single eigenvalue ...
(defun kth-eigenvalue-of-tridiagonal-matrix
X       (diag
X        off-diag
X        k
X        &key
X        (squared-off-diag (map 'vector #'(lambda (x) (* x x))  off-diag))
X        (macheps single-float-epsilon))
X  (let* ((n (length diag))
X         (n-1 (1- n))
X         (n-2 (1- n-1)))
X    ;;; Set to zero all elements of squared-off-diag that correspond to
X    ;;; small elements of off-diag (see do-loop 40 of tridib) ...
X    (let ((previous-abs (aref diag 0)))
X      (dotimes (i n-1)
X        (if (<= (aref off-diag i) (* macheps (+ previous-abs
X                                                (setf previous-abs
X                                                      (aref diag (1+ i))))))
X          (setf (aref squared-off-diag i) 0.0))))
X   ;;; Use Equation (6) of Barth, Martin and Wilkinson to find
X   ;;; upper and lower bounds for all eigenvalues ...
X    (let* ((abs-off-diag-behind (abs (aref off-diag 0)))
X           (abs-off-diag-ahead  (abs (aref off-diag n-2)))
X           sum-of-abs
X           (lower-bound-all-eigenvalues
X            (min (- (aref diag 0) abs-off-diag-behind)
X                 (- (aref diag n-1) abs-off-diag-ahead)))
X           (upper-bound-all-eigenvalues
X            (max (+ (aref diag 0) abs-off-diag-behind)
X                 (+ (aref diag n-1) abs-off-diag-ahead)))
X           (i 0))
X      (dotimes (j n-2)
X        (setf abs-off-diag-ahead (abs (aref off-diag (incf i)))
X              sum-of-abs (+ abs-off-diag-behind abs-off-diag-ahead))
X        (if (> lower-bound-all-eigenvalues (- (aref diag i) sum-of-abs))
X          (setf lower-bound-all-eigenvalues (- (aref diag i) sum-of-abs)))
X        (if (< upper-bound-all-eigenvalues (+ (aref diag i) sum-of-abs))
X          (setf upper-bound-all-eigenvalues (+ (aref diag i) sum-of-abs)))
X        (setf abs-off-diag-behind abs-off-diag-ahead))
X      ;;; Expand upper and lower bounds a little (evidently to avoid
X      ;;; numerical problems -- see code following do-loop 40 of tridib) ...
X      (let ((delta (* N macheps (max (abs lower-bound-all-eigenvalues)
X                                     (abs upper-bound-all-eigenvalues)))))
X        (decf lower-bound-all-eigenvalues delta)
X        (incf upper-bound-all-eigenvalues delta))
X      ;;; use bisection to isolate kth eigenvalue ...
X      (do* ((upper-target-count (- n k))
X            (lower-target-count (1- upper-target-count))
X            current-count
X            (U upper-bound-all-eigenvalues)
X            (L lower-bound-all-eigenvalues)
X            (mid-point (/ (+ U L) 2.0) (/ (+ U L) 2.0)))
X           ((<= (abs (- U L)) (* macheps (+ (abs U) (abs L))))
X            mid-point)
X        (setf current-count (sturm-sequence-count
X                             mid-point
X                             diag
X                             off-diag
X                             squared-off-diag
X                             macheps))
X        (if (<= current-count lower-target-count)
X          (setf L mid-point)
X          (setf U mid-point))))))
X
;-------------------------------------------------------------------------------
;;; not currently used -- superceded by fast-tridiag-eigenvalue->dpss! ...
(defun tridiag-eigenvalue->dpss!
X       (eigenvalue
X        order-of-taper
X        diag
X        off-diag
X        &key
X        (eps (* 10 single-float-epsilon))
X        (maximum-number-of-iterations 25)
X        (b (generate-initial-guess-at-dpss (length diag) order-of-taper)))
X  (x+b! diag (- eigenvalue))
X  (let* ((N (length diag))
X         (b-old (copy-seq b))
X         (end-sum (max 2 (truncate N (1+ order-of-taper)))))
X    (dotimes (i maximum-number-of-iterations (values b nil))
X      (symmetric-tridiagonal-solve! (copy-seq diag)
X                                       (copy-seq off-diag)
X                                       b)
X       (a*x! (* (/ (sqrt (sum-of-squares b)))
X               (if (plusp (sum b :end end-sum))
X                 1.0
X                 -1.0))
X            b)
X      (if (< (compare-seqs b b-old) eps) (return (values b (1+ i))))
X      (replace b-old b))))
X
;-------------------------------------------------------------------------------
(defun toeplitz
X       (r
X        y
X        &optional
X        (x (make-array (length y))))
X  "given
X   [1] r (required)
X       ==> vector of length 2*N-1 representing
X           the Toeplitz matrix R(i,j) with 0 <= i,j <= N-1,
X           via the 2*N-1 values
X           r(0,N-1), r(0,N-2),..., r(0,1),
X           r(0,0), r(1,0),..., r(N-1,0)
X   [2] y (required)
X       ==> right-hand side vector of length N
X   [3] x (optional; vector of same size as y)
X       <== solution vector of length N
solves the Toeplitz system sum_{j=0}^{N-1} r_{N-1+i-j} x_j = y_j
for i = 0, 1, ..., N-1 and
returns
X   [1] x (a vector of length N)
X       or
X       nil (in case of a failure of the Levinson method
X       due to a singular principal minor
---
Note: this function is essentially a Lisp version
X      of the Fortran routine toeplz from Numerical Recipes"
X  (let* ((n (array-total-size y))
X         (index-of-0 (1- n))
X         (g (make-array n))
X         (h (make-array n))
X         m1 sxn sd sgn shn sgd k pp qq pt1 pt2 qt1 qt2)
X    (block main
X      (if (zerop (aref r index-of-0)) (return-from main nil))
X      (setf (aref x 0) (/ (aref y 0)
X                          (aref r index-of-0)))
X      (if (= n 1) (return-from main x))
X      (setf (aref g 0) (/ (aref r (1- index-of-0))
X                          (aref r index-of-0))
X            (aref h 0) (/ (aref r (1+ index-of-0))
X                          (aref r index-of-0)))
X      (dotimes (m n)
X        (setf m1 (1+ m)
X              sxn (- (aref y m1))
X              sd (- (aref r index-of-0)))
X        (dotimes (j (1+ m))
X          (incf sxn (* (aref r (+ index-of-0 (- m1 j)))
X                       (aref x j)))
X          (incf sd (* (aref r (+ index-of-0 (- m1 j)))
X                      (aref g (- m j)))))
X        (if (zerop sd) (return-from main nil))
X        (setf (aref x m1) (/ sxn sd))
X        (dotimes (j (1+ m))
X          (decf (aref x j) (* (aref x m1) (aref g (- m j)))))
X        (if (= (1+ m1) n) (return-from main x))   ;usual exit point
X        (setf sgn (- (aref r (- index-of-0 (1+ m1))))
X              shn (- (aref r (+ index-of-0 (1+ m1))))
X              sgd (- (aref r index-of-0)))
X        (dotimes (j (1+ m))
X          (incf sgn (* (aref r (+ index-of-0 (- j m1)))
X                       (aref g j)))
X          (incf shn (* (aref r (+ index-of-0 (- m1 j)))
X                       (aref h j)))
X          (incf sgd (* (aref r (+ index-of-0 (- j m1)))
X                       (aref h (- m j)))))
X        (if (or (zerop sd) (zerop sgd)) (return nil))
X        (setf (aref g m1) (/ sgn sgd)
X              (aref h m1) (/ shn sd)
X              k m
X              pp (aref g m1)
X              qq (aref h m1))
X        (dotimes (j (truncate (+ m 2) 2))
X          (setf pt1 (aref g j)
X                pt2 (aref g k)
X                qt1 (aref h j)
X                qt2 (aref h k)
X                (aref g j) (- pt1 (* pp qt2))
X                (aref g k) (- pt2 (* pp qt1))
X                (aref h j) (- qt1 (* qq pt2))
X                (aref h k) (- qt2 (* qq pt1)))
X          (decf k))))
X    ;;; NOTE: this function will exit from this point (with a value of nil
X    ;;;       since this is the value of dotimes) only if there was some sort of
X    ;;;       gross screw-up
X    ))
X
#|
;;; [ 3 2 0 ] [ x_0 ]   [  7 ]
;;; [ 2 3 2 ] [ x_1 ] = [ 14 ]
;;; [ 1 2 3 ] [ x_2 ]   [ 14 ]
(let* ((r-test #(0 2 3 2 1))
X       (y-test-1 #(7 14 14))
X       (y-test-2 #(3 2 1))
X       (y-test-3 #(2 3 2))
X       (y-test-4 #(0 2 3)))
X  (print (toeplitz r-test y-test-1))
X  (print (toeplitz r-test y-test-2))
X  (print (toeplitz r-test y-test-3))
X  (print (toeplitz r-test y-test-4)))
;==> #(1 2 3) 
;    #(1 0 0) 
;    #(0 1 0) 
;    #(0 0 1) 
;    #(0 0 1)
|#
X
;-------------------------------------------------------------------------------
(defun Householder-reduction-of-real-sym-matrix!
X       (real-sym-matrix &key (prepare-for-eigenvectors-p t))
X  "given
X   [1] real-sym-matrix (required)
X       ==> real, symmetric matrix (n by n)
X   [2] prepare-for-eigenvectors-p (keyword; t)
X       ==> if true, crunch real-sym-matrix for use with
X           eigenvalues-and-vectors-of-sym-tridiag-matrix
returns
X   [1] n-dimensional vector with diagonal elements of
X       associated tridiagonal matrix
X   [2] (n-1)-dimensional vector with off-diagonal elements of
X       associated tridiagonal matrix
---
See tred2 in Numerical Recipes"
X  (let* ((n (car (array-dimensions real-sym-matrix)))
X         (diag (make-array n))
X         (off-diag (make-array n))
X         i-local l-local f g h hh scale)
X    (if (> n 1)
X      (dotimes (i-index (1- n))
X        (setf i-local (- n i-index 1)
X              l-local (- i-local 1)
X              h 0.0
X              scale 0.0)
X        (cond ((plusp l-local)
X               ;;; BEGIN LOOP 11
X               (dotimes (k (1+ l-local))
X                 (incf scale (abs (aref real-sym-matrix i-local k))))
X               ;;; END LOOP 11
X               (cond ((zerop scale)
X                      (setf (aref off-diag i-local)
X                            (aref real-sym-matrix i-local l-local)))
X                     (t
X                      ;;; BEGIN LOOP 12
X                      (dotimes (k (1+ l-local))
X                        (setf (aref real-sym-matrix i-local k)
X                              (/ (aref real-sym-matrix i-local k) scale))
X                        (incf h (expt (aref real-sym-matrix i-local k) 2)))
X                      ;;; END LOOP 12
X                      (setf f (aref real-sym-matrix i-local l-local)
X                            g (- (sign (sqrt h) f))
X                            (aref off-diag i-local) (* scale g)
X                            h (- h (* f g))
X                            (aref real-sym-matrix i-local l-local) (- f g)
X                            f 0.0)
X                      ;;; BEGIN LOOP 15
X                      (dotimes (j (1+ l-local))
X                        (if prepare-for-eigenvectors-p
X                          (setf (aref real-sym-matrix j i-local)
X                                (/ (aref real-sym-matrix i-local j) h)))
X                        (setf g 0.0)
X                        ;;; BEGIN LOOP 13
X                        (dotimes (k (1+ j))
X                          (incf g (* (aref real-sym-matrix j k)
X                                     (aref real-sym-matrix i-local k))))
X                        ;;; END LOOP 13
X                        (if (> l-local j)
X                          ;;; BEGIN LOOP 14
X                          (dotimes (k (- l-local j))
X                            (incf g (* (aref real-sym-matrix (+ k j 1) j)
X                                       (aref real-sym-matrix i-local (+ k j 1)))))
X                          ;;; END LOOP 14
X                          )
X                        (setf (aref off-diag j) (/ g h)
X                              f (+ f (* (aref off-diag j)
X                                        (aref real-sym-matrix i-local j)))))
X                      ;;; END LOOP 15
X                      (setf hh (/ f (+ h h)))
X                      ;;; BEGIN LOOP 17
X                      (dotimes (j (1+ l-local))
X                        (setf f (aref real-sym-matrix i-local j)
X                              g (- (aref off-diag j) (* hh f))
X                              (aref off-diag j) g)
X                        ;;; BEGIN LOOP 16
X                        (dotimes (k (1+ j))
X                          (decf (aref real-sym-matrix j k)
X                                (+ (* f (aref off-diag k))
X                                   (* g (aref real-sym-matrix i-local k)))))
X                        ;;; END LOOP 16
X                        )
X                      ;;; END LOOP 17
X                      )))
X              (t
X               (setf (aref off-diag i-local)
X                     (aref real-sym-matrix i-local l-local))))
X        (setf (aref diag i-local) h)))
X    ;;; END OF ``IF'' CLAUSE ...
X    (if prepare-for-eigenvectors-p
X      (setf (aref diag 0) 0.0))
X    (setf (aref off-diag 0) 0.0)
X    ;;; BEGIN LOOP 23
X    (dotimes (i-local n)
X      (when prepare-for-eigenvectors-p
X        (setf l-local (1- i-local))
X        (when (/= (aref diag i-local) 0.0)
X          (dotimes (j (1+ l-local))
X            (setf g 0.0)
X            (dotimes (k (1+ l-local))
X              (incf g (* (aref real-sym-matrix i-local k)
X                         (aref real-sym-matrix k j))))
X            (dotimes (k (1+ l-local))
X              (decf (aref real-sym-matrix k j)
X                    (* g (aref real-sym-matrix k i-local)))))))
X      (setf (aref diag i-local) (aref real-sym-matrix i-local i-local))
X      (when prepare-for-eigenvectors-p
X        (setf (aref real-sym-matrix i-local i-local) 1.0)
X        (when (>= l-local 0)
X          (dotimes (j (1+ l-local))
X            (setf (aref real-sym-matrix i-local j) 0.0
X                  (aref real-sym-matrix j i-local) 0.0)))))
X    ;;; END LOOP 23
X    (values diag (make-array (1- n)
X                             :displaced-to off-diag
X                             :displaced-index-offset 1))))
X
;-------------------------------------------------------------------------------
(defun eigenvalues-and-vectors-of-sym-tridiag-matrix
X       (diag
X        off-diag
X        &key
X        (maximum-number-of-iterations 30)
X        (return-eigenvectors-p t)
X        (z (if return-eigenvectors-p
X             (let* ((n (length diag))
X                    (temp (make-array `(,n ,n)
X                                      :initial-element 0.0)))
X               (dotimes (i n temp)
X                 (setf (aref temp i i) 1.0))))))
X  "given
X   [1] diag (required)
X       ==> sequence of diagonal elements
X           (length of, say, n)
X   [2] off-diag (required)
X       ==> sequence of off-diagonal elements
X           (length of n-1)
X   [3] maximum-number-of-iterations (keyword; 30)
X       ==> maximum number of iterations
X   [4] return-eigenvectors-p (keyword; t)
X       ==> if true, returns eigenvalues and eigenvectors;
X           if nil, just do the eigenvalues (this is faster)
X   [5] z (keyword; n by n diagonal matrix)
X       ==> usually bound to output from
X           Householder-reduction-of-real-sym-matrix
returns
X   [1] vector with n eigenvalues
X       (NOT necessarily ordered by size)
X   [2] if return-eigenvectors-p is true, an n by n array
X       whose jth column is the eigenvector corresponding
X       to the jth eigenvalue; if return-eigenvectors-p is nil,
X       that's what you get!
---
This is based upon tqli in Numerical Recipes"
X  (assert (= (length diag) (1+ (length off-diag))))
X  (let* ((n (length diag))
X         (d (copy-seq diag))
X         (e (make-array n :initial-element 0.0))
X         dd m m-local g s c p i-local f b r)
X    (cond
X     ((<= n 1) d)
X     (t
X      (dotimes (i (1- n))
X        (setf (elt e i) (elt off-diag i)))
X      ;;; LOOP 15
X      ;;; routine returns d --- vector of eigenvalues
X      (dotimes (l n (values d z))
X        (do ((iteration-count 0 (+ iteration-count 1)))
X            (nil)    ;do forever
X          (setf m (dotimes (mm (- n l 1) (1- n))
X                    (setf m-local (+ mm l)
X                          dd (+ (abs (elt d m-local))
X                                (abs (elt d (1+ m-local)))))
X                    (if (= (+ (abs (elt e  m-local)) dd)
X                           dd)
X                      (return m-local))))
X          (cond
X           ((= m l) (return))
X           (t
X            (assert (< iteration-count maximum-number-of-iterations))
X            (setf g (/ (- (elt d (1+ l)) (elt d l))
X                       (* 2.0 (elt e l)))
X                  r (sqrt (+ 1.0 (* g g)))
X                  g (+ (/ (elt e l) (+ g (sign r g)))
X                       (elt d m)
X                       (- (elt d l)))
X                  s 1.0
X                  c 1.0
X                  p 0.0)
X            ;;; LOOP 14
X            (dotimes (i (- m l))
X              (setf i-local (- m 1 i)
X                    f (* s (elt e i-local))
X                    b (* c (elt e i-local)))
X              (if (>= (abs f) (abs g))
X                (setf c (/ g f)
X                      r (sqrt (+ 1.0 (* c c)))
X                      (elt e (1+ i-local)) (* f r)
X                      s (/ r)
X                      c (* c s))
X                (setf s (/ f g)
X                      r (sqrt (+ 1.0 (* s s)))
X                      (elt e (1+ i-local)) (* g r)
X                      c (/ r)
X                      s (* s c)))
X              (setf g (- (elt d (1+ i-local)) p)
X                    r (+ (* s (- (elt d i-local) g))
X                         (* 2.0 c b))
X                    p (* s r)
X                    (elt d (1+ i-local)) (+ g p)
X                    g (- (* c r) b))
X              (when return-eigenvectors-p
X                (dotimes (k n)
X                  (setf f (aref z k (1+ i-local))
X                        (aref z k (1+ i-local)) (+ (* s (aref z k i-local))
X                                                   (* c f))
X                        (aref z k i-local) (- (* c (aref z k i-local))
X                                              (* s f)))))
X              )
X            (setf (elt d l) (-  (elt d l) p)
X                  (elt e l) g
X                  (elt e m) 0.0)))))))))
X
#|
;;; The following set of abscissa points and weights
;;; are used as defaults in dpss-tapers-Thomson-approx and
;;; were computed using the following Lisp form:
(multiple-value-setq (*abscissas-32-point* *weights-32-point*)
X  (gauss-legendre-quadrature -1.0 1.0 32))
|#
(defvar *abscissas-32-point* (vector -0.9972638618494816D0
X                                     -0.9856115115452684D0
X                                     -0.9647622555875064D0
X                                     -0.9349060759377397D0
X                                     -0.8963211557660521D0
X                                     -0.84936761373257D0
X                                     -0.7944837959679424D0
X                                     -0.7321821187402897D0
X                                     -0.6630442669302152D0
X                                     -0.5877157572407623D0
X                                     -0.5068999089322294D0
X                                     -0.4213512761306353D0
X                                     -0.3318686022821277D0
X                                     -0.2392873622521371D0
X                                     -0.1444719615827965D0
X                                     -4.830766568773831D-2 
X                                     4.830766568773831D-2 
X                                     0.1444719615827965D0
X                                     0.2392873622521371D0
X                                     0.3318686022821277D0
X                                     0.4213512761306353D0
X                                     0.5068999089322294D0
X                                     0.5877157572407623D0
X                                     0.6630442669302152D0
X                                     0.7321821187402897D0
X                                     0.7944837959679424D0
X                                     0.84936761373257D0
X                                     0.8963211557660521D0
X                                     0.9349060759377397D0
X                                     0.9647622555875064D0
X                                     0.9856115115452684D0
X                                     0.9972638618494816D0))
X
(defvar *weights-32-point* (vector 7.018610009469521D-3 
X                                   1.627439473090571D-2 
X                                   2.539206530926214D-2 
X                                   3.427386291302141D-2 
X                                   4.283589802221816D-2 
X                                   0.0509980592623745D0
X                                   5.868409347853513D-2 
X                                   6.582222277636181D-2 
X                                   7.234579410884862D-2 
X                                   7.819389578707044D-2 
X                                   8.331192422694672D-2 
X                                   8.765209300440374D-2 
X                                   0.0911738786957639D0
X                                   9.384439908080441D-2 
X                                   9.563872007927485D-2 
X                                   9.654008851472785D-2 
X                                   9.654008851472785D-2 
X                                   9.563872007927485D-2 
X                                   9.384439908080441D-2 
X                                   0.0911738786957639D0
X                                   8.765209300440374D-2 
X                                   8.331192422694672D-2 
X                                   7.819389578707044D-2 
X                                   7.234579410884862D-2 
X                                   6.582222277636181D-2 
X                                   5.868409347853513D-2 
X                                   0.0509980592623745D0
X                                   4.283589802221816D-2 
X                                   3.427386291302141D-2 
X                                   2.539206530926214D-2 
X                                   1.627439473090571D-2 
X                                   7.018610009469521D-3))
X
;-------------------------------------------------------------------------------
;;; This and the following two functions are used to generate Greenhall's
;;; trig prolate data tapers.
(defun trig-prolate-sign-cosine-coefficients (M k)
X  (case M
X    (2 (case k
X         (0 '(0.8202108 0.4041691 0.0165649))
X         (1 '(0.0 -0.7007932 -0.0942808))
X         (otherwise
X          (error "can't handle k = ~A for M = 2" k ))))
X    (3 (case k
X         (0 '(0.7499700 0.4596063 0.0867984 0.0007513))
X         (1 '(0.0 -0.6507499 -0.2765560 -0.0064282))
X         (2 '(0.4969513 -0.3050683 -0.5312499 -0.0350227))
X         (3 '(0.0 -0.2731233 0.6397174 0.1271430))
X         (otherwise
X          (error "can't handle k = ~A for M = 3" k ))))
X    (4 (case k
X         (0 '(0.6996910 0.4830013 0.1473918 0.0141997 0.0000368))
X         (1 '(0.0 -0.5927723 -0.3805986 -0.0613650 -0.0003329))
X         (2 '(0.4783016 -0.1666510 -0.5724443 -0.1736202 -0.0022015))
X         (3 '(0.0 -0.3540569 0.4929565 0.3626279 0.0117722))
X         (4 '(0.3862293 -0.3223025 0.0856254 0.5584413 0.0484379))
X         (otherwise
X          (error "can't handle k = ~A for M = 4" k ))))
X    (5 (case k
X         (0 '(0.6632850 0.4915713 0.1927963 0.0347859 0.0019243 0.0000018))
X         (1 '(0.0 -0.5401300 -0.4383060 -0.1266343 -0.0105462 -0.0000191))
X         (2 '(0.4560698 -0.0704481 -0.5519198 -0.2915206 -0.0379143 -0.0001319))
X         (3 '(0.0 -0.3866087 0.3363930 0.4760267 0.1037856 0.0007467))
X         (4 '(0.3821638 -0.2527019 -0.1138304 0.5457777 0.2286313 0.0037712))
X         (5 '(0.0 -0.2216043 0.3885522 -0.3657298 -0.4072901 -0.0165910))
X         (6 '(0.3246026 -0.2957322 0.1964585 0.0266965 -0.5631039 -0.0588589))
X         (otherwise
X          (error "can't handle k = ~A for M = 5" k ))))
X    (otherwise
X     (error "can't handle M = ~A" M))))
X
;-------------------------------------------------------------------------------
(defun produce-upsilon-func (coefficients)
X  (if (zerop (car coefficients))
X    #'(lambda (t-index)
X        (let ((sum 0.0))
X          (dotimes (j (1- (length coefficients)) sum)
X            (incf sum (* 2.0 (elt coefficients (1+ j))
X                         (sin (* 2.0 pi (1+ j) t-index)))))))
X    #'(lambda (t-index)
X        (let ((sum (car coefficients)))
X          (dotimes (j (1- (length coefficients)) sum)
X            (incf sum (* 2.0 (elt coefficients (1+ j))
X                         (cos (* 2.0 pi (1+ j) t-index)))))))))
X
;-------------------------------------------------------------------------------
(defun generate-tri-prolate-taper (N generating-function)
X  (let ((the-taper (make-array N))
X        (factor (float (/ (sqrt N))))
X        (factor-top (float (1- N)))
X        (factor-bot (float (* 2 N))))
X    (dotimes (i N the-taper)
X      (setf (aref the-taper i)
X            (* factor
X               (funcall generating-function
X                        (/ (- (* 2 i) factor-top)
X                           factor-bot)))))))
SHAR_EOF
chmod 0644 multitaper.lisp ||
echo 'restore of multitaper.lisp failed'
Wc_c="`wc -c < 'multitaper.lisp'`"
test 90780 -eq "$Wc_c" ||
	echo 'multitaper.lisp: original size 90780, current size' "$Wc_c"
fi
# ============= nonparametric.lisp ==============
if test -f 'nonparametric.lisp' -a X"$1" != X"-c"; then
	echo 'x - skipping nonparametric.lisp (File already exists)'
else
echo 'x - extracting nonparametric.lisp (Text)'
sed 's/^X//' << 'SHAR_EOF' > 'nonparametric.lisp' &&
;;;-*- Mode: LISP; Package: :SAPA; Syntax: COMMON-LISP -*-
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;
;  nonparametric.lisp
;
;  a collection of Lisp functions for nonparametric spectral estimation ...
;  Note:  before compiling and loading nonparametric.lisp,
;         you should compile and load (in the order listed)
;            sapa-package.lisp, utilities.lisp, basic-math.lisp,
;            basic-statistics.lisp, dft-and-fft.lisp, tapers.lisp,
;            filtering.lisp, random.lisp and acvs.lisp
;
;  6/3/93
;
;  SAPA, Version 1.0; Copyright 1993, Donald B. Percival, All Rights Reserved
;
;  Use and copying of this software and preparation of derivative works
;  based upon this software are permitted.  Any distribution of this
;  software or derivative works must comply with all applicable United
;  States export control laws.
; 
;  This software is made available AS IS, and no warranty -- about the
;  software, its performance, or its conformity to any
;  specification -- is given or implied. 
;
;  Comments about this software can be addressed to dbp@apl.washington.edu
;-------------------------------------------------------------------------------
;;; (compile-file "ccl:SAPA;nonparametric.lisp")
;;; (load "ccl:SAPA;nonparametric.fasl")
;-------------------------------------------------------------------------------
(if (not (find-package :SAPA))
X  (defpackage :SAPA (:USE :COMMON-LISP)))
X
(in-package :SAPA)
X
(export '(;;; functions to compute nonparametric spectral estimates ...
X          periodogram
X          direct-spectral-estimate
X          lag-window-spectral-estimate
X          wosa-spectral-estimate
X
X          ;;; function to correct a spectral estimate for prewhitening ...
X          postcolor-spectral-estimate
X
X          ;;; functions to compute spectral and smoothing windows ...
X          spectral-window-for-direct-spectral-estimate
X          spectral-window-for-lag-window-spectral-estimate
X          smoothing-window-for-lag-window-spectral-estimate
X
X          ;;; functions to bandwidths, degrees of freedom, etc. ...
X          Grenander-smoothing-window-bandwidth
X          Jenkins-smoothing-window-bandwidth
X          equivalent-degrees-of-freedom
X          bandwidth&confidence-intervals-for-sdf-dB
X          
X          ;;; functions useful for objectively choosing lag window parameter ...
X          time-series-bandwidth
X          sample-cepstrum
X          cepstrum->I_m          
X
X          ;;; functions to compute the cumulative periodogram test statistic ...
X          cumulative-periodogram
X          quantile-of-Kolmogorov-test-statistic
X
X          ;;; functions to convert from one- to two-sided representation ...
X          one-sided-freq->two-sided-freq
X          one-sided-sdf->two-sided-sdf
X
X          ;;; functions concerning specific lag windows ...
X          bartlett-lag-window
X          bartlett-m->bandwidth
X          bartlett-bandwidth->m
X          bartlett-N-m->degrees-of-freedom
X          daniell-lag-window
X          daniell-m->bandwidth
X          daniell-bandwidth->m
X          daniell-N-m->degrees-of-freedom
X          parzen-lag-window
X          parzen-m->bandwidth
X          parzen-bandwidth->m
X          parzen-N-m->degrees-of-freedom
X          papoulis-lag-window
X          papoulis-m->bandwidth
X          papoulis-bandwidth->m
X          papoulis-N-m->degrees-of-freedom
X          ))
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  periodogram
;;;                 direct-spectral-estimate
;;;                 lag-window-spectral-estimate
;;;                 wosa-spectral-estimate
;;;  compute four common nonparametric spectral density estimates.
;;;  The first two functions take a time series as input and return
;;;  a direct spectral estimate (the periodogram is a special case of such 
;;;  an estimate).  The function lag-window-spectral-estimate takes an estimate
;;;  of the autocovariance sequence (acvs) and returns a lag window spectral
;;;  estimate.  Note that there is an option in direct-spectral-estimate for
;;;  returning the acvs estimate corresponding to a direct spectral estimate
;;;  (the acvs estimate corresponding to the periodogram is just the usual
;;;  biased estimator of the acvs, which can be computed using the Lisp
;;;  function acvs).  The function wosa-spectral-estimate computes
;;;  a nonparametric estimate using Welch's overlapped segment averaging (wosa)
;;;  technique.   Chapter 6 of the SAPA book is devoted to nonparametric
;;;  spectral estimation.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun periodogram
X       (time-series
X        &key
X        (center-data t)
X        (start 0)
X        (end (length time-series))
X        (N-nonzero-freqs :half-next-power-of-2)
X        (return-est-for-0-freq-p nil)
X        (sampling-time 1.0)
X        (scratch-dft (make-array (get-N-dft N-nonzero-freqs (- end start))))
X        (sdf-transformation #'convert-to-dB)
X        (result-sdf (make-array (get-N-freqs N-nonzero-freqs (- end start)
X                                             return-est-for-0-freq-p)))
X        (return-frequencies-p t)
X        (freq-transformation nil)
X        (result-freq (if return-frequencies-p
X                       (make-array (get-N-freqs N-nonzero-freqs (- end start)
X                                                return-est-for-0-freq-p)))))
X  "given
X   [1] time-series (required)
X       ==> a vector of time series values
X   [2] center-data (keyword; t)
X       ==> if t, subtract sample mean from time series;
X           if a number, subtracts that number from time series;
X           if nil, time-series is not centered
X   [3] start (keyword; 0)
X       ==> start index of vector to be used
X   [4] end (keyword; length of time-series)
X       ==> 1 + end index of vector to be used
X   [5] N-nonzero-freqs (keyword; :half-next-power-of-2)
X       ==> specifies at how many nonzero frequencies
X           periodogram is to be computed -- choices are:
X           :half-next-power-of-2
X            ==> 1/2 * power of two >= sample size;
X           :next-power-of-2
X            ==> power of two >= sample size;
X           :twice-next-power-of-2
X            ==> 2 * power of two >= sample size;
X           :Fourier
X            ==> just at Fourier frequencies
X            -- or --
X           any power of 2 >= 1/2 * [power of two >= sample size]
X   [6] return-est-for-0-freq-p (keyword; nil)
X       ==> if t, periodogram is computed at zero frequency;
X           otherwise, it is not computed.
X   [7] sampling-time (keyword; 1.0)
X       ==> sampling time (called delta t in the SAPA book)
X   [8] scratch-dft (keyword; vector of correct length)
X       ==> vector in which the in-place dft is done
X   [9] sdf-transformation (keyword; #'convert-to-dB)
X       ==> a function of one argument or nil;
X           if bound to a function, the function is used
X           to transform all elements of result-sdf
X  [10] result-sdf (keyword; vector of correct length)
X       <== vector into which periodogram is placed;
X           it must be exactly of the length dictated
X           by N-nonzero-freqs and return-est-for-0-freq-p
X  [11] return-frequencies-p (keyword; t)
X       ==> if t, the frequencies associated with the periodogram
X           are computed and returned in result-freq
X  [12] freq-transformation (keyword; nil)
X       ==> a function of one argument or nil;
X           if bound to a function, the function is used
X           to transform all elements of result-freq
X           (ignored unless return-frequencies-p is true)
X  [13] result-freq (keyword; nil or vector of correct length)
X       <== not used if return-frequencies-p nil; otherwise,
X           vector of length dictated by
X           N-nonzero-freqs and return-est-for-0-freq-p
X           into which the frequencies associated with the values
X           in result-sdf are placed
returns
X   [1] result-sdf, a vector holding
X       the properly transformed periodogram
X   [2] result-freq (if return-frequencies-p is t),
X       a vector holding the properly transformed
X       frequencies associated with values in result-sdf
X        -- or --
X       nil (if return-frequencies-p is nil)
X   [3] the length of the vector result-sdf
---
Note: see Section 6.3 of the SAPA book"
X  (let* ((sample-size (- end start))
X         (N-freqs (get-N-freqs N-nonzero-freqs sample-size
X                               return-est-for-0-freq-p))
X         (N-dft (get-N-dft N-nonzero-freqs sample-size))
X         (offset (if return-est-for-0-freq-p 0 1))
X         (fiddle-factor-sdf (/ sampling-time sample-size))
X         (fiddle-factor-freq (/ (* N-dft sampling-time))))
X    ;;; put centered series into first part of scratch-dft
X    (center&taper-time-series time-series
X                              :center-data center-data
X                              :start start
X                              :end end
X                              :result scratch-dft)
X    ;;; zero the rest of scratch-dft ...
X    (fill scratch-dft 0.0 :start sample-size :end N-dft)
X    (dft! scratch-dft :N N-dft)
X    (dotimes (i N-freqs)
X      (setf (aref result-sdf i)
X            (* fiddle-factor-sdf
X               (expt (abs (aref scratch-dft (+ i offset))) 2))))
X    (if sdf-transformation
X      (transform-a-sequence! sdf-transformation result-sdf))
X    (when return-frequencies-p
X      (dotimes (i N-freqs)
X        (setf (aref result-freq i)
X              (* fiddle-factor-freq (float (+ i offset)))))
X      (if freq-transformation
X        (transform-a-sequence! freq-transformation result-freq)))
X    (values result-sdf result-freq N-freqs)))
X
#|
;;; Here we test periodogram via the Parseval result stated
;;; in Exercise [6.6c] of the SAPA book ...
(dolist (N '(5 8 9 22 32 63 64 65) (values))
X  (let* ((time-series (x+b! (ranorms N) 100))
X         (sigma^2 (sample-variance time-series))
X         (sampling-time 0.25)
X         (the-periodogram (periodogram
X                           time-series
X                           :N-nonzero-freqs :Fourier
X                           :sdf-transformation nil
X                           :sampling-time sampling-time))
X         (N-f (length the-periodogram))
X         (Parseval-sum (if (evenp N)
X                         (/ (+ (* 2 (sum the-periodogram :end (1- N-f)))
X                               (svref the-periodogram (1- N-f)))
X                            (* N sampling-time))
X                         (/ (* 2 (sum the-periodogram))
X                            (* N sampling-time)))))
X    (format t "~& N = ~2D, N-f = ~2D: ~11,8F ~11,8F ~F"
X            N N-f sigma^2 Parseval-sum (/ Parseval-sum sigma^2))))
;==>
X N =  5, N-f =  2:  0.51511055  0.51511055 1.0000000000000004
X N =  8, N-f =  4:  0.68129255  0.68129255 0.9999999999999999
X N =  9, N-f =  4:  0.97793576  0.97793576 0.9999999999999989
X N = 22, N-f = 11:  0.92517636  0.92517636 1.000000000000006
X N = 32, N-f = 16:  0.74719322  0.74719322 0.999999999999999
X N = 63, N-f = 31:  0.96915207  0.96915207 1.0000000000000078
X N = 64, N-f = 32:  1.12571507  1.12571507 1.0000000000000024
X N = 65, N-f = 32:  0.64325949  0.64325949 1.000000000000028
;;; Note: the important thing here is that the last column of numbers
;;;       be close to unity -- since we are using random numbers,
;;;       these three columns of numbers will change each time
;;;       the above Lisp form is evaluated.
X
;;; For this example we use the 20 point time series of Section 6.16
;;; of the SAPA book:
(let ((20-pt-ts #(71.0  63.0  70.0  88.0  99.0  90.0 110.0 135.0 128.0 154.0
X                  156.0 141.0 131.0 132.0 141.0 104.0 136.0 146.0 124.0 129.0)))
X  (multiple-value-bind (the-periodogram freqs N-f)
X                       (periodogram 
X                        20-pt-ts
X                        :N-nonzero-freqs :Fourier
X                        :sampling-time 0.25)
X    (dotimes (i N-f)
X      (format t "~&~6,4F: ~8,4F"
X              (svref freqs i)
X              (svref the-periodogram i)))
X    (values)))
;==>
0.2000:  30.8586
0.4000:  24.9217
0.6000:  21.0531
0.8000:  19.1324
1.0000:  12.5527
1.2000:  18.2452
1.4000:  12.7906
1.6000:  19.9132
1.8000:  11.6759
2.0000:   5.0515
|#
X
;-------------------------------------------------------------------------------
(defun direct-spectral-estimate
X       (time-series
X        &key
X        (center-data t)
X        (start 0)
X        (end (length time-series))
X        (N-nonzero-freqs :half-next-power-of-2)
X        (return-est-for-0-freq-p nil)
X        (sampling-time 1.0)
X        (scratch-dft (make-array (get-N-dft N-nonzero-freqs (- end start))))
X        (data-taper nil)
X        (data-taper-parameters)
X        (recenter-after-tapering-p t)
X        (restore-power-option-p t)
X        (return-acvs-p nil)
X        (result-acvs (make-array (- end start)))
X        (sdf-transformation #'convert-to-dB)
X        (result-sdf (make-array (get-N-freqs N-nonzero-freqs (- end start)
X                                             return-est-for-0-freq-p)))
X        (return-frequencies-p t)
X        (freq-transformation nil)
X        (result-freq (if return-frequencies-p
X                       (make-array (get-N-freqs N-nonzero-freqs (- end start)
X                                                return-est-for-0-freq-p)))))
X  "given
X   [1] time-series (required)
X       ==> a vector of time series values
X   [2] center-data (keyword; t)
X       ==> if t, subtract sample mean from time series;
X           if a number, subtracts that number from time series;
X           if nil, time-series is not centered
X   [3] start (keyword; 0)
X       ==> start index of vector to be used
X   [4] end (keyword; length of time-series)
X       ==> 1 + end index of vector to be used
X   [5] N-nonzero-freqs (keyword; :half-next-power-of-2)
X       ==> specifies at how many nonzero frequencies
X           direct spectral estimate is to be computed -- choices are:
X           :half-next-power-of-2
X            ==> 1/2 * power of two >= sample size;
X           :next-power-of-2
X            ==> power of two >= sample size;
X           :twice-next-power-of-2
X            ==> 2 * power of two >= sample size;
X           :Fourier
X            ==> just at Fourier frequencies
X            -- or --
X           any power of 2 >= 1/2 * [power of two >= sample size]
X   [6] return-est-for-0-freq-p (keyword; nil)
X       ==> if t, sdf is computed at zero frequency;
X           otherwise, it is not computed.
X   [7] sampling-time (keyword; 1.0)
X       ==> sampling time (called delta t in the SAPA book)
X   [8] scratch-dft (keyword; vector of correct length)
X       ==> vector in which the in-place dft is done
X   [9] data-taper (keyword; nil)
X       ==> nil or a tapering function
X  [10] data-taper-parameters (keyword)
X       ==> parameters for tapering function (not used
X           if data-taper is nil)
X  [11] recenter-after-tapering-p (keyword; t)
X       ==> if t and data-taper is a function,
X           centers tapered series by subtracting
X           off its sample mean
X  [12] restore-power-option-p (keyword; t)
X       ==> if t and data-taper is a function,
X           normalizes tapered series to have same
X           sum of squares as before tapering
X  [13] return-acvs-p (keyword; nil)
X       ==> if t, computes acvs corresponding
X           to direct spectral estimate
X  [14] result-acvs (keyword; vector of correct length)
X       <== vector into which acvs values are placed
X           (not used if return-acvs-estimate-p is nil)
X  [15] sdf-transformation (keyword; #'convert-to-dB)
X       ==> a function of one argument or nil;
X           if bound to a function, the function is used
X           to transform all elements of result-sdf
X  [16] result-sdf (keyword; vector of correct length)
X       <== vector into which direct spectral estimates are placed;
X           it must be exactly of the length dictated
X           by N-nonzero-freqs and return-est-for-0-freq-p
X  [17] return-frequencies-p (keyword; t)
X       ==> if t, the frequencies associated with the spectral estimate
X           are computed and returned in result-freq
X  [18] freq-transformation (keyword; nil)
X       ==> a function of one argument or nil;
X           if bound to a function, the function is used
X           to transform all elements of result-freq
X           (ignored unless return-frequencies-p is true)
X  [19] result-freq (keyword; nil or vector of correct length)
X       <== not used if return-frequencies-p nil; otherwise,
X           vector of length dictated by
X           N-nonzero-freqs and return-est-for-0-freq-p
X           into which the frequencies associated with the values
X           in result-sdf are placed
returns
X   [1] result-sdf, a vector holding
X       the properly transformed sdf
X   [2] result-freq (if return-frequencies-p is t),
X       where result-freq is a vector holding
X       the properly transformed frequencies
X       associated with values in result-sdf
X        -- or --
X       nil (if return-frequencies-p is nil)
X   [3] the length of the vector result-sdf
X   [4] C_h, the variance inflation factor due to the data taper
X   [5] result-acvs (if return-acvs-p is t),
X       a vector holding the acvs corresponding
X       to the direct spectral estimate
X        -- or --
X       nil (if return-acvs-p is nil)
---
Note: see Section 6.3 of the SAPA book"
X  (let* ((sample-size (- end start))
X         (N-freqs (get-N-freqs N-nonzero-freqs sample-size
X                               return-est-for-0-freq-p))
X         (N-dft (get-N-dft N-nonzero-freqs sample-size))
X         (offset (if return-est-for-0-freq-p 0 1))
X         (fiddle-factor-sdf (/ sampling-time sample-size))
X         (fiddle-factor-freq (/ (* N-dft sampling-time))))
X    ;;; put centered & tapered series into first part of scratch-dft
X    (multiple-value-bind (junk more-junk C_h)
X                         (center&taper-time-series
X                          time-series
X                          :center-data center-data
X                          :start start
X                          :end end
X                          :data-taper data-taper
X                          :data-taper-parameters data-taper-parameters
X                          :recenter-after-tapering-p recenter-after-tapering-p
X                          :restore-power-option-p restore-power-option-p
X                          :result scratch-dft)
X      (declare (ignore junk more-junk))
X      (if return-acvs-p
X        (acvs scratch-dft
X              :end sample-size
X              :center-data-p nil
X              :result result-acvs))
X      ;;; zero the rest of scratch-dft ...
X      (fill scratch-dft 0.0 :start sample-size :end N-dft)
X      (dft! scratch-dft :N N-dft)
X      (dotimes (i N-freqs)
X        (setf (aref result-sdf i)
X              (* fiddle-factor-sdf
X                 (expt (abs (aref scratch-dft (+ i offset))) 2))))
X      (if sdf-transformation
X        (transform-a-sequence! sdf-transformation result-sdf))
X      (when return-frequencies-p
X        (dotimes (i N-freqs)
X          (setf (aref result-freq i)
X                (* fiddle-factor-freq (float (+ i offset)))))
X        (if freq-transformation
X          (transform-a-sequence! freq-transformation result-freq)))
X      (values result-sdf result-freq N-freqs C_h result-acvs))))
X
#|
;;; Here we test the direct spectral estimate using the Parseval result
;;; analogous to the one we used to check the periodogram ... 
(dolist (N '(5 8 9 22 32 63 64 65) (values))
X  (let* ((time-series (x+b! (ranorms N) 100))
X         (sigma^2 (sample-variance time-series))
X         (sampling-time 0.25)
X         (the-sdf-est (direct-spectral-estimate
X                       time-series
X                       :data-taper #'dpss-data-taper!
X                       :data-taper-parameters 4.0
X                       :N-nonzero-freqs :Fourier
X                       :sdf-transformation nil
X                       :sampling-time sampling-time))
X         (N-f (length the-sdf-est))
X         (Parseval-sum (if (evenp N)
X                         (/ (+ (* 2 (sum the-sdf-est  :end (1- N-f)))
X                               (svref the-sdf-est  (1- N-f)))
X                            (* N sampling-time))
X                         (/ (* 2 (sum the-sdf-est ))
X                            (* N sampling-time)))))
X    (format t "~& N = ~2D, N-f = ~2D: ~11,8F ~11,8F ~F"
X            N N-f sigma^2 Parseval-sum (/ Parseval-sum sigma^2))))
;==>
X N =  5, N-f =  2:  0.42023448  0.42023448 1.0000000000000004
X N =  8, N-f =  4:  0.66547264  0.66547264 1.0000000000000004
X N =  9, N-f =  4:  0.46234261  0.46234261 0.9999999999999981
X N = 22, N-f = 11:  0.58047776  0.58047776 1.0000000000000089
X N = 32, N-f = 16:  1.46881228  1.46881228 0.9999999999999989
X N = 63, N-f = 31:  1.03414687  1.03414687 1.0000000000000056
X N = 64, N-f = 32:  0.53556515  0.53556515 1.000000000000004
X N = 65, N-f = 32:  1.12152795  1.12152795 1.0000000000000302
;;; Note: again, the important thing is that the last column of numbers
;;;       be close to unity -- since we are using random numbers,
;;;       these three columns of numbers will change each time
;;;       the above Lisp form is evaluated.
|#
X
;-------------------------------------------------------------------------------
(defun lag-window-spectral-estimate
X       (acvs
X        lag-window-function
X        &key
X        (max-lag (1- (length acvs)))
X        (N-ts (length acvs))
X        (N-nonzero-freqs :half-next-power-of-2)
X        (return-est-for-0-freq-p nil)
X        (sampling-time 1.0)
X        (scratch-dft (make-array (get-N-dft N-nonzero-freqs (1+ max-lag))))
X        (C_h 1.0)
X        (sdf-transformation #'convert-to-dB)
X        (result-sdf (make-array (get-N-freqs N-nonzero-freqs (1+ max-lag)
X                                             return-est-for-0-freq-p)))
X        (return-frequencies-p t)
X        (freq-transformation nil)
X        (result-freq (if return-frequencies-p
X                       (make-array (get-N-freqs N-nonzero-freqs (1+ max-lag)
X                                                return-est-for-0-freq-p)))))
X  "given
X   [1] acvs (required)
X       ==> vector containing autocovariance sequence
X   [2] lag-window-function (required)
X       ==> function of a single variable that computes the value
X           of the lag window for a given lag
X   [3] max-lag (keyword; (1- (length acvs)))
X       ==> maximum lag in acvs to be used
X   [4] N-ts (keyword; length of acvs)
X       ==> length of the time series from which acvs was constructed;
X           this is needed to compute equivalent degrees of freedom
X   [5] N-nonzero-freqs (keyword; :half-next-power-of-2)
X       ==> specifies at how many nonzero frequencies
X           direct spectral estimate is to be computed -- choices are:
X           :half-next-power-of-2
X            ==> 1/2 * power of two >= sample size;
X           :next-power-of-2
X            ==> power of two >= sample size;
X           :twice-next-power-of-2
X            ==> 2 * power of two >= sample size;
X           :Fourier
X            ==> just at Fourier frequencies
X            -- or --
X           any power of 2 >= 1/2 * [power of two >= sample size]
X   [6] return-est-for-0-freq-p (keyword; nil)
X       ==> if t, sdf is computed at zero frequency;
X           otherwise, it is not computed.
X   [7] sampling-time (keyword; 1.0)
X       ==> sampling time (called delta t in the SAPA book)
X   [8] scratch-dft (keyword; vector of correct length)
X       ==> vector in which the in-place dft is done
X   [9] C_h (keyword; 1.0)
X       ==> variance inflation factor due to tapering
X  [10] sdf-transformation (keyword; #'convert-to-dB)
X       ==> a function of one argument or nil;
X           if bound to a function, the function is used
X           to transform all elements of result-sdf
X  [11] result-sdf (keyword; vector of correct length)
X       <== vector into which lag window spectral estimates are placed;
X           it must be exactly of the length dictated
X           by N-nonzero-freqs and return-est-for-0-freq-p
X  [12] return-frequencies-p (keyword; t)
X       ==> if t, the frequencies associated with the spectral estimate
X           are computed and returned in result-freq
X  [13] freq-transformation (keyword; nil)
X       ==> a function of one argument or nil;
X           if bound to a function, the function is used
X           to transform all elements of result-freq
X           (ignored unless return-frequencies-p is true)
X  [14] result-freq (keyword; nil or vector of correct length)
X       <== not used if return-frequencies-p nil; otherwise,
X           vector of length dictated by
X           N-nonzero-freqs and return-est-for-0-freq-p
X           into which the frequencies associated with the values
X           in result-sdf are placed
returns
X   [1] result-sdf, a vector holding
X       the properly transformed sdf
X   [2] result-freq (if return-frequencies-p is t),
X       where result-freq is a vector holding
X       the properly transformed frequencies
X       associated with values in  result-sdf
X        -- or --
X       nil (if return-frequencies-p is nil)
X   [3] the length of the vector result-sdf
X   [4] the equivalent degrees of freedom
X   [5] the smoothing window bandwidth
---
Note: see Section 6.7 of the SAPA book"
X  ;;; Note: in what follows, we assume that the lag window
X  ;;;       at lag 0 is unity (see Equation (240a) of the SAPA book)
X  (let* ((N-acvs (1+ max-lag))
X         (N-freqs (get-N-freqs N-nonzero-freqs N-acvs
X                               return-est-for-0-freq-p))
X         (N-dft (get-N-dft N-nonzero-freqs N-acvs))
X         (offset (if return-est-for-0-freq-p 0 1))
X         (fiddle-factor-freq (/ (* N-dft sampling-time)))
X         (B_W-bot 1.0))
X    (cond
X     ((<= N-dft (* 2 max-lag))
X      (setf (aref scratch-dft 0) (aref acvs 0))
X      (let ((tau 1))
X        (dotimes (i max-lag)
X          (let ((lag-window-value (funcall lag-window-function tau)))
X            (incf B_W-bot (* 2 (expt lag-window-value 2)))
X            (setf (aref scratch-dft tau)
X                  (* 2 lag-window-value (aref acvs tau))))
X          (incf tau)))
X      (fill scratch-dft 0.0 :start max-lag :end N-dft))
X     (t
X      (let ((tau 1)
X            (tau-backward (1- N-dft)))
X        (dotimes (i max-lag)
X          (let ((lag-window-value (funcall lag-window-function tau)))
X            (incf B_W-bot (* 2 (expt lag-window-value 2)))
X            (setf (aref scratch-dft tau)
X                  (setf (aref scratch-dft tau-backward)
X                        (* lag-window-value (aref acvs tau))))
X            (incf tau)
X            (decf tau-backward)))
X        (fill scratch-dft 0.0 :start tau :end (1+ tau-backward)))))
X    (dft! scratch-dft :N N-dft)
X    (dotimes (i N-freqs)
X      (setf (aref result-sdf i)
X            (* sampling-time
X               (realpart (aref scratch-dft (+ i offset))))))
X    (if sdf-transformation
X      (transform-a-sequence! sdf-transformation result-sdf))
X    (when return-frequencies-p
X      (dotimes (i N-freqs)
X        (setf (aref result-freq i)
X              (* fiddle-factor-freq (float (+ i offset)))))
X      (if freq-transformation
X        (transform-a-sequence! freq-transformation result-freq)))
X    (let ((B_W (/ (* sampling-time B_W-bot))))
X      (values result-sdf
X              result-freq
X              N-freqs
X              (/ (* 2 N-ts B_W sampling-time) C_h)  ;equivalent dof
X              B_W                                   ;smoothing window bandwidth
X              ))))
X
#|
;;; For this example we use the 20 point time series of Section 6.16
;;; of the SAPA book:
(let* ((20-pt-ts #(71.0  63.0  70.0  88.0  99.0  90.0 110.0 135.0 128.0 154.0
X                   156.0 141.0 131.0 132.0 141.0 104.0 136.0 146.0 124.0 129.0))
X       (sampling-time 0.25)
X       (the-acvs (acvs 20-pt-ts))
X       (Parzen-15-sdf-est (lag-window-spectral-estimate
X                           the-acvs 
X                           #'(lambda (tau)
X                               (parzen-lag-window
X                                tau 15))
X                           :N-nonzero-freqs :Fourier
X                           :sampling-time sampling-time))
X       (Parzen-10-sdf-est (lag-window-spectral-estimate
X                           the-acvs 
X                           #'(lambda (tau)
X                               (parzen-lag-window
X                                tau 10))
X                           :N-nonzero-freqs :Fourier
X                           :sampling-time sampling-time)))
X  (multiple-value-bind (Parzen-5-sdf-est freqs N-f)
X                       (lag-window-spectral-estimate
X                        the-acvs 
X                        #'(lambda (tau)
X                            (parzen-lag-window
X                             tau 5))
X                        :N-nonzero-freqs :Fourier
X                        :sampling-time sampling-time)
X    (dotimes (i N-f)
X      (format t "~&~6,4F: ~8,4F ~8,4F ~8,4F"
X              (svref freqs i)
X              (svref Parzen-5-sdf-est i)
X              (svref Parzen-10-sdf-est i)
X              (svref Parzen-15-sdf-est i)))
X    (values)))
;==>
0.2000:  27.0227  28.2599  28.8918
0.4000:  26.0828  26.1869  25.8978
0.6000:  24.5508  22.3966  20.3750
0.8000:  22.5325  18.8636  18.0633
1.0000:  20.2961  17.8869  17.7836
1.2000:  18.2975  17.5004  17.8190
1.4000:  16.8731  16.8297  16.3072
1.6000:  15.9210  15.8348  16.4136
1.8000:  15.2629  14.0008  13.5807
2.0000:  15.0010  12.3206   9.7972
|#
X
;-------------------------------------------------------------------------------
(defun wosa-spectral-estimate
X       (time-series
X        block-size
X        &key
X        (proportion-of-overlap 0.5)
X        (oversampling-factor 1)
X        (center-data t)
X        (start 0)
X        (end (length time-series))
X        (return-est-for-0-freq-p t)
X        (sampling-time 1.0)
X        (scratch-dft (make-array (* oversampling-factor block-size)))
X        (data-taper #'Hanning-data-taper!)
X        (data-taper-parameters nil)
X        (restore-power-option-p t)
X        (sdf-transformation #'convert-to-dB)
X        (result-sdf (make-array (wosa-get-N-freqs block-size
X                                                  oversampling-factor
X                                                  return-est-for-0-freq-p)))
X        (return-sdf-estimates-for-each-block-p t)
X        (return-frequencies-p t)
X        (freq-transformation nil)
X        (result-freq (if return-frequencies-p
X                       (make-array (wosa-get-N-freqs block-size
X                                                     oversampling-factor
X                                                     return-est-for-0-freq-p)))))
X  "given
X   [1] time-series (required)
X       ==> a vector of real-valued numbers
X   [2] block-size (required)
X       ==> a power of two
X   [3] proportion-of-overlap (keyword; 0.5)
X       ==> number greater than 0 and less than 1
X   [4] oversampling-factor (keyword; 1)
X       ==> a factor that controls the number of frequencies
X           at which the wosa spectral estimate is computed;
X           this factor should be an integer power of two
X           such as 1, 2, 4, etc; for example,
X           1 yields Fourier frequencies for block-size;
X           2 yields grid twice as fine as Fourier frequencies;
X           4 yields griid 4 times as fine as Fourier frequencies;
X           etc.
X   [5] center-data (keyword; t)
X       ==> if t, subtract sample mean from time series;
X           if a number, subtracts that number from time series;
X           if nil, time-series is not centered
X   [6] start (keyword; 0)
X       ==> start index of vector to be used
X   [7] end (keyword; length of time-series)
X       ==> 1 + end index of vector to be used
X   [8] return-est-for-0-freq-p (keyword; nil)
X       ==> if t, sdf is computed at zero frequency;
X           otherwise, it is not computed.
X   [9] sampling-time (keyword; 1.0)
X       ==> sampling time (called delta t in the SAPA book)
X  [10] scratch-dft (keyword; vector of correct length)
X       ==> vector in which the in-place dft is done
X  [11] data-taper (keyword; #'Hanning-data-taper!)
X       ==> a tapering function or nil
X  [12] data-taper-parameters (keyword; nil)
X       ==> parameters for tapering function (not used
X           if data-taper is nil); the default of nil
X           is appropriate for the Hanning data taper
X           because it does not have any parameters
X  [13] restore-power-option-p (keyword; t)
X       ==> if t and data-taper is non-nil,
X           normalizes tapered series to have same
X           sum of squares as before tapering
X  [14] sdf-transformation (keyword; #'convert-to-dB)
X       ==> a function of one argument or nil;
X           if bound to a function, the function is used
X           to transform all elements of result-sdf
X  [15] result-sdf (keyword; vector of correct length)
X       <== vector into which wosa sdf estimate is placed;
X           it must be EXACTLY of the length dictated
X           by block-size, oversampling-factor and return-est-for-0-freq-p
X  [16] return-sdf-estimates-for-each-block-p (keyword; t)
X       ==> if t, individual spectra for each block are returned in a list;
X           note that these spectra are untransformed
X           (i.e., the option sdf-transformation applies only
X           to the final wosa estimate)
X  [17] return-frequencies-p (keyword; t)
X       ==> if t, the frequencies associated with the spectral estimate
X           are computed and returned in result-freq
X  [18] freq-transformation (keyword; nil)
X       ==> a function of one argument or nil;
X           if bound to a function, the function is used
X           to transform all elements of result-freq
X           (ignored unless return-frequencies-p is true)
X  [19] result-freq (keyword; nil or vector of correct length)
X       <== not used if return-frequencies-p nil; otherwise,
X           vector of length dictated by
X           block-size, oversampling-factor and return-est-for-0-freq-p
X           into which the frequencies associated with the values
X           in result-sdf are placed
returns
X   [1] wosa spectral estimate
X   [2] associated frequencies
X   [3] equivalent degrees of freedom
X   [4] list of individual direct spectral estimates"
X  (let* ((centered-time-series (center&taper-time-series
X                                time-series
X                                :center-data center-data
X                                :start start
X                                :end end))
X         (N-dft (* oversampling-factor block-size))
X         (N-freqs (wosa-get-N-freqs block-size
X                                    oversampling-factor
X                                    return-est-for-0-freq-p))
X         (offset-freq (if return-est-for-0-freq-p 0 1))
X         (sample-size (- end start))
X         (number-of-blocks (calculate-number-of-blocks
X                            sample-size block-size
X                            proportion-of-overlap))
X         (list-of-individual-sdfs '())
X         (vector-with-data-taper (if data-taper
X                                   (center&taper-time-series
X                                    (make-array block-size
X                                                :initial-element 1.0)
X                                    :center-data nil
X                                    :data-taper data-taper
X                                    :data-taper-parameters data-taper-parameters
X                                    :recenter-after-tapering-p nil
X                                    :restore-power-option-p nil)
X                                   (make-array block-size
X                                               :initial-element 1.0)))
X         (fiddle-factor-sdf (/ sampling-time block-size))
X         (fiddle-factor-freq (/ (* N-dft sampling-time)))
X         offset-block)
X    (fill result-sdf 0.0)
X    (dotimes (k number-of-blocks)
X      (setf offset-block
X            (get-offset-to-kth-block sample-size block-size number-of-blocks k))
X      (replace scratch-dft
X               centered-time-series
X               :start2 offset-block
X               :end2 (+ offset-block block-size))
X      (fill scratch-dft 0.0 :start block-size :end N-dft)
X      (when data-taper
X        (let ((sum-of-squares-before (if restore-power-option-p
X                                       (sum-of-squares
X                                        scratch-dft
X                                        :end block-size))))
X          (map-into scratch-dft #'* scratch-dft vector-with-data-taper)
X          (if restore-power-option-p
X            (let ((mult-factor (sqrt (/ sum-of-squares-before
X                                        (sum-of-squares
X                                         scratch-dft
X                                         :end block-size)))))
X              (dotimes (i block-size)
X                (multf (aref scratch-dft i) mult-factor))))))
X      (dft! scratch-dft :N N-dft)
X      (dotimes (i N-freqs)
X        (setf (aref scratch-dft i)
X              (* fiddle-factor-sdf
X                 (expt (abs (aref scratch-dft (+ i offset-freq))) 2))))
X      (if return-sdf-estimates-for-each-block-p
X        (push (subseq scratch-dft 0 N-freqs) list-of-individual-sdfs))
X      ;;; Note: here we make use of the assumption that result-sdf is
X      ;;;       of exactly the correct size
X      (map-into result-sdf #'+ result-sdf scratch-dft))
X    (a*x! (/ number-of-blocks) result-sdf)
X    (if sdf-transformation
X      (transform-a-sequence! sdf-transformation result-sdf))
X    (when return-frequencies-p
X      (dotimes (i N-freqs)
X        (setf (aref result-freq i)
X              (* fiddle-factor-freq (float (+ i offset-freq)))))
X      (if freq-transformation
X        (transform-a-sequence! freq-transformation result-freq)))
X    (values result-sdf
X            result-freq
X            (equivalent-dof-for-wosa sample-size
X                                     block-size
X                                     number-of-blocks
X                                     vector-with-data-taper)
X            (reverse list-of-individual-sdfs)
X            )))
X
#|
;;; For this example we use the 20 point time series of Section 6.16
;;; of the SAPA book:
(let* ((20-pt-ts #(71.0  63.0  70.0  88.0  99.0  90.0 110.0 135.0 128.0 154.0
X                   156.0 141.0 131.0 132.0 141.0 104.0 136.0 146.0 124.0 129.0))
X       (sampling-time 0.25)
X       (wosa-4 (wosa-spectral-estimate
X                20-pt-ts
X                4
X                :oversampling-factor 4
X                :sampling-time sampling-time)))
X  (multiple-value-bind (wosa-8 freqs)
X                       (wosa-spectral-estimate
X                        20-pt-ts
X                        8
X                        :oversampling-factor 2
X                        :sampling-time sampling-time)
X    (dotimes (i (length freqs))
X      (format t "~&~6,4F: ~8,4F ~8,4F"
X              (svref freqs i)
X              (svref wosa-4 i)
X              (svref wosa-8 i))))
X  (values))
;==>
0.0000:  26.9058  28.3938
0.2500:  26.4416  27.3650
0.5000:  25.0489  23.7518
0.7500:  22.7732  17.1140
1.0000:  19.9187  16.5887
1.2500:  17.4236  17.5460
1.5000:  16.1819  16.2470
1.7500:  15.7165  14.1435
2.0000:  15.5530  12.5478
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The function   postcolor-spectral-estimate
;;;  takes a spectral estimate for a prewhitened time series and corrects
;;;  it for the effect of prewhitening.  One caveat here is that the spectral
;;;  estimate is assumed to be unconverted (e.g., not expressed in decibel
;;;  units), so care must be taken to insure that this is true by setting
;;;  the sdf-transformation keyword to nil in whatever function is used
;;;  to compute the spectral estimate for the prewhitened time series.
;;;  Sections 6.5 and 9.10 of the SAPA book discuss prewhitening and
;;;  postcoloring.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun postcolor-spectral-estimate
X       (precolored-sdf-estimate
X        prewhitening-filter
X        sample-size
X        &key
X        (N-nonzero-freqs :half-next-power-of-2)
X        (includes-est-for-0-freq-p nil)
X        (scratch-dft (make-array (get-N-dft N-nonzero-freqs sample-size)))
X        (sdf-transformation #'convert-to-dB)
X        (result-sdf (make-array (length precolored-sdf-estimate))))
X  "given
X   [1] precolored-sdf-estimate (required)
X       ==> vector containing sdf estimate to be postcolored;
X           note that this estimate is assumed to be untransformed
X           (i.e., not expressed in dB, etc)
X   [2] prewhitening-filter (required)
X       ==> vector with coefficients of prewhitening filter
X   [3] sample-size (required)
X       ==> length of the time series from which precolored-sdf-estimate
X           was constructed; this is needed to get corresponding
X           grid of frequencies
X   [4] N-nonzero-freqs (keyword; :half-next-power-of-2)
X       ==> specifies at how many nonzero frequencies
X           precolored-sdf-estimate was computed -- choices are:
X           :half-next-power-of-2
X            ==> 1/2 * power of two >= sample size;
X           :next-power-of-2
X            ==> power of two >= sample size;
X           :twice-next-power-of-2
X            ==> 2 * power of two >= sample size;
X           :Fourier
X            ==> just at Fourier frequencies
X            -- or --
X           any power of 2 >= 1/2 * [power of two >= sample size]
X   [5] includes-est-for-0-freq-p (keyword; nil)
X       ==> if t, first element of precolored-sdf-estimate
X           corresponds to zero frequency
X   [6] scratch-dft (keyword; vector of correct length)
X       ==> vector in which the in-place dft is done
X   [7] sdf-transformation (keyword; #'convert-to-dB)
X       ==> a function of one argument or nil;
X           if bound to a function, the function is used
X           to transform all elements of result-sdf
X   [8] result-sdf (keyword; vector of correct length)
X       <== vector into which postcolored sdf estimate is placed;
X           it must be exactly the same length as precolored-sdf-estimate
returns
X   [1] result-sdf, a vector holding
X       the postcolored sdf estimate
---
Note: see Equation (438) of the SAPA book"
X  (let ((N-freqs (get-N-freqs N-nonzero-freqs sample-size
X                              includes-est-for-0-freq-p))
X        (N-dft (get-N-dft N-nonzero-freqs sample-size))
X        (offset (if includes-est-for-0-freq-p 0 1)))
X    ;;; put prewhitening filter into first part of scratch-dft
X    (copy-vector prewhitening-filter scratch-dft)
X    ;;; zero the rest of scratch-dft ...
X    (fill scratch-dft 0.0 :start (length prewhitening-filter) :end N-dft)
X    (dft! scratch-dft :N N-dft)
X    (dotimes (i N-freqs)
X      (setf (aref result-sdf i)
X            (/ (aref precolored-sdf-estimate i)
X               (expt (abs (aref scratch-dft (+ i offset))) 2))))
X    (if sdf-transformation
X      (transform-a-sequence! sdf-transformation result-sdf))
X    (values result-sdf)))
X
#|
;;; For this example we use the 20 point time series of Section 6.16
;;; of the SAPA book:
(let* ((20-pt-ts #(71.0  63.0  70.0  88.0  99.0  90.0 110.0 135.0 128.0 154.0
X                   156.0 141.0 131.0 132.0 141.0 104.0 136.0 146.0 124.0 129.0))
X       (N-ts (length 20-pt-ts))
X       (sampling-time 0.25)
X       (pgram-20-pt-ts-dB (periodogram 20-pt-ts :sampling-time sampling-time)))
X  (multiple-value-bind (pgram-differenced-20-pt-ts freqs N-f)
X                       (periodogram (difference 20-pt-ts)
X                                    :sampling-time sampling-time
X                                    :sdf-transformation nil)
X    (let ((pgram-differenced-20-pt-ts-dB
X           (postcolor-spectral-estimate pgram-differenced-20-pt-ts
X                                        #(-1 1)
X                                        (1- N-ts))))
X      (dotimes (i N-f)
X        (format t "~&~6,4F: ~8,4F ~8,4F"
X                (svref freqs i)
X                (svref pgram-20-pt-ts-dB i)
X                (svref pgram-differenced-20-pt-ts-dB i)))
X      (values))))
;==>
0.1250:  31.2344  30.7002
0.2500:  29.0872  26.2775
0.3750:  26.0741  20.3138
0.5000:  -8.9025  14.7938
0.6250:  21.0714  15.2171
0.7500:  16.4107   4.0106
0.8750:  18.5221  15.2630
1.0000:  12.5527   6.9111
1.1250:  20.5795  20.3538
1.2500:  16.3219  20.0302
1.3750:  15.2798   7.3147
1.5000:  14.1776  14.0639
1.6250:  19.3791  16.7867
1.7500:   8.9029  13.9247
1.8750:  10.9060   2.5244
2.0000:   5.0515   4.4350
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  spectral-window-for-direct-spectral-estimate
;;;                 spectral-window-for-lag-window-spectral-estimate
;;;                 smoothing-window-for-lag-window-spectral-estimate
;;;  compute the spectral windows associated with direct and lag window
;;;  spectral estimates and the smoothing window associated with
;;;  a lag window spectral estimate.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun spectral-window-for-direct-spectral-estimate
X       (sample-size
X        &key
X        (N-nonzero-freqs :twice-next-power-of-2)
X        (sampling-time 1.0)
X        (scratch-dft (make-array (get-N-dft N-nonzero-freqs sample-size)))
X        (data-taper nil)
X        (data-taper-parameters)
X        (spec-wind-transformation #'careful-convert-to-dB)
X        (result-spec-wind
X         (make-array (get-N-freqs N-nonzero-freqs sample-size t)))
X        (return-frequencies-p t)
X        (freq-transformation nil)
X        (result-freq
X         (if return-frequencies-p
X           (make-array (get-N-freqs N-nonzero-freqs sample-size t)))))
X  "given
X   [1] sample-size (required)
X       ==> sample size for which spectral window is to be computed
X   [2] N-nonzero-freqs (keyword; :twice-next-power-of-2)
X       ==> specifies at how many nonzero frequencies
X           spectral window is to be computed -- choices are:
X           :half-next-power-of-2
X            ==> 1/2 * power of two >= sample size;
X           :next-power-of-2
X            ==> power of two >= sample size;
X           :twice-next-power-of-2
X            ==> 2 * power of two >= sample size;
X           :Fourier
X            ==> just at Fourier frequencies
X            -- or --
X           any power of 2 >= 1/2 * [power of two >= sample size]
X   [3] sampling-time (keyword; 1.0)
X       ==> sampling time (called delta t in the SAPA book)
X   [4] scratch-dft (keyword; vector of correct length)
X       ==> vector in which the in-place dft is done
X   [5] data-taper (keyword; nil)
X       ==> nil or a tapering function
X   [6] data-taper-parameters (keyword)
X       ==> parameters for tapering function (not used
X           if data-taper is nil)
X   [7] spec-wind-transformation (keyword; #'convert-to-dB)
X       ==> a function of one argument or nil;
X           if bound to a function, the function is used
X           to transform all elements of result-spec-wind
X   [8] result-spec-wind (keyword; vector of correct length)
X       <== vector into which spectral window values are placed;
X           it must be exactly of the length dictated by N-nonzero-freqs
X   [9] return-frequencies-p (keyword; t)
X       ==> if t, the frequencies associated with the spectral window
X           are computed and returned in result-freq
X  [10] freq-transformation (keyword; nil)
X       ==> a function of one argument or nil;
X           if bound to a function, the function is used
X           to transform all elements of result-freq
X  [11] result-freq (keyword; nil or vector of correct length)
X       <== not used if return-frequencies-p nil; otherwise,
X           vector of N-nonzero-freqs +1 into which
X           the frequencies associated with the values
X           in result-spec-wind are placed
returns
X   [1] result-spec-wind, a vector holding
X       the properly transformed spectral window
X   [2] result-freq (if return-frequencies-p is t),
X       where result-freq is a vector holding
X       the properly transformed frequencies
X       associated with values in result-spec-wind
X        -- or --
X       nil (if return-frequencies-p is nil)
X   [3] the length of the vector result-spec-wind
---
Note: see equation below Equation (207a) of the SAPA book"
X  (let* ((N-freqs (get-N-freqs N-nonzero-freqs sample-size t))
X         (N-dft (get-N-dft N-nonzero-freqs sample-size))
X         (fiddle-factor-freq (/ (* N-dft sampling-time))))
X    ;;; put data taper into first part of scratch-dft
X    (fill scratch-dft (float (/ (sqrt sample-size))) :end sample-size)
X    (if data-taper
X      (funcall data-taper (make-array
X                           sample-size
X                           :displaced-to scratch-dft)
X               :taper-parameter data-taper-parameters
X               :normalization :N))
X    ;;; zero the rest of scratch-dft ...
X    (fill scratch-dft 0.0 :start sample-size :end N-dft)
X    (dft! scratch-dft :N N-dft)
X    ;;; Note: here is where we use the assumption that
X    ;;;       result-spec-wind is exactly of the size
X    ;;;       N-nonzero-freqs +1
X    ;;; Note: Allegro Common Lisp does not properly handle
X    ;;;       map-into, and it will not allow a redefinition
X    ;;;       of that function, so we must resort to this hack
X    ;;;       (for details, see hacks.lisp).
X    #-allegro
X    (map-into result-spec-wind
X                   (if spec-wind-transformation
X                     #'(lambda (x)
X                         (funcall spec-wind-transformation
X                                  (* sampling-time
X                                     (expt (abs x) 2))))
X                     #'(lambda (x)
X                         (* sampling-time
X                            (expt (abs x) 2))))
X                   scratch-dft)
X    #+allegro
X    (sapa-map-into result-spec-wind
X                   (if spec-wind-transformation
X                     #'(lambda (x)
X                         (funcall spec-wind-transformation
X                                  (* sampling-time
X                                     (expt (abs x) 2))))
X                     #'(lambda (x)
X                         (* sampling-time
X                            (expt (abs x) 2))))
X                   scratch-dft)
X    (when return-frequencies-p
X      (map-into result-freq
X                (if freq-transformation
X                  #'(lambda (k)
X                      (funcall freq-transformation
X                               (* k fiddle-factor-freq)))
X                  #'(lambda (k)
X                      (* k fiddle-factor-freq)))
X                (iota 0 (1- N-freqs))))
X    (values result-spec-wind result-freq N-freqs)))
X
#|
;;; yields Fejer's kernel (in dB) for N = 4
;;; (cf. top plot of Figure 200 of the SAPA book)
(spectral-window-for-direct-spectral-estimate 4)
;==> #(6.020599913279622 5.16438561858617 2.322606875058723 -3.9256792411212302
X      -100.0 -7.427827503856934 -5.332906831698538 -8.862378656807074 -100.0)
;    #(0.0 0.0625 0.125 0.1875 0.25 0.3125 0.375 0.4375 0.5)
;    9
X
;;; spectral window (in dB) for N = 4, NW=1 dpss data taper
(spectral-window-for-direct-spectral-estimate
X 4
X :data-taper #'dpss-data-taper!
X :data-taper-parameters  1.0)
;==> #(5.774834924698505 5.089606295677452 2.904268902242502
X      -1.3159448282420236 -9.584646039813563 -25.51680903895323
X     -13.166120020034743 -15.171623681945857 -100.0)
;    #(0.0 0.0625 0.125 0.1875 0.25 0.3125 0.375 0.4375 0.5)
;    9
|#
X
;-------------------------------------------------------------------------------
(defun spectral-window-for-lag-window-spectral-estimate
X       (sample-size
X        lag-window-function
X        &key
X        (N-nonzero-freqs :twice-next-power-of-2)
X        (sampling-time 1.0)
X        (scratch-dft (make-array (get-N-dft N-nonzero-freqs sample-size)))
X        (data-taper nil)
X        (data-taper-parameters)
X        (spec-wind-transformation #'careful-convert-to-dB)
X        (result-spec-wind
X         (make-array (get-N-freqs N-nonzero-freqs sample-size
X                                  t)))
X        (return-frequencies-p t)
X        (freq-transformation nil)
X        (result-freq
X         (if return-frequencies-p
X           (make-array (get-N-freqs N-nonzero-freqs sample-size t)))))
X  "given
X   [1] sample-size (required)
X       ==> sample size for which spectral window is to be computed
X   [2] lag-window-function (required)
X       ==> function of a single variable that computes the value
X           of the lag window for a given lag
X   [3] N-nonzero-freqs (keyword; :twice-next-power-of-2)
X       ==> specifies at how many nonzero frequencies
X           spectral window is to be computed -- choices are:
X           :half-next-power-of-2
X            ==> 1/2 * power of two >= sample size;
X           :next-power-of-2
X            ==> power of two >= sample size;
X           :twice-next-power-of-2
X            ==> 2 * power of two >= sample size;
X           :Fourier
X            ==> just at Fourier frequencies
X            -- or --
X           any power of 2 >= 1/2 * [power of two >= sample size]
X   [4] sampling-time (keyword; 1.0)
X       ==> sampling time (called delta t in the SAPA book)
X   [5] scratch-dft (keyword; vector of correct length)
X       ==> vector in which the in-place dft is done
X   [6] data-taper (keyword; nil)
X       ==> nil or a tapering function
X   [7] data-taper-parameters (keyword)
X       ==> parameters for tapering function (not used
X           if data-taper is nil)
X   [8] spec-wind-transformation (keyword; #'convert-to-dB)
X       ==> a function of one argument or nil;
X           if bound to a function, the function is used
X           to transform all elements of result-spec-wind
X   [9] result-spec-wind (keyword; vector of correct length)
X       <== vector into which spectral window values are placed;
X           it must be exactly of the length dictated by N-nonzero-freqs
X  [10] return-frequencies-p (keyword; t)
X       ==> if t, the frequencies associated with the spectral window
X           are computed and returned in result-freq
X  [11] freq-transformation (keyword; nil)
X       ==> a function of one argument or nil;
X           if bound to a function, the function is used
X           to transform all elements of result-freq
X           (ignored unless return-frequencies-p is true)
X  [12] result-freq (keyword; nil or vector of correct length)
X       <== not used if return-frequencies-p nil; otherwise,
X           vector of N-nonzero-freqs +1 into which
X           the frequencies associated with the values
X           in result-spec-wind are placed
returns
X   [1] result-spec-wind, a vector holding
X       the properly transformed spectral window
X   [2] result-freq (if return-frequencies-p is t),
X       where result-freq is a vector holding
X       the properly transformed frequencies
X       associated with values in result-spec-wind
X        -- or --
X       nil (if return-frequencies-p is nil)
X   [3] the length of the vector result-spec-wind
---
Note: see equation below Equation (244a) of the SAPA book"
X  ;;; put data taper into first part of scratch-dft
X  (fill scratch-dft 1.0 :end sample-size)
X  (if data-taper
X    (funcall data-taper (make-array
X                         sample-size
X                         :displaced-to scratch-dft)
X             :taper-parameter data-taper-parameters
X             :normalization :N))
X  (acvs scratch-dft :end sample-size :center-data-p nil :result scratch-dft)
X  (lag-window-spectral-estimate
X   scratch-dft
X   lag-window-function
X   :max-lag (1- sample-size)
X   :N-ts sample-size
X   :N-nonzero-freqs N-nonzero-freqs
X   :return-est-for-0-freq-p t
X   :sampling-time sampling-time
X   :scratch-dft scratch-dft
X   :sdf-transformation spec-wind-transformation
X   :result-sdf result-spec-wind
X   :return-frequencies-p return-frequencies-p
X   :freq-transformation freq-transformation
X   :result-freq result-freq)
X  (values result-spec-wind result-freq
X          (get-N-freqs N-nonzero-freqs sample-size t)))
X
#|
(multiple-value-bind (the-spec-wind-5 freqs N-f)
X                     (spectral-window-for-lag-window-spectral-estimate 
X                      20
X                      #'(lambda (tau)
X                          (parzen-lag-window tau 5))
X                      :N-nonzero-freqs :Fourier)
X  (let ((the-spec-wind-10 (spectral-window-for-lag-window-spectral-estimate 
X                             20
X                             #'(lambda (tau)
X                                 (parzen-lag-window tau 10))
X                             :N-nonzero-freqs :Fourier))
X        (the-spec-wind-15 (spectral-window-for-lag-window-spectral-estimate 
X                             20
X                             #'(lambda (tau)
X                                 (parzen-lag-window tau 15))
X                             :N-nonzero-freqs :Fourier)))
X    (dotimes (i N-f)
X      (format t "~&~6,4F: ~8,4F ~8,4F ~8,4F"
X              (svref freqs i)
X              (svref the-spec-wind-5  i)
X              (svref the-spec-wind-10 i)
X              (svref the-spec-wind-15 i)))
X    (values)))
;==>
0.0000:   5.4920   8.2174   9.6800
0.0500:   5.0695   6.6304   6.4004
0.1000:   3.7838   1.7435  -2.6115
0.1500:   1.5802  -5.7033  -8.2604
0.2000:  -1.6189 -10.3761 -10.2127
0.2500:  -5.8104 -11.9382 -12.6393
0.3000: -10.3570 -12.6772 -13.8851
0.3500: -13.1849 -14.3975 -14.7067
0.4000: -13.9394 -15.3618 -15.4563
0.4500: -14.2726 -15.5651 -15.7582
0.5000: -14.4370 -15.3760 -15.8875
|#
X
;-------------------------------------------------------------------------------
(defun smoothing-window-for-lag-window-spectral-estimate
X       (sample-size
X        lag-window-function
X        &key
X        (N-nonzero-freqs :twice-next-power-of-2)
X        (sampling-time 1.0)
X        (scratch-dft (make-array (get-N-dft N-nonzero-freqs sample-size)))
X        (smooth-wind-transformation #'careful-convert-to-dB)
X        (result-smooth-wind
X         (make-array (get-N-freqs N-nonzero-freqs sample-size
X                                  t)))
X        (return-frequencies-p t)
X        (freq-transformation nil)
X        (result-freq
X         (if return-frequencies-p
X           (make-array (get-N-freqs N-nonzero-freqs sample-size t)))))
X  "given
X   [1] sample-size (required)
X       ==> sample size for which smoothing window is to be computed
X   [2] lag-window-function (required)
X       ==> function of a single variable that computes the value
X           of the lag window for a given lag
X   [3] N-nonzero-freqs (keyword; :twice-next-power-of-2)
X       ==> specifies at how many nonzero frequencies
X           smoothing window is to be computed -- choices are:
X           :half-next-power-of-2
X            ==> 1/2 * power of two >= sample size;
X           :next-power-of-2
X            ==> power of two >= sample size;
X           :twice-next-power-of-2
X            ==> 2 * power of two >= sample size;
X           :Fourier
X            ==> just at Fourier frequencies
X            -- or --
X           any power of 2 >= 1/2 * [power of two >= sample size]
X   [4] sampling-time (keyword; 1.0)
X       ==> sampling time (called delta t in the SAPA book)
X   [5] scratch-dft (keyword; vector of correct length)
X       ==> vector in which the in-place dft is done
X   [6] smooth-wind-transformation (keyword; #'convert-to-dB)
X       ==> a function of one argument or nil;
X           if bound to a function, the function is used
X           to transform all elements of result-smooth-wind
X   [7] result-smooth-wind (keyword; vector of correct length)
X       <== vector into which smoothing window values are placed;
X           it must be exactly of the length dictated by N-nonzero-freqs
X   [8] return-frequencies-p (keyword; t)
X       ==> if t, the frequencies associated with the smoothing window
X           are computed and returned in result-freq
X   [9] freq-transformation (keyword; nil)
X       ==> a function of one argument or nil;
X           if bound to a function, the function is used
X           to transform all elements of result-freq
X           (ignored unless return-frequencies-p is true)
X  [10] result-freq (keyword; nil or vector of correct length)
X       <== not used if return-frequencies-p nil; otherwise,
X           vector of N-nonzero-freqs +1 into which
X           the frequencies associated with the values
X           in result-smooth-wind are placed
returns
X   [1] result-smooth-wind, a vector holding
X       the properly transformed smoothing window
X   [2] result-freq (if return-frequencies-p is t),
X       where result-freq is a vector holding
X       the properly transformed frequencies
X       associated with values in result-smooth-wind
X        -- or --
X       nil (if return-frequencies-p is nil)
X   [3] the length of the vector result-smooth-wind
---
Note: see equation below Equation (237c) of the SAPA book"
X  (fill scratch-dft 1.0 :end sample-size)
X  (lag-window-spectral-estimate
X   scratch-dft
X   lag-window-function
X   :max-lag (1- sample-size)
X   :N-ts sample-size
X   :N-nonzero-freqs N-nonzero-freqs
X   :return-est-for-0-freq-p t
X   :sampling-time sampling-time
X   :scratch-dft scratch-dft
X   :sdf-transformation smooth-wind-transformation
X   :result-sdf result-smooth-wind
X   :return-frequencies-p return-frequencies-p
X   :freq-transformation freq-transformation
X   :result-freq result-freq)
X  (values result-smooth-wind result-freq
X          (get-N-freqs N-nonzero-freqs sample-size t)))
X
#|
(multiple-value-bind (the-smooth-wind-5 freqs N-f)
X                     (smoothing-window-for-lag-window-spectral-estimate 
X                      20
X                      #'(lambda (tau)
X                          (parzen-lag-window tau 5))
X                      :N-nonzero-freqs :Fourier)
X  (let ((the-smooth-wind-10 (smoothing-window-for-lag-window-spectral-estimate 
X                             20
X                             #'(lambda (tau)
X                                 (parzen-lag-window tau 10))
X                             :N-nonzero-freqs :Fourier))
X        (the-smooth-wind-15 (smoothing-window-for-lag-window-spectral-estimate 
X                             20
X                             #'(lambda (tau)
X                                 (parzen-lag-window tau 15))
X                             :N-nonzero-freqs :Fourier)))
X    (dotimes (i N-f)
X      (format t "~&~6,4F: ~8,4F ~8,4F ~8,4F"
X              (svref freqs i)
X              (svref the-smooth-wind-5  i)
X              (svref the-smooth-wind-10 i)
X              (svref the-smooth-wind-15 i)))
X    (values)))
;==>
0.0000:   5.7426   8.7506  10.5116
0.0500:   5.2934   6.9265   6.2889
0.1000:   3.9171   0.9068 -10.3965
0.1500:   1.5182 -12.1526 -27.9226
0.2000:  -2.1035 -149.5459 -16.4171
0.2500:  -7.3518 -20.9691 -35.2827
0.3000: -15.0060 -18.0163 -29.3197
0.3500: -22.4607 -26.4962 -27.4921
0.4000: -21.0266 -148.5768 -35.3403
0.4500: -20.6780 -29.5773 -31.7854
0.5000: -20.9691 -23.9794 -35.2827
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  Grenander-smoothing-window-bandwidth
;;;                 Jenkins-smoothing-window-bandwidth
;;;                 equivalent-degrees-of-freedom
;;;                 bandwidth&confidence-intervals-for-sdf-dB
;;;  compute the two measures of the smoothing window bandwidths for a lag
;;;  window estimator, the equivalent degrees of freedom for a lag window
;;;  estimator, and a point-wise confidence interval for the true sdf at
;;;  given frequency.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun Grenander-smoothing-window-bandwidth
X       (lag-window-function
X        max-lag
X        &key
X        (sampling-time 1.0))
X  "given a lag window function, a maximum lag and the sampling time,
returns Grenander's measure of smoothing window bandwidth
---
Note: see Equations (241c) and (241b) of the SAPA book;
X      unfortunately, because of the square root operation,
X      this measure can be complex-valued for certain lag windows"
X  (let ((sum 0.0)
X        (tau 1))
X    (dotimes (i max-lag (/ (sqrt (1+ (* (/ 12 (* pi pi)) sum)))
X                           sampling-time))
X      (incf sum (/ (* (expt -1 tau) (funcall lag-window-function tau))
X                   (* tau tau)))
X      (incf tau))))
X
#|
(dolist (m '(5 10 13) (values))
X  (print (Grenander-smoothing-window-bandwidth
X          #'(lambda (tau)
X              (parzen-lag-window tau m))
X          20
X          :sampling-time 1/4)))
;==>
1.4444046775094332 
0.7450539900333778 
0.5761648410683915 
|#
X
;-------------------------------------------------------------------------------
(defun Jenkins-smoothing-window-bandwidth
X       (lag-window-function
X        max-lag
X        &key
X        (sampling-time 1.0))
X  "given a lag window function, a maximum lag and the sampling time,
returns Jenkins' measure of smoothing window bandwidth
---
Note: see Equation (242c) of the SAPA book"
X  (let ((sum 0.0)
X        (tau 1))
X    (dotimes (i max-lag (/ (* sampling-time (+ 1.0 (* 2.0 sum)))))
X      (incf sum (expt (funcall lag-window-function tau) 2))
X      (incf tau))))
X
#|
(dolist (m '(5 10 13) (values))
X  (print (Jenkins-smoothing-window-bandwidth
X          #'(lambda (tau)
X              (parzen-lag-window tau m))
X          20
X          :sampling-time 1/4)))
;==>
1.4822720265623148 
0.7417022065640646 
0.5705456589192851 
|#
X
;-------------------------------------------------------------------------------
(defun equivalent-degrees-of-freedom
X       (lag-window-function
X        max-lag
X        &key
X        (sampling-time 1.0)
X        (sample-size (1+ max-lag))
X        (C_h 1.0))
X  "given a lag window function, a maximum lag, the sampling time,
the sample size and the variance inflation factor C_h,
returns the corresponding equivalent degrees of freedom
---
Note: see Equation (255a) of the SAPA book"
X  (/ (* 2 sample-size sampling-time
X        (Jenkins-smoothing-window-bandwidth
X         lag-window-function max-lag :sampling-time sampling-time))
X     C_h))
X
#|
;;; two examples discussed on page 302 of the SAPA book ...
(equivalent-degrees-of-freedom
X #'(lambda (tau) (parzen-lag-window tau 150))
X 150
X :sampling-time 1/4
X :sample-size 1024
X :C_h 1.96)
;==> 12.917060857770348
X
(equivalent-degrees-of-freedom
X #'(lambda (tau) (parzen-lag-window tau 55))
X 55
X :sampling-time 1/4
X :sample-size 1024
X :C_h 1.96)
;==> 35.22834596647285
|#
X
;-------------------------------------------------------------------------------
(defun bandwidth&confidence-intervals-for-sdf-dB
X       (sdf-dB
X        freq
X        sample-size
X        &key
X        (confidence-level 0.95)
X        (lag-window-function nil)
X        (max-lag (1- sample-size))
X        (sampling-time 1.0)
X        (C_h 1.0))
X  "given
X   [1] sdf-dB (required)
X       ==> an sdf estimate (in decibels) for a single frequency
X   [2] freq (required)
X       ==> the frequency associated with sdf-dB
X   [3] sampling-size (required)
X       ==> the sample size of the time series
X           from which sdf-dB was computed
X   [4] confidence-level (keyword; 0.95)
X       ==> the level of the confidence interval
X           (e.g., 0.95 yields a 95% confidence interval)
X   [5] lag-window-function (keyword; nil)
X       ==> the lag window function used to create
X           sdf-dB (nil means a lag window function was NOT used)
X   [6] max-lag (keyword; sample-size - 1)
X       ==> the maximum lag used in conjunction with lag-window-function
X           to create sdf-dB (not used in lag-window-function is nil)
X   [7] sampling-time (keyword; 1.0)
X       ==> sampling time (called delta t in the SAPA book)
X   [8] C_h (keyword; 1.0)
X       ==> the variance inflation factor due to tapering
X           (see Table 248 of the SAPA book); the default value
X           of 1.0 is appropriate for a rectangular data taper
returns
X   [1] left-hand size of interval describing the bandwidth of
X       the spectral estimate (interval is centered on freq)
X   [2] right-hand size of bandwidth interval
X   [3] lower limit (in decibels) of confidence interval for
X       true sdf based upon sdf-dB
X   [4] upper limit (in decibels) of confidence interval
---
Note: see Sections 6.7 and 6.10 of the SAPA book"
X  (let* ((half-B_W (* 0.5 (if lag-window-function
X                            (Jenkins-smoothing-window-bandwidth
X                             lag-window-function max-lag
X                             :sampling-time sampling-time)
X                            (/ C_h (* sampling-time sample-size)))))
X         (nu (if lag-window-function
X               (equivalent-degrees-of-freedom
X                lag-window-function
X                max-lag
X                :sampling-time sampling-time
X                :sample-size sample-size
X                :C_h C_h)
X               2))
X         (acc-level (if (>= nu 150) :quick :accurate))
X         (for-upper (/ (- 1.0 confidence-level) 2))
X         (lower-dB (convert-to-dB (/ nu
X                                     (quantile-of-chi-square-distribution
X                                      nu (+ confidence-level for-upper)
X                                      :accuracy-level acc-level))))
X         (upper-dB (convert-to-dB (/ nu
X                                     (quantile-of-chi-square-distribution
X                                      nu for-upper
X                                      :accuracy-level acc-level)))))
X    (values (- freq half-B_W)
X            (+ freq half-B_W)
X            (+ sdf-dB lower-dB)
X            (+ sdf-dB upper-dB))))
X
#|
;;; example discussed on page 302 of the SAPA book ...
(bandwidth&confidence-intervals-for-sdf-dB
X 6.3
X 1.0
X 1024
X :lag-window-function #'(lambda (tau) (parzen-lag-window tau 150))
X :max-lag 150
X :sampling-time 1/4
X :C_h 1.96)
;==>  0.9752759382019239
;     1.024724061798076
;     3.4987218467878627
;    10.458588688308273
;;; Note that (- 1.024724061798076 0.9752759382019239) ==> 0.0494,
;;; in agreement with the result stated on page 301 of the SAPA book
X
;;; example discussed on page 299 of the SAPA book ...
(bandwidth&confidence-intervals-for-sdf-dB
X 0
X 1.0
X 1024
X :sampling-time 1/4
X :C_h 1.96)
;==> 0.996171875
;    1.003828125
;   -5.668944635031615
;   15.965738982351615
;;; Note that (- 1.003828125 0.996171875) ==> 0.008,
;;; in agreement with the result stated on page 299 of the SAPA book
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  time-series-bandwidth
;;;                 sample-cepstrum
;;;                 cepstrum->I_m
;;;  can be used to objectively determine an appropriate value for
;;;  the lag window parameter m.  The function time-series-bandwidth
;;;  computes an estimate of the spectral bandwidth for a time series
;;;  (see the discussion in Section 6.14 of the SAPA book).  The functions
;;;  sample-cepstrum and cepstrum->I_m can be used to compute an approximation
;;;  to the mean integrated squared error (see the discussion in Section 6.15
;;;  of the SAPA book).
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun time-series-bandwidth (acvs &key (sampling-time 1.0))
X  "given
X   [1] acvs (required)
X       ==> a vector of length N with values of
X           the acvs (or acs) from lag 0 to N-1
X   [2] sampling-time (keyword; 1.0)
X       ==> the sampling time
returns
X   [1] unbiased estimate of trace bandwidth
X       (\tilde B_T in Equation (280))
X   [2] biased estimate of trace bandwidth
X       (\hat B_T in equation just above Equation (28)))
---
Note: see Section 6.14 of the SAPA book"
X  (let ((N (length acvs))
X        (tau 1)
X        (sum 0.0)
X        (s_0^2 (* (elt acvs 0) (elt acvs 0))))
X    (dotimes (i (1- N) )
X      (incf sum (* (elt acvs tau) (elt acvs tau)
X                   (- 1.0 (/ tau N))))
X      (incf tau))
X    (let ((hat-B_T (/ s_0^2 (* 2.0 sampling-time
X                               (+ s_0^2 (* 2.0 sum))))))
X      (values (- (* 5/3 hat-B_T)
X                 (/ (* N sampling-time)))
X              hat-B_T))))
X
#|
;;; For this example we use the 20 point time series of Section 6.16
;;; of the SAPA book:
(let* ((20-pt-ts #(71.0  63.0  70.0  88.0  99.0  90.0 110.0 135.0 128.0 154.0
X                   156.0 141.0 131.0 132.0 141.0 104.0 136.0 146.0 124.0 129.0))
X       (sampling-time 0.25)
X       (the-acvs (acvs 20-pt-ts)))
X  (time-series-bandwidth the-acvs :sampling-time sampling-time))
;==> 0.7889256471629704
;    0.5933553882977822A
|#
X
;-------------------------------------------------------------------------------
(defun sample-cepstrum
X       (log-spectrum
X        &key
X        (sampling-time 1.0)
X        (zero-frequency-p nil)
X        (Nyquist-frequency-p t)
X        (result (make-array (if zero-frequency-p
X                              (length log-spectrum)
X                              (1+ (length log-spectrum))))))
X  "given
X   [1] log-spectrum (required)
X       ==> a vector containing log spectrum
X   [2] sampling-time (keyword; 1.0)
X       ==> the sample time (i.e., delta t)
X   [3] zero-frequency-p (keyword; nil)
X       ==> if t, the first element of log-spectrum
X           is associated with 0 frequency;
X           if nil, first element goes with
X           lowest nonzero frequency
X   [4] Nyquist-frequency-p (keyword; nil)
X       ==> if t, the last element of log-spectrum
X           is associated with Nyquist frequency;
X           if nil, last element goes with
X           a frequency just less than Nyquist
X   [5] result (keyword; vector of appropriate length)
X       <== vector to hold sample cepstrum
returns
X   [1] result, i.e., the sample cepstrum from lag 0 to lag N_U
X       where N_U = length of log-spectrum - 1 if zero-frequency-p is t
X       or        = length of log-spectrum     if zero-frequency-p is nil
---
Note: see Equation 282 of the SAPA book"
X  (let* ((N_U (if zero-frequency-p (1- (length log-spectrum))
X                  (length log-spectrum)))
X         (N_L (if Nyquist-frequency-p (1- N_U) N_U))
X         (N-prime (+ N_U N_L 1))
X         (for-DFT (make-array N-prime))
X         (k (if zero-frequency-p 1 0))
X         (i+1 1))
X    #+mcl(declare (dynamic-extent for-DFT))
X    ;;; NOTE: if the log spectrum is not supplied for 0 frequency,
X    ;;;       this sets it equal to that of the closest non-zero frequency.
X    (setf (aref for-DFT 0) (elt log-spectrum 0))
X    (dotimes (i N_U)
X      (setf (aref for-DFT i+1) (elt log-spectrum k))
X      (incf i+1)
X      (incf k))
X    (setf k (if zero-frequency-p N_L (1- N_L)))
X    (dotimes (i N_L)
X      (setf (aref for-DFT i+1) (elt log-spectrum k))
X      (incf i+1)
X      (decf k))
X    (inverse-dft! for-DFT :sampling-time sampling-time)
X    (dotimes (i (1+ N_U) result)
X      (setf (aref result i) (realpart (aref for-DFT i))))))
X
#|
;;; For this example we use the 20 point time series of Section 6.16
;;; of the SAPA book:
(let* ((20-pt-ts #(71.0  63.0  70.0  88.0  99.0  90.0 110.0 135.0 128.0 154.0
X                   156.0 141.0 131.0 132.0 141.0 104.0 136.0 146.0 124.0 129.0))
X       (sampling-time 0.25)
X       (pgram-20-pt-ts-log (periodogram 20-pt-ts
X                                          :sampling-time sampling-time
X                                          :N-nonzero-freqs :Fourier
X                                          :sdf-transformation #'log))
X       (pgram-19-pt-ts-log (periodogram 20-pt-ts
X                                          :end 19
X                                          :sampling-time sampling-time
X                                          :N-nonzero-freqs :Fourier
X                                          :sdf-transformation #'log)))
X  (let ((cepstrum-20 (sample-cepstrum pgram-20-pt-ts-log))
X        (cepstrum-19 (sample-cepstrum pgram-19-pt-ts-log
X                                      :Nyquist-frequency-p nil)))
X    (dotimes (i (length cepstrum-20))
X      (format t "~&~2D: ~8,4F"
X              i
X              (svref cepstrum-20 i)))
X    (format t "~&...")
X    (dotimes (i (length cepstrum-19))
X      (format t "~&~2D: ~8,4F"
X              i
X              (svref cepstrum-19 i)))
X    (values)))
;==>
X 0:   4.3542
X 1:   0.9286
X 2:   0.2987
X 3:   0.3236
X 4:  -0.1946
X 5:   0.2022
X 6:  -0.1170
X 7:   0.1663
X 8:  -0.2263
X 9:  -0.1352
10:   0.2587
...
X 0:   4.4931
X 1:   0.9671
X 2:   0.1895
X 3:   0.0822
X 4:   0.2625
X 5:   0.0137
X 6:  -0.0898
X 7:   0.1136
X 8:  -0.1625
X 9:  -0.0512
|#
X
;-------------------------------------------------------------------------------
(defun cepstrum->I_m
X       (cepstrum
X        &key
X        (sampling-time 1.0)
X        (N-prime (- (* 2 (length cepstrum)) 2))
X        (lag-window-function #'parzen-lag-window)
X        (lambda-factor 1.0)
X        (result (make-array (length cepstrum))))
X  "given
X   [1] cepstrum (required)
X       ==> a vector of length N_U+1 with values of the cepstrum
X           from lag 0 to N_U
X   [2] sampling-time (keyword; 1.0)
X       ==> the sample time (i.e., delta t)
X   [3] N-prime (keyword; 2 * N_U)
X       ==> effective sample size, as discussed in Section 6.15
X           of the SAPA book (default value ok if N-prime is even,
X           but must be replaced if N-prime is odd)
X   [4] lag-window-function (keyword; #'parzen-lag-window)
X       ==> lag window function of 2 parameters,
X           tau (the lag) and m (the window parameter)
X   [5] lambda-factor (keyword; 1.0)
X       ==> lambda, as discussed in Section 6.15 of the SAPA book
X   [6] result (keyword; vector of same length as cepstrum)
X       <== vector to hold I_m
returns
X   [1] the value m at which I_m is minimized
X   [2] result, i.e., I_m for m = 0 to N_U
X   [3] the minimum value of I_m for m = 0 to N_U
---
Note: see Equation (283b) of the SAPA book"
X  (let* ((N_U (1- (length cepstrum)))
X         (N_L (- N-prime N_U 1))
X         (a-constant (/ (* lambda-factor lambda-factor
X                           +euler-constant+ +euler-constant+)
X                        sampling-time))
X         (another-constant (/ (* lambda-factor lambda-factor pi pi)
X                              (* 6.0 N-prime sampling-time)))
X         (yet-another-constant (/ another-constant
X                                  sampling-time))
X         first-sum second-sum wtm tau)
X    ;;; loop over all possible lag window parameters
X    (dotimes (m (1+ N_U))
X      (setf first-sum 0.0
X            second-sum 0.0
X            tau 1)
X      (dotimes (i N_L)
X        (setf wtm (funcall lag-window-function tau m))
X        (incf first-sum (* (- (expt (aref cepstrum tau) 2)
X                              yet-another-constant)
X                           (expt (- 1.0 wtm) 2)))
X        (incf second-sum (* wtm wtm))
X        (incf tau))
X      ;;; double these ...
X      (incf first-sum first-sum)
X      (incf second-sum second-sum)
X      ;;; add in the odd term if need be ...
X      (when (not (= N_L N_U))
X        (setf wtm (funcall lag-window-function tau m))
X        (incf first-sum (* (- (expt (aref cepstrum tau) 2)
X                              yet-another-constant)
X                           (expt (- 1.0 wtm) 2)))
X        (incf second-sum (* wtm wtm)))
X      (setf (aref result m) (+ a-constant
X                            (* sampling-time first-sum)
X                            (* another-constant second-sum))))
X    (multiple-value-bind (min-I_m min-m)
X                         (min-of-seq result)
X      (values min-m result min-I_m))))
X
#|
;;; For this example we use the 20 point time series of Section 6.16
;;; of the SAPA book:
(let* ((20-pt-ts #(71.0  63.0  70.0  88.0  99.0  90.0 110.0 135.0 128.0 154.0
X                   156.0 141.0 131.0 132.0 141.0 104.0 136.0 146.0 124.0 129.0))
X       (sampling-time 0.25)
X       (pgram-20-pt-ts-log (periodogram 20-pt-ts
X                                          :sampling-time sampling-time
X                                          :N-nonzero-freqs :Fourier
X                                          :sdf-transformation #'log))
X       (pgram-19-pt-ts-log (periodogram 20-pt-ts
X                                          :end 19
X                                          :sampling-time sampling-time
X                                          :N-nonzero-freqs :Fourier
X                                          :sdf-transformation #'log)))
X  (let ((cepstrum-20 (sample-cepstrum pgram-20-pt-ts-log))
X        (cepstrum-19 (sample-cepstrum pgram-19-pt-ts-log
X                                      :Nyquist-frequency-p nil)))
X    (multiple-value-bind (m-min I_m-20 I_m-min)
X                         (cepstrum->I_m cepstrum-20)
X      (print m-min)
X      (print I_m-min)
X      (dotimes (i (length I_m-20))
X        (format t "~&~2D: ~8,4F"
X                i
X                (svref I_m-20 i))))
X      (format t "~&...")
X      (multiple-value-bind (m-min I_m-19 I_m-min)
X                           (cepstrum->I_m cepstrum-19 :N-prime 19)
X        (print m-min)
X        (print I_m-min)
X        (dotimes (i (length I_m-19))
X          (format t "~&~2D: ~8,4F"
X                  i
X                  (svref I_m-19 i))))
X    (values)))
;==>
5 
-0.05104128407495598 
X 0:   1.3292
X 1:   1.3292
X 2:   0.6569
X 3:   0.1269
X 4:  -0.0210
X 5:  -0.0510
X 6:  -0.0336
X 7:   0.0100
X 8:   0.0715
X 9:   0.1449
10:   0.2252
...
4 
-0.4396526809044076 
X 0:   0.9690
X 1:   0.9690
X 2:   0.2372
X 3:  -0.3244
X 4:  -0.4397
X 5:  -0.4118
X 6:  -0.3428
X 7:  -0.2564
X 8:  -0.1608
X 9:  -0.0586
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  cumulative-periodogram
;;;                 quantile-of-Kolmogorov-test-statistic
;;;  can be used to compute and evaluate the normalized cumulative periodogram,
;;;  a useful way of evaluating the hypothesis that a time series can be
;;;  regarded as white noise.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun cumulative-periodogram
X       (time-series
X        &key
X        (sampling-time 1.0)
X        (center-data t)
X        (start 0)
X        (end (length time-series))
X        (cumulative-periodogram-test-p t)
X        (significance-level 0.95)
X        (scratch-dft (make-array (- end start)))
X        (return-frequencies-p t)
X        (result-freq (if return-frequencies-p
X                       (make-array (truncate (1- (- end start)) 2))))
X        (result (make-array (truncate (1- (- end start)) 2))))
X  "given
X   [1] time-series (required)
X       ==> a vector containing time series values
X   [2] sampling-time (keyword; 1.0)
X       ==> the sample time (i.e., delta t)
X   [3] center-data (keyword; t)
X       ==> if t, subtract sample mean from time series;
X           if a number, subtracts that number from time series;
X           if nil, time-series is not centered
X   [4] start (keyword; 0)
X       ==> start index of vector to be used
X   [5] end (keyword; length of time-series)
X       ==> 1 + end index of vector to be used
X   [6] cumulative-periodogram-test-p (keyword; t)
X       ==> if t, returns values associated with
X           Kolmogorov test statistic
X   [7] significance-level (keyword; 0.95)
X       ==> level of significance at which Kolmogorov test
X           is performed
X   [8] scratch-dft (keyword; vector of size (- end start))
X       ==> scratch vector used for in-place dft calculations
X   [9] return-frequencies-p (keyword; t)
X       ==> if t, the frequencies associated with the cumulative
X           periodogram are computed and returned in result-freq
X  [10] result-freq (keyword; nil or vector of correct length)
X       <== not used if return-frequencies-p nil; otherwise,
X           vector of correct length into which the frequencies
X           associated with the values in result are placed
X  [11] result (keyword; vector of correct length)
X       <== vector to hold normalized cumulative periodogram
returns
X   [1] vector with normalized cumulative periodogram 
X   [2] vector with associated frequencies
X and -- if cumulative-periodogram-test-p is true --
X   [3] Kolmogorov test statistic
X   [4] index of maximum deviation
X   [5] frequency of maximum deviation
X   [6] quantile of Kolmogorov test statistic (under null hypothesis)
X   [7] either :reject or :fail-to-reject, depending on whether
X       or not we reject or fail to reject the null hypothesis
X   [8] slope of upper and lower lines
X   [9] intercept of upper line
X  [10] intercept of lower line
---
Note: see Section 6.6 of the SAPA book"
X  (let* ((N (- end start))
X         ;;; M represents the number of Fourier frequencies
X         ;;; that are greater than zero and less than the Nyquist frequency
X         (M (truncate (/ (1- N) 2)))
X         ;;; M-1 is the effective number of samples for
X         ;;; the Kolmogorov test statistic
X         (M-1 (1- M))
X         (sampling-time-over-N (/ sampling-time N)))
X    ;;; Here we center the time series ...
X    (center&taper-time-series time-series
X                              :center-data center-data
X                              :start start
X                              :end end
X                              :result scratch-dft)
X    (dft! scratch-dft :N N)
X    (dotimes (i M)
X      (setf (aref scratch-dft i)
X            (* (expt (abs (aref scratch-dft (1+ i)))
X                     2)
X               sampling-time-over-N)))
X    (cumulative-sums scratch-dft :end M :result result)
X    (let ((factor (/ (elt result (1- M)))))
X      (a*x! factor result)
X      (if return-frequencies-p
X        (let ((1-over-N-sampling-time (/ (* N sampling-time))))
X          (dotimes (i M)
X            (setf (aref result-freq  i)
X                  (* (1+ i) 1-over-N-sampling-time)))))
X      (if cumulative-periodogram-test-p
X        (let* ((D+ (- (/ M-1) (elt result 0)))
X               (D- (elt result 0))
X               (D (max D+ D-))
X               (j-D 0)
X               (j 0)
X               (k 1)
X               D+current D-current KS-quantile)
X          (dotimes (i (1- M-1))
X            (incf j)
X            (incf k)
X            (setf D+current (- (/ k M-1) (elt result j))
X                  D-current (- (elt result j) (/ j M-1)))
X            (if (< D (max D+current D-current))
X              (setf j-D j
X                    D (max D+current D-current))))
X          (setf KS-quantile (quantile-of-Kolmogorov-test-statistic
X                             M-1
X                             :significance-level significance-level))
X          (values result
X                  result-freq
X                  D
X                  j-D
X                  (aref result-freq j-D)
X                  KS-quantile
X                  (if (> D KS-quantile) :reject :fail-to-reject)
X                  ;;; slope of upper/lowe lines
X                  (/ (* N sampling-time) M-1)
X                  ;;; intecept of upper line
X                  (- KS-quantile (/ M-1))
X                  ;;; intecept of lower line
X                  (- KS-quantile)))
X        (values result result-freq)))))
X
#|
;;; For this example we use the 20 point time series of Section 6.16
;;; of the SAPA book:
(let* ((20-pt-ts #(71.0  63.0  70.0  88.0  99.0  90.0 110.0 135.0 128.0 154.0
X                   156.0 141.0 131.0 132.0 141.0 104.0 136.0 146.0 124.0 129.0))
X       (fd (difference 20-pt-ts))
X       (sampling-time 0.25))
X  (let* ((cum-per-results (multiple-value-list
X                           (cumulative-periodogram
X                            20-pt-ts
X                            :sampling-time sampling-time)))
X         (cum-per (elt cum-per-results 0))
X         (cum-per-freqs (elt cum-per-results 1))
X         (KS-results (subseq cum-per-results 2)))
X    (print KS-results)
X    (dotimes (i (length cum-per))
X      (format t "~&~6,4F: ~8,4F"
X              (svref cum-per-freqs i)
X              (svref cum-per i))))
X  (format t "~&...")
X  (let* ((cum-per-results (multiple-value-list
X                           (cumulative-periodogram
X                            fd
X                            :sampling-time sampling-time)))
X         (cum-per (elt cum-per-results 0))
X         (cum-per-freqs (elt cum-per-results 1))
X         (KS-results (subseq cum-per-results 2)))
X    (print KS-results)
X    (dotimes (i (length cum-per))
X      (format t "~&~6,4F: ~8,4F"
X              (svref cum-per-freqs i)
X              (svref cum-per i))))
X    (values))
;==>
(0.6571803654140053 1 0.4 0.454 :reject 0.625 0.329 -0.454) 
0.2000:   0.6233
0.4000:   0.7822
0.6000:   0.8474
0.8000:   0.8893
1.0000:   0.8985
1.2000:   0.9326
1.4000:   0.9423
1.6000:   0.9925
1.8000:   1.0000
...
(0.2987185709446501 4 1.0526315789473684 0.454 :fail-to-reject 0.59375 0.329 -0.454) 
0.2105:   0.0977
0.4211:   0.1381
0.6316:   0.1800
0.8421:   0.2368
1.0526:   0.3263
1.2632:   0.6725
1.4737:   0.7619
1.6842:   0.9916
1.8947:   1.0000
|#
X
;-------------------------------------------------------------------------------
(defun quantile-of-Kolmogorov-test-statistic
X       (sample-size
X        &key
X        (significance-level 0.95))
X  "given
X   [1] sample-size (required)
X       ==> a positive integer
X   [2] significance-level (keyword; 0.95)
X       ==> 0.99, 0.98, 0.95, 0.90, or 0.80;
X           0.75 is also ok if sample-size > 40
returns
X   [1] quantile of two-sided Kolmogorov test statistic
---
Note: see Table 14, page 462, of Conover (2nd edition, 1980);
X      Stephens (1974); and Diggle (1990, page 55);
X      significance-level is currently limited to one of
X      these values: 0.99, 0.98, 0.95, 0.90, or 0.80;
X      if sample-size is greater than 40, 0.75 is also ok"
X  (cond
X   ((> sample-size 40) (let ((sqrt-sample-size (sqrt sample-size)))
X                         (/ (case significance-level
X                              (0.99 1.628)
X                              (0.98 1.52)
X                              (0.95 1.358)
X                              (0.90 1.224)
X                              (0.80 1.07)
X                              (0.75 1.02))
X                            (+ sqrt-sample-size
X                               0.12
X                               (/ 0.11 sqrt-sample-size)))))
X   (t
X    (nth (position significance-level
X                   '(0.80 0.90 0.95 0.98 0.99) :test #'=)
X         (nth (1- sample-size)
X              '((0.900 0.950 0.975 0.990 0.995)   ;;;  1 = sample size
X                (0.684 0.776 0.842 0.900 0.929)   ;;;  2
X                (0.565 0.636 0.708 0.785 0.829)   ;;;  3
X                (0.493 0.565 0.624 0.689 0.734)   ;;;  4
X                (0.447 0.509 0.563 0.627 0.669)   ;;;  5
X                (0.410 0.468 0.519 0.577 0.617)   ;;;  6
X                (0.381 0.436 0.483 0.538 0.576)   ;;;  7
X                (0.358 0.410 0.454 0.507 0.542)   ;;;  8
X                (0.339 0.387 0.430 0.480 0.513)   ;;;  9
X                (0.323 0.369 0.409 0.457 0.489)   ;;; 10
X                (0.308 0.352 0.391 0.437 0.468)   ;;; 11
X                (0.296 0.338 0.375 0.419 0.449)   ;;; 12
X                (0.285 0.325 0.361 0.404 0.432)   ;;; 13
X                (0.275 0.314 0.349 0.390 0.418)   ;;; 14
X                (0.266 0.304 0.338 0.377 0.404)   ;;; 15
X                (0.258 0.295 0.327 0.366 0.392)   ;;; 16
X                (0.250 0.286 0.318 0.355 0.381)   ;;; 17
X                (0.244 0.279 0.309 0.346 0.371)   ;;; 18
X                (0.237 0.271 0.301 0.337 0.361)   ;;; 19
X                (0.232 0.265 0.294 0.329 0.352)   ;;; 20
X                (0.226 0.259 0.287 0.321 0.344)   ;;; 21
X                (0.221 0.253 0.281 0.314 0.337)   ;;; 22
X                (0.216 0.247 0.275 0.307 0.330)   ;;; 23
X                (0.212 0.242 0.269 0.301 0.323)   ;;; 24
X                (0.208 0.238 0.264 0.295 0.317)   ;;; 25
X                (0.204 0.233 0.259 0.290 0.311)   ;;; 26
X                (0.200 0.229 0.254 0.284 0.305)   ;;; 27
X                (0.197 0.225 0.250 0.279 0.300)   ;;; 28
X                (0.193 0.221 0.246 0.275 0.295)   ;;; 29
X                (0.190 0.218 0.242 0.270 0.290)   ;;; 30
X                (0.187 0.214 0.238 0.266 0.285)   ;;; 31
X                (0.184 0.211 0.234 0.262 0.281)   ;;; 32
X                (0.182 0.208 0.231 0.258 0.277)   ;;; 33
X                (0.179 0.205 0.227 0.254 0.273)   ;;; 34
X                (0.177 0.202 0.224 0.251 0.269)   ;;; 35
X                (0.174 0.199 0.221 0.247 0.265)   ;;; 36
X                (0.172 0.196 0.218 0.244 0.262)   ;;; 37
X                (0.170 0.194 0.215 0.241 0.258)   ;;; 38
X                (0.168 0.191 0.213 0.238 0.255)   ;;; 39
X                (0.165 0.189 0.210 0.235 0.252)   ;;; 40
X                ))))))
X
#|
(quantile-of-Kolmogorov-test-statistic 17)
;==> 0.318
(quantile-of-Kolmogorov-test-statistic 17 :significance-level 0.9)
;==> 0.286
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  one-sided-freq->two-sided-freq
;;;                 one-sided-sdf->two-sided-sdf
;;;  allow conversion of a one-sided sdf estimate (or spectral/smoothing window)
;;;  into its two-sided representation (i.e., involving negative as well as
;;;  positive frequencies).
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun one-sided-freq->two-sided-freq
X       (one-sided-freq
X        &key
X        (result (make-array (1- (* 2 (length one-sided-freq))))))
X  "given a vector of N frequencies
0, f_1, ..., f_{N-1},
returns the vector of 2N-1 frequencies
-f_{N-1}, ..., -f_1, 0, f_1, ..., f_{N-1}"
X  (let* ((n (length one-sided-freq))
X         (j-up (1- n))
X         (j-down (1- n)))
X    (setf (aref result j-up) (aref one-sided-freq 0))
X    (dotimes (i (1- n) result)
X      (incf j-up)
X      (decf j-down)
X      (setf (aref result j-down)
X            (- (setf (aref result j-up)
X                     (aref one-sided-freq (1+ i))))))))
X
#|
(one-sided-freq->two-sided-freq #(0 1 2 3 4))
;==> #(-4 -3 -2 -1 0 1 2 3 4)
|#
X
;-------------------------------------------------------------------------------
(defun one-sided-sdf->two-sided-sdf
X       (one-sided-sdf
X        &key
X        (result (make-array (1- (* 2 (length one-sided-sdf))))))
X  "given a one-sided sdf
S(0), S(f_1), ..., S_(f_{N-1}) of length N,
returns the two-sized sdf
S_(-f_{N-1}), ..., S(-f_1), S(0), S(f_1), ..., S_(f_{N-1})
of length 2N-1, where S_(-f_k) = S_(f_k)"
X  (let* ((n (length one-sided-sdf))
X         (j-up (1- n))
X         (j-down (1- n)))
X    (dotimes (i n result)
X      (setf (aref result j-up)
X            (setf (aref result j-down)
X                  (aref one-sided-sdf i)))
X      (incf j-up)
X      (decf j-down))))
X
#|
(one-sided-sdf->two-sided-sdf #(1 7 2 5))
;==> #(5 2 7 1 7 2 5)
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  bartlett-lag-window
;;;                 bartlett-m->bandwidth
;;;                 bartlett-bandwidth->m
;;;                 bartlett-N-m->degrees-of-freedom
;;;                 daniell-lag-window
;;;                 daniell-m->bandwidth
;;;                 daniell-bandwidth->m
;;;                 daniell-N-m->degrees-of-freedom
;;;                 parzen-lag-window
;;;                 parzen-m->bandwidth
;;;                 parzen-bandwidth->m
;;;                 parzen-N-m->degrees-of-freedom
;;;                 papoulis-lag-window
;;;                 papoulis-m->bandwidth
;;;                 papoulis-bandwidth->m
;;;                 papoulis-N-m->degrees-of-freedom
;;;  implement four of the lag windows discussed in Section 6.10
;;;  of the SAPA book and allow computation of the smoothing window bandwidth
;;;  and equivalent degrees of freedom associated with these windows.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;
;;; Bartlett lag window ...
;;;
(defun bartlett-lag-window (tau m)
X  "given the lag tau and window parameter m,
returns the value of the Bartlett lag window
---
Note: see Equation (260) of the SAPA book
X      or Priestley, page 439, Equation (6.2.65)"
X  (assert (not (minusp m)))
X  (if (zerop tau)
X    1.0
X    (let ((abs-tau (abs tau)))
X      (if (< abs-tau m)
X        (- 1.0 (float (/ abs-tau m)))
X        0.0))))
X
#|
(map 'vector #'(lambda (lag) (bartlett-lag-window lag 4))
X     (iota -5 5))
;==> #(0.0 0.0 0.25 0.5 0.75 1.0 0.75 0.5 0.25 0.0 0.0)
|#
X
;-------------------------------------------------------------------------------
(defun bartlett-m->bandwidth (m &key (sampling-time 1.0))
X  "given window parameter m and sampling time,
returns bandwidth B_W for the Bartlett smoothing window
---
Note: see Table 269 of the SAPA book"
X  (/ 1.5 (* m sampling-time)))
X
#|
(map 'vector #'bartlett-m->bandwidth #(1 10 100 1000))
;==> #(1.5 0.15 0.015 0.0015)
|#
X
;-------------------------------------------------------------------------------
(defun bartlett-bandwidth->m (B_W &key (sampling-time 1.0))
X  "given desired smoothing window bandwidth B_W and sampling time,
returns 
X   [1] window parameter m required to approximately achieve B_W
X       using the Bartlett lag window
X   [2] actual B_W achieved
---
Note: see Table 269 of the SAPA book"
X  (let ((m (max 1 (round (/ 1.5 (* B_W sampling-time))))))
X    (values m (bartlett-m->bandwidth m))))
X
#|
(map 'vector #'bartlett-bandwidth->m #(0.1 0.01 0.001 0.0001))
;==> #(15 150 1500 15000)
|#
X
;-------------------------------------------------------------------------------
(defun bartlett-N-m->degrees-of-freedom (N m &key (C_h 1.0))
X  "given sample size N, window parameter m and variance inflation factor C_h,
X  returns equivalent degrees of freedom nu for Bartlett lag window
---
Note: see Table 269 of the SAPA book"
X  (/ (* 3.0 N) (* m C_h)))
X
#|
(bartlett-N-m->degrees-of-freedom 512 20)
;==> 76.8
|#
X
;-------------------------------------------------------------------------------
;;;
;;; Daniell lag window ...
;;;
(defun daniell-lag-window (tau m)
X  "given the lag tau and window parameter m,
returns the value of the Daniell lag window
---
Note: see equation between Equations (264a) and (264b) of the SAPA book
X      or Priestley, page 441, Equation (6.2.73)"
X  (assert (plusp m))
X  (cond
X   ((zerop tau) 1.0)
X   (t
X    (let ((ratio (/ (* pi tau) m)))
X      (/ (sin ratio) ratio)))))
X
#|
(map 'vector #'(lambda (lag) (daniell-lag-window lag 4))
X     (iota -5 5))
;==> #(-0.18006326323142122 3.898043091051478E-17 0.3001054387190354
X        0.6366197723675814 0.9003163161571061
X        1.0
X        0.9003163161571061 0.6366197723675814
X        0.3001054387190354 3.898043091051478E-17 -0.18006326323142122)
|#
X
;-------------------------------------------------------------------------------
(defun daniell-m->bandwidth (m &key (sampling-time 1.0))
X  "given window parameter m and sampling time,
returns bandwidth B_W for the Daniell smoothing window
---
Note: see Table 269 of the SAPA book"
X  (/ (* m sampling-time)))
X
#|
(map 'vector #'daniell-m->bandwidth #(1 10 100 1000))
;==> #(1.0 0.1 0.01 0.001)
|#
X
;-------------------------------------------------------------------------------
(defun daniell-bandwidth->m (B_W &key (sampling-time 1.0))
X  "given desired smoothing window bandwidth B_W and sampling time,
returns 
X   [1] window parameter m required to approximately achieve B_W
X       using the Daneill lag window
X   [2] actual B_W achieved (in fact, this is always equal to B_W,
X       but we return it anyway for consistency with other lag windows)
---
Note: see Table 269 of the SAPA book"
X  (values (/ (* B_W sampling-time)) B_W))
X
#|
(map 'vector #'daniell-bandwidth->m #(0.1 0.01 0.001 0.0001))
;==> #(10.0 100.0 1000.0 10000.0)
|#
X
;-------------------------------------------------------------------------------
(defun daniell-N-m->degrees-of-freedom (N m &key (C_h 1.0))
X  "given sample size N, window parameter m and variance inflation factor C_h,
X  returns equivalent degrees of freedom nu for Daniell lag window
---
Note: see Table 269 of the SAPA book"
X  (/ (* 2.0 N) (* m C_h)))
X
#|
(daniell-N-m->degrees-of-freedom 512 20)
;==> 51.2
|#
X
;-------------------------------------------------------------------------------
;;;
;;; Parzen lag window ...
;;;
(defun parzen-lag-window (tau m)
X  "given the lag tau and window parameter m,
returns the value of the Parzen lag window
---
Note: see the equation on page 265 of the SAPA book
X      or Priestley, page 443, Equation (6.2.82)"
X  (assert (not (minusp m)))
X  (cond
X   ((zerop tau) 1.0)
X   ((zerop m) 0.0)
X   (t
X    (let* ((abs-tau (abs tau))
X           (ratio (float (/ abs-tau m))))
X      (cond
X       ((<= abs-tau (/ m 2))
X        (1+ (* 6 (- (expt ratio 3) (* ratio ratio)))))
X       ((<= abs-tau m)
X        (* 2.0 (expt (- 1.0 ratio) 3)))
X       (t
X        0.0))))))
X
#|
(map 'vector #'(lambda (lag) (parzen-lag-window lag 4))
X     (iota -5 5))
;==> #(0.0 0.0 0.03125 0.25 0.71875 1.0 0.71875 0.25 0.03125 0.0 0.0)
|#
X
;-------------------------------------------------------------------------------
(defun parzen-m->bandwidth (m &key (sampling-time 1.0))
X  "given window parameter m and sampling time,
returns bandwidth B_W for the Parzen smoothing window
---
Note: see Table 269 of the SAPA book"
X  (/ 1.85 (* m sampling-time)))
X
#|
(map 'vector #'parzen-m->bandwidth #(1 10 100 1000))
;==> #(1.85 0.185 0.018500000000000003 0.00185)
|#
X
;-------------------------------------------------------------------------------
(defun parzen-bandwidth->m (B_W &key (sampling-time 1.0))
X  "given desired smoothing window bandwidth B_W and sampling time,
returns
X   [1] window parameter m required to approximately achieve B_W
X       using the Parzen lag window
X   [2] actual B_W achieved
---
Note: see Table 269 of the SAPA book"
X  (let ((m (max 1 (round (/ 1.85 (* B_W sampling-time))))))
X    (values m (parzen-m->bandwidth m))))
X
#|
(map 'vector #'parzen-bandwidth->m #(0.1 0.01 0.001 0.0001))
;==> #(18 185 1850 18500)
|#
X
;-------------------------------------------------------------------------------
(defun parzen-N-m->degrees-of-freedom (N m &key (C_h 1.0))
X  "given sample size N, window parameter m and variance inflation factor C_h,
X  returns equivalent degrees of freedom nu for Parzen lag window
---
Note: see Table 269 of the SAPA book"
X  (/ (* 3.71 N) (* m C_h)))
X
#|
(parzen-N-m->degrees-of-freedom 512 20)
;==> 94.976
|#
X
;-------------------------------------------------------------------------------
;;;
;;; Papoulis lag window ...
;;;
(defun papoulis-lag-window (tau m)
X  "given the lag tau and window parameter m,
returns the value of the Papoulis lag window
---
Note: see equation near bottom of page 266 of the SAPA book"
X  (assert (not (minusp m)))
X  (cond
X   ((zerop tau) 1.0)
X   ((zerop m) 0.0)
X   ((< (abs tau) m)
X    (let* ((ratio (float (/ tau m)))
X           (prod (* pi ratio)))
X      (+ (/ (abs (sin prod)) pi)
X         (* (- 1.0 (abs ratio))
X            (cos prod)))))
X   (t
X    0.0)))
X
#|
(map 'vector #'(lambda (lag) (papoulis-lag-window lag 4))
X     (iota -5 5))
;==> #(0.0 0.0 0.048302383742639676 0.31830988618379075 0.7554091649291872
X       1.0
X       0.7554091649291872 0.31830988618379075 0.048302383742639676 0.0 0.0)
|#
X
;-------------------------------------------------------------------------------
(defun papoulis-m->bandwidth (m &key (sampling-time 1.0))
X  "given window parameter m and sampling time,
returns bandwidth B_W for the Papoulis smoothing window
---
Note: see Table 269 of the SAPA book"
X  (/ 1.70 (* m sampling-time)))
X
#|
(map 'vector #'papoulis-m->bandwidth #(1 10 100 1000))
;==> #(1.7 0.16999999999999998 0.017 0.0017)
|#
X
;-------------------------------------------------------------------------------
(defun papoulis-bandwidth->m (B_W &key (sampling-time 1.0))
X  "given desired smoothing window bandwidth B_W and sampling time,
returns
X   [1] window parameter m required to approximately achieve B_W
X       using the Papoulis lag window
X   [2] actual B_W achieved
---
Note: see Table 269 of the SAPA book"
X  (let ((m (max 1 (round (/ 1.70 (* B_W sampling-time))))))
X    (values m (papoulis-m->bandwidth m))))
X
#|
(map 'vector #'papoulis-bandwidth->m #(0.1 0.01 0.001 0.0001))
;==> #(17 170 1700 17000)
|#
X
;-------------------------------------------------------------------------------
(defun papoulis-N-m->degrees-of-freedom (N m &key (C_h 1.0))
X  "given sample size N, window parameter m and variance inflation factor C_h,
X  returns equivalent degrees of freedom nu for Papoulis lag window
---
Note: see Table 269 of the SAPA book"
X  (/ (* 3.41 N) (* m C_h)))
X
#|
(papoulis-N-m->degrees-of-freedom 512 20)
;==> 87.296
X
;-------------------------------------------------------------------------------
;;; Pathology checks.
(bartlett-lag-window 0 0)  ;1.0
(bartlett-lag-window 1 0)  ;0.0
(daniell-lag-window 0 0)   ;assert error
(daniell-lag-window 1 0)   ;assert error
(parzen-lag-window 0 0)    ;1.0
(parzen-lag-window 1 0)    ;0.0
(papoulis-lag-window 0 0)  ;1.0
(papoulis-lag-window 1 0)  ;0.0
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  Everything below here consists of internal symbols in the SAPA package
;;;  and should be regarded as "dirty laundry" ...
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun get-N-dft (N-nonzero-freqs sample-size)
X  (case N-nonzero-freqs
X    (:half-next-power-of-2
X     (next-power-of-2 sample-size))
X    (:next-power-of-2
X     (* 2 (next-power-of-2 sample-size)))
X    (:twice-next-power-of-2
X     (* 4 (next-power-of-2 sample-size)))
X    (:Fourier sample-size)
X    (otherwise
X     (if (and (integerp N-nonzero-freqs)
X              (power-of-2 N-nonzero-freqs)
X              (>= (* 2 N-nonzero-freqs) sample-size))
X       (* 2 N-nonzero-freqs)
X       (error "N-nonzero-freqs (~A) is should be either
[a] a power of 2 greater than or equal to ~F (half the sample size)
or
[b] one of the following keywords: half-next-power-of-2,
X    next-power-of-2, twice-next-power-of-2 or Fourier"
X              N-nonzero-freqs
X              (float (/ sample-size 2))
X              )))))
X
#|
(get-N-dft :half-next-power-of-2 77)  ;==> 128
(get-N-dft :next-power-of-2 77)       ;==> 256
(get-N-dft :twice-next-power-of-2 77) ;==> 512
(get-N-dft :Fourier 76)               ;==> 76
(get-N-dft :Fourier 77)               ;==> 77
(get-N-dft 64 77)                     ;==> 128
(get-N-dft 32 77)                     ;==> error
(get-N-dft :foobar 78)                ;==> error
|#
X
;-------------------------------------------------------------------------------
(defun get-N-freqs
X       (N-nonzero-freqs
X        sample-size
X        return-est-for-0-freq-p)
X  (let ((temp (truncate (get-N-dft N-nonzero-freqs sample-size) 2)))
X    (if return-est-for-0-freq-p
X      (1+ temp)
X      temp)))
X
#|
(get-N-freqs :half-next-power-of-2 77 t)     ;==> 65
(get-N-freqs :half-next-power-of-2 77 nil)   ;==> 64
(get-N-freqs :next-power-of-2 77 t)          ;==> 129
(get-N-freqs :next-power-of-2 77 nil)        ;==> 128
(get-N-freqs :twice-next-power-of-2 77 t)    ;==> 257
(get-N-freqs :twice-next-power-of-2 77 nil)  ;==> 256
(get-N-freqs :Fourier 76 nil)                ;==> 38 
(get-N-freqs :Fourier 77 nil)                ;==> 38
(get-N-freqs :Fourier 78 nil)                ;==> 39
(get-N-freqs :Fourier 76 t)                  ;==> 39 
(get-N-freqs :Fourier 77 t)                  ;==> 39
(get-N-freqs :Fourier 78 t)                  ;==> 40
(get-N-freqs :foobar 78 t)                   ;==> error
(get-N-freqs 64 78 t)                        ;==> 65
(get-N-freqs 63 78 t)                        ;==> error
(get-N-freqs 32 78 t)                        ;==> error
|#
X
;-------------------------------------------------------------------------------
;;; The next 5 functions support wosa-spectral-estimate
(defun wosa-get-N-freqs
X       (block-size
X        oversampling-factor
X        return-est-for-0-freq-p)
X  (let ((temp (* oversampling-factor (/ block-size 2))))
X    (if return-est-for-0-freq-p
X      (1+ temp)
X      temp)))
X
#|
(wosa-get-N-freqs 256 1 t)         ;==> 129
(wosa-get-N-freqs 256 1 nil)       ;==> 128
(wosa-get-N-freqs 256 2 t)         ;==> 257
(wosa-get-N-freqs 256 2 nil)       ;==> 256
|#
X
;-------------------------------------------------------------------------------
(defun calculate-number-of-blocks (sample-size block-size proportion-of-overlap)
X  (+ 1 (truncate
X        (/ (float (- sample-size block-size))
X           (* (float block-size) (- 1.0 proportion-of-overlap))))))
X
;-------------------------------------------------------------------------------
(defun get-offset-to-kth-block (sample-size block-size number-of-blocks k)
X  (if (>  number-of-blocks 1)
X    (truncate (* k (- sample-size block-size))
X              (1- number-of-blocks))
X    0))
X
#|
(calculate-number-of-blocks 1024 256 0.5)  ;==> 7
(dotimes (k 7 (values))
X  (print (get-offset-to-kth-block 1024 256 7 k)))
;==> 0 128 256 384 512 640 768
X
(calculate-number-of-blocks 1024 128 0.5)  ;==> 15
(dotimes (k 15 (values))
X  (print (get-offset-to-kth-block 1024 128 15 k)))
;==> 0 64 128 192 256 320 384 448 512 576 640 704 768 832 896
|#
X
;-------------------------------------------------------------------------------
(defun equivalent-dof-for-wosa
X       (N N_S N_B vector-with-data-taper)
X  "given the sample size N, the block size N_S, the number of blocks N_B,
and the values of the data taper,
returns the equivalent degrees of freedom for wosa
using Equation (292b) of the SAPA book"
X  (let ((acs-taper (acvs vector-with-data-taper
X                         :center-data-p nil :acs-p t))
X        (n-shift (round (if (= N_B 1)
X                          0
X                          (/ (- N N_S) (1- N_B))))))
X    (transform-a-sequence! #'(lambda (x) (expt x 2)) acs-taper)
X    (let ((sum 0.0)
X          (m 0))
X      (dotimes (m-1 (1- N_B))
X        (incf m)
X        (if (< (* m n-shift) N_S)
X          (incf sum (* (- N_B m)
X                       (aref acs-taper (* m n-shift))))
X          (return)))
X      (values (/ (* 2 N_B) (1+ (/ (* 2.0 sum) N_B)))))))
X
#|
(dolist (N_S '(64 16 8) (values))
X  (print
X   (equivalent-dof-for-wosa
X    1024
X    N_S
X    (calculate-number-of-blocks 1024 N_S 0.5)
X    (Hanning-data-taper! (make-array N_S :initial-element 1.0)))))
;==>  58.45100403500567
;    233.79011528097186
;    453.18366286133556
|#
X
;-------------------------------------------------------------------------------
(defun equivalent-dof-for-wosa-standard-case (N_S)
X  "given number of blocks N_S,
returns the equivalent degrees of freedom for wosa
with 50% overlap and Hanning data taper
using Equation (294) of the SAPA book"
X  (float (/ (* 36 N_S N_S)
X            (1- (* 19 N_S)))))
#|
(dolist (N_S '(64 16 8) (values))
X  (print (equivalent-dof-for-wosa-standard-case
X          (calculate-number-of-blocks 1024 N_S 0.5))))
;==>  58.83673469387755
;    240.73134328358208
;    483.2576383154418
;;; Note: comparison of these results with the above
;;;       indicates that Equation (294) agrees with
;;;       Equation (292b) to within about 10%.
|#
SHAR_EOF
chmod 0644 nonparametric.lisp ||
echo 'restore of nonparametric.lisp failed'
Wc_c="`wc -c < 'nonparametric.lisp'`"
test 111652 -eq "$Wc_c" ||
	echo 'nonparametric.lisp: original size 111652, current size' "$Wc_c"
fi
# ============= parametric.lisp ==============
if test -f 'parametric.lisp' -a X"$1" != X"-c"; then
	echo 'x - skipping parametric.lisp (File already exists)'
else
echo 'x - extracting parametric.lisp (Text)'
sed 's/^X//' << 'SHAR_EOF' > 'parametric.lisp' &&
;;;-*- Mode: LISP; Package: :SAPA; Syntax: COMMON-LISP -*-
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;
;  parametric.lisp
;
;  a collection of Lisp functions for parametric spectral estimation ...
;  Note:  before compiling and loading parametric.lisp,
;         you should compile and load (in the order listed)
;            sapa-package.lisp, utilities.lisp, basic-math.lisp,
;            matrix.lisp, basic-statistics.lisp, dft-and-fft.lisp
;            filtering.lisp, random.lisp and acvs.lisp
;
;  6/3/93
;
;  SAPA, Version 1.0; Copyright 1993, Donald B. Percival, All Rights Reserved
;
;  Use and copying of this software and preparation of derivative works
;  based upon this software are permitted.  Any distribution of this
;  software or derivative works must comply with all applicable United
;  States export control laws.
; 
;  This software is made available AS IS, and no warranty -- about the
;  software, its performance, or its conformity to any
;  specification -- is given or implied. 
;
;  Comments about this software can be addressed to dbp@apl.washington.edu
;-------------------------------------------------------------------------------
;;; (compile-file "ccl:SAPA;parametric.lisp")
;;; (load "ccl:SAPA;parametric.fasl")
;-------------------------------------------------------------------------------
(if (not (find-package :SAPA))
X  (defpackage :SAPA (:USE :COMMON-LISP)))
X
(in-package :SAPA)
X
(export '(;;; functions for fitting AR models to data ...
X          yule-walker-algorithm-given-data
X          yule-walker-algorithm-given-acvs
X          burg-algorithm
X          fast-forward-backward-ls
X          forward-backward-ls
X          forward-ls
X          two-step-Burg-algorithm
X
X          ;;; functions for converting from one AR quantity to another ...
X          ar-coeffs->acvs
X          ar-coeffs->variance
X          ar-coeffs->reflection-coeffs
X          ar-coeffs->prewhitening-filter
X          reflection-coeffs->ar-coeffs
X          reflection-coeffs->variance
X
X          ;;; functions for computing AR sdf estimates ...
X          ar-coeffs->sdf
X          ar-coeffs->sdf-at-single-freq
X          sdf->sdf-with-accurate-peaks
X          integrate-sdf
X
X          ;;; functions for searching for peaks in an sdf estimate ...
X          find-peak-of-ar-sdf-using-quadratic-approx
X          find-peak-or-valley-of-sdf-using-bisection+Newton-Raphson
X
X          ;;; functions for generating observed innovations ...
X          generate-forward-innovations
X          generate-backward-innovations
X
X          ;;; function for regression analysis with AR errors ...
X          regression-with-ar-errors
X          ))
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  yule-walker-algorithm-given-data
;;;                 yule-walker-algorithm-given-acvs
;;;                 burg-algorithm
;;;                 fast-forward-backward-ls
;;;                 forward-backward-ls
;;;                 forward-ls
;;;                 two-step-Burg-algorithm
;;;  all estimate the parameters of an AR model for a given time series
;;;  and model order p.  Each function returns at least a vector of estimated
;;;  AR coefficients and an estimate of the associated innovations
;;;  variance, from which an AR sdf estimate can be computed or a
;;;  prewhitening filter be created.  There are also other intermediate
;;;  computations returned by these functions that might be of future
;;;  interest -- what these are depend on the nature of the particular
;;;  estimation scheme.  Chapter 9 of the SAPA book is devoted to
;;;  parametric spectral estimation.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun yule-walker-algorithm-given-data
X       (time-series
X        p
X        &key
X        (start 0)
X        (end (length time-series))
X        (center-data t)  ;t, nil or value to be subtracted off ...
X        (AR-coeffs (make-array p))
X        (approximate-MSEs (make-array (1+ p)))
X        (reflection-coeffs (make-array p))
X        (forward-prediction-errors (make-array (- end start p)))
X        (backward-prediction-errors (make-array (- end start p))))
X  "given
X   [1] time-series (required)
X       ==> a vector containing a real-valued time series x_t
X   [2] p (required)
X       ==> autoregressive model order;
X           p should be an integer > 0 and < end - start
X   [3] start (keyword; 0)
X       ==> start index of vector to be used
X   [4] end (keyword; length of the-seq)
X       ==> 1 + end index of vector to be used
X   [5] center-data (keyword; t)
X       ==> if t, subtract sample mean from time series;
X           if a number, subtracts that number from time series;
X           if nil, the effect of this function is to return
X           a copy of the relevant portion of time-series
X   [6] AR-coeffs (keyword; a vector of length p)
X       <== estimates of AR coefficients phi_{1,p}, ..., phi_{p,p}
X   [7] approximate-MSEs (keyword; a vector of length p+1)
X       <== estimates of innovations variance (mean square errors)
X           for models of order 0, 1, ..., p (the innovations variance
X           for the zeroth order model is the sample variance)
X   [8] reflection-coeffs (keyword; a vector of length p)
X       <== estimates of reflection coefficients phi_{1,1}, ..., phi_{p,p}
X   [9] forward-prediction-errors (keyword; a vector of length end - start - p)
X       <== computed forward prediction errors
X  [10] backward-prediction-errors (keyword; a vector of length end - start - p)
X       <== computed backward prediction errors
uses the Yule-Walker method to estimate the phi's in the model
X   x_t = phi_{1,p} x_{t-1} + ... + phi_{p,p} x_{t-p} + e_t
and returns
X   [1] AR-coeffs, a vector containing phi_{1,p}, ..., phi_{p,p}
X   [2] estimate of the innovations variance, i.e.,
X       the variance of the e_t's; this number is also given in
X       (elt approximate-MSEs p)
X   [3] approximate-MSEs
X   [4] reflection-coeffs
X   [5] forward-prediction-errors
X   [6] backward-prediction-errors
---
Note: see p. 420 of the SAPA book"
X  (let ((N-ts (- end start)))
X    (assert (and (integerp p) (> p 0) (< p N-ts)))
X    (let ((forward-j (make-array (+ N-ts p) :initial-element 0.0)))
X      (if center-data
X        (replace forward-j (center&taper-time-series 
X                            time-series
X                            :center-data center-data
X                            :start start
X                            :end end))
X        (replace forward-j time-series :start2 start :end2 end))
X      (let ((scratch (make-array (+ N-ts p)))
X            (backward-j (circular-shift-sequence forward-j)))
X        ;;; zeros at end don't matter ...
X        (setf (elt approximate-MSEs 0)
X              (float (/ (sum-of-squares forward-j)
X                        N-ts)))
X        (dotimes (j p)
X          (setf (elt reflection-coeffs j) (/ (dot-product backward-j forward-j)
X                                             (dot-product backward-j backward-j))
X                (elt approximate-MSEs (1+ j)) (* (elt approximate-MSEs j)
X                                                 (- 1.0
X                                                    (expt
X                                                     (elt reflection-coeffs j)
X                                                     2))))
X          (a*x+y! (- (elt reflection-coeffs j))
X                  backward-j
X                  forward-j
X                  :result scratch)
X          (a*x+y! (- (elt reflection-coeffs j))
X                  forward-j
X                  backward-j
X                  :result backward-j)
X          (circular-shift-sequence backward-j :result backward-j)
X          (copy-vector scratch forward-j))
X        ;;; At this point we have the reflection coefficients and
X        ;;; approximate mse's, from which we can now generate
X        ;;; the AR coefficients.
X        (reflection-coeffs->ar-coeffs
X         reflection-coeffs
X         :p p
X         :scratch scratch
X         :result AR-coeffs)
X        (values AR-coeffs
X                (elt approximate-MSEs p)
X                approximate-MSEs
X                reflection-coeffs
X                ;;; forward prediction errors ...
X                (copy-vector forward-j forward-prediction-errors
X                             :start p)
X                ;;; backward prediction errors ...
X                (copy-vector backward-j backward-prediction-errors
X                             :start (1+ p)))))))
X
#|
(yule-walker-algorithm-given-data
X #(71 63 70 88 99 90 110 135 128 154 156 141 131 132 141 104 136 146 124 129)
X 3)
;==> #(0.7564278899100625 -0.046223877768903764 0.03592061284674202)
;    351.6342314211137
;    #(782.64 352.21671265971577 352.088527878192 351.6342314211137)
;    #(0.7415951139732189 -0.01907713943936502 0.03592061284674202)
;    #(5.606819467197429 3.602049495972566 -13.13807178312505 13.531670850282099 22.591971410846554 -2.0709627659062573 29.661217150751934 10.772509547539435 -4.287081120362002 -3.7821509501885506 4.016928556685038 12.33707108178707 -31.065349921167176 29.302576292578536 13.063314822386415 -13.692737312779835 7.261455431834833)
;    #(-6.385268577444424 -19.24336074828875 -25.027314595591257 -16.48244854145692 1.3518647019431667 -21.36965186210829 -21.537852188257478 9.887122637410975 -16.148745552011714 8.004246630102024 20.852505588217184 13.139722849466036 4.1283725348055595 -4.53922156293019 33.568568323879624 -26.384631892925015 -3.14543916717523)     
X
(yule-walker-algorithm-given-data
X #(71 63 70 88 99 90 110 135 128 154 156 141 131 132 141 104 136 146 124 129)
X 3
X :start 3
X :end 9)
;==> #(0.5519841575716381 -0.2118319667773691 -0.20557341091030015)
;    214.6294761698192
;    #(322.8888888888889 253.32586368772996 224.1000395950692 214.6294761698192)
;    #(0.46415462261986695 -0.3396593834915565 -0.20557341091030015)
;    #(5.629285177048483 19.944421844632704 1.5312965426964307)
;    #(-18.722444902065973 6.621387124383649 -9.561510733986982)
X
(yule-walker-algorithm-given-data
X #(88 99 90 110 135 128)
X 3)
;==> same results as above
|#
X
;-------------------------------------------------------------------------------
(defun yule-walker-algorithm-given-acvs
X       (acvs
X        p
X        &key
X        (AR-coeffs (make-array p))
X        (approximate-MSEs (make-array (1+ p)))
X        (reflection-coeffs (make-array p)))
X  "given
X   [1] acvs (required)
X       ==> a vector containing values of the acvs
X           for a real-valued time series x_t
X           from lag 0 up to (and including) p
X   [2] p (required)
X       ==> autoregressive model order;
X           p should be an integer > 0 and < end - start
X   [3] AR-coeffs (keyword; a vector of length p)
X       <== estimates of AR coefficients phi_{1,p}, ..., phi_{p,p}
X   [4] approximate-MSEs (keyword; a vector of length p+1)
X       <== estimates of innovations variance (mean square errors)
X           for models of order 0, 1, ..., p (the innovations variance
X           for the zeroth order model is the sample variance)
X   [5] reflection-coeffs (keyword; a vector of length p)
X       <== estimates of reflection coefficients phi_{1,1}, ..., phi_{p,p}
uses the Yule-Walker method to estimate the phi's in the model
X   x_t = phi_{1,p} x_{t-1} + ... + phi_{p,p} x_{t-p} + e_t
and returns
X   [1] AR-coeffs, a vector containing phi_{1,p}, ..., phi_{p,p}
X   [2] estimate of the innovations variance, i.e.,
X       the variance of the e_t's; this number is also given in
X       (elt approximate-MSEs p)
X   [3] approximate-MSEs
X   [4] reflection-coeffs
---
Note: see Sections 9.3 and 9.4 of the SAPA book"
X  (assert (and (integerp p) (> p 0) (< p (length acvs))))
X  (setf (elt approximate-MSEs 0) (elt acvs 0))
X  ;;;
X  ;;; Main computational loop - once for each AR order
X  ;;;
X  (let ((scratch (make-array p))
X        top k)
X    (dotimes (k-minus-1 p)
X      (setf k (1+ k-minus-1))
X      (setf top (elt acvs k))
X      (dotimes (j k-minus-1)
X        (decf top (* (elt AR-coeffs j)
X                     (elt acvs (- k-minus-1 j)))))
X      ;;; Calculate k-th order reflection coefficient
X      (setf (elt reflection-coeffs k-minus-1)
X            (setf (elt AR-coeffs k-minus-1)
X                  (/ top (elt approximate-MSEs k-minus-1))))
X      ;;; Update mean square error using approximate method.
X      (setf (elt approximate-MSEs k)
X            (* (elt approximate-MSEs k-minus-1)
X               (- 1.0
X                  (expt (abs (elt AR-coeffs k-minus-1)) 2))))
X      ;;; Use Levinson-Durbin recursions to generate
X      ;;; the remaining (k)-th order AR coefficients
X      (dotimes (j k-minus-1)
X        (setf (elt scratch j) (elt AR-coeffs j)))
X      (dotimes (j k-minus-1)
X        (setf (elt AR-coeffs j)
X              (- (elt scratch j)
X                 (* (elt AR-coeffs k-minus-1)
X                    (elt scratch (- k-minus-1 j 1))))))))
X  (values AR-coeffs
X          (elt approximate-MSEs p)
X          approximate-MSEs
X          reflection-coeffs))
X
#|
(yule-walker-algorithm-given-acvs
X (acvs
X  #(71 63 70 88 99 90 110 135 128 154 156 141 131 132 141 104 136 146 124 129))
X 3)
;==> #(0.7564278899100633 -0.046223877768905255 0.03592061284674314)
;    351.634231421114
;    #(782.6400000000014 352.21671265971605 352.0885278781923 351.634231421114)
;    #(0.7415951139732192 -0.019077139439365634 0.03592061284674314)
;;; ...good agreement with results from yule-walker-algorithm-given-data
X
(yule-walker-algorithm-given-acvs
X (acvs #(88 99 90 110 135 128))
X 3)
;==> #(0.5519841575716382 -0.21183196677736926 -0.2055734109103001)
;    214.62947616981913
;    #(322.88888888888886 253.32586368772988 224.1000395950691 214.62947616981913)
;    #(0.464154622619867 -0.33965938349155667 -0.2055734109103001)
;;; ...good agreement with results from yule-walker-algorithm-given-data
|#
X
;-------------------------------------------------------------------------------
(defun burg-algorithm
X       (time-series
X        p
X        &key
X        (start 0)
X        (end (length time-series))
X        (center-data t)  ;t, nil or value to be subtracted off ...
X        (AR-coeffs (make-array p))
X        (approximate-MSEs (make-array (1+ p)))
X        (reflection-coeffs (make-array p))
X        (forward-prediction-errors (make-array (- end start p)))
X        (backward-prediction-errors (make-array (- end start p)))
X        (exact-MSEs (make-array (1+ p))))
X  "given
X   [1] time-series (required)
X       ==> a vector containing a real-valued or
X           complex-valued time series x_t
X   [2] p (required)
X       ==> autoregressive model order;
X           p should be an integer > 0 and < end - start
X   [3] start (keyword; 0)
X       ==> start index of vector to be used
X   [4] end (keyword; length of the-seq)
X       ==> 1 + end index of vector to be used
X   [5] center-data (keyword; t)
X       ==> if t, subtract sample mean from time series;
X           if a number, subtracts that number from time series;
X           if nil, the effect of this function is to return
X           a copy of the relevant portion of time-series
X   [6] AR-coeffs (keyword; a vector of length p)
X       <== estimates of AR coefficients phi_{1,p}, ..., phi_{p,p}
X   [7] approximate-MSEs (keyword; a vector of length p+1)
X       <== estimates of innovations variance (mean square errors)
X           for models of order 0, 1, ..., p (the innovations variance
X           for the zeroth order model is the sample variance)
X   [8] reflection-coeffs (keyword; a vector of length p)
X       <== estimates of reflection coefficients phi_{1,1}, ..., phi_{p,p}
X   [9] forward-prediction-errors (keyword; a vector of length end - start - p)
X       <== computed forward prediction errors
X  [10] backward-prediction-errors (keyword; a vector of length end - start - p)
X       <== computed backward prediction errors
X  [11] exact-MSEs (keyword; a vector of length p+1)
X       <== another set of estimates of innovations variance
X           for models of order 0, 1, ..., p; these estimates
X           are based on Equation (419a) of the SAPA book
uses Burg's algorithm to estimate the phi's in the model
X   x_t = phi_{1,p} x_{t-1} + ... + phi_{p,p} x_{t-p} + e_t
and returns
X   [1] AR-coeffs, a vector containing phi_{1,p}, ..., phi_{p,p}
X   [2] estimate of the innovations variance, i.e.,
X       the variance of the e_t's; this number is also given in
X       (elt approximate-MSEs p)
X   [3] approximate-MSEs
X   [4] reflection-coeffs
X   [5] forward-prediction-errors
X   [6] backward-prediction-errors
X   [7] exact-MSEs
---
Note: see Section 9.5 of the SAPA book"
X  (let ((N-ts (- end start)))
X    (assert (and (integerp p) (> p 0) (< p N-ts)))
X    (let* ((forward-scratch (center&taper-time-series
X                             time-series
X                             :center-data center-data
X                             :start start
X                             :end end))
X           (backward-scratch (copy-seq forward-scratch))
X           (sh2-old (float (/ (sum-of-squares forward-scratch)
X                              N-ts)))
X           (mseb-old sh2-old)
X           (scratch (make-array p))
X           top bottom k)
X      ;(declare (dynamic-extent forward-scratch backward-scratch))
X      ;;; Set zeroth order mean square prediction error
X      (setf (aref approximate-MSEs 0)
X            (setf (aref exact-MSEs 0) sh2-old))
X      ;;;
X      ;;; Main computational loop - once for each AR order
X      ;;;
X      (dotimes (k-minus-1 p)
X        (setf k (1+ k-minus-1)
X              top 0.0
X              bottom 0.0)
X        (dotimes
X          (j (- N-ts k))
X          (incf top 
X                (* (aref forward-scratch (1+ j))
X                   (conjugate (aref backward-scratch j))))
X          (incf bottom
X                (+ (expt (abs (aref forward-scratch (1+ j))) 2)
X                   (expt (abs (aref backward-scratch j)) 2))))
X        ;;; Calculate k-th order reflection coefficient
X        (setf (aref reflection-coeffs k-minus-1)
X              (setf (aref AR-coeffs k-minus-1)
X                    (/ (* 2.0 top) bottom)))
X        ;;; Update mean square error using:
X        ;;; (1) approximate method
X        (setf (aref approximate-MSEs k)
X              (setf sh2-old
X                    (* sh2-old
X                       (- 1.0
X                          (expt (abs (aref AR-coeffs k-minus-1)) 2)))))
X        ;;; (2) exact method
X        (setf (aref exact-MSEs k)
X              (setf mseb-old
X                    (* (/ (- (* mseb-old (1+ (- N-ts k)))
X                             (/ (+ (expt (abs (aref forward-scratch 0)) 2)
X                                   (expt (abs (aref backward-scratch (- N-ts k))) 2))
X                                2.0))
X                          (- N-ts k))
X                       (- 1.0
X                          (expt (abs (aref AR-coeffs k-minus-1)) 2)))))
X        ;;; Update prediction errors
X        (dotimes
X          (j (- N-ts k))
X          (setf (aref forward-scratch j)
X                (- (aref forward-scratch (1+ j))
X                   (* (aref AR-coeffs k-minus-1)
X                      (aref backward-scratch j))))
X          (decf (aref backward-scratch j)
X                (* (conjugate (aref AR-coeffs k-minus-1))
X                   (aref forward-scratch (1+ j)))))
X        ;;; Use Levinson-Durbin recursions to generate
X        ;;; the remaining (k)-th order AR coefficients
X        (cond
X         ((> k 1)
X          (dotimes
X            (j k-minus-1)
X            (setf (aref scratch j) (aref AR-coeffs j)))
X          (dotimes
X            (j k-minus-1)
X            (setf (aref AR-coeffs j)
X                  (- (aref scratch j)
X                     (* (aref AR-coeffs k-minus-1)
X                        (conjugate (aref scratch (- k-minus-1 j 1))))))))))
X      (values AR-coeffs
X              (elt approximate-MSEs p)
X              approximate-MSEs
X              reflection-coeffs
X              (copy-vector forward-scratch forward-prediction-errors
X                           :end (- N-ts p))
X              (copy-vector backward-scratch backward-prediction-errors
X                           :end (- N-ts p))
X              exact-MSEs))))
X
#|
(burg-algorithm
X #(71 63 70 88 99 90 110 135 128 154 156 141 131 132 141 104 136 146 124 129)
X 3)
;==> #(0.6990553262846771 -0.005832588617464066 0.15582353261696208)
;    271.77271305830277
;    #(782.64 281.6806619464492 278.53583421681344 271.77271305830277)
;    #(0.8000556894184591 0.10566226444248338 0.15582353261696208)
;    #(10.648141558530694 10.352562066664447 -7.322824655671381 16.228008168577496 25.480349486540202 2.523029895325208 30.445761242415312 10.333906323267405 -4.82179229692915 -8.375709073465082 -0.7842917051141942 9.795680071680922 -31.93174995009409 25.82996688737918 11.841978873870683 -11.196460846389641 4.254729174305107)
;    #(-4.066643091642674 -18.56910263930764 -22.68552824408707 -15.544100783114937 -2.031539389627614 -23.776066471565816 -25.344689597046013 4.388697925766872 -18.43772239114538 7.634913453192931 19.906593928671484 10.500567986583063 5.619476664683157 -4.874180094467894 28.619274487654387 -27.264052349707473 -3.1620402252232624)
;    #(782.64 274.8400818422606 278.4733264361111 272.60343018395645)
X
(burg-algorithm
X #(71 63 70 88 99 90 110 135 128 154 156 141 131 132 141 104 136 146 124 129)
X 3
X :start 3
X :end 9)
;==> #(0.7110790464963852 -0.17650146677307085 -0.2576098242725891)
;    168.87383595887832
;    #(322.8888888888889 212.4012089094338 180.87736848843608 168.87383595887832)
;    #(0.5849656667871341 -0.3852485990185058 -0.2576098242725891)
;    #(7.817702402342421 19.841316338455563 -3.7241189069459537)
;    #(-16.50310608308572 10.866880277657888 -9.745432752851167)
;    #(322.8888888888889 202.24161908203828 190.11187075352709 159.01307057103358)
X
(burg-algorithm
X #(88 99 90 110 135 128)
X 3)
;==> same results as above
|#
X
;-------------------------------------------------------------------------------
;;;  adapted from Fortran subroutine MODCOVAR, , Marple, 1987 --
;;;  this works with both real and complex-valued time series
(defun fast-forward-backward-ls
X       (time-series
X        p
X        &key
X        (start 0)
X        (end (length time-series))
X        (center-data t)
X        (AR-coeffs (make-array p))
X        (pseudo-reflection-coeffs (make-array p)))
X  "given
X   [1] time-series (required)
X       ==> a vector containing a real-valued or
X           complex-valued time series x_t
X   [2] p (required)
X       ==> autoregressive model order;
X           p should be an integer > 0 and < end - start
X   [3] start (keyword; 0)
X       ==> start index of vector to be used
X   [4] end (keyword; length of the-seq)
X       ==> 1 + end index of vector to be used
X   [5] center-data (keyword; t)
X       ==> if t, subtract sample mean from time series;
X           if a number, subtracts that number from time series;
X           if nil, the effect of this function is to return
X           a copy of the relevant portion of time-series
X   [6] AR-coeffs (keyword; a vector of length p)
X       <== estimates of AR coefficients phi_{1,p}, ..., phi_{p,p}
X   [7] pseudo-reflection-coeffs (keyword; a vector of length p)
X       <== sequence of reflection coefficients computed
X           stepwise by the algorithm, but these do not correspond
X           to the estimated phi_{k,p} (even if the estimated
X           AR model is stationary)
uses the forward/backward least squares method
to estimate the phi's in the model
X   x_t = phi_{1,p} x_{t-1} + ... + phi_{p,p} x_{t-p} + e_t
and returns
X   [1] AR-coeffs, a vector containing phi_{1,p}, ..., phi_{p,p}
X   [2] estimate of the innovations variance, i.e.,
X       the variance of the e_t's
X   [3] pseudo-reflection-coeffs
---
Note: see Section 9.7 of the SAPA book;
X      this function is an adaptation of
X      the Fortran routine modcovar, pp. 258-60,
X      ``Digital Spectral Analysis with Applications''
X      by Marple, 1987"
X  (let ((N-ts (- end start)))
X    (assert (and (integerp p) (> p 0) (< p N-ts)))
X    (let* ((centered-ts (center&taper-time-series
X                         time-series
X                         :center-data center-data
X                         :start start
X                         :end end))
X           (N-ts-1 (1- N-ts))
X           (r1 (let ((temp 0.0))
X                 (dotimes (k (- N-ts 2) temp)
X                   (incf temp
X                         (* 2.0 (expt (abs (aref centered-ts (1+ k))) 2))))))
X           (r2 (expt (abs (aref centered-ts 0)) 2))
X           (r3 (expt (abs (aref centered-ts N-ts-1)) 2))
X           (r4 (/ (+ r1 (* 2.0 (+ r2 r3)))))
X           r5
X           (estimated-mse (+ r1 r2 r3))
X           (delta (- 1.0 (* r2 r4)))
X           (gamma (- 1.0 (* r3 r4)))
X           (lambda (* (conjugate (* (aref centered-ts 0) (aref centered-ts N-ts-1)))
X                      r4))
X           (c (make-array (1+ p)))
X           (d (make-array (1+ p)))
X           (r (make-array p)))
X      (setf (aref c 0) (* (aref centered-ts N-ts-1) r4)
X            (aref d 0) (* (conjugate (aref centered-ts 0)) r4))
X      (cond
X       ((zerop p) (values
X                   (flip-sign-of-coefficients! AR-coeffs)
X                   (/ (+ (* 0.5 r1) r2 r3) N-ts)))
X       (t
X        ;;;
X        ;;; Main loop
X        ;;;
X        (do ((m 1 (+ m 1))
X             (save1 0.0 0.0)
X             save2 save3 save4
X             m-k-1 m-k-2 theta psi xi ef eb c1 c2 c3 c4)
X            (nil)
X          (dotimes (k-m (- N-ts m))
X            (incf save1 (* (aref centered-ts (+ k-m m))
X                           (conjugate (aref centered-ts k-m)))))
X          (setf save1 (* 2.0 save1)
X                (aref r (1- m)) (conjugate save1)
X                theta (* (aref centered-ts N-ts-1) (aref d 0))
X                psi (* (aref centered-ts N-ts-1) (aref c 0))
X                xi (* (conjugate (aref centered-ts 0)) (aref d 0)))
X          (if (> m 1) (dotimes (k (1- m))
X                        (incf theta (* (aref centered-ts (- N-ts k 2))
X                                       (aref d (1+ k))))
X                        (incf psi (* (aref centered-ts (- N-ts k 2))
X                                     (aref c (1+ k))))
X                        (incf xi (* (conjugate (aref centered-ts (1+ k)))
X                                    (aref d (1+ k))))
X                        (decf (aref r k) (+ (* (aref centered-ts (- N-ts m))
X                                               (conjugate
X                                                (aref centered-ts (+ k 1 (- N-ts m)))))
X                                            (* (conjugate
X                                                (aref centered-ts (1- m)))
X                                               (aref centered-ts (- m k 2)))))
X                        (incf save1 (* (conjugate (aref r k))
X                                       (aref AR-coeffs (- m k 2))))))
X          ;;; Note: Marple (1987) claims that c1 is less than unity
X          ;;;       in magnitude (true) and that it is in fact
X          ;;;       part of the reflection coefficient sequence
X          ;;;       for the forward/backward least squares method
X          ;;;       (false -- else f/b ls would be truly order
X          ;;;       recursive as, for example, Burg is).
X          (setf c1 (/ (- save1) estimated-mse)
X                (aref AR-coeffs (1- m)) c1
X                (aref pseudo-reflection-coeffs (1- m)) c1
X                estimated-mse (* estimated-mse
X                                 (- 1.0 (expt (abs c1) 2))))
X          (if (> m 1) (dotimes (k (truncate m 2))
X                        (setf m-k-2 (- m k 2)
X                              save1 (aref AR-coeffs k)
X                              (aref AR-coeffs k)
X                              (+ save1
X                                 (* c1 (conjugate
X                                        (aref AR-coeffs m-k-2)))))
X                        (if (/= k m-k-2)
X                          (incf (aref AR-coeffs m-k-2)
X                                (* c1 (conjugate save1))))))
X          (if (= m p) (return (values
X                               (flip-sign-of-coefficients! AR-coeffs)
X                               (* 0.5 (/ estimated-mse (float (- N-ts m))))
X                               pseudo-reflection-coeffs)))
X          ;;;
X          ;;; Time update of c, d vectors and gamma, delta, lambda scalars
X          ;;;
X          (setf r1 (/ (- (* delta gamma)
X                         (expt (abs lambda) 2)))
X                c1 (* (+ (* theta (conjugate lambda)) (* psi delta))
X                      r1)
X                c2  (* (+ (* psi lambda) (* theta gamma))
X                       r1)
X                c3 (* (+ (* xi (conjugate lambda)) (* theta delta))
X                      r1)
X                c4 (* (+ (* theta lambda) (* xi gamma))
X                      r1))
X          (dotimes (k (1+ (truncate (1- m) 2)))
X            (setf m-k-1 (- m k 1)
X                  save1 (conjugate (aref c k))
X                  save2 (conjugate (aref d k))
X                  save3 (conjugate (aref c m-k-1))
X                  save4 (conjugate (aref d m-k-1)))
X            (incf (aref c k) (+ (* c1 save3) (* c2 save4)))
X            (incf (aref d k) (+ (* c3 save3) (* c4 save4)))
X            (when (/= k m-k-1)
X              (incf (aref c m-k-1) (+ (* c1 save1) (* c2 save2)))
X              (incf (aref d m-k-1) (+ (* c3 save1) (* c4 save2)))))
X          (setf r2 (expt (abs psi) 2)
X                r3 (expt (abs theta) 2)
X                r4 (expt (abs xi) 2)
X                r5 (- gamma (* (+ (* r2 delta)
X                                  (* r3 gamma)
X                                  (* 2.0 (realpart
X                                          (* psi lambda (conjugate theta)))))
X                               r1))
X                r2 (- delta (* (+ (* r3 delta)
X                                  (* r4 gamma)
X                                  (* 2.0 (realpart
X                                          (* theta lambda (conjugate xi)))))
X                               r1))
X                gamma r5
X                delta r2)
X          (incf lambda (+ (* c3 (conjugate psi)) (* c4 (conjugate theta))))
X          (if  (<= estimated-mse 0.0)
X            (error "estimated mse <= 0.0 (~F)" estimated-mse))
X          (if (not (and (plusp delta) (<= delta 1.0)
X                        (plusp gamma) (<= gamma 1.0)))
X            (error "delta' and gamma' are not in the range 0 to 1 (~F,~F)"
X                   delta gamma))
X          ;;;
X          ;;; Time update of A vector; order updates of c, d vectors
X          ;;; and gamma, delta, lambda scalars
X          ;;;
X          (setf r1 (/ estimated-mse)
X                r2 (/ (- (* delta gamma)
X                         (expt (abs lambda) 2)))
X                ef (aref centered-ts m)
X                eb (aref centered-ts (- N-ts m 1)))
X          (dotimes (k m)
X            (incf ef (* (aref AR-coeffs k)
X                        (aref centered-ts (- m k 1))))
X            (incf eb (* (conjugate (aref AR-coeffs k))
X                        (aref centered-ts (+ (- N-ts m) k)))))
X          (setf c1 (* eb r1)
X                c2 (* (conjugate ef) r1)
X                c3 (* (+ (* (conjugate eb) delta) (* ef lambda))
X                      r2)
X                c4 (* (+ (* ef gamma) (conjugate (* eb lambda)))
X                      r2))
X          (dotimes (k m)
X            (setf m-k-1 (- m k 1)
X                  save1 (aref AR-coeffs m-k-1)
X                  (aref AR-coeffs m-k-1) (+ save1
X                                                   (* c3 (aref c m-k-1))
X                                                   (* c4 (aref d m-k-1)))
X                  (aref c (1+ m-k-1)) (+ (aref c m-k-1) (* c1 save1))
X                  (aref d (1+ m-k-1)) (+ (aref d m-k-1) (* c2 save1))))
X          (setf (aref c 0) c1
X                (aref d 0) c2
X                r3 (expt (abs eb) 2)
X                r4 (expt (abs ef) 2))
X          (decf estimated-mse (* r2 (+ (* r3 delta) (* r4 gamma)
X                                       (* 2.0 (realpart (* ef eb lambda))))))
X          (decf delta (* r4 r1))
X          (decf gamma (* r3 r1))
X          (incf lambda (* (conjugate (* ef eb)) r1))
X          (if  (<= estimated-mse 0.0)
X            (error "estimated mse' <= 0.0 (~F)" estimated-mse))
X          (if (not (and (plusp delta) (<= delta 1.0)
X                        (plusp gamma) (<= gamma 1.0)))
X            (error "delta' and gamma' are not in the range 0 to 1 (~F,~F)"
X                   delta gamma))))))))
X
#|
(fast-forward-backward-ls
X #(71 63 70 88 99 90 110 135 128 154 156 141 131 132 141 104 136 146 124 129)
X 3)
;==> #(0.6736988680045488 -0.017857067378710748 0.1515182334143321)
;    271.5712900083261
;    #(-0.8000556894184591 -0.10369838100465911 -0.1515182334143321)
X
(fast-forward-backward-ls
X #(71 63 70 88 99 90 110 135 128 154 156 141 131 132 141 104 136 146 124 129)
X 3
X :start 3
X :end 9)
;==> #(0.5310692482787671 -0.33739103356904027 -0.34747481540511216)
;    137.56127670899096
;    #(-0.5849656667871341 0.406220671240301 0.34747481540511216)
X
(fast-forward-backward-ls
X #(88 99 90 110 135 128)
X 3)
;==> same results as above
|#
X
;-------------------------------------------------------------------------------
;;; This implements the forward/backward least squares algorithm
;;; a straight-forward manner, but it is considerably slower than
;;; fast-forward-backward-ls
(defun forward-backward-ls
X       (time-series
X        p
X        &key
X        (start 0)
X        (end (length time-series))
X        (center-data t)
X        (AR-coeffs (make-array p)))
X  "given
X   [1] time-series (required)
X       ==> a vector containing a real-valued or
X           complex-valued time series x_t
X   [2] p (required)
X       ==> autoregressive model order;
X           p should be an integer > 0 and < end - start
X   [3] start (keyword; 0)
X       ==> start index of vector to be used
X   [4] end (keyword; length of the-seq)
X       ==> 1 + end index of vector to be used
X   [5] center-data (keyword; t)
X       ==> if t, subtract sample mean from time series;
X           if a number, subtracts that number from time series;
X           if nil, the effect of this function is to return
X           a copy of the relevant portion of time-series
X   [6] AR-coeffs (keyword; a vector of length p)
X       <== estimates of AR coefficients phi_{1,p}, ..., phi_{p,p}
uses the forward/backward least squares method
to estimate the phi's in the model
X   x_t = phi_{1,p} x_{t-1} + ... + phi_{p,p} x_{t-p} + e_t
and returns
X   [1] AR-coeffs, a vector containing phi_{1,p}, ..., phi_{p,p}
X   [2] estimate of the innovations variance, i.e.,
X       the variance of the e_t's
---
Note: see Section 9.7 of the SAPA book;
X      this function is an adaptation of
X      the Fortran routine covmcov, pp. 262-4,
X      ``Modern Spectrum Estimation: Theory and Application''
X      by Kay, 1988"
X  (let ((N-ts (- end start)))
X    (assert (and (integerp p) (> p 0) (< p N-ts)))
X    (let* ((centered-ts (center&taper-time-series
X                         time-series
X                         :center-data center-data
X                         :start start
X                         :end end))
X           (N-ts-p (- N-ts p))
X           (p-1 (1- p))
X           (CC (make-array (/ (* p (1+ p)) 2)))
X           (L 0)
X           ;;; if zero is set to 0 rather than 0.0
X           ;;; and if centered-ts contains integers,
X           ;;; this algorithm will use rational arithmetic
X           (zero 0.0))
X      ;(declare (dynamic-extent CC))
X      ;;; Compute autocorrelation estimates and insert them in matrix CC
X      ;;; (stored in linearized form); see Kay's Equations (7.8) and (7.22)
X      (dotimes (k p)
X        (dotimes (j (1+ k))
X          (setf (svref CC L) zero)
X          (dotimes (i N-ts-p)
X            (incf (svref CC L)
X                  (+ (* (conjugate (aref centered-ts (+ p-1 (- i j))))
X                        (aref centered-ts (+ p-1 (- i k))))
X                     (* (aref centered-ts (1+ (+ i j)))
X                        (conjugate (aref centered-ts (1+ (+ i k))))
X                        ))))
X          (incf L)))
X      ;;; Compute right-hand vector of covariance equations;
X      ;;; see Kay's Equation (7.8)
X      (dotimes (j p)
X        (setf (aref AR-coeffs j) zero)
X        (dotimes (i N-ts-p)
X          (decf (aref AR-coeffs j)
X                (+ (* (conjugate (aref centered-ts (+ p-1 (- i j))))
X                      (aref centered-ts (+ p i)))
X                   (* (aref centered-ts (1+ (+ i j)))
X                      (conjugate (aref centered-ts i)))))))
X      ;;; Compute AR filter parameter estimate
X      (cholesky! CC AR-coeffs)
X      ;;; Compute estimate of innovations variance (Kay's Equation (7.9))
X      (let ((sum zero)
X            C)
X        (dotimes (k (1+ p) (values
X                            (flip-sign-of-coefficients! AR-coeffs)
X                            (/ (realpart sum) (* 2 N-ts-p))))
X          (setf C zero)
X          (dotimes (i N-ts-p)
X            (incf C (+ (* (conjugate (aref centered-ts (+ p i)))
X                          (aref centered-ts (- (+ p i) k)))
X                       (* (aref centered-ts i)
X                          (conjugate (aref centered-ts (+ i k)))))))
X          (incf sum (if (zerop k)
X                      C
X                      (* C (aref AR-coeffs (1- k))))))))))
X
#|
(forward-backward-ls
X #(71 63 70 88 99 90 110 135 128 154 156 141 131 132 141 104 136 146 124 129)
X 3)
;==> #(0.6736988680045498 -0.017857067378712857 0.15151823341433324)
;    271.571290008326
;;; agrees well with fast-forward-backward-ls ...
X
(forward-backward-ls
X #(71 63 70 88 99 90 110 135 128 154 156 141 131 132 141 104 136 146 124 129)
X 3
X :start 3
X :end 9)
;==> #(0.5310692482787671 -0.33739103356903977 -0.3474748154051125)
;    137.5612767089909
;;; agrees well with fast-forward-backward-ls ...
X
(forward-backward-ls
X #(88 99 90 110 135 128)
X 3)
;==> same results as above
|#
X
;-------------------------------------------------------------------------------
(defun forward-ls
X       (time-series
X        p
X        &key
X        (start 0)
X        (end (length time-series))
X        (center-data t)
X        (AR-coeffs (make-array p)))
X  "given
X   [1] time-series (required)
X       ==> a vector containing a real-valued or
X           complex-valued time series x_t
X   [2] p (required)
X       ==> autoregressive model order;
X           p should be an integer > 0 and < end - start
X   [3] start (keyword; 0)
X       ==> start index of vector to be used
X   [4] end (keyword; length of the-seq)
X       ==> 1 + end index of vector to be used
X   [5] center-data (keyword; t)
X       ==> if t, subtract sample mean from time series;
X           if a number, subtracts that number from time series;
X           if nil, the effect of this function is to return
X           a copy of the relevant portion of time-series
X   [6] AR-coeffs (keyword; a vector of length p)
X       <== estimates of AR coefficients phi_{1,p}, ..., phi_{p,p}
uses the forward least squares method
to estimate the phi's in the model
X   x_t = phi_{1,p} x_{t-1} + ... + phi_{p,p} x_{t-p} + e_t
and returns
X   [1] AR-coeffs, a vector containing phi_{1,p}, ..., phi_{p,p}
X   [2] estimate of the innovations variance, i.e.,
X       the variance of the e_t's
---
Note: see Section 9.7 of the SAPA book;
X      this function is an adaptation of
X      the Fortran routine covmcov, pp. 262-4,
X      ``Modern Spectrum Estimation: Theory and Application''
X      by Kay, 1988"
X  (let ((N-ts (- end start)))
X    (assert (and (integerp p) (> p 0) (< p N-ts)))
X    (let* ((centered-ts (center&taper-time-series
X                         time-series
X                         :center-data center-data
X                         :start start
X                         :end end))
X           (N-ts-p (- N-ts p))
X           (p-1 (1- p))
X           (CC (make-array (/ (* p (1+ p)) 2)))
X           (L 0))
X      ;(declare (dynamic-extent CC))
X      ;;; Compute autocorrelation estimates and insert them in matrix CC
X      ;;; (stored in linearized form); see Kay's Equation (7.8)
X      (dotimes (k p)
X        (dotimes (j (1+ k))
X          (setf (svref CC L) 0.0)
X          (dotimes (i N-ts-p)
X            (incf (svref CC L) (* (conjugate (aref centered-ts (+ p-1 (- i j))))
X                                  (aref centered-ts (+ p-1 (- i k))))))
X          (incf L)))
X      ;;; Compute right-hand vector of covariance equations;
X      ;;; see Kay's Equation (7.8)
X      (dotimes (j p)
X        (setf (aref AR-coeffs j) 0.0)
X        (dotimes (i N-ts-p)
X          (decf (aref AR-coeffs j)
X                (* (conjugate (aref centered-ts (+ p-1 (- i j))))
X                   (aref centered-ts (+ p i))))))
X      ;;; Compute AR filter parameter estimate
X      (cholesky! CC AR-coeffs)
X      ;;; Compute estimate of innovations variance (Kay's Equation (7.9))
X      (let ((sum 0.0)
X            C)
X        (dotimes (k (1+ p) (values
X                            (flip-sign-of-coefficients! AR-coeffs)
X                            (/ (realpart sum) N-ts-p)))
X          (setf C 0.0)
X          (dotimes (i N-ts-p)
X            (incf C (* (conjugate (aref centered-ts (+ p i)))
X                       (aref centered-ts (- (+ p i) k)))))
X          (incf sum (if (zerop k)
X                      C
X                      (* C (aref AR-coeffs (1- k))))))))))
X
#|
(forward-ls
X #(71 63 70 88 99 90 110 135 128 154 156 141 131 132 141 104 136 146 124 129)
X 3)
;==> #(0.5055761990420168 -0.016273987434470755 0.17082594780884994)
;    235.7359874382611
X
(forward-ls
X #(71 63 70 88 99 90 110 135 128 154 156 141 131 132 141 104 136 146 124 129)
X 3
X :start 3
X :end 9)
;==> #(0.7287366959740863 -1.3206848681166132 -0.13280888477556643)
;    -1.4210854715202004E-14  ;;; opps -- should be positive!
X
(forward-ls
X #(88 99 90 110 135 128)
X 3)
;==> same results as above
|#
X
;-------------------------------------------------------------------------------
(defun two-step-Burg-algorithm
X       (time-series
X        p
X        &key
X        (start 0)
X        (end (length time-series))
X        (center-data t)  ;t, nil or value to be subtracted off ...
X        (AR-coeffs (make-array p))
X        (approximate-MSEs (make-array (1+ p)))
X        (reflection-coeffs (make-array p))
X        (forward-prediction-errors (make-array (- end start p)))
X        (backward-prediction-errors (make-array (- end start p)))
X        (exact-MSEs (make-array (1+ p))))
X  "given
X   [1] time-series (required)
X       ==> a vector containing a real-valued or
X           complex-valued time series x_t
X   [2] p (required)
X       ==> autoregressive model order;
X           p should be an EVEN integer > 0 and < end - start
X   [3] start (keyword; 0)
X       ==> start index of vector to be used
X   [4] end (keyword; length of the-seq)
X       ==> 1 + end index of vector to be used
X   [5] center-data (keyword; t)
X       ==> if t, subtract sample mean from time series;
X           if a number, subtracts that number from time series;
X           if nil, the effect of this function is to return
X           a copy of the relevant portion of time-series
X   [6] AR-coeffs (keyword; a vector of length p)
X       <== estimates of AR coefficients phi_{1,p}, ..., phi_{p,p}
X   [7] approximate-MSEs (keyword; a vector of length p+1)
X       <== estimates of innovations variance (mean square errors)
X           for models of order 0, 1, ..., p (the innovations variance
X           for the zeroth order model is the sample variance)
X   [8] reflection-coeffs (keyword; a vector of length p)
X       <== estimates of reflection coefficients phi_{1,1}, ..., phi_{p,p}
X   [9] forward-prediction-errors (keyword; a vector of length end - start - p)
X       <== computed forward prediction errors
X  [10] backward-prediction-errors (keyword; a vector of length end - start - p)
X       <== computed backward prediction errors
X  [11] exact-MSEs (keyword; a vector of length p+1)
X       <== another set of estimates of innovations variance
X           for models of order 0, 1, ..., p; these estimates
X           are based on Equation (419a) of the SAPA book
uses the two step Burg algorithm to estimate the phi's in the model
X   x_t = phi_{1,p} x_{t-1} + ... + phi_{p,p} x_{t-p} + e_t
and returns
X   [1] AR-coeffs, a vector containing phi_{1,p}, ..., phi_{p,p}
X   [2] estimate of the innovations variance, i.e.,
X       the variance of the e_t's; this number is also given in
X       (elt approximate-MSEs p)
X   [3] approximate-MSEs
X   [4] reflection-coeffs
X   [5] forward-prediction-errors
X   [6] backward-prediction-errors
X   [7] exact-MSEs
---
Note: see Section 9.5 of the SAPA book"
X  (let ((N-ts (- end start)))
X    (assert (and (evenp p) (integerp p) (> p 0) (< p N-ts)))
X    (let* ((forward-scratch (center&taper-time-series
X                             time-series
X                             :center-data center-data
X                             :start start
X                             :end end))
X           (backward-scratch (copy-seq forward-scratch)))
X      ;(declare (dynamic-extent forward-scratch backward-scratch))
X      ;;; Set zeroth order mean square prediction error
X      (setf (aref approximate-MSEs 0)
X            (setf (aref exact-MSEs 0)
X                  (float (/ (sum-of-squares forward-scratch)
X                            N-ts))))
X      (dotimes (i (/ p 2))
X        (multiple-value-bind (ss ss-m-1 ref-coeff-m-1)
X                             (one-iteration-of-two-step-Burg
X                              (* i 2) 1
X                              N-ts
X                              forward-scratch
X                              backward-scratch
X                              AR-coeffs
X                              1.0)
X          (let* ((mse-index-k (* (1+ i) 2))
X                 (mse-index-k-minus-1 (1- mse-index-k))
X                 (rc-index-k mse-index-k-minus-1))
X            (setf (aref exact-MSEs mse-index-k-minus-1) ss-m-1)
X            (setf (aref exact-MSEs mse-index-k) ss)
X            (setf (aref approximate-MSEs mse-index-k-minus-1)
X                  (* (aref approximate-MSEs (1- mse-index-k-minus-1))
X                     (- 1.0 (expt ref-coeff-m-1 2))))
X            (setf (aref approximate-MSEs mse-index-k)
X                  (* (aref approximate-MSEs mse-index-k-minus-1)
X                     (- 1.0 (expt (aref AR-coeffs rc-index-k) 2))))
X            (setf (aref reflection-coeffs (1- rc-index-k)) ref-coeff-m-1)
X            (setf (aref reflection-coeffs rc-index-k)
X                  (aref AR-coeffs rc-index-k)))))
X      (values AR-coeffs
X              (elt approximate-MSEs p)
X              approximate-MSEs
X              reflection-coeffs
X              (copy-vector forward-scratch forward-prediction-errors
X                           :start p)
X              (copy-vector backward-scratch backward-prediction-errors
X                           :end (- N-ts p))
X              exact-MSEs))))
X
#|
(two-step-burg-algorithm
X #(71 63 70 88 99 90 110 135 128 154 156 141 131 132 141 104 136 146 124 129)
X 4)
;==> #(0.7475191644352183 -0.01583855174510637 0.46440844342490983 -0.4428652096133231)
;    245.25050973984762
;    #(782.64 317.1313317877149 313.7211062417939 305.0871233404041 245.25050973984762)
;    #(0.771228137382831 0.10369838100465757 0.16589516290564163 -0.4428652096133231)
;    #(7.541189677934279 -16.19020798032216 5.452393054435857 18.222543695391344 1.9025289159837546 20.25717110563984 -0.04210392524786144 -1.8028505644380866 -15.733061990760522 3.0902299405116476 19.03618232922566 -26.674545308896946 28.233150201353947 9.98969974364716 2.190241048379951 -7.4530347618842505)
;    #(-0.9808164276349656 -23.02263641146349 -16.2667771190106 -4.848573521669076 -1.5143975597659964 -10.303462500921928 -20.291200583839142 2.7814408815304876 -21.08490573505456 8.26942730376083 24.845207695933038 -3.2294508544868155 17.520376061131927 0.7742590731821393 23.552182767386412 -24.778733173674873)
;    #(782.64 275.4746826295362 277.97870057136765 271.85461785540076 229.6941141226679)
X
;;; reconstruction of forward prediction errors ...
(filter-time-series 
X (center&taper-time-series
X  #(71 63 70 88 99 90 110 135 128 154 156 141 131 132 141 104 136 146 124 129))
X (ar-coeffs->prewhitening-filter
X  #(0.7475191644352183 -0.01583855174510637 0.46440844342490983 -0.4428652096133231))
X )
;==> #(7.541189677934277 -16.190207980322164 5.452393054435857 18.222543695391344 1.9025289159837548 20.25717110563984 -0.04210392524786499 -1.8028505644380886 -15.73306199076052 3.09022994051165 19.03618232922566 -26.674545308896946 28.233150201353947 9.989699743647158 2.1902410483799546 -7.453034761884252)
;    16
;;; good agreement with what two-step-burg-algorithm gave ...
X
;;; reconstruction of backward prediction errors ...
(filter-time-series 
X (center&taper-time-series
X  #(71 63 70 88 99 90 110 135 128 154 156 141 131 132 141 104 136 146 124 129))
X (reverse
X  (ar-coeffs->prewhitening-filter
X   #(0.7475191644352183 -0.01583855174510637 0.46440844342490983 -0.4428652096133231)))
X )
;==> #(-0.9808164276349629 -23.022636411463495 -16.266777119010605 -4.8485735216690795 -1.514397559765996 -10.303462500921928 -20.291200583839146 2.7814408815304876 -21.084905735054566 8.26942730376083 24.845207695933038 -3.2294508544868163 17.520376061131927 0.7742590731821384 23.55218276738641 -24.778733173674876)
;    16
;;; good agreement with what two-step-burg-algorithm gave ...
X
;;; second test ...
(two-step-burg-algorithm
X #(71 63 70 88 99 90 110 135 128 154 156 141 131 132 141 104 136 146 124 129)
X 4
X :start 3
X :end 9)
;==> #(0.7203613053679895 -0.836626599020122 0.6819467895912388 -0.9407536937974277)
;    23.9213100745261
;    #(322.8888888888889 249.5010816657653 208.3296023555775 208.04307310163574 23.9213100745261)
;    #(0.47674418604651164 -0.4062206712403009 0.0370859144408096 -0.9407536937974277)
;    #(-2.6359115620116818 5.573399522283869)
;    #(-4.997928279988322 5.5839098523630515)
;    #(322.8888888888889 205.84237425635473 183.62625741893694 169.17806158288212 23.542537082241278)
X
(two-step-burg-algorithm
X #(88 99 90 110 135 128)
X 4)
;==> same results as above
X
;;; reconstruction of exact MSE ...
(/ (+ (sum-of-squares #(-2.6359115620116818 5.573399522283869))
X      (sum-of-squares #(-4.997928279988322 5.5839098523630515)))
X   4)
;==> 23.542537082241275
;;; agrees with last element of exact-MSEs above (23.542537082241278) ...
X
;;; reconstruction of approximate MSE ...
(* (sample-variance
X    #(88 99 90 110 135 128))
X   (- 1.0 (expt 0.47674418604651164 2))
X   (- 1.0 (expt -0.4062206712403009 2))
X   (- 1.0 (expt 0.0370859144408096 2))
X   (- 1.0 (expt -0.9407536937974277 2)))
;==> 23.9213100745261
;;; agrees prefectly with above computation (23.9213100745261) ...
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  ar-coeffs->acvs
;;;                 ar-coeffs->variance
;;;                 ar-coeffs->reflection-coeffs
;;;                 ar-coeffs->prewhitening-filter
;;;                 reflection-coeffs->ar-coeffs
;;;                 reflection-coeffs->variance
;;;  take the AR coefficients or the reflection coefficients and convert 
;;;  them to some other quantity of interest.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun ar-coeffs->acvs
X       (AR-coeffs
X        variance
X        max-lag
X        &key
X        (process-variance-p t)
X        (list-of-lower-order-phi nil)
X        (list-of-lower-order-pev nil)
X        (result (make-array (1+ max-lag))))
X  "given
X   [1] AR-coeffs (required)
X       ==> vector of length p with real or
X           complex-valued AR coefficients phi_{j,p}
X           for an AR(p) model of order p:
X           x_t = phi_{1,p} * x_{t-1} + phi_{2,p} * x_{t-2}
X                 + ... + phi_{p,p} * x_{t-p} + e_t
X   [2] variance (required)
X       ==> process variance (if process-variance-p is true)
X           or innovations variance (if nil)
X   [3] max-lag (required)
X       ==> acvs to be computed out to this lag
X   [4] process-variance-p (keyword; t)
X       ==> if true, required variance parameter
X            is taken to be process variance;
X           otherwise, it is taken to be the innovations
X           variance (both variance and process-variance-p
X           are ignored if the next two keyword parameters
X           are supplied)
X   [5] list-of-lower-order-phi (keyword; nil)
X       ==> if supplied, bypasses call to step-down-Levinson-Durbin-recursions,
X           in which case list-of-lower-order-pev must be supplied also
X   [6] list-of-lower-order-pev (keyword; nil)
X       ==> if supplied, bypasses call to step-down-Levinson-Durbin-recursions
X   [7] result (keyword; vector of size max-lag + 1)
X       <== vector to hold acvs for lags 0, 1, ..., max-lag
returns
X   [1] result, i.e., the acvs
---
Note: see Section 9.4 of the SAPA book"
X  (if (not list-of-lower-order-pev)
X    (multiple-value-setq
X      (list-of-lower-order-phi list-of-lower-order-pev)
X      (step-down-Levinson-Durbin-recursions
X       AR-coeffs variance
X       :process-variance-p process-variance-p)))
X  (let ((p (length AR-coeffs))
X        (current-lag 0)
X        current-phi-vector)
X    (setf (aref result 0) (nth 0 list-of-lower-order-pev))
X    (when (plusp max-lag)
X      (dotimes (k (min p max-lag))
X        (incf current-lag)   ; this is k+1
X        (setf (aref result current-lag)
X              (*  (aref (nth k list-of-lower-order-phi) k)   ;reflection coeff
X                  (nth k list-of-lower-order-pev)))
X        (when (plusp k)
X          (setf current-phi-vector
X                (nth (1- k) list-of-lower-order-phi))
X          (dotimes (j (length current-phi-vector))
X            (incf (aref result current-lag)
X                  (* (aref current-phi-vector j)
X                     (aref result (- current-lag j 1)))))))
X      (when (> max-lag p)
X        (setf current-phi-vector (nth (1- p) list-of-lower-order-phi))
X        (dotimes (k (- max-lag p))
X          (incf current-lag)
X          (setf (aref result current-lag) 0.0)
X          (dotimes (j p)
X            (incf (aref result current-lag)
X                  (* (aref current-phi-vector j)
X                     (aref result (- current-lag j 1))))))))
X    (values result)))
X
#|
;;; coefficients are for the AR(4) process
;;; described by Equation (46a) of the SAPA book ...
(ar-coeffs->acvs
X #(2.7607 -3.8106  2.6535 -0.9238)
X 1.0 6 :process-variance-p nil)
;==> #(761.7172900314962 545.7530764979967 27.14232320481949 -487.6647246019976 -705.2431860566808 -520.8142016231673 -69.504506541503)
;;; agrees out to lag 6 with lower plot of Figure 148
|#
X
;-------------------------------------------------------------------------------
(defun ar-coeffs->variance (ar-coeffs innovations-variance)
X  "given
X   [1] AR-coeffs (required)
X       ==> a sequence of AR coefficients
X           phi_{1,p}, ..., phi_{p,p}
X   [2] innovations-variance (required)
X       ==> innovations variance of the process
returns
X   [1] variance of the process"
X  (svref (ar-coeffs->acvs AR-coeffs innovations-variance 0
X                          :process-variance-p nil)
X         0))
X
#|
;;; coefficients are for the AR(4) process
;;; described by Equation (46a) of the SAPA book ...
(ar-coeffs->variance
X #(2.7607 -3.8106  2.6535 -0.9238)
X 1.0)
;==> 761.7172900314962
|#
X
;-------------------------------------------------------------------------------
(defun ar-coeffs->reflection-coeffs
X       (AR-coeffs
X        &key
X        (list-of-lower-order-phi
X         (step-down-Levinson-Durbin-recursions
X          AR-coeffs 1.0))
X        (result (make-array (length AR-coeffs))))
X  "given
X   [1] AR-coeffs (required)
X       ==> vector of length p with real or
X           complex-valued AR coefficients phi_{j,p}
X           for an AR(p) model of order p:
X           x_t = phi_{1,p} * x_{t-1} + phi_{2,p} * x_{t-2}
X                 + ... + phi_{p,p} * x_{t-p} + e_t
X   [2] list-of-lower-order-phi (keyword; calls step-down-Levinson-Durbin-recursions)
X       ==> must contain the list that is returned by
X           the function step-down-Levinson-Durbin-recursions;
X           this keyword parameter is useful if the result
X           of calling that function is already available
X   [3] results (keyword; vector of length p)
X       <== vector of size p into which
X            the reflection coefficients phi_{j,j}
X            of orders j = 1, 2, ..., p are placed
returns
X   [1] results, i.e., the reflection coefficients
---
Note: see Section 9.4 of the SAPA book"
X  (map-into result #'(lambda (a-seq)
X                        (elt a-seq (1- (length a-seq))))
X            list-of-lower-order-phi))
X
#|
;;; test case uses results of first test for two-step-burg-algorithm ...
(ar-coeffs->reflection-coeffs
X #(0.7475191644352183 -0.01583855174510637 0.46440844342490983 -0.4428652096133231)
X )
;==> #(0.7712281373828311 0.10369838100465756 0.16589516290564163 -0.4428652096133231)
;;; agrees with what two-step-burg-algorithm returned ...
|#
X
;-------------------------------------------------------------------------------
(defun ar-coeffs->prewhitening-filter
X       (AR-coeffs
X        &key
X        (result (make-array (1+ (length AR-coeffs)))))
X  "given
X   [1] AR-coeffs (required)
X       ==> a sequence of AR coefficients
X           phi_{1,p}, ..., phi_{p,p}
X   [2] result (keyword; vector of size (1+ (length AR-coeffs)))
X       <== sequence to hold coefficients for prewhitening filter
X           1, -phi_{1,p}, ..., -phi_{p,p}
returns
X   [1] result, the prewhitening filter
---
Note: this filter can be used to compute
X      forward prediction errors; to obtain
X      the filter needed to compute
X      backward prediction errors,
X      reverse the elements in result; i.e.,
X      evaluate (reverse result)"
X  (setf (elt result 0) 1.0)
X  (dotimes (i (length AR-coeffs) result)
X    (setf (elt result (1+ i)) (* -1 (elt AR-coeffs i)))))
X
#|
;;; forward prewhitening filter ...
(ar-coeffs->prewhitening-filter #(0.75 -0.5))
;==> #(1.0 -0.75 0.5)
X
;;; backward prewhitening filter ...
(reverse (ar-coeffs->prewhitening-filter #(0.75 -0.5)))
;==> #(0.5 -0.75 1.0)
|#
X
;-------------------------------------------------------------------------------
;;; Note: inclusion of p as a keyword is useful with routines that
;;;       map from the reflection coefficients because truncation
;;;       of this sequence at a certain point yields a valid lower
;;;       order model (NOT always the case with the AR coefficients)
(defun reflection-coeffs->ar-coeffs
X       (reflection-coeffs
X        &key
X        (p (length reflection-coeffs))
X        (scratch (make-array p))
X        (result (make-array p)))
X  "given
X   [1] reflection-coeffs (required)
X       ==> vector of length p with real-valued 
X           reflection coefficients
X           phi_{1,1}, phi_{2,2}, ... , phi_{p,p}
X   [2] p (keyword; length of reflection-coeffs)
X       ==> AR model order
X   [3] scratch (keyword; vector of length p)
X       ==> vector of size p used for intermediate results
X   [4] results (keyword; vector of length p)
X       <== vector of size p into which are placed
X           the AR coefficients phi_{j,p} for model
X           x_t = phi_{1,p} * x_{t-1} + phi_{2,p} * x_{t-2}
X                 + ... + phi_{p,p} * x_{t-p} + e_t
returns
X   [1] results, i.e., the AR coefficients
---
Note: see Section 9.4 of the SAPA book"
X  (when (plusp p)
X    (setf (aref result 0) (aref reflection-coeffs 0))
X    (if (> p 1)
X      (dotimes (k-minus-1 p result)
X        (setf (aref result k-minus-1)
X              (aref reflection-coeffs k-minus-1))
X        ;;; Use Levinson-Durbin recursions to generate
X        ;;; the remaining (k)-th order AR coefficients
X        (copy-vector result scratch :end k-minus-1)
X        (dotimes (j k-minus-1)
X          (setf (aref scratch j) (aref result j)))
X        (dotimes (j k-minus-1)
X          (setf (aref result j)
X                (- (aref scratch j)
X                   (* (aref result k-minus-1)
X                      (aref scratch (- k-minus-1 j 1)))))))
X      result)))
X
#|
;;; test case uses results of first test for two-step-burg-algorithm ...
(reflection-coeffs->ar-coeffs
X #(0.771228137382831 0.10369838100465757 0.16589516290564163 -0.4428652096133231)
X )
;==> #(0.7475191644352183 -0.01583855174510637 0.46440844342490983 -0.4428652096133231)
;;; agrees with what two-step-burg-algorithm returned ...
|#
X
;-------------------------------------------------------------------------------
(defun reflection-coeffs->variance
X       (reflection-coeffs
X        innovations-variance
X        &key
X        (p (length reflection-coeffs)))
X  "given
X   [1] reflection-coeffs (required)
X       ==> vector of length p with real-valued 
X           reflection coefficients
X           phi_{1,1}, phi_{2,2}, ... , phi_{p,p}
X   [2] innovations-variance (required)
X       ==> innovations variance of the process
X   [3] p (keyword; length of reflection-coeffs)
X       ==> AR model order
returns
X   [1] variance of the process
---
Note: this is a recursive application of
X      Equation (404c) of the SAPA book"
X  (dotimes (k p innovations-variance)
X    (divf innovations-variance
X          (- 1.0 (expt (abs (aref reflection-coeffs k)) 2)))))
X
#|
;;; test case uses results of first test for two-step-burg-algorithm ...
(reflection-coeffs->variance
X #(0.771228137382831 0.10369838100465757 0.16589516290564163 -0.4428652096133231)
X 245.25050973984762)
;==> 782.6399999999999
(sample-variance
X #(71 63 70 88 99 90 110 135 128 154 156 141 131 132 141 104 136 146 124 129))
;==> 19566/25
(float 19566/25)
;==> 782.64
;;; good agreement ...
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  ar-coeffs->sdf
;;;                 ar-coeffs->sdf-at-single-freq
;;;                 form-spectrum-with-accurate-peaks
;;;                 integrate-sdf
;;;  compute or manipulate an AR sdf estimate.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun ar-coeffs->sdf
X       (AR-coeffs
X        innovations-variance
X        &key
X        (N-nonzero-freqs 256)
X        (sampling-time 1.0)
X        (return-est-for-0-freq-p t)
X        (sdf-transformation #'convert-to-dB)
X        (result-sdf (make-array
X                     (if return-est-for-0-freq-p
X                       (1+ N-nonzero-freqs)
X                       N-nonzero-freqs)))
X        (return-frequencies-p nil)
X        (freq-transformation nil)
X        (result-freq (if return-frequencies-p
X                       (make-array
X                        (if return-est-for-0-freq-p
X                          (1+ N-nonzero-freqs)
X                          N-nonzero-freqs)))))
X  "given
X   [1] AR-coeffs (required)
X       ==> vector with real-valued AR coefficients phi_{j,p}
X           for an AR(p) model of order p:
X           x_t = phi_{1,p} * x_{t-1} + phi_{2,p} * x_{t-2}
X                 + ... + phi_{p,p} * x_{t-p} + e_t
X   [2] innovations-variance (required)
X       ==> innovations variance of the process
X   [3] N-nonzero-freqs (keyword; 256)
X       ==> the number of nonzero frequencies at which the AR sdf
X           is to be computed (need NOT be a power of 2);
X           the first nonzero frequency (and also spacing
X           between frequencies) is given by
X           1/(2 * N-nonzero-freqs * sampling time)
X           = Nyquist frequency/N-nonzero-freqs;
X           the last nonzero frequency is the Nyquist frequency
X   [4] sampling-time (keyword; 1.0)
X       ==> sampling time (called delta t in the SAPA book)
X   [5] return-est-for-0-freq-p (keyword; t)
X       ==> if t, sdf is computed at zero frequency;
X           otherwise, it is not computed.
X   [6] sdf-transformation (keyword; #'convert-to-dB)
X       ==> a function of one argument or nil;
X           if bound to a function, the function is used
X           to transform all elements of result-sdf
X   [7] result-sdf (keyword; vector of correct length)
X       <== vector of length N-nonzero-freqs (if return-est-for-0-freq-p is nil)
X           or N-nonzero-freqs +1 (if return-est-for-0-freq-p is t)
X           into which the properly transformed sdf is placed
X   [8] return-frequencies-p (keyword; nil)
X       ==> if t, the frequencies associated with the transfer function
X           are computed and returned in result-freq
X   [9] freq-transformation (keyword; nil)
X       ==> a function of one argument or nil;
X           if bound to a function, the function is used
X           to transform all elements of result-freq
X           (ignored unless return-frequencies-p is true)
X  [10] result-freq (keyword; nil or vector of correct length)
X       <== not used if return-frequencies-p nil; otherwise,
X           vector of length N-nonzero-freqs (if return-est-for-0-freq-p is nil)
X           or N-nonzero-freqs +1 (if return-est-for-0-freq-p is t)
X           into which the frequencies associated with the values
X           in result-sdf are placed
returns
X   [1] result-sdf, a vector holding
X       the properly transformed sdf
X   [2] nil (if return-frequencies-p is nil) or
X       result-freq (if return-frequencies-p is t),
X       where result-freq is a vector holding
X       the properly transformed frequencies
X       associated with values in  result-sdf
X   [3] the length of the vector result-sdf
---
Note: see Equation (392b) of the SAPA book"
X  (let* ((N-dft (* 2 N-nonzero-freqs))
X         (scratch (make-array N-dft))
X         (offset (if return-est-for-0-freq-p 0 1))
X         (N-freqs (if return-est-for-0-freq-p
X                    (1+ N-nonzero-freqs)
X                    N-nonzero-freqs))
X         (p (length AR-coeffs))
X         (num-constant (* innovations-variance sampling-time)))
X    ;(declare (dynamic-extent scratch))
X    (ar-coeffs->prewhitening-filter AR-coeffs :result scratch)
X    (fill scratch 0.0 :start (1+ p))
X    (dft! scratch)
X    (dotimes (i N-freqs)
X      (setf (aref result-sdf i)
X            (/ num-constant
X               (expt (abs (aref scratch (+ i offset))) 2))))
X    ;;; convert sdf if needed
X    (if sdf-transformation
X      (transform-a-sequence! sdf-transformation result-sdf))
X    ;;; create vector of associated frequencies ...
X    (when return-frequencies-p
X      (let ((freq-constant (* 2 N-nonzero-freqs sampling-time)))
X        (if freq-transformation
X          (dotimes (i N-freqs)
X            (setf (aref result-freq i)
X                  (funcall freq-transformation
X                           (/ (float (+ i offset)) freq-constant))))
X          (dotimes (i N-freqs)
X            (setf (aref result-freq i)
X                  (/ (float (+ i offset)) freq-constant))))))
X    (values result-sdf result-freq N-freqs)))
X
#|
(multiple-value-bind (sdf freq)
X                     (ar-coeffs->sdf
X                      #(0.6990553262846771
X                        -0.005832588617464066
X                        0.15582353261696208)
X                      271.77271305830277
X                      :return-frequencies-p t
X                      :sampling-time 0.25)
X  (let ((count 0))
X    (dolist (i '(0 1 2 3 4  126 127 128 129 130 252 253 254 255 256) (values))
X      (format t "~&~7,5F: ~10,5F" (svref freq i) (svref sdf i))
X      (if (zerop (mod (incf count) 5))
X        (format t "~&...")))))
;==>
0.00000:   34.74458
0.00781:   34.69756
0.01562:   34.55951
0.02344:   34.33891
0.03125:   34.04806
...
0.98437:   17.27619
0.99219:   17.25707
1.00000:   17.23760
1.00781:   17.21776
1.01562:   17.19750
...
1.96875:   12.93563
1.97656:   12.93224
1.98437:   12.92982
1.99219:   12.92836
2.00000:   12.92788
X
(multiple-value-bind (sdf freq N-freq)
X                     (ar-coeffs->sdf
X                      #(0.6990553262846771
X                        -0.005832588617464066
X                        0.15582353261696208)
X                      271.77271305830277
X                      :N-nonzero-freqs 5
X                      :return-frequencies-p t
X                      :sampling-time 0.25)
X  (dotimes (i N-freq (values))
X    (format t "~&~7,5F: ~10,5F" (svref freq i) (svref sdf i))
X    ))
;==>
0.00000:   34.74458
0.40000:   20.98899
0.80000:   17.73602
1.20000:   16.53083
1.60000:   14.12722
2.00000:   12.92788
|#
X
;-------------------------------------------------------------------------------
(defun ar-coeffs->sdf-at-single-freq
X       (AR-coeffs
X        innovations-variance
X        freq
X        &key
X        (sampling-time 1.0)
X        (sdf-transformation #'convert-to-dB))
X  "given
X   [1] AR-coeffs (required)
X       ==> vector with real-valued AR coefficients phi_{j,p}
X           for an AR(p) model of order p:
X           x_t = phi_{1,p} * x_{t-1} + phi_{2,p} * x_{t-2}
X                 + ... + phi_{p,p} * x_{t-p} + e_t
X   [2] innovations-variance (required)
X       ==> innovations variance of the process
X   [3] freq (required)
X       ==> the single frequency at which the sdf is
X           to be computed
X   [4] sampling-time (keyword; 1.0)
X       ==> sampling time (called delta t in the SAPA book)
X   [6] sdf-transformation (keyword; #'convert-to-dB)
X       ==> a function of one argument or nil;
X           if bound to a function, the function is used
X           to transform the value of the sdf at freq
returns
X   [1] properly transformed AR sdf evaluated at freq
X   [2] nil (if return-frequencies-p is nil) or
X       result-freq (if return-frequencies-p is t),
X       where result-freq is a vector holding
X       the properly transformed frequencies
X       associated with values in  result-sdf
X   [3] the length of the vector result-sdf
---
Note: see Equation (392b) of the SAPA book"
X  (let ((omega (* 2 pi freq sampling-time))
X        (complex-sum 1.0))
X    (dotimes (k (length AR-coeffs))
X      (decf complex-sum (* (aref AR-coeffs k)
X                           (exp (complex 0.0 (- (* omega (1+ k))))))))
X    (let ((result (/ (* innovations-variance sampling-time)
X                     (expt (abs complex-sum) 2))))
X      (if sdf-transformation
X        (funcall sdf-transformation result)
X        result))))
X
#|
(ar-coeffs->sdf-at-single-freq
X #(0.6990553262846771 -0.005832588617464066 0.15582353261696208)
X 271.77271305830277
X (* 4/256 2.0)        ;==> 0.03125 = freq
X :sampling-time 0.25)
;==> 34.048055842320565
|#
X
;-------------------------------------------------------------------------------
(defun sdf->sdf-with-accurate-peaks
X       (sdf
X        freqs
X        AR-coeffs
X        innovations-variance
X        &key
X        (sampling-time 1.0)
X        (sdf-transformation #'convert-to-dB)
X        (search-method :quadratic-approx)  ; or :bisection+Newton-Raphson
X        (dB-tolerance 0.1)
X        (accuracy (* 10.0 single-float-epsilon))
X        (maximum-number-of-iterations 20))
X  "given
X   [1] sdf (required)
X       ==> vector with sdf for a real-valued process
X           computed over a grid of frequencies (the grid
X           need not be uniform) -- note that the sdf values
X           can be either untransformed or transformed via
X           an order preserving transformation (such as decibels)
X   [2] freqs (required)
X       ==> vector with frequencies corresponding to values
X           in sdf (must be same size as vector with sdf values)
X   [3] AR-coeffs (required)
X       ==> vector with real-valued AR coefficients phi_{j,p}
X           for an AR(p) model of order p:
X           x_t = phi_{1,p} * x_{t-1} + phi_{2,p} * x_{t-2}
X                 + ... + phi_{p,p} * x_{t-p} + e_t
X   [4] innovations-variance (required)
X       ==> innovations variance of the process
X   [5] sampling-time (keyword; 1.0)
X       ==> sampling time (called delta t in the SAPA book)
X   [6] sdf-transformation (keyword; #'convert-to-dB)
X       ==> a function of one argument or nil;
X           if bound to a function, the function is used
X           to transform any newly computed values of the sdf
X   [7] search-method (keyword; :quadratic-approx)
X       ==> selects between two methods to search
X           for the location of peaks;
X           if set to :quadratic-approx, a quadratic approximation is used;
X           if set to :bisection+Newton-Raphson, an interval search
X           using a combination of bisection and Newton-Raphson is used
X   [8] dB-tolerance (keyword; 0.1)
X       ==> convergence criterion (used if search-method
X           is :quadratic-approx)
X   [9] accuracy (keyword; 10 * single-float-epsilon)
X       ==> convergence criterion (used if search-method
X           is :bisection+Newton-Raphson)
X  [10] maximum-number-of-iterations (keyword; 20)
X       ==> controls the number of iterations (used if search-method
X           is :bisection+Newton-Raphson)
returns
X   [1] a vector with values of sdf merged
X       with an additional set of values
X       chosen to represent peaks accurately
X   [2] the corresponding augmented vector of frequencies
X   [3] the length of the two vector returned
---
Note: see pages 524--5 of the SAPA book"
X  (let ((n (length sdf))
X        (additional-freqs '())
X        (additional-sdf '())
X        (ersatz-acs (if (eq search-method :bisection+Newton-Raphson)
X                      (acvs (ar-coeffs->prewhitening-filter AR-coeffs)
X                            :center-data-p nil :acs-p t))))
X    ;;; search for local peaks ...
X    (dotimes (i (- n 2))
X      (if (and (> (aref sdf (+ i 1)) (aref sdf i))
X               (> (aref sdf (+ i 1)) (aref sdf (+ i 2))))
X        ;;; a local peak has been found, so we attempt to
X        ;;; more accurately determine its location ...
X        (let* ((left-freq (aref freqs i))
X               (mid-freq (aref freqs (+ i 1)))
X               (right-freq (aref freqs (+ i 2)))
X               (a-new-freq
X                (case search-method
X                  (:quadratic-approx
X                   (find-peak-of-ar-sdf-using-quadratic-approx
X                    left-freq mid-freq right-freq
X                    AR-coeffs
X                    :sampling-time sampling-time
X                    :dB-tolerance dB-tolerance))
X                  (:bisection+Newton-Raphson
X                   (find-peak-or-valley-of-sdf-using-bisection+Newton-Raphson
X                    left-freq
X                    right-freq
X                    ersatz-acs
X                    :sampling-time sampling-time
X                    :accuracy accuracy 
X                    :maximum-number-of-iterations maximum-number-of-iterations))
X                  (otherwise
X                   (error ":search-method (~A) must be set to :quadratic-approx or :bisection+Newton-Raphson"
X                          search-method)))))
X          ;;; use a-new-freq only if it is distinct from all
X          (if (not (member a-new-freq `(,left-freq ,mid-freq ,right-freq)
X                           :test #'=))
X            (let ((a-symmetric-freq
X                   (find-symmetric-frequency
X                    left-freq mid-freq right-freq a-new-freq)))
X              (cond
X               (a-symmetric-freq
X                (push (min a-new-freq a-symmetric-freq) additional-freqs)
X                (push (ar-coeffs->sdf-at-single-freq
X                       AR-coeffs
X                       innovations-variance
X                       (car additional-freqs)
X                       :sampling-time sampling-time
X                       :sdf-transformation sdf-transformation)
X                      additional-sdf)
X                (push (max a-new-freq a-symmetric-freq) additional-freqs)
X                (push (ar-coeffs->sdf-at-single-freq
X                       AR-coeffs
X                       innovations-variance
X                       (car additional-freqs)
X                       :sampling-time sampling-time
X                       :sdf-transformation sdf-transformation)
X                      additional-sdf))
X               (t
X                (push a-new-freq additional-freqs)
X                (push (ar-coeffs->sdf-at-single-freq
X                       AR-coeffs
X                       innovations-variance
X                       a-new-freq
X                       :sampling-time sampling-time
X                       :sdf-transformation sdf-transformation)
X                      additional-sdf))))))))
X    ;;; all local peaks have now been found,
X    ;;; so we now merge the new material with the old ...
X    (cond
X     ((plusp (length additional-freqs))
X      (setf additional-freqs (reverse additional-freqs)
X            additional-sdf (reverse additional-sdf))
X      (let* ((n-additional (length additional-freqs))
X             (n-total (+ n n-additional))
X             (combined-freqs (make-array n-total))
X             (combined-sdf (make-array n-total))
X             (k 0)
X             (k-add 0))
X        (dotimes (i n-total (values combined-sdf combined-freqs n-total))
X          (cond
X           ((or (= k-add n-additional)
X                (< (aref freqs k) (elt additional-freqs k-add)))
X            (setf (aref combined-freqs i) (aref freqs k)
X                  (aref combined-sdf i) (aref sdf k))
X            (incf k))
X           (t
X            (setf (aref combined-freqs i) (elt additional-freqs k-add)
X                  (aref combined-sdf i) (elt additional-sdf k-add))
X            (incf k-add))))))
X     (t
X      (values sdf freqs n)))))
X
#|
;;; coefficients are for the AR(4) process
;;; described by Equation (46a) and shown in Figure 148
;;; of the SAPA book ...
(multiple-value-bind (sdf freqs N-freq)
X                     (ar-coeffs->sdf
X                      #(2.7607 -3.8106  2.6535 -0.9238)
X                      1.0
X                      :N-nonzero-freqs 50
X                      :return-frequencies-p t)
X  (dotimes (i N-freq)
X    (if (<= 0.0999 (svref freqs i) 0.1401)
X      (format t "~&~7,5F: ~10,5F"
X              (svref freqs i)
X              (svref sdf i))))
X  (format t "~&finding peaks using quadratic approximation ...")
X  (multiple-value-bind (extended-sdf extended-freq N-freq)
X                       (sdf->sdf-with-accurate-peaks
X                        sdf
X                        freqs
X                        #(2.7607 -3.8106  2.6535 -0.9238)
X                        1.0)
X    (dotimes (i N-freq (values))
X      (if (<= 0.0999 (svref extended-freq i) 0.1401)
X        (format t "~&~7,5F: ~10,5F"
X                (svref extended-freq i)
X                (svref extended-sdf i)))))
X  (format t "~&finding peaks using bisection + Newton-Raphson ...")
X  (multiple-value-bind (extended-sdf extended-freq N-freq)
X                       (sdf->sdf-with-accurate-peaks
X                        sdf
X                        freqs
X                        #(2.7607 -3.8106  2.6535 -0.9238)
X                        1.0
X                        :search-method :bisection+Newton-Raphson)
X    (dotimes (i N-freq (values))
X      (if (<= 0.0999 (svref extended-freq i) 0.1401)
X        (format t "~&~7,5F: ~10,5F"
X                (svref extended-freq i)
X                (svref extended-sdf i))))))
;==>
0.10000:   31.48832
0.11000:   43.67693
0.12000:   36.12925
0.13000:   35.59567
0.14000:   42.14003
finding peaks using quadratic approximation ...
0.10000:   31.48832
0.11000:   43.67693
0.11066:   43.61783
0.11133:   43.24097
0.12000:   36.12925
0.13000:   35.59567
0.13903:   42.05075
0.13951:   42.18959
0.14000:   42.14003
finding peaks using bisection + Newton-Raphson ...
0.10000:   31.48832
0.11000:   43.67693
0.11022:   43.69726
0.11044:   43.67740
0.12000:   36.12925
0.13000:   35.59567
0.13928:   42.14315
0.13964:   42.19592
0.14000:   42.14003
X
;;; Let's see if bisection + Newton-Raphson did a decent job on first peak ...
(dolist (a-freq '(0.11021 0.11022 0.11023) (values))
X  (print (ar-coeffs->sdf-at-single-freq
X          #(2.7607 -3.8106  2.6535 -0.9238)
X          1.0
X          a-freq)))
;==>
43.69721694178745 
43.6972567318079 
43.697213122598804
;;; ... on second peak ...
(dolist (a-freq '(0.13963 0.13964 0.13965) (values))
X  (print (ar-coeffs->sdf-at-single-freq
X          #(2.7607 -3.8106  2.6535 -0.9238)
X          1.0
X          a-freq)))
;==>
42.195895846811524 
42.19592320450518 
42.19586698000828 
;;; ... results look quite good
|#
X
;-------------------------------------------------------------------------------
;;; Because display of AR sdf's is problematic due to the possibility
;;; of missing very sharp peaks, it is useful to compare the variance
;;; of a displayed sdf with the sample variance for a time series.
;;; If the sdf is evaluated over an irregular grid of frequencies from
;;; f = 0 to f = Nyquist, the following function may be used to compute
;;; its corresponding variance.  It uses a trapazoid approximation.
(defun integrate-sdf
X       (sdf
X        freqs
X        &key
X        (fiddle-factor 2.0)
X        (return-integrated-sdf-p nil)
X        (result (when return-integrated-sdf-p
X                  (make-array (length sdf)))))
X  "given
X   [1] sdf (required)
X       ==> vector with sdf for a real or complex-valued process
X           computed over a grid of frequencies (the grid
X           need not be uniform) -- note that the sdf values
X           must be untransformed (e.g., not in decibels)
X   [2] freqs (required)
X       ==> vector with frequencies over which
X           sdf has been computed
X           (must be same size as vector with sdf values)
X   [3] fiddle-factor (keyword; 2.0)
X       ==> values returned are all multiplied by this factor;
X           the default of 2 is useful for real-valued processes
X           since these have symmetric two-sided sdf's.
X   [4] return-integrated-sdf-p (keyword; nil)
X       ==> if t, integrated sdf as a function of frequency
X           is computed and returned in result
X   [5] result (keyword; vector of length sdf)
X       <== vector with integrated sdf (if return-integrated-sdf-p is t)
X           or nil (if return-integrated-sdf-p is nil)
returns
X   [1] the integral of the sdf * fiddle-factor
X   [2] result, i.e., the integrated sdf
X       if return-integrated-sdf-p is t" 
X  (let ((sum 0.0)
X        delta-f
X        the-increment)
X    (if return-integrated-sdf-p (setf (svref result 0) 0.0))
X    (dotimes (i (1- (length sdf)) (values (* fiddle-factor sum) result))
X      (setf delta-f (- (aref freqs (1+ i)) (aref freqs i))
X            ;;; area of "supporting" rectangle + area of triangle on top of it
X            the-increment (+ (* delta-f (min (aref sdf (1+ i))
X                                             (aref sdf i)))
X                             (* 0.5 delta-f (abs (- (aref sdf (1+ i))
X                                                    (aref sdf i))))))
X      (incf sum the-increment)
X      (if return-integrated-sdf-p
X        (setf (svref result (1+ i))
X              (+ (* fiddle-factor the-increment) (svref result i)))))))
X
#|
;;; coefficients are for the AR(4) process
;;; described by Equation (46a) and shown in Figure 148
;;; of the SAPA book ...
(ar-coeffs->variance
X #(2.7607 -3.8106  2.6535 -0.9238)
X 1.0)
;==> 761.7172900314962
(dolist (n-freq '(50 128 512) (values))
X  (format t "~&number of nonzero frequencies = ~D:" n-freq)
X  (multiple-value-bind (sdf freqs)
X                       (ar-coeffs->sdf
X                        #(2.7607 -3.8106  2.6535 -0.9238)
X                        1.0
X                        :N-nonzero-freqs n-freq
X                        :return-frequencies-p t
X                        :sdf-transformation nil)
X    (print (integrate-sdf sdf freqs))
X    (multiple-value-bind (extended-sdf extended-freq)
X                         (sdf->sdf-with-accurate-peaks
X                          sdf
X                          freqs
X                          #(2.7607 -3.8106  2.6535 -0.9238)
X                          1.0
X                          :sdf-transformation nil)
X      (print (integrate-sdf extended-sdf extended-freq)))
X    (multiple-value-bind (extended-sdf extended-freq)
X                         (sdf->sdf-with-accurate-peaks
X                          sdf
X                          freqs
X                          #(2.7607 -3.8106  2.6535 -0.9238)
X                          1.0
X                          :sdf-transformation nil
X                          :search-method :bisection+Newton-Raphson)
X      (print (integrate-sdf extended-sdf extended-freq)))))
;==>
number of nonzero frequencies = 50:
1005.8010282299119 
1019.4741312708759 
1023.8046237341628 
number of nonzero frequencies = 128:
765.713758711881 
796.38238484366 
796.5233656088724 
number of nonzero frequencies = 512:
761.7172884723057 
762.1844965261544 
762.155417339912 
;;; Note that the equally spaced grid here did better than
;;; either of the unequally spaced grids (created to get a better
;;; representation of the two spectral peaks) -- this doesn't
;;; always happen!
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  find-peak-of-ar-sdf-using-quadratic-approx
;;;                 find-peak-or-valley-of-sdf-using-bisection+Newton-Raphson
;;;  can be used to determine the frequencies at which an sdf has a peak
;;;  or a valley.  For details, see pages 479--80 and 524--5 of the SAPA book.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun find-peak-of-ar-sdf-using-quadratic-approx
X       (f_1
X        f_2
X        f_3
X        AR-coeffs
X        &key
X        (sampling-time 1.0)
X        (dB-tolerance 0.1))
X  "given
X   [1] f_1 (required)
X       ==> frequency f_1 such that f_1 < f_2 and S(f_1) < S(f_2),
X           where S(.) is the sdf associated with AR-coeffs
X   [2] f_2 (required)
X       ==> middle frequency
X   [3] f_3 (required)
X       ==> frequency f_3 such that f_3 > f_2 and S(f_3) < S(f_2)
X   [4] AR-coeffs (required)
X       ==> vector with real-valued AR coefficients phi_{j,p}
X           for an AR(p) model of order p:
X           x_t = phi_{1,p} * x_{t-1} + phi_{2,p} * x_{t-2}
X                 + ... + phi_{p,p} * x_{t-p} + e_t
X   [5] sampling-time (keyword; 1.0)
X       ==> sampling time (called delta t in the SAPA book)
X   [6] dB-tolerance (keyword; 0.1)
X       ==> convergence criterion
returns
X   [1] estimate of frequency of peak/valley in sdf
---
Note: for details on this method, see pages 20-21 of Adby and Dempster,
X      ``Introduction to Optimization Methods''; note in particular that
X      we search the inverse of the AR portion of the spectrum instead
X      of the spectrum itself in the belief that the quadratic approximation
X      is better with this reformulation"
X   (let* ((S_1 (/ (ar-coeffs->sdf-at-single-freq
X                   AR-coeffs 1.0 f_1
X                   :sampling-time sampling-time
X                   :sdf-transformation nil)))
X          (S_2 (/ (ar-coeffs->sdf-at-single-freq
X                   AR-coeffs 1.0 f_2
X                   :sampling-time sampling-time
X                   :sdf-transformation nil)))
X          (S_3 (/ (ar-coeffs->sdf-at-single-freq
X                   AR-coeffs 1.0 f_3
X                   :sampling-time sampling-time
X                   :sdf-transformation nil)))
X          (a_23 (- f_2 f_3))
X          (a_31 (- f_3 f_1))
X          (a_12 (- f_1 f_2))
X          (b_23 (- (* f_2 f_2) (* f_3 f_3)))
X          (b_31 (- (* f_3 f_3) (* f_1 f_1)))
X          (b_12 (- (* f_1 f_1) (* f_2 f_2)))
X          (f_4 (/ (+ (* b_23 S_1) (* b_31 S_2) (* b_12 S_3))
X                  (* 2.0 (+ (* a_23 S_1) (* a_31 S_2) (* a_12 S_3)))))
X          (S_4 (/ (ar-coeffs->sdf-at-single-freq
X                   AR-coeffs 1.0 f_4
X                   :sampling-time sampling-time
X                   :sdf-transformation nil))))
X    (cond
X     ((< (abs (convert-to-dB (/ S_2 S_4))) dB-tolerance)
X      (values f_4))
X     ((and (< f_1 f_4 f_2) (< S_4 S_1) (< S_4 S_2))
X      (find-peak-of-ar-sdf-using-quadratic-approx
X       f_1 f_4 f_2
X       AR-coeffs
X       :sampling-time sampling-time
X       :dB-tolerance dB-tolerance))
X     ((and (< f_2 f_4 f_3) (< S_4 S_2) (< S_4 S_3))
X      (find-peak-of-ar-sdf-using-quadratic-approx
X       f_2 f_4 f_3
X       AR-coeffs
X       :sampling-time sampling-time
X       :dB-tolerance dB-tolerance))
X     ((and (< f_2 f_4 f_3) (< S_2 S_4))
X      (find-peak-of-ar-sdf-using-quadratic-approx
X       f_1 f_2 f_4
X       AR-coeffs
X       :sampling-time sampling-time
X       :dB-tolerance dB-tolerance))
X     ((and (< f_1 f_4 f_2) (< S_2 S_4))
X      (find-peak-of-ar-sdf-using-quadratic-approx
X       f_4 f_2 f_3
X       AR-coeffs
X       :sampling-time sampling-time
X       :dB-tolerance dB-tolerance))
X     (t
X      (format t "~&find-peak-of-ar-sdf-using-quadratic-approx: non-convergence at f = ~F, ~F, ~F"
X              f_1 f_2 f_3)
X      (values f_2)))))
X
#|
;;; coefficients are for the AR(4) process
;;; described by Equation (46a) and shown in Figure 148
;;; of the SAPA book ...
(find-peak-of-ar-sdf-using-quadratic-approx
X 0.11000
X 0.11025
X 0.11050
X #(2.7607 -3.8106  2.6535 -0.9238)
X )
;==> 0.11022137202732714
|#
X
;-------------------------------------------------------------------------------
(defun find-peak-or-valley-of-sdf-using-bisection+Newton-Raphson
X       (f-left
X        f-right
X        ersatz-acvs
X        &key
X        (sampling-time 1.0)
X        (N (length ersatz-acvs))
X        (accuracy (* 10.0 single-float-epsilon))
X        (maximum-number-of-iterations 20))
X  "given
X   [1] f-left (required)
X       ==> left bracket  for frequency of peak/valley
X   [2] f-right (required)
X       ==> right bracket for frequency of peak/valley
X   [3] ersatz-acvs (required)
X       ==> vector of length N with acvs (or acs) s_0 to s_{N-1}
X           corresponding to sdf (only elements 1 to N-1 are used)
X   [4] sampling-time (keyword; 1.0)
X       ==> positive number
X   [5] N (keyword; length of ersatz-acvs)
X       ==> positive integer
X   [6] accuracy (keyword; 10.0 * single-float-epsilon)
X       ==> passed to bisection-with-Newton-Raphson
X   [7] maximum-number-of-iterations (keyword; 20)
X       ==> passed to bisection-with-Newton-Raphson
returns
X   [1] estimate of frequency of peak/valley in sdf
X   [2] number of iterations required
---
Note: searches for peak using a combination of
bisection and Newton-Raphson;
see pages 479--80 and 524--5 of the SAPA book"
X  (let ((N-1 (1- N)))
X    (multiple-value-bind (omega-answer iteration-count)
X                         (bisection-with-Newton-Raphson
X                          #'(lambda (omega)
X                              (let ((first-deriv 0.0)
X                                    (tau 0))
X                                (dotimes (tau-1 N-1 (* -2.0 first-deriv))
X                                  (incf tau)
X                                  (incf first-deriv
X                                        (* tau
X                                           (elt ersatz-acvs tau)
X                                           (sin (* tau omega)))))))
X                          #'(lambda (omega)
X                              (let ((second-deriv 0.0)
X                                    (tau 0))
X                                (dotimes (tau-1 N-1 (* -2.0 second-deriv))
X                                  (incf tau)
X                                  (incf second-deriv
X                                        (* tau tau
X                                           (elt ersatz-acvs tau)
X                                           (cos (* tau omega)))))))
X                          (* 2 pi f-left sampling-time)
X                          (* 2 pi f-right sampling-time)
X                          :accuracy accuracy
X                          :maximum-number-of-iterations
X                          maximum-number-of-iterations)
X      (values (/ omega-answer (* 2 pi sampling-time)) iteration-count))))
X
#|
;;; coefficients are for the AR(4) process
;;; described by Equation (46a) and shown in Figure 148
;;; of the SAPA book ...
(find-peak-or-valley-of-sdf-using-bisection+Newton-Raphson
X 0.1100
X 0.1105
X (acvs (ar-coeffs->prewhitening-filter #(2.7607 -3.8106  2.6535 -0.9238))
X                            :center-data-p nil :acs-p t)
X )
;==> 0.1102197683614135
;    6
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  generate-forward-innovations
;;;                 generate-backward-innovations
;;;  compute the observed forward and backward prediction errors for
;;;  a given time series and AR prewhitening filter.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun generate-forward-innovations
X       (time-series
X        AR-coeffs
X        &key
X        (start 0)
X        (end (length time-series))
X        (center-data t)
X        (result (make-array (- (- end start) (length AR-coeffs)))))
X  "given
X   [1] time-series (required)
X       ==> a vector containing a real-valued 
X           time series x_t
X   [2] AR-coeffs (required)
X       ==> vector with real-valued AR coefficients phi_{j,p}
X           for an AR(p) model of order p:
X           x_t = phi_{1,p} * x_{t-1} + phi_{2,p} * x_{t-2}
X                 + ... + phi_{p,p} * x_{t-p} + e_t
X   [3] start (keyword; 0)
X       ==> start index of vector to be used
X   [4] end (keyword; length of the-seq)
X       ==> 1 + end index of vector to be used
X   [5] center-data (keyword; t)
X       ==> if t, subtract sample mean from time series;
X           if a number, subtracts that number from time series
X           (the vector time-series is not altered)
X   [6] result (keyword; vector of appropriate length)
X       <== vector to contain forward innovations
X                            p
X           f_t = x_{t+p} - SUM phi_{j,p} x_{t-j},  t = 1, ..., N - p
X                           j=1
returns
X   [1] result, i.e., the fprward innovations"
X  (filter-time-series
X   (center&taper-time-series time-series
X                       :center-data center-data
X                       :start start
X                       :end end)
X   (ar-coeffs->prewhitening-filter AR-coeffs)
X   :result result))
X
#|
(generate-forward-innovations
X #(71 63 70 88 99 90 110 135 128 154 156 141 131 132 141 104 136 146 124 129)
X #(0.7475191644352183 -0.01583855174510637 0.46440844342490983 -0.4428652096133231)
X )
;==> #(7.541189677934277 -16.190207980322164 5.452393054435857 18.222543695391344 1.9025289159837548 20.25717110563984 -0.04210392524786499 -1.8028505644380886 -15.73306199076052 3.09022994051165 19.03618232922566 -26.674545308896946 28.233150201353947 9.989699743647158 2.1902410483799546 -7.453034761884252)
;    16
;;; good agreement with what two-step-burg-algorithm gave ...
|#
X
;-------------------------------------------------------------------------------
(defun generate-backward-innovations
X       (time-series
X        AR-coeffs
X        &key
X        (start 0)
X        (end (length time-series))
X        (center-data t)
X        (result (make-array (- (- end start) (length AR-coeffs)))))
X  "given
X   [1] time-series (required)
X       ==> a vector containing a real-valued 
X           time series x_t
X   [2] AR-coeffs (required)
X       ==> vector with real-valued AR coefficients phi_{j,p}
X           for an AR(p) model of order p:
X           x_t = phi_{1,p} * x_{t-1} + phi_{2,p} * x_{t-2}
X                 + ... + phi_{p,p} * x_{t-p} + e_t
X   [3] start (keyword; 0)
X       ==> start index of vector to be used
X   [4] end (keyword; length of the-seq)
X       ==> 1 + end index of vector to be used
X   [5] center-data (keyword; t)
X       ==> if t, subtract sample mean from time series;
X           if a number, subtracts that number from time series
X           (the vector time-series is not altered)
X   [6] result (keyword; vector of appropriate length)
X       <== vector to contain backward innovations
X                        p
X           b_t = x_t - SUM phi_{j,p} x_{t+j},  t = 1, ..., N-p
X                       j=1
returns
X   [1] result, i.e., the backward innovations"
X  (filter-time-series
X   (center&taper-time-series time-series
X                       :center-data center-data
X                       :start start
X                       :end end)
X   (reverse (ar-coeffs->prewhitening-filter AR-coeffs))
X   :result result))
X
#|
(generate-backward-innovations
X #(71 63 70 88 99 90 110 135 128 154 156 141 131 132 141 104 136 146 124 129)
X #(0.7475191644352183 -0.01583855174510637 0.46440844342490983 -0.4428652096133231)
X )
;==> #(-0.9808164276349629 -23.022636411463495 -16.266777119010605 -4.8485735216690795 -1.514397559765996 -10.303462500921928 -20.291200583839146 2.7814408815304876 -21.084905735054566 8.26942730376083 24.845207695933038 -3.2294508544868163 17.520376061131927 0.7742590731821384 23.55218276738641 -24.778733173674876)
;    16
;;; good agreement with what two-step-burg-algorithm gave ...
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The function   regression-with-ar-errors
;;;  does linear regression with AR errors.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun regression-with-ar-errors
X       (time-series
X        ind-vars
X        p
X        &key
X        (AR-coeff-estimator #'burg-algorithm)
X        (maximum-number-of-iterations 5)
X        (eps 1.0e-7)
X        (trace-AR-parameter-estimates-p nil)
X        (trace-LS-parameter-estimates-p nil))
X  "given
X   [1] time-series (required)
X       ==> a vector with dependent variables
X   [2] ind-vars (required)
X       ==> a list of vector(s) with independent variable(s)
X   [3] p (required)
X       ==> autoregressive order
X   [4] AR-coeff-estimator (keyword; #'burg-algorithm)
X       ==> AR estimation procedure
X   [5] maximum-number-of-iterations (keyword; 5)
X       ==> ... like I said ...
X   [6] eps (keyword; 1.0e-7)
X       ==> a small number -- used to test for convergence
X           in estimates of both regression and AR coefficients
X   [7] trace-AR-parameter-estimates-p (keyword; nil)
X       ==> if t, results of each iteration
X           on AR parameters are printed
X   [8] trace-LS-parameter-estimates-p (keyword; nil)
X       ==> if t, results of each iteration
X           on least squares parameters are printed
returns
X   [1] estimates of least squares coefficients
X   [2] estimates of AR coefficients
X   [3] variance of residuals
X   [4] either number of iterations required for convergence
X       or nil if maximum number of iterations completed
X       but convergence not yet achieved
---
Note: this function is based on Newton, TIMESLAB, page 241;
X      it has NOT been thoroughly tested (particularly with
X      regard to the convergence criterion),
X      so a heavy dose of the usual ``caveat emptor''
X      is appropriate here!!!"
X  ;;; First, we obtain the ordinary least squares estimates;
X  ;;; the least squares residuals are pointed to by transformed-ts
X  ;;; and should resemble a realization of an AR(p) process.
X  (flet ((trace-parameter-estimates
X             (estimates iter-tag label)
X           (format t "~&~A, iteration ~D:" label (1+ iter-tag))
X           (dotimes (j (length estimates))
X             (format t "~&   ~D:  ~F" j (elt estimates j)))))
X    (multiple-value-bind (ols-estimates transformed-ts)
X                         (ordinary-least-squares-Cholesky
X                          time-series
X                          ind-vars
X                          :compute-residuals-p t)
X      (when trace-LS-parameter-estimates-p
X        (trace-parameter-estimates
X         ols-estimates -1 "LS estimates"))
X      (let (previous-ols-estimates
X            (transformed-ind-vars '())
X            (n (length time-series))
X            (n-ind (length ind-vars))
X            (current-whitening-filter (make-array (1+ p)))
X            current-filter-length
X            current-pev-factor
X            an-ind-var-to-adjust an-ind-var
X            innovations-variance
X            AR-coeffs
X            previous-AR-coeffs)
X        ;;; Create space for transformed independent variables ...
X        (dotimes (i n-ind)
X          (push (make-array n) transformed-ind-vars))
X        ;;; Loop until there is little change in LS and AR parameter estimates ...
X        (dotimes (iter maximum-number-of-iterations
X                       (values ols-estimates
X                               AR-coeffs
X                               (ar-coeffs->variance
X                                AR-coeffs
X                                innovations-variance)
X                               nil))
X          (setf previous-ols-estimates ols-estimates
X                previous-AR-coeffs AR-coeffs)
X          ;;; Fit the AR model to the current set of LS residuals
X          (multiple-value-setq (AR-coeffs
X                                innovations-variance)
X            (funcall AR-coeff-estimator transformed-ts p :center-data nil))
X          (when trace-AR-parameter-estimates-p
X            (trace-parameter-estimates
X             AR-coeffs iter "AR estimates"))
X          (multiple-value-bind (list-of-pe-coeff list-of-pev)
X                               (step-down-Levinson-Durbin-recursions
X                                AR-coeffs
X                                innovations-variance
X                                :process-variance-p nil)
X            ;;; Transform the original time series and vectors of independent
X            ;;; variables in accordance with fitted AR(p) model ...
X            (dotimes (i n)
X              ;;; Normalized prediction error filter is constant
X              ;;; for i greater than or equal to p ...
X              (when (<= i p)
X                (setf current-pev-factor (/ (sqrt (elt list-of-pev i)))
X                      current-filter-length (1+ i)
X                      (aref current-whitening-filter 0) current-pev-factor)
X                (dotimes (j i)
X                  (setf (aref current-whitening-filter (1+ j))
X                        (* (- (aref (elt list-of-pe-coeff (1- i)) j))
X                           current-pev-factor))))
X              (setf (aref transformed-ts i) 0.0)
X              (dotimes (j current-filter-length)
X                (incf (aref transformed-ts i)
X                      (* (aref time-series (- i j))
X                         (aref current-whitening-filter j))))
X              (dotimes (k n-ind)
X                (setf an-ind-var-to-adjust (elt transformed-ind-vars k)
X                      an-ind-var (elt ind-vars k))
X                (setf (aref an-ind-var-to-adjust i) 0.0)
X                (dotimes (j current-filter-length)
X                  (incf (aref an-ind-var-to-adjust i)
X                        (* (aref an-ind-var (- i j))
X                           (aref current-whitening-filter j))))))
X            ;;; Apply ordinary least squares to obtain new least squares
X            ;;; estimates for the regression parameters ...
X            (setf ols-estimates (ordinary-least-squares-Cholesky
X                                 transformed-ts
X                                 transformed-ind-vars))
X            (when trace-LS-parameter-estimates-p
X              (trace-parameter-estimates
X               ols-estimates iter "LS estimates"))
X            ;;; convergence test (needs to be tuned up ...)
X            (if (and previous-AR-coeffs
X                     (< (compare-seqs
X                         previous-ols-estimates
X                         ols-estimates)
X                        eps)
X                     (< (compare-seqs
X                         previous-AR-coeffs
X                         AR-coeffs)
X                        eps))
X              (return (values ols-estimates
X                              AR-coeffs
X                              (ar-coeffs->variance
X                                AR-coeffs
X                                innovations-variance)
X                              (1+ iter))))
X            ;;; Compute the residuals with respect to the original
X            ;;; dependent and independent variables ...
X            (copy-vector time-series transformed-ts)
X            (dotimes (i n)
X              (dotimes (j n-ind)
X                (decf (aref transformed-ts i)
X                      (* (aref ols-estimates j)
X                         (aref (elt ind-vars j) i)))))))))))
X
#|
;;; coefficients are for the AR(2) process
;;; described by Equation (45) of the SAPA book ...
(regression-with-ar-errors
X (x+y! (sample-from-a-function #'(lambda (x) (+ pi (* 7 x)))
X                               :N 500)
X       (generate-ar-time-series
X        #(0.75 -0.5)
X        1.0
X        500
X        :process-variance-p t))
X (list (make-array 500 :initial-element 1.0)
X       (iota 0 499 :type 'vector))
X 2
X :trace-AR-parameter-estimates-p t
X :trace-LS-parameter-estimates-p t)
;==>
LS estimates, iteration 0:
X   0:  3.1724362189599566
X   1:  6.999820737703846
AR estimates, iteration 1:
X   0:  0.6959987596867
X   1:  -0.5249062346802079
LS estimates, iteration 1:
X   0:  3.174976778340132
X   1:  6.999811113561525
AR estimates, iteration 2:
X   0:  0.6959967546900704
X   1:  -0.52491167609844
LS estimates, iteration 2:
X   0:  3.1749768748282645
X   1:  6.999811113149029
AR estimates, iteration 3:
X   0:  0.6959967546455309
X   1:  -0.5249116762847136
LS estimates, iteration 3:
X   0:  3.174976874824981
X   1:  6.99981111314904
#(3.174976874824981 6.99981111314904)      ;truth is pi   and  7
#(0.6959967546455309 -0.5249116762847136)  ;truth is 0.75 and -0.5
0.9537407718258482                         ;truth is 1.0
3
;;; Note: don't expect agreement with the above numbers
;;;       (due to implementation dependence on random number
;;;       generation)
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  Everything below here consists of internal symbols in the SAPA package
;;;  and should be regarded as "dirty laundry" ...
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;; handles the different conventions for AR parameters
(defun flip-sign-of-coefficients! (a-sequence)
X  (map-into a-sequence #'(lambda (x) (* -1 x)) a-sequence))
X
;-------------------------------------------------------------------------------
(defun one-iteration-of-two-step-Burg (k step N Fk Bk phi limit)
X  ;;; Input:
X  ;;;  k	Order of the AR model fit so far (integer)
X  ;;;  step	Factor relating AR coefficients to time lag; i.e. j (integer)
X  ;;;  N	Number of data points in the original sequence (integer)
X  ;;;  Fk	kth order forward differences (real*8 vector of size N; 1:N)
X  ;;;  Bk	kth order backward differences (real*8 vector of size N; 1:N)
X  ;;;  phi	The Autoregresive model coefficients for original order
X  ;;;           (real vector of size k; 1:p)
X  ;;;  limit	limit on the magnitude of the piviot elements (real)
X  ;;;
X  ;;;  Output:
X  ;;;  Fk	forward differences that correspond to the new order (updated)
X  ;;;  Bk	Backward differences that correspond to the new order (updated)
X  ;;;  phi	The Autoregresive model coefficients for new order (updated)
X  ;;;
X  ;;; Model: X(j)= phi(1)*X(j-1*step) + ... + X(j-k*step) + epsilon(j)
X  (let ((f (make-array 3 :initial-element 0.0))
X        (g (make-array 3 :initial-element 0.0))
X        (f-prime (make-array 2))
X        (g-prime (make-array 2))
X        (f*f (make-array 5))
X        (f*g (make-array 5))
X        (g*g (make-array 5))
X        (f-prime*f*f (make-array 6))
X        (g-prime*f*g (make-array 6))
X        (f-prime*g*g (make-array 6))       
X        (A (make-array 6))   ;A is f'*f**2 - 2*g'*f*g + f'*g**2
X        (A-flip (make-array 6))
X        (roots (make-array 6))
X        (new-order (+ k 2))
X        (degree-of-A 0)
X        j F0 F1 B0 B1 a-small-number)
X    (dotimes (j-times (1+ (- N (* step (1+ new-order)))))
X      (setf j (+ j-times (1- (* step (1+ new-order)))))
X      (setf F0 (aref Fk j))  ; FO = F_k-2
X      (setf F1 (aref Fk (- j step)))  ;F1= S_1*step F_k-2
X      (setf B0 (aref Bk (- j (* new-order step))))  ;B0= S_k*step B_k-2
X      (setf B1 (aref Bk (- j (* (1- new-order) step))))  ;B1= S_(k-1)*step B_k-2
X      
X      (incf (aref f 0) (+ (* F0 F0) (* B0 B0)))
X      (decf (aref f 1) (+ (* 2 F0 B1) (* 2 B0 F1)))
X      (incf (aref f 2) (+ (* F1 F1) (* B1 B1)))
X
X      (incf (aref g 0) (* 2 F0 B0))
X      (decf (aref g 1) (+ (* 2 F0 F1) (* 2 B0 B1)))
X      (incf (aref g 2) (* 2 F1 B1)))  ;end of 100
X
X    (cond
X     ((= new-order 2)
X      ;;; In this case we must solve f-prime = 0.
X      (setf (aref A 0) (aref f 1))
X      (setf (aref A 1) (* 2 (aref f 2)))
X      (setf degree-of-A 1))
X     (t
X      ;;; In this case we must solve Equation 7:
X      ;;; f'*f**2 - 2*g'*f*g + f'*g**2 = 0
X      (setf (aref f-prime 0) (aref f 1))
X      (setf (aref f-prime 1) (* 2 (aref f 2)))
X      (setf (aref g-prime 0) (aref g 1))
X      (setf (aref g-prime 1) (* 2 (aref g 2)))
X      (multiply-2-polynomials f f :product-polynomial f*f)
X      (multiply-2-polynomials f-prime f*f :product-polynomial f-prime*f*f)
X      (multiply-2-polynomials f g :product-polynomial f*g)
X      (multiply-2-polynomials g-prime f*g :product-polynomial g-prime*f*g)
X      (multiply-2-polynomials g g :product-polynomial g*g)
X      (multiply-2-polynomials f-prime g*g :product-polynomial f-prime*g*g)
X      ;;; A= f'*f**2 - 2*g'*f*g + f'*g**2
X      (setf degree-of-A 0)
X      (dotimes (j 6)
X        (setf (aref A j)
X              (- (+ (aref f-prime*f*f j) (aref f-prime*g*g j))
X                 (* 2 (aref g-prime*f*g j)))))
X      (setf a-small-number (* (absolute-sum-of-vector-elements A) 1.0e-5))
X      (dotimes (j 6)
X        (if (> (abs (aref A j)) a-small-number)
X          (setf degree-of-A j)))))
X
X    ;;; Convention conflict: must flip A for use with zeros-of-polynomial
X    (dotimes (j (1+ degree-of-A))
X      (setf (aref A-flip (- degree-of-A j)) (aref A j)))
X    ;;; HACK: need to examine condition code and exit if roots
X    ;;;       were not properly found
X    (zeros-of-polynomial A-flip
X                         :degree-of-polynomial degree-of-A
X                         :the-roots roots)
X
X    ;;; Check the roots for minimum of h= f - g**2/f
X    ;;; Start with hmin as value on boundary; i.e. +/- one
X    (let* ((pmin limit)
X           (hmin (h-func f g pmin))
X           (ptry (- limit))
X           (htry (h-func f g ptry))
X           (sum-of-squares-new-order-minus-1 0.0)
X           (sum-of-squares 0.0)
X           reflection-coefficient-new-order-minus-1)
X      (cond
X       ((< htry hmin)
X        (setf pmin ptry)
X        (setf hmin htry)))
X
X      (dotimes (j degree-of-A)
X        (setf ptry (realpart (aref roots j)))
X        (cond
X         ((<= (abs ptry) limit)
X          (setf htry (h-func f g ptry))
X          (cond
X           ((< htry hmin)
X            (setf pmin ptry)
X            (setf hmin htry))))))   ;end of 220
X           
X      ;;; Set new pivot elements
X      (setf reflection-coefficient-new-order-minus-1
X            (setf (aref phi (- new-order 2)) pmin))
X      (setf (aref phi (1- new-order))
X            (/ (eval-poly-2 g pmin) (eval-poly-2 f pmin)))
X
X     ;;; Set the new differences
X      (dotimes (j-times (1+ (- N (* step new-order))))
X        (setf j (+ j-times (1- (* step new-order))))
X        (setf F0 (aref Fk j))
X        (setf B0 (aref Bk (- j (* (1- new-order) step))))
X	(setf (aref Fk j) (- F0 (* (aref phi (- new-order 2)) B0)))
X        (setf (aref Bk (- j (* (1- new-order) step)))
X              (- B0 (* (aref phi (- new-order 2)) F0)))
X        (incf sum-of-squares-new-order-minus-1
X              (+ (expt (aref Fk j) 2)
X                 (expt (aref Bk (- j (* (1- new-order) step))) 2))))
X
X      (dotimes (j-times (1+ (- N (* step (1+ new-order)))))
X        (setf j (+ j-times (1- (* step (1+ new-order)))))
X        (setf F0 (aref Fk j))
X        (setf B0 (aref Bk (- j (* new-order step))))
X	(setf (aref Fk j) (- F0 (* (aref phi (1- new-order)) B0)))
X        (setf (aref Bk (- j (* new-order step)))
X              (- B0 (* (aref phi (1- new-order)) F0)))
X        (incf sum-of-squares
X              (+ (expt (aref Fk j) 2)
X                 (expt (aref Bk (- j (* new-order step))) 2))))
X
X      ;;; Set New autoregressive coefficients
X      (let (phij phikj)
X        (dotimes (j (truncate (1- new-order) 2))
X          (setf	phij (aref phi j))
X          (setf phikj (aref phi (- new-order j 3)))
X          (setf (aref phi j) (- phij (* (aref phi (- new-order 2)) phikj)))
X          (setf (aref phi (- new-order j 3))
X                (- phikj (* (aref phi (- new-order 2)) phij))))
X
X        (dotimes (j (truncate new-order 2))
X          (setf	phij (aref phi j))
X          (setf phikj (aref phi (- new-order j 2)))
X          (setf (aref phi j) (- phij (* (aref phi (1- new-order)) phikj)))
X          (setf (aref phi (- new-order j 2))
X                (- phikj (* (aref phi (1- new-order)) phij))))
X
X        (values
X         ;;; rms of forward and backward prediction errors
X         (/ sum-of-squares (* 2.0 (- (1+ N) (* (1+ new-order) step))))
X         (/ sum-of-squares-new-order-minus-1
X            (* 2.0 (- (1+ N) (* new-order step))))
X         reflection-coefficient-new-order-minus-1
X         new-order
X         )))))
X
;-------------------------------------------------------------------------------
(defun h-func (f g x)
X  ;;; h(p)= fpoly(p) - gpoly(p)**2/fpoly(p)
X  (let ((temp (eval-poly-2 f x)))
X    (- temp
X       (/ (expt (eval-poly-2 g x) 2)
X          temp))))
X
;-------------------------------------------------------------------------------
(defun absolute-sum-of-vector-elements (vector)
X  #-allegro (reduce #'+ vector :key #'(lambda (x) (abs x)))
X  #+allegro (let ((the-sum 0.0))
X              (dotimes (i (length vector) the-sum)
X                (incf the-sum (abs (aref vector i)))))
X  )
X
;-------------------------------------------------------------------------------
(defun eval-poly-2 (polynomial x)
X  (+ (aref polynomial 0)
X     (* (aref polynomial 1) x)
X     (* (aref polynomial 2) x x)))
X
;-------------------------------------------------------------------------------
;;; Given equally spaced frequencies f_1 < f_2 < f_3 and a distinct third
;;; frequency f_N satisfying f_1 < f_N < f_3, this routine finds another
;;; frequency f_A satisfying f_1 < f_A < f_3 such that, in the set of four
;;; frequencies f_1, f_2, f_3 and f_A, two of them can be written as
;;; f_N - delta and f_N + delta for some delta > 0; if in fact f_A is
;;; equal to f_1, f_2 or f_3, nil is returned  (see the examples below).
(defun find-symmetric-frequency (f_1 f_2 f_3 f_N)
X  ;;; assumes that f_1, f_2, and f_3 are equally spaced
X  (let ((delta-f (- f_2 f_1))
X        (f_S (- (* 2.0 f_N) f_1))
X        (f_trans (/ (+ f_1 f_2) 2.0)))
X    (if (> f_N f_trans)
X      (decf f_S (* (ceiling (/ (- f_N f_trans) delta-f)) delta-f)))
X    (if (or (= f_S f_1) (= f_S f_2) (= f_S f_3) (= f_S f_N))
X      nil f_S)))
X
#|
(map 'list #'(lambda (x)
X               (find-symmetric-frequency 1.0 2.0 3.0 x))
X     '(1.1 1.46 1.5 1.56 2.1 2.46 2.5 2.56))
;==> (1.2000000000000002 1.92 nil 1.12 2.2 2.92 nil 2.12)
|#
SHAR_EOF
chmod 0644 parametric.lisp ||
echo 'restore of parametric.lisp failed'
Wc_c="`wc -c < 'parametric.lisp'`"
test 115199 -eq "$Wc_c" ||
	echo 'parametric.lisp: original size 115199, current size' "$Wc_c"
fi
# ============= random.lisp ==============
if test -f 'random.lisp' -a X"$1" != X"-c"; then
	echo 'x - skipping random.lisp (File already exists)'
else
echo 'x - extracting random.lisp (Text)'
sed 's/^X//' << 'SHAR_EOF' > 'random.lisp' &&
;;;-*- Mode: LISP; Package: :SAPA; Syntax: COMMON-LISP -*-
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;
;  random.lisp
;
;  a collection of Lisp functions to simulate stationary random processes ...
;  Note:  before compiling and loading random.lisp,
;         you should compile and load (in the order listed)
;            sapa-package.lisp, utilities.lisp and dft-and-fft.lisp
;
;  6/3/93
;
;  SAPA, Version 1.0; Copyright 1993, Donald B. Percival, All Rights Reserved
;
;  Use and copying of this software and preparation of derivative works
;  based upon this software are permitted.  Any distribution of this
;  software or derivative works must comply with all applicable United
;  States export control laws.
; 
;  This software is made available AS IS, and no warranty -- about the
;  software, its performance, or its conformity to any
;  specification -- is given or implied. 
;
;  Comments about this software can be addressed to dbp@apl.washington.edu
;-------------------------------------------------------------------------------
;;; (compile-file "ccl:SAPA;random.lisp")
;;; (load "ccl:SAPA;random.fasl")
;-------------------------------------------------------------------------------
(if (not (find-package :SAPA))
X  (defpackage :SAPA (:USE :COMMON-LISP)))
X
(in-package :SAPA)
X
(export '(;;; functions to generate normally distributed random deviates ...
X          ranorm
X          ranorms
X
X          ;;; white noise with different distributions ...
X          generate-white-noise
X          
X          ;;; normally distributed moving average and autoregressive processes
X          generate-ma-time-series
X          generate-ar-time-series
X          step-down-Levinson-Durbin-recursions
X
X          ;;; simulate a stationary process using frequency domain techniques
X          simulate-time-series-from-sdf
X          acvs-for-time-series-simulated-from-sdf
X          ))
X
;-------------------------------------------------------------------------------
(defconstant +sapa-sqrt-8-over-e+ (sqrt (/ 8.0d0 (exp 1.0d0))))
(defconstant +sapa-4-time-exp-of-1-over-4+ (* 4.0d0 (exp 0.25d0)))
(defconstant +sapa-4-time-exp-of-minus-1-point-35+ (* 4.0d0 (exp (- 1.35d0))))
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  ranorm
;;;                 ranorms
;;;  generate one or more uncorrelated normally distributed deviates.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun ranorm
X       ()
X  "returns a random deviate from a normal distribution
with zero mean and unit variance"
X  (let ((u (random 1.0)))
X    (cond
X     ((= u 0.0) (ranorm))   ;bad choice!
X     (t
X      (let* ((x (/ (* (- (random 1.0) 0.5)
X                      +sapa-sqrt-8-over-e+) u))
X             (xs (* x x)))
X        (cond
X         ((<= xs (- 5.0 (* u +sapa-4-time-exp-of-1-over-4+)))
X          x)   ;done
X         ((>= xs (+ 1.4 (/ +sapa-4-time-exp-of-minus-1-point-35+ u)))
X          (ranorm))   ;do it again
X         ((<= xs (- (* 4.0 (log u))))
X          x)
X         (t
X          (ranorm))))))))
X
#|
;;; each evaluation of ranorm produces a different random deviate,
;;; so don't expect to get the same result!
(ranorm)  ;==> -0.2719851788703296
|#
X
;-------------------------------------------------------------------------------
(defun ranorms
X       (n
X        &key
X        (result (make-array n)))
X  "given:
X   [1] n (required)
X       ==> number of normal random deviates
X           to be generated
X   [2] result (keyword; vector of length n)
X       <== vector to hold random deviates
returns:
X   [1] result, a vector with n normal random deviates
X       (i.e., a realization of length n of a white
X       noise process from a normal distribution)"
X  (dotimes (i n result)
X    (setf (svref result i) (ranorm))))
X
#|
(ranorms 3)
;==> #(-1.1128679745343053 1.4875994154684462 -0.06612547414267966)
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The function  generate-white-noise
;;;  generate a sequence of white noise with a specified distribution.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun generate-white-noise
X       (n
X        &key
X        (distribution :normal)
X        (result (make-array n)))
X  "given
X   [1] n (required)
X       ==> a sample size
X   [2] distribution (keyword; :normal)
X       ==> either a keyword or a function with no arguments
X           that returns a random deviate from a particular
X           distribution. Choices are
X            :binary :cauchy :chi-square-2
X            :double-exponential :exponential
X            :extreme-value :Gaussian (same as :normal)
X            :logistic :lognormal :normal :uniform
X   [3] result (keyword; vector of length n)
X       <== a sequence in which results are to be stored
return
X   [1] result, containing n samples from a white noise process
X       with a distribution specified by the keyword distribution"
X  (let ((generating-function
X         (if (keywordp distribution)
X           (case distribution
X             (:binary #'(lambda ()
X                          (if (<= (random 1.0) 0.5)
X                            0 1)))
X             (:cauchy #'(lambda ()
X                          (tan (* pi (- (random 1.0) 0.5)))))
X             (:chi-square-2 #'(lambda ()
X                                (* 2 (- (log (- 1.0 (random 1.0)))))))
X             (:double-exponential #'(lambda ()
X                                      (let ((U (random 1.0)))
X                                        (if (<= U 0.5)
X                                          (log (* 2 U))
X                                          (- (log (* 2 (- 1.0 U))))))))
X             (:exponential #'(lambda ()
X                               (- (log (- 1.0 (random 1.0))))))
X             (:extreme-value #'(lambda ()
X                                 (log (- (log (- 1.0 (random 1.0)))))))
X             (:Gaussian #'ranorm)
X             (:logistic #'(lambda ()
X                            (let ((U (random 1.0)))
X                              (log (/ U (- 1.0 U))))))
X             (:lognormal #'(lambda ()
X                             (exp (ranorm))))
X             (:normal #'ranorm)
X             (:uniform #'(lambda () (random 1.0)))
X             (otherwise
X              (error "~&generate-white-noise: unknown distribution keyword (~A)"
X                     distribution)))
X           distribution)))
X    (dotimes (i n result)
X      (setf (elt result i) (funcall generating-function)))))
X
#|
(generate-white-noise 10 :distribution :binary)
;==> #(1 0 1 0 1 0 0 0 0 0)
(generate-white-noise 3 :distribution :cauchy)
;==> #(-27.51099440214296 -0.6075812207412569 0.36550059169755994)
(generate-white-noise 3 :distribution :lognormal)
;==> #(0.5822708648379962 3.833558848750849 0.3592735989502123)
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  generate-ma-time-series
;;;                 generate-ar-time-series
;;;                 step-down-Levinson-Durbin-recursions
;;;  can be used to generate realizations of specified lengths from
;;;  normally distributed moving average (MA) or autoregressive processes (AR)
;;;  (these processes are discussed in Chapter 2 and 9 of the SAPA book).
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun generate-ma-time-series
X       (coeffs
X        variance
X        sample-size
X        &key
X        (process-variance-p t)
X        (result (make-array sample-size)))
X  "given
X   [1] coeffs (required)
X       ==> sequence of length q with MA coefficients:
X           x_t = e_t - coeffs_0*e_{t-1}
X                 - coeffs_1*e_{t-2}
X                 - ... - coeffs_{q-1}*e_{t-q}
X       (see Equation (43a) in the SAPA book)
X   [2] variance (required)
X       ==> process variance or innovations variance
X           (see keyword process-variance-p)
X   [3] sample-size (required)
X       ==> length of generated time series
X   [4] process-variance-p (keyword; t)
X       ==> if t, variance is taken to be process variance
X           if nil, variance is taken to be innovations variance
X   [5] result (keyword; vector of length sample-size)
X       <== vector with simulated series
generates realization of zero mean normally distributed MA(q) process and
returns
X   [1] vector of length sample-size with realization"
X  (let ((q+1 (1+ (length coeffs)))
X        (sd (sqrt (if process-variance-p
X                    (/ variance (1+ (sum-of-squares coeffs)))
X                    variance)))
X        (scratch-1 (make-array (1+ (length coeffs))))
X        (scratch-2 (make-array (1+ (length coeffs)))))
X    #+mcl(declare (dynamic-extent scratch-1 scratch-2))
X    (generate-white-noise q+1 :result scratch-1)
X    (a*x! sd scratch-1)
X    (dotimes (i (1- q+1))
X      (setf (aref scratch-2 (1+ i)) (- (elt coeffs i))))
X    (setf (aref scratch-2 0) 1.0)
X    (dotimes (i sample-size result)
X      (setf (aref result i)
X            (dot-product scratch-2 scratch-1))
X      (circular-shift-sequence scratch-1 :result scratch-1)
X      (setf (aref scratch-1 0) (* sd (ranorm))))))
X
#|
;;; See top plot of Figure 44 of the SAPA book.
(generate-ma-time-series #(1.0) 1.0 4 :process-variance-p nil)
;==> #(-0.555869655695869 -0.6755180985339049 0.3536345814801179 0.7421943411372343)
;;; See bottom plot of Figure 44 of the SAPA book.
(generate-ma-time-series #(-1.0) 1.0 4 :process-variance-p nil)
#(0.2105828208289749 1.1611109339530934 1.7492356283253514 0.07032022771512858)
|#
X
;-------------------------------------------------------------------------------
(defun generate-ar-time-series
X       (coeffs
X        variance
X        sample-size
X        &key
X        (process-variance-p t)
X        (list-of-lower-order-phi nil)
X        (list-of-lower-order-pev nil)
X        (result (make-array
X                 sample-size
X                 :initial-element 0.0)
X                result-supplied-p))
X  "given
X   [1] coeffs (required)
X       ==> sequence of length p with AR coefficients:
X           x_t = coeffs_0*x_{t-1} + coeffs_1*x_{t-2}
X                 + ... + coeffs_{p-1}*x_{t-p} + e_t
X           (see Equation (392a) in the SAPA book;
X           the coefficients can be real or complex-valued)
X   [2] variance (required)
X       ==> process variance or innovations variance
X           (see keyword process-variance-p)
X   [3] sample-size (required)
X       ==> length of generated time series
X   [4] process-variance-p (keyword; t)
X       ==> if t, variance is taken to be process variance
X           if nil, variance is taken to be innovations variance
X   [5] list-of-lower-order-phi (keyword; nil)
X       ==> to bypass call to step-down-Levinson-Durbin-recursions
X           (useful for multiple realizations)
X   [6] list-of-lower-order-pev (keyword; nil)
X       ==> to bypass call to step-down-Levinson-Durbin-recursions
X   [7] result (keyword; vector of length sample-size)
X       <== vector with simulated series
generates realization of zero mean normally distributed AR(p) process and
returns
X   [1] result, vector of length sample-size with realization
X   [2] list-of-lower-order-phi (can be used on subsequent calls
X       to this function to bypass call to
X       step-down-Levinson-Durbin-recursions)
X   [3] list-of-lower-order-pev (can be also used on subsequent calls)
---
Note: this function generates the proper stationary initial conditions"
X  (if result-supplied-p (fill result 0.0 :end sample-size))
X  (if (not list-of-lower-order-pev)
X    (multiple-value-setq (list-of-lower-order-phi list-of-lower-order-pev)
X      (step-down-Levinson-Durbin-recursions
X       coeffs variance
X       :process-variance-p process-variance-p)))
X  (let ((p (length coeffs))
X        (sd (sqrt (nth 0 list-of-lower-order-pev))))
X    (dotimes (i sample-size (values result
X                                    list-of-lower-order-phi
X                                    list-of-lower-order-pev))
X      (cond
X       ((< i p)
X        (dotimes (j i)
X          (incf (aref result i)
X                (* (aref result (- i j 1))
X                   (aref (nth (1- i) list-of-lower-order-phi) j))))
X        (incf (aref result i) (* sd (ranorm)))
X        (setf sd (sqrt (nth (1+ i) list-of-lower-order-pev))))
X       (t
X        (dotimes (j p)
X          (incf (aref result i) (* (aref result (- i j 1))
X                                   (elt coeffs j))))
X        (incf (aref result i) (* sd (ranorm))))))))
X
#|
;;; AR(2) process described by Equation (45) of the SAPA book:
(generate-ar-time-series #(0.75 -0.5) 1.0 4 :process-variance-p nil)
;==> #(-0.5263377203675295 -1.5334622586392335 -1.628199801947855 1.0935751830627396)
;    (#(0.5) #(0.75 -0.5))
;    (1.7777777777777777 1.3333333333333333 1.0)
X
;;; AR(4) process described by Equation (46a) of the SAPA book:
(generate-ar-time-series #(2.7607 -3.8106  2.6535 -0.9238)
X                         1.0 4 :process-variance-p nil)
;==> #(-6.887266139675448 -31.09831681568327 -34.70711181816937 -17.248058702547237)
X     (#(0.7164772070165697)
X      #(1.4197722065109775 -0.9816013581547822)
X      #(2.1105749802378733 -1.9807672315209488 0.7037508332562511)
X      #(2.7607 -3.8106 2.6535 -0.9238))
;    (761.7172900314962 370.6976500615111 13.515181723106767 6.821582066770185 1.0)
|#
X
;-------------------------------------------------------------------------------
(defun step-down-Levinson-Durbin-recursions
X       (coeffs
X        variance
X        &key
X        (process-variance-p t))
X  "given
X   [1] coeffs (required)
X       ==> sequence of length p with AR coefficients:
X           x_t = coeffs_0*x_{t-1} + coeffs_1*x_{t-2}
X                 + ... + coeffs_{p-1}*x_{t-p} + e_t
X           (see Equation (392a) in the SAPA book;
X           the coefficients can be real or complex-valued)
X   [2] variance (required)
X       ==> process variance or innovations variance
X           (see keyword process-variance-p)
X   [3] process-variance-p (keyword; t)
X       ==> if t, variance is taken to be process variance
X           if nil, variance is taken to be innovations variance
computes best linear prediction coefficients of orders 1, 2, ..., p-1
and prediction error variances of orders 0 (process variance), 1, ..., p and
returns
X   [1] list of vectors with best linear prediction coefficients
X       going from order 1 to order p;
X   [2] list of prediction error variances going from order 0 to order p
---
Note: see item [4] of the Comments and Extensions
X      to Section 9.4 of the SAPA book.  The values returned by
X      this function can be used to set the keyword parameters
X      list-of-lower-order-phi and list-of-lower-order-pev in
X      the function generate-ar-time-series"
X  (let* ((p (length coeffs))
X         (p-1 (1- p))
X         (list-of-lower-order-phi `(,coeffs))
X         (list-of-lower-order-pev `(,variance))
X         k-1 phi-k phi-k-1 den)
X    ;;; generate lower order coefficients via step-down Levinson-Durbin
X    (dotimes (i p-1)
X      (setf k-1 (- p-1 i)
X            phi-k (car list-of-lower-order-phi)
X            phi-k-1 (car (push (make-array k-1) list-of-lower-order-phi))
X            den (- 1.0 (realpart
X                        (* (aref phi-k k-1) (conjugate (aref phi-k k-1))))))
X      (dotimes (j k-1)
X        (setf (aref phi-k-1 j) (/ (+ (aref phi-k j)
X                                     (* (aref phi-k k-1)
X                                        (conjugate (aref phi-k (- k-1 j 1)))))
X                                  den))))
X    (values (if (plusp p) list-of-lower-order-phi)
X            (cond
X             (process-variance-p
X              (dotimes (i p (reverse list-of-lower-order-pev))
X                (push (* (car list-of-lower-order-pev)
X                         (- 1.0
X                            (realpart
X                             (* (aref (nth i list-of-lower-order-phi) i)
X                                (conjugate
X                                 (aref (nth i list-of-lower-order-phi) i))))))
X                      list-of-lower-order-pev)))
X             (t
X              (dotimes (i p list-of-lower-order-pev)
X                (push (/ (car list-of-lower-order-pev)
X                         (- 1.0
X                            (realpart
X                             (* (aref (nth (- p-1 i) list-of-lower-order-phi)
X                                      (- p-1 i))
X                                (conjugate
X                                 (aref (nth (- p-1 i) list-of-lower-order-phi)
X                                       (- p-1 i)))))))
X                      list-of-lower-order-pev)))))))
X
#|
;;; AR(2) process described by Equation (45) of the SAPA book:
(step-down-Levinson-Durbin-recursions #(0.75 -0.5)
X                                      1.0 :process-variance-p nil)
;==> (#(0.5) #(0.75 -0.5))
;    (1.7777777777777777 1.3333333333333333 1.0)
X
;;; AR(4) process described by Equation (46a) of the SAPA book:
(step-down-Levinson-Durbin-recursions #(2.7607 -3.8106  2.6535 -0.9238)
X                                      1.0 :process-variance-p nil)
;==> (#(0.7164772070165697)
X      #(1.4197722065109775 -0.9816013581547822)
X      #(2.1105749802378733 -1.9807672315209488 0.7037508332562511)
X      #(2.7607 -3.8106 2.6535 -0.9238))
;    (761.7172900314962 370.6976500615111 13.515181723106767 6.821582066770185 1.0)
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  simulate-time-series-from-sdf
;;;                 acvs-for-time-series-simulated-from-sdf
;;;  can be used to generate realizations of specified lengths
;;;  from normally distributed stationary processes with a specified
;;;  spectral density function or autocovariance sequence.  These functions
;;;  use frequency domain techniques.  The techniques are discussed
;;;  in the paper ``Simulating Gaussian Random Processes
;;;  with Specified Spectra'' by Percival in "Computing Science and Statistics",
;;;  vol 24, 534--8 (contains Proc. of 24th Symp. on the Interface of Computer
;;;  Science and Statistics).
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun simulate-time-series-from-sdf
X       (n-series
X        sdf
X        &key
X        (sampling-time 1.0)
X        (n-total
X         (* 4 (next-power-of-2 n-series)))
X        (result (make-array n-series)))
X  "given
X   [1] n-series (required)
X       ==> length of time series to be simulated
X   [2] sdf (required)
X       ==> spectral density function (sdf) defined for
X           0 <= f <= 1/(2.0 sampling-time) -- the
X           sdf is assumed to be two-sided and symmetric
X           about f = 0
X   [3] sampling-time (keyword; 1.0)
X       ==> the assumed sampling time (delta t)
X   [4] n-total (keyword; 4 * next power of 2 for n-series)
X       ==> a power of 2 controlling degree of accuracy
X           of the approximation (the larger, the better --
X           see Percival, 1992, for details)
X   [5] result (keyword; vector of length n-series)
X       <== a vector to contain simulated time series
returns
X   [1] a vector of length n-series generated from sdf
---
Note: the method used here is an approximate frequency domain
X      technique; in particular, it will NOT simulate the DC component
X      correctly, so beware of using this function in studies
X      where the process mean is of important"
X  ;;; This assertion fails if N is not a power of two 
X  (assert n-total
X          ()
X          "n-total = ~D not a power of 2"
X          n-total)
X  (let* ((n-sampling-time (float (* n-total sampling-time)))
X         (vector-for-fft (make-array n-total))
X         (j 0)
X         (n-j n-total)
X         (n-over-2 (/ n-total 2))
X         (1-over-sqrt-n-sampling-time (/ (sqrt (* n-total sampling-time)))))
X    (setf (svref vector-for-fft 0)
X          (* (sqrt (funcall sdf 0.0))
X             (ranorm))
X          (svref vector-for-fft n-over-2)
X          (* (sqrt (funcall sdf (/ 0.5 sampling-time)))
X             (ranorm)))
X    (dotimes (i (1- n-over-2))
X      (incf j) (decf n-j)
X      (setf (svref vector-for-fft n-j)
X            (conjugate
X             (setf (svref vector-for-fft j)
X                   (* (sqrt (* 0.5 (funcall sdf (/ j n-sampling-time))))
X                      (complex (ranorm) (ranorm)))))))
X    (fft! vector-for-fft)
X    (dotimes (i n-series result)
X      (setf (svref result i)
X            (* (realpart (svref vector-for-fft i))
X               1-over-sqrt-n-sampling-time)))))
X
#|
(simulate-time-series-from-sdf
X 4
X #'(lambda (f) (1+ f)))
;==> #(-0.775118374075842 -0.13513773684982933 -1.4345372479995293 0.31204493073149375)
|#
X
;-------------------------------------------------------------------------------
(defun acvs-for-time-series-simulated-from-sdf
X       (n-series
X        sdf
X        &key
X        (sampling-time 1.0)
X        (n-total (* 4 (next-power-of-2 n-series)))
X        (result (make-array n-series)))
X  "given
X   [1] n-series (required)
X       ==> length of time series simulated using
X           simulate-time-series-from-sdf
X           (must be a power of 2)
X   [2] sdf (required)
X       ==> spectral density function (sdf) defined for
X           0 <= f <= 1/(2.0 sampling-time) -- the
X           sdf is assumed to be two-sided and symmetric
X           about f = 0
X   [3] sampling-time (keyword; 1.0)
X       ==> the assumed sampling time (delta t)
X   [4] n-total (keyword; 4 * next power of 2 for n-series)
X       ==> a power of 2 controlling degree of accuracy
X           of the approximation (the larger, the better --
X           see Percival, 1992, for details)
X   [5] result (keyword; vector of length n-series)
X       <== a vector to contain simulated time series
returns
X   [1] a vector with the theoretical acvs from lag 0
X       to n-series - 1 for a time series generated via a call
X       to simulate-time-series-from-sdf with n-series,
X       sdf, sampling-time and n-total set to the same values"
X  (assert (power-of-2 n-total))
X  (let* ((sigma^2 (make-array n-total))
X         (n-total-sampling-time (* n-total sampling-time))
X         (j 1)
X         (N-j (1- n-total)))
X    (setf (aref sigma^2 0)
X          (/ (funcall sdf 0.0) n-total-sampling-time)
X          (aref sigma^2 (/ n-total 2))
X          (/ (funcall sdf (/ 0.5 sampling-time)) n-total-sampling-time))
X    (dotimes (i (- (/ n-total 2) 1))
X      (setf (aref sigma^2 j)
X            (/ (funcall sdf (/ j n-total-sampling-time))
X               n-total-sampling-time)
X            (aref sigma^2 N-j) (aref sigma^2 j))
X      (incf j)
X      (decf N-j))
X    (fft! sigma^2)
X    (dotimes (i n-series result)
X      (setf (aref result i) (realpart (aref sigma^2 i))))))
X
#|
(acvs-for-time-series-simulated-from-sdf
X 4
X #'(lambda (f) (1+ f)))
;==> #(1.25 -0.10263336862925071 0.0 -0.012655581284545123)
|#
SHAR_EOF
chmod 0644 random.lisp ||
echo 'restore of random.lisp failed'
Wc_c="`wc -c < 'random.lisp'`"
test 23412 -eq "$Wc_c" ||
	echo 'random.lisp: original size 23412, current size' "$Wc_c"
fi
# ============= sapa-package.lisp ==============
if test -f 'sapa-package.lisp' -a X"$1" != X"-c"; then
	echo 'x - skipping sapa-package.lisp (File already exists)'
else
echo 'x - extracting sapa-package.lisp (Text)'
sed 's/^X//' << 'SHAR_EOF' > 'sapa-package.lisp' &&
;;;-*- Mode: LISP; Package: :CL-USER; Syntax: COMMON-LISP -*-
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;
;  sapa-package.lisp
;
;  In some implementations of Lisp, the first line of a Lisp file (the so-called
;  mode line) is parsed automatically when opened by an editor, and,
;  if the package indicated by the mode line does not yet exist,
;  all sorts of strange things can happen.   To circumvent this annoyance,
;  we have created this short file that does nothing more than define
;  the SAPA package from within the CL-USER package, which -- according to
;  Steele2 -- is ALWAYS supposed to exist.  The idea is that you should
;  compile and load this file PRIOR to attempting to do anything with
;  any of the other files in sapaclisp contribution to StatLib.
;
;  6/3/93
;
;  SAPA, Version 1.0; Copyright 1993, Donald B. Percival, All Rights Reserved
;
;  Use and copying of this software and preparation of derivative works
;  based upon this software are permitted.  Any distribution of this
;  software or derivative works must comply with all applicable United
;  States export control laws.
; 
;  This software is made available AS IS, and no warranty -- about the
;  software, its performance, or its conformity to any
;  specification -- is given or implied.
;
;  Comments about this software can be addressed to dbp@apl.washington.edu
;-------------------------------------------------------------------------------
;;; (compile-file "home:SAPA;sapa-package.lisp")
;;; (load "home:SAPA;sapa-package.fasl")
;-------------------------------------------------------------------------------
(in-package :CL-USER)
X
(if (not (find-package :SAPA))
X  (defpackage :SAPA (:USE :COMMON-LISP)))
SHAR_EOF
chmod 0644 sapa-package.lisp ||
echo 'restore of sapa-package.lisp failed'
Wc_c="`wc -c < 'sapa-package.lisp'`"
test 1822 -eq "$Wc_c" ||
	echo 'sapa-package.lisp: original size 1822, current size' "$Wc_c"
fi
# ============= tapers.lisp ==============
if test -f 'tapers.lisp' -a X"$1" != X"-c"; then
	echo 'x - skipping tapers.lisp (File already exists)'
else
echo 'x - extracting tapers.lisp (Text)'
sed 's/^X//' << 'SHAR_EOF' > 'tapers.lisp' &&
;;;-*- Mode: LISP; Package: :SAPA; Syntax: COMMON-LISP -*-
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;
;  tapers.lisp
;
;  a collection of Lisp functions to taper a time series ...
;  Note:  before compiling and loading tapers.lisp,
;         you should compile and load (in the order listed)
;            sapa-package.lisp, utilities.lisp and basic-statistics.lisp
;
;  6/3/93
;
;  SAPA, Version 1.0; Copyright 1993, Donald B. Percival, All Rights Reserved
;
;  Use and copying of this software and preparation of derivative works
;  based upon this software are permitted.  Any distribution of this
;  software or derivative works must comply with all applicable United
;  States export control laws.
; 
;  This software is made available AS IS, and no warranty -- about the
;  software, its performance, or its conformity to any
;  specification -- is given or implied. 
;
;  Comments about this software can be addressed to dbp@apl.washington.edu
;-------------------------------------------------------------------------------
;;; (compile-file "ccl:SAPA;tapers.lisp")
;;; (load "ccl:SAPA;tapers.fasl")
;-------------------------------------------------------------------------------
(if (not (find-package :SAPA))
X  (defpackage :SAPA (:USE :COMMON-LISP)))
X
(in-package :SAPA)
X
(export '(;;; functions to compute data tapers discussed in the SAPA book ...
X          cosine-data-taper!
X          Hanning-data-taper!
X          dpss-data-taper!
X
X          ;;; function to taper a time series using a supplied vector ...
X          supplied-data-taper!
X          
X          ;;; function to center and taper a time series ...
X          center&taper-time-series
X          ))
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  cosine-data-taper!
;;;                 Hanning-data-taper!
;;;                 dpss-data-taper!
;;;  implement tapering of a time series (see Section 6.4 of the SAPA book).
;;;  In all cases, the time series that is supplied to these functions
;;;  are altered (hence the ``!'' is attached to all three names)
;;;  UNLESS the user supplies a value for the keyword parameter result.
;;;  All these functions return two values, namely, a vector containing
;;;  the tapered time series (same as the vector used as input to the
;;;  function unless the keyword parameter result is specified) and
;;;  C_h, the variance inflation factor due to tapering (see Equation (251b)
;;;  in the SAPA book.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun cosine-data-taper!
X       (a-time-series
X        &key
X        (taper-parameter 1.0)
X        (normalization :unity)
X        (result a-time-series))
X  "given
X   [1] a-time-series (required)
X       <=> a vector containing a time series;
X           this vector is modified UNLESS keyword result
X           is bound to a different vector
X   [2] taper-parameter (keyword; 1.0, i.e., Hanning data taper)
X       ==> a number >= 0.0 and <= 1.0 that specifies the degree
X           of tapering (parameter p of Equation (209)).
X   [3] normalization (keyword; :unity)
X       ==> if :unity, taper normalized such that
X           its sum of squares is equal to unity (Equation (208a));
X           if :N, taper normalized such that
X           its sum of squares is equal to the number of points
X           in the time series
X   [4] result (keyword; a-time-series)
X       <=> a vector to contain time series multiplied
X           by cosine data taper
returns
X   [1] a vector containing the tapered time series
X   [2] C_h, the variance inflation factor
X       computed using Equation (251b) in the SAPA book"
X  (let* ((n (length a-time-series))
X         (n-i-1 (1- n))
X         (sum-of-squares 0.0)
X         (C_h-top 0.0)
X         ;;; function ``1-'' converts low-cutoff to 0 based indexing
X         (low-cutoff (1- (/ (floor (* taper-parameter n)) 2)))
X         ;;; function ``1+'' converts low-cutoff to SAPA convention;
X         ;;; function ``1-'' converts high-cutoff to 0 based indexing
X         (high-cutoff (1- (- (1+ N) (1+ low-cutoff))))
X         (two-pi/denom (/ (* 2 pi) (1+ (floor (* taper-parameter n))))))
X    (dotimes (i (truncate n 2))
X      (let ((factor
X             (cond ((<= 0 i low-cutoff)
X                    (* 0.5 (- 1.0 (cos (* two-pi/denom (1+ i))))))
X                   ((< low-cutoff i high-cutoff)
X                    1.0)
X                   (t
X                    (error
X                     "fatal error in cosine-data-taper!")))))
X        (incf sum-of-squares (* 2 (expt factor 2)))
X        (incf C_h-top (* 2 (expt factor 4)))
X        (setf (aref result i) (* factor (aref a-time-series i))
X              (aref result n-i-1) (* factor (aref a-time-series n-i-1)))
X        (decf n-i-1)))
X    (if (oddp n)
X      (let ((factor
X             (cond ((<= 0 n-i-1 low-cutoff)
X                    (* 0.5 (- 1.0 (cos (* two-pi/denom (1+ n-i-1))))))
X                   ((< low-cutoff n-i-1 high-cutoff)
X                    1.0)
X                   (t
X                    (error
X                     "fatal error in cosine-data-taper!")))))
X        (incf sum-of-squares (expt factor 2))
X        (incf C_h-top (expt factor 4))
X        (setf (aref result n-i-1) (* factor (aref a-time-series n-i-1)))))
X    (let ((factor (sqrt (/ (if (eq normalization :unity) 1 N)
X                           sum-of-squares))))
X      (values (a*x! factor result)
X              (/ (* N C_h-top) (expt sum-of-squares 2))  ; C_h
X              ))))
X
#|
(dolist (N (list 100 1000))
X  (dolist (p (list 0.0 0.2 0.5 1.0))
X    (multiple-value-bind (taper C_h)
X                         (cosine-data-taper!
X                          (make-array N  :initial-element 1.0)
X                          :taper-parameter p)
X      (format t "~&;;; p = ~3,1F, N = ~D,   ~8,6F,   ~8,6F"
X              p N C_h (sum-of-squares taper)))))
X
;;; p = 0.0, N = 100,    1.000000,   1.000000
;;; p = 0.2, N = 100,    1.110360,   1.000000
;;; p = 0.5, N = 100,    1.338254,   1.000000
;;; p = 1.0, N = 100,    1.925193,   1.000000
X
;;; p = 0.0, N = 1000,   1.000000,   1.000000
;;; p = 0.2, N = 1000,   1.115727,   1.000000
;;; p = 0.5, N = 1000,   1.346217,   1.000000
;;; p = 1.0, N = 1000,   1.942502,   1.000000
X
;;; good agreement with Table 248 of SAPA
X
(dolist (N (list 3 33 99))
X  (dolist (p (list 0.0 0.2 0.5 1.0))
X    (multiple-value-bind (taper C_h)
X                         (cosine-data-taper!
X                          (make-array N  :initial-element 1.0)
X                          :taper-parameter p)
X      (format t "~&;;; p = ~3,1F, N = ~D,   ~8,6F,   ~8,6F"
X              p N C_h (sum-of-squares taper)))))
X
;;; p = 0.0, N = 3,    1.000000,   1.000000
;;; p = 0.2, N = 3,    1.000000,   1.000000
;;; p = 0.5, N = 3,    1.000000,   1.000000
;;; p = 1.0, N = 3,    1.500000,   1.000000
X
;;; p = 0.0, N = 33,   1.000000,   1.000000
;;; p = 0.2, N = 33,   1.087192,   1.000000
;;; p = 0.5, N = 33,   1.307487,   1.000000
;;; p = 1.0, N = 33,   1.887255,   1.000000
X
;;; p = 0.0, N = 99,   1.000000,   1.000000
;;; p = 0.2, N = 99,   1.105163,   1.000000
;;; p = 0.5, N = 99,   1.333636,   1.000000
;;; p = 1.0, N = 99,   1.925000,   1.000000
|#
X
;-------------------------------------------------------------------------------
(defun Hanning-data-taper!
X       (a-time-series
X        &key
X        taper-parameter
X        (normalization :unity)  ;or :N
X        (result a-time-series))
X  "calls cosine-data-taper! with taper-parameter set to 1.0"
X  (declare (ignore taper-parameter))
X  (cosine-data-taper!
X   a-time-series
X   :taper-parameter 1.0
X   :normalization normalization
X   :result result))
X
#|
(dolist (N (list 100 1000))
X  (multiple-value-bind (taper C_h)
X                       (Hanning-data-taper!
X                        (make-array N  :initial-element 1.0))
X    (format t "~&;;; p = ~3,1F, N = ~D,   ~8,6F,   ~8,6F"
X            1.0 N C_h (sum-of-squares taper))))
X
;;; p = 1.0, N = 100,    1.925193,   1.000000
;;; p = 1.0, N = 1000,   1.942502,   1.000000
X
(dolist (N (list 3 33 99))
X  (multiple-value-bind (taper C_h)
X                       (Hanning-data-taper!
X                        (make-array N  :initial-element 1.0))
X    (format t "~&;;; p = ~3,1F, N = ~D,   ~8,6F,   ~8,6F"
X            1.0 N C_h (sum-of-squares taper))))
X
;;; p = 1.0, N = 3,    1.500000,   1.000000
;;; p = 1.0, N = 33,   1.887255,   1.000000
;;; p = 1.0, N = 99,   1.925000,   1.000000
|#
X
;-------------------------------------------------------------------------------
(defun dpss-data-taper!
X       (a-time-series
X        &key
X        (taper-parameter 4.0)   ;NW
X        (normalization :unity)  ;or :N
X        (result a-time-series))
X  "given
X   [1] a-time-series (required)
X       <=> a vector containing a time series;
X           this vector is modified UNLESS keyword result
X           is bound to a different vector
X   [2] taper-parameter (keyword; 4.0)
X       ==> a number > 0.0 that specifies NW,
X           the product of the sample size and
X           the half-width of the concentration interval
X           (see SAPA, page 211)
X   [3] normalization (keyword; :unity)
X       ==> if :unity, taper normalized such that
X           its sum of squares is equal to unity (Equation (208a));
X           if :N, taper normalized such that
X           its sum of squares is equal to the number of points
X           in the time series
X   [4] result (keyword; a-time-series)
X       <=> a vector to contain time series multiplied
X           by cosine data taper
returns
X   [1] a vector containing the tapered time series
X   [2] C_h, the variance inflation factor
X       computed using Equation (251b) in the SAPA book"
X  (let* ((N (length a-time-series))
X         (N-i-1 (1- N))
X         (sum-of-squares 0.0)
X         (C_h-top 0.0)
X         (W (/ taper-parameter N))
X         (beta-pi (* pi W (1- N))))
X    (dotimes (i (truncate N 2))
X      (let ((factor
X             (bessi0-nr
X              (* beta-pi
X                 (sqrt (- 1.0 (expt (1- (float (/ (1+ (* 2 i)) N)))
X                                    2)))))))
X        (incf sum-of-squares (* 2 (expt factor 2)))
X        (incf C_h-top (* 2 (expt factor 4)))
X        (setf (aref result i) (* factor (aref a-time-series i))
X              (aref result N-i-1) (* factor (aref a-time-series N-i-1)))
X        (decf N-i-1)))
X    (if (oddp N)
X      (let ((factor
X             (bessi0-nr
X              (* beta-pi
X                 (sqrt (- 1.0 (expt (1- (float (/ (1+ (* 2  n-i-1)) N)))
X                                    2)))))))
X        (incf sum-of-squares (expt factor 2))
X        (incf C_h-top (expt factor 4))
X        (setf (aref result n-i-1) (* factor (aref a-time-series n-i-1)))))
X    (let ((factor (sqrt (/ (if (eq normalization :unity) 1 N)
X                           sum-of-squares))))
X      (values (a*x! factor result)
X              (/ (* N C_h-top) (expt sum-of-squares 2))  ; C_h
X              ))))
X
#|
(dolist (N (list 100 1000))
X  (dolist (NW (list 1.0 2.0 4.0 8.0))
X    (multiple-value-bind (taper C_h)
X                         (dpss-data-taper!
X                          (make-array N :initial-element 1.0)
X                          :taper-parameter NW)
X      (format t "~&;;; NW = ~3,1F, N = ~D,   ~8,6F,   ~8,6F"
X              NW N C_h (sum-of-squares taper)))))
X
;;; NW = 1.0, N = 100,    1.404251,   1.000000
;;; NW = 2.0, N = 100,    1.995994,   1.000000
;;; NW = 4.0, N = 100,    2.820333,   1.000000
;;; NW = 8.0, N = 100,    3.984642,   1.000000
X
;;; NW = 1.0, N = 1000,   1.410677,   1.000000
;;; NW = 2.0, N = 1000,   2.005059,   1.000000
;;; NW = 4.0, N = 1000,   2.833080,   1.000000
;;; NW = 8.0, N = 1000,   4.002674,   1.000000
X
;;; These are slightly off from Table 248 of SAPA
;;; (evidently computed using a less accurate approximation
;;; to 0th order dpss):
;;;  NW = 1,              1.34    [should be 1.41, a 5% error]
;;;  NW = 2,              1.96    [should be 2.01, a 3% error]
;;;  NW = 4,              2.80    [should be 2.83, a 1% error]
;;;  NW = 8,              3.94    [should be 4.00, a 1% error]
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The function   supplied-data-taper!
;;;  tapers a time series using a taper supplied in the form of a vector.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun supplied-data-taper!
X       (a-time-series
X        &key
X        (taper-parameter (Hanning-data-taper!
X                          (make-array (length a-time-series)
X                                      :initial-element 1.0)))
X        (normalization :unity)  ;or :N
X        (result a-time-series))
X  "given
X   [1] a-time-series (required)
X       <=> a vector containing a time series;
X           this vector is modified UNLESS keyword result
X           is bound to a different vector
X   [2] taper-parameter (keyword; vector with Hanning taper)
X       ==> a vector containing the supplied data taper
X           (must be the same length as a-time-series);
X           unmodified upon return
X   [3] normalization (keyword; :unity)
X       ==> if :unity, taper normalized such that
X           its sum of squares is equal to unity (Equation (208a));
X           if :N, taper normalized such that
X           its sum of squares is equal to the number of points
X           in the time series
X   [4] result (keyword; a-time-series)
X       <=> a vector to contain time series multiplied
X           by the supplied data taper
returns
X   [1] a vector containing the tapered time series
X   [2] C_h, the variance inflation factor
X       computed using Equation (251b) in the SAPA book"
X  (let ((N (length a-time-series))
X        (sum-of-squares (sum-of-squares taper-parameter)))
X    (x*y! a-time-series taper-parameter :result result)
X    (let ((factor (sqrt (/ (if (eq normalization :unity) 1 N)
X                           sum-of-squares))))
X      (values (a*x! factor result)
X              (/ (* N
X                    #-allegro
X                    (reduce #'+ taper-parameter
X                            :key #'(lambda (x) (expt x 4)))
X                    #+allegro
X                    (let ((SSSS 0.0)
X                          (j start))
X                      (dotimes (i N SSSS)
X                        (incf SSSS (expt (aref taper-parameter i) 4))))
X                    )
X                 (expt sum-of-squares 2))  ; C_h
X              ))))
X
#|
(supplied-data-taper! #(5 4 3 2 1)
X                      :taper-parameter #(0 1 2 1 0))
;==>
#(0.0 1.632993161855452 2.449489742783178 0.816496580927726 0.0)
5/2
X
(supplied-data-taper! #(5 4 3 2 1)
X                      :taper-parameter #(0 2 4 2 0))
;==>
#(0.0 1.632993161855452 2.449489742783178 0.816496580927726 0.0)
5/2
X
(supplied-data-taper! #(5 4 3 2 1)
X                      :taper-parameter #(0 1 2 1 0)
X                      :normalization :N)
;==>
#(0.0 3.6514837167011076 5.477225575051661 1.8257418583505538 0.0)
5/2
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The function   center&taper-time-series
;;;  centers and/or tapers a time series.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun center&taper-time-series
X       (time-series
X        &key
X        (center-data t)  ;t, nil or value to be subtracted off ...
X        (start 0)
X        (end (length time-series))
X        (data-taper nil)
X        (data-taper-parameters)
X        (recenter-after-tapering-p t)
X        (restore-power-option-p t)
X        (result (make-array (- end start))))
X  "given
X   [1] time-series (required)
X       ==> a sequence of time series values
X   [2] center-data (keyword; t)
X       ==> if t, subtract sample mean from time series;
X           if a number, subtracts that number from time series;
X           if nil, the effect of this function is to return
X           a copy of the relevant portion of time-series
X   [3] start (keyword; 0)
X       ==> start index of sequence to be used
X   [4] end (keyword; length of time-series)
X       ==> 1 + end index of sequence to be used
X   [5] data-taper (keyword; nil)
X       ==> nil or a tapering function
X   [6] data-taper-parameters (keyword)
X       ==> parameters for tapering function (not used
X           if data-taper is nil)
X   [7] recenter-after-tapering-p (keyword; t)
X       ==> if t and data-taper is a function,
X           centers tapered series by subtracting
X           off its sample mean
X   [8] restore-power-option-p (keyword; t)
X       ==> if t and data-taper is a function,
X           normalizes tapered series to have same
X           sum of squares as before tapering
X   [9] result (keyword; vector of size (- end start))
X       <== sequence to hold centered time series
returns
X   [1] result, the sequence containing the centered and/or
X       tapered time series (the contents of time-series
X       are unaltered unless it is also bound to result)
X   [2] the number used to center the time series
X       (this is nil if center-data is nil)
X   [3] C_h, the variance inflation factor due to the data taper
---
Note: see also center&prewhiten&taper-time-series
X      in filtering.lisp"
X  ;;; transfer relevant part of time-series to lower part of result
X  (if (not (and (eq time-series result) (zerop start)))
X    (copy-vector time-series result :start start :end end))
X  (let ((N (- end start))
X        (center-factor nil)
X        (C_h 1.0))
X    ;;; center the time series if required ...
X    (when center-data
X      (setf center-factor (if (numberp center-data)
X                            center-data
X                            (sample-mean result :end N)))
X      (dotimes (i N)
X        (decf (elt result i) center-factor)))
X    ;;; taper the time series if required ...
X    (when data-taper
X      (let ((sum-of-squares-before (if restore-power-option-p
X                                     (sum-of-squares result :end N))))
X        ;;; Note: a more elegant way to do this is to use the macro
X        ;;;       nth-value (Steele2, p. 184), but I have elected not
X        ;;;       to use it because [1] some versions of Common Lisp don't
X        ;;;       seem to have everything in Steele2 and [2] I don't
X        ;;;       want to get into defining any more macros than absolutely
X        ;;;       necessary (multf and divf are more than enough!).
X        (setf C_h (multiple-value-bind (junk non-junk)
X                                       (funcall data-taper
X                                                (make-array
X                                                 N
X                                                 :displaced-to result)
X                                                :taper-parameter
X                                                data-taper-parameters
X                                                :normalization :N)
X                    (declare (ignore junk))
X                    non-junk))
X        ;;; Note: Although the sample mean has already been subtracted
X        ;;; off, tapering reintroduces a nonzero mean in the tapered
X        ;;; series.  Here we remove it.  This amounts to estimating
X        ;;; the process mean using the alternative estimator discussed
X        ;;; at the end of Section 6.4 of SAPA.
X        (if recenter-after-tapering-p
X          (let ((recenter-factor (sample-mean result :end N)))
X            (dotimes (i N)
X              (decf (elt result i) recenter-factor))))
X        (if restore-power-option-p
X          (let ((mult-factor (sqrt (/ sum-of-squares-before
X                                      (sum-of-squares result :end N)))))
X            (dotimes (i N)
X              (multf (elt result i) mult-factor))))))
X    (values result center-factor C_h)))
X
#|
(center&taper-time-series #(71 63 70 88 99 90 110))
;==> #(-94/7 -150/7 -101/7 25/7 102/7 39/7 179/7)
;    591/7
;    1.0
X
(sum (center&taper-time-series #(71 63 70 88 99 90 110)))
;==> 0
X
(center&taper-time-series #(71 63 70 88 99 90 110) :start 1 :end 3)
;==> #(-7/2 7/2)
;    133/2
;    1.0
X
(let ((a #(71 63 70 88 99 90 110)))
X  (center&taper-time-series a :start 1 :end 3 :result a)
X  a)
;==> #(-7/2 7/2 70 88 99 90 110)
X
(sample-mean-and-variance #(71 63 70 88 99 90 110))
;==> 591/7
;    12304/49
(float 12304/49)
;==> 251.10204081632654
X
(center&taper-time-series
X #(71 63 70 88 99 90 110)
X :data-taper #'dpss-data-taper!
X :data-taper-parameters 1.0)
;==> #(-7.8468508813076685 -22.4213972133101 -19.89677763285509
X        6.044565561735114 21.153297157294492 6.4937737920259515
X       16.473389216417303)
;    591/7
;    1.3015152611329361
(sample-mean-and-variance (center&taper-time-series
X                           #(71 63 70 88 99 90 110)
X                           :data-taper #'dpss-data-taper!
X                           :data-taper-parameters 1.0))
;==> 0.0
;    251.10204081632654
X
(center&taper-time-series
X #(71 63 70 88 99 90 110)
X :data-taper #'dpss-data-taper!
X :data-taper-parameters 1.0
X :recenter-after-tapering-p nil)
;==> #(-8.369372277293243 -22.935860652629373 -20.412636880641728
X        5.514363897068667 20.61474219629164 5.963323769683637
X       15.937421676973303)
;    591/7
;    1.3015152611329361
(sample-mean-and-variance (center&taper-time-series
X                           #(71 63 70 88 99 90 110)
X                           :data-taper #'dpss-data-taper!
X                           :data-taper-parameters 1.0
X                           :recenter-after-tapering-p nil))
;==> -0.5268597529352993
;    250.82445961706344
(+ 250.82445961706344 (expt -0.5268597529352993 2))
;==> 251.10204081632648
X
(center&taper-time-series
X #(71 63 70 88 99 90 110)
X :data-taper #'dpss-data-taper!
X :data-taper-parameters 1.0
X :recenter-after-tapering-p nil
X :restore-power-option-p nil)
;==> #(-7.201263387970674 -19.73471463773794 -17.563655881256334
X        4.744727026618859 17.737553464239067 5.1310258058118645
X       13.713044111135646)
;    591/7
;    1.3015152611329361
(sample-mean-and-variance (center&taper-time-series
X                           #(71 63 70 88 99 90 110)
X                           :data-taper #'dpss-data-taper!
X                           :data-taper-parameters 1.0
X                           :recenter-after-tapering-p nil
X                           :restore-power-option-p nil))
;==> -0.453326214165644
;    185.695553371949
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  Everything below here consists of internal symbols in the SAPA package
;;;  and should be regarded as "dirty laundry" ...
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;; The Lisp function bessi0-nr is based on the Fortran function bessi0,
;;; Numerical Recipes, Second Edition, pages 230--1 (the coefficients are
;;; in fact from Abramowitz and Stegun).
(defun bessi0-nr (x)
X  ;;; translation of BESSI0, Numerical Recipes, p. 177
X  (cond
X   ((<= (abs x) 3.75)
X    (let ((y (expt (/ x 3.75) 2)))
X      (+ 1.0
X         (* y
X            (+ 3.5156229
X               (* y
X                  (+ 3.0899424
X                     (* y
X                        (+ 1.2067492
X                           (* y
X                              (+ 0.2659732
X                                 (* y
X                                    (+ 0.0360768
X                                       (* y 0.0045813))))))))))))))
X   (t
X    (let* ((ax (abs x))
X           (y (/ 3.75 ax)))
X      (* (/ (exp ax) (sqrt ax))
X         (+ 0.39894228
X            (* y
X               (+ 0.01328592
X                  (* y
X                     (+ 0.00225319
X                        (* y
X                           (+ -0.00157565
X                              (* y
X                                 (+ 0.00916281
X                                    (* y
X                                       (+ -0.02057706
X                                          (* y
X                                             (+ 0.02635537
X                                                (* y
X                                                   (+ -0.01647633
X                                                      (* y 0.00392377)
X                                                      ))))))))))))))))))))
X
#|
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;; This version is from Kaiser (1974), Proc. IEEE Int. Symp.
;;; on Circuits and Systems, pp. 20-23.  We use the routine from
;;; Numerical Recipes instead, but the difference between
;;; the two bessel functions seems to be quite small (-2e-7 and 9e-8).
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun mod-bessel-first-kind-of-0-order (x)
X  (do* ((d 0.0 (+ d 2.0))
X        (ds 1.0 (/ (* ds x x) (* d d)))
X        (s 1.0 (+ s ds)))
X       ((not (plusp (- ds (* 0.2e-8 s)))) s)))
|#
SHAR_EOF
chmod 0644 tapers.lisp ||
echo 'restore of tapers.lisp failed'
Wc_c="`wc -c < 'tapers.lisp'`"
test 25572 -eq "$Wc_c" ||
	echo 'tapers.lisp: original size 25572, current size' "$Wc_c"
fi
# ============= utilities.lisp ==============
if test -f 'utilities.lisp' -a X"$1" != X"-c"; then
	echo 'x - skipping utilities.lisp (File already exists)'
else
echo 'x - extracting utilities.lisp (Text)'
sed 's/^X//' << 'SHAR_EOF' > 'utilities.lisp' &&
;;;-*- Mode: LISP; Package: :SAPA; Syntax: COMMON-LISP -*-
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;
;  utilities.lisp
;
;  a collection of Lisp functions for certain basic operations ...
;  Note:  before compiling and loading utilities.lisp,
;         you should compile and load
;            sapa-package.lisp
;
;  6/3/93
;
;  SAPA, Version 1.0; Copyright 1993, Donald B. Percival, All Rights Reserved
;
;  Use and copying of this software and preparation of derivative works
;  based upon this software are permitted.  Any distribution of this
;  software or derivative works must comply with all applicable United
;  States export control laws.
; 
;  This software is made available AS IS, and no warranty -- about the
;  software, its performance, or its conformity to any
;  specification -- is given or implied. 
;
;  Comments about this software can be addressed to dbp@apl.washington.edu
;-------------------------------------------------------------------------------
;;; (compile-file "ccl:SAPA;utilities.lisp")
;;; (load "ccl:SAPA;utilities.fasl")
;-------------------------------------------------------------------------------
(if (not (find-package :SAPA))
X  (defpackage :SAPA (:USE :COMMON-LISP)))
X
(in-package :SAPA)
X
(export '(;;; two macros like incf for multiplication and division ...
X          multf
X          divf
X
X          ;;; functions dealing with logs and decibels ...
X          convert-to-dB
X          careful-convert-to-dB
X          convert-from-dB
X          log10
X          power-of-2
X          next-power-of-2
X          
X          ;;; functions to manipulate sequence(s) and return sequence(s) ...
X          a*x
X          a*x!
X          a*x+b
X          a*x+b!
X          a*x+y
X          a*x+y!
X          x+y
X          x+y!
X          x+b
X          x+b!
X          x-y
X          x-y!
X          x*y
X          x*y!
X          add-sequences
X          multiply-sequences
X          circular-shift-sequence
X          copy-vector
X          cumulative-sums
X          difference
X          transform-a-sequence
X          transform-a-sequence!
X          linear-interpolation!
X
X          ;;; functions to manipulate sequence(s) and return scalar(s) ...
X          dot-product
X          euclidean-norm
X          sum-of-squares
X          sum
X          min-and-max-of-seq
X          min-of-seq
X          max-of-seq
X          binary-search
X          compare-seqs
X
X          ;;; functions to generate sequence(s) ...
X          iota
X          sample-from-a-function
X
X          ;;; miscellaneous functions ...
X          sign
X          sampling-time->Nyquist-frequency
X          Nyquist-frequency->sampling-time
X          ))
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The macros  multf
;;;              divf
;;;  act like incf and decf in Common Lisp, but use multiplication
;;;  and division rather than addition and subtraction:
;;;    (multf x a)   is  somewhat like  (setf x (* x a))
;;;    (divf  x a)         ...          (setf x (/ x a))
;;;  These are the ONLY macros used in SAPA.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(define-modify-macro multf (&optional (delta 1)) * "Like <incf>")
X
(define-modify-macro divf (&optional (delta 1)) / "Like <incf>")
X
;-------------------------------------------------------------------------------
;;; The following constants are used with the next collection of functions:
(defconstant +sapa-10-over-log-10+ (/ 10.0d0 (log 10.0)))
(defconstant +sapa-2-pi+ (* 2 pi))
(defconstant +sapa-minus-2-pi+ (* -2 pi))
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  convert-to-dB
;;;                 convert-from-dB
;;;                 log10
;;;                 power-of-2
;;;                 next-power-of-2
;;;  all take a single number and perform some sort of a ``log or exponent''
;;;  operation on it:
;;;  convert-to-dB converts the number to decibels;
;;;  convert-from-dB converts the number back from its decibel formulation;
;;;  log10 returns the log base 10 of the number
;;;  power-of-2 tests the number to see if it is a power of 2; and
;;;  next-power-of-2 returns a power of 2 greater than or equal to the number.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun convert-to-dB (num)
X  "given
X   [1] num (required)
X       ==> a number
returns
X   [1] 10 log_10(num)
---
Note: num is intended to be a positive number,
but this is not checked"
X  (* (log num) +sapa-10-over-log-10+))
X
;;; (convert-to-dB 0.0001)  ;==> -39.999999999999986
X
;-------------------------------------------------------------------------------
(defun careful-convert-to-dB (num &optional (zero-mapping -100.0))
X  "given
X   [1] num (required)
X       ==> a number
X   [1] zero-mapping (required)
X       ==> a number
returns
X   [1] 10 log_10(num) if num > 0 or
X       zero-mapping   if num <= 0"
X  (if (plusp num)
X    (* (log num) +sapa-10-over-log-10+)
X    zero-mapping))
X
#|
(careful-convert-to-dB 0.0001)
;==> -39.999999999999986
(careful-convert-to-dB 0.0)
;==> -100.0
(careful-convert-to-dB 0.0 -1000.0)
;==> -1000.0
(careful-convert-to-dB 0.0 nil)
;==> nil
|#
X
;-------------------------------------------------------------------------------
(defun convert-from-dB (num)
X  "given
X   [1] num (required)
X       ==> a number
returns
X   [1] 10^(x/10)"
X  (exp (/ num +sapa-10-over-log-10+)))
X
;;; (convert-from-dB 20.0)  ;==> 100.00000000000004
X
;-------------------------------------------------------------------------------
(defun log10 (num)
X  "given
X   [1] num (required)
X       ==> a number
returns
X   [1] log_10(num)
---
Note: shorthand for (log num 10)"
X  (log num 10))
X
;;; (log10 100.0)  ;==> 2.0
X
;-------------------------------------------------------------------------------
(defun power-of-2 (n)
X  "given
X   [1] n (required)
X       ==> an integer
returns
X   [1] nil if n is not a power of 2 
X           -- or --
X       m   if n = 2^m for a positive integer m"
X  (let ((answer+1 (power-of-2+1 n)))
X    (if (numberp answer+1)
X      (1- answer+1))))
X
;;; (power-of-2 1023)  ;==> nil
;;; (power-of-2 32)    ;==> 5
X
;-------------------------------------------------------------------------------
(defun next-power-of-2 (n)
X  "given
X   [1] n (required)
X       ==> an integer
returns
X   [1] an integer power of 2 greater than or
X       equal to n"
X  (expt 2 (ceiling (log n 2))))
X
;;; (next-power-of-2 31)   ;==> 32
;;; (next-power-of-2 32)   ;==> 32
;;; (next-power-of-2 33)   ;==> 64
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  a*x
;;;                 a*x!
;;;                 a*x+b
;;;                 a*x+b!
;;;                 a*x+y
;;;                 a*x+y!
;;;                 x+y
;;;                 x+y!
;;;                 x+b
;;;                 x+b!
;;;                 x-y
;;;                 x-y!
;;;                 x*y
;;;                 x*y!
;;;                 add-sequences
;;;                 multiply-sequences
;;;                 circular-shift-sequence
;;;                 copy-vector
;;;                 cumulative-sums
;;;                 difference
;;;                 transform-a-sequence
;;;                 transform-a-sequence!
;;;                 linear-interpolation!
;;;  all take sequence(s) of numbers as input and return sequence(s)
;;;  of numbers as output.  Function names that end with a ``!''
;;;  indicate that the function will modify an input sequence
;;;  UNLESS the function has a keyword parameter with the name ``result''
;;;  with a value supplied by the user. 
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun a*x (a x)
X  "given
X   [1] a (required)
X       ==> a number
X   [2] x (required)
X       ==> a sequence of numbers
returns
X   [1] a vector of numbers obtained by multiplying
X       each element of x by a
----
Note: corresponds to SSCAL in LINPACK's bla."
X  (map 'vector #'(lambda (y) (* a y)) x))
X
#|
(setf *a* (vector 10.0 20.0 30.0))
(a*x pi *a*)  ;==> #(31.41592653589793 62.83185307179586 94.24777960769379)
*a*  ;==> #(10.0 20.0 30.0)
|#
X
;-------------------------------------------------------------------------------
(defun a*x! (a x &key (result x))
X  "given
X   [1] a (required)
X       ==> a number
X   [2] x (required)
X       ==> a sequence of numbers
X   [2] result (keyword; x)
X       <== a sequence of numbers
returns
X   [1] a sequence of numbers obtained by multiplying
X       each element of x by a
----
Note: corresponds to SSCAL in LINPACK's bla."
X  (map-into result #'(lambda (y) (* a y)) x))
X
#|
(setf *a* (vector 10.0 20.0 30.0))
(a*x! pi *a*)  ;==> #(31.41592653589793 62.83185307179586 94.24777960769379)
*a*  ;==> #(31.41592653589793 62.83185307179586 94.24777960769379)
(setf *a* (vector 10.0 20.0 30.0))
(setf *b* `(,nil ,nil ,nil))
(a*x! pi *a* :result *b*)
*a*  ;==> #(10.0 20.0 30.0)
*b*  ;==> (31.41592653589793 62.83185307179586 94.24777960769379)
|#
X
;-------------------------------------------------------------------------------
(defun a*x+b (a x b)
X  "given
X   [1] a (required)
X       ==> a number
X   [2] x (required)
X       ==> a sequence of numbers
X   [3] b (required)
X       ==> a number
returns
X   [1] a vector of numbers obtained by multiplying
X       each element of x by a and adding b"
X  (map 'vector #'(lambda (y) (+ b (* a y))) x))
X
#|
(setf *x* (vector  1.0  2.0  3.0))
(a*x+b pi *x* 10.0)
;==> #(13.141592653589793 16.283185307179586 19.42477796076938)
*x*  ;==> #(1.0 2.0 3.0)
|#
X
;-------------------------------------------------------------------------------
(defun a*x+b! (a x b &key (result x))
X  "given
X   [1] a (required)
X       ==> a number
X   [2] x (required)
X       ==> a sequence of numbers
X   [3] b (required)
X       ==> a number
X   [4] result (keyword; x)
X       <== a sequence of numbers
returns
X   [1] a sequence of numbers obtained by multiplying
X       each element of x by a and adding b"
X  (map-into result #'(lambda (y) (+ b (* a y))) x))
X
#|
(setf *x* (vector  1.0  2.0  3.0))
(a*x+b! pi *x* 10.0)
;==> #(13.141592653589793 16.283185307179586 19.42477796076938)
*x*  ;==> #(13.141592653589793 16.283185307179586 19.42477796076938)
|#
X
;-------------------------------------------------------------------------------
(defun a*x+y (a x y)
X  "given
X   [1] a (required)
X       ==> a number
X   [2] x (required)
X       ==> a sequence of numbers
X   [3] y (required)
X       ==> a sequence of numbers
returns
X   [1] a vector of numbers obtained by multiplying
X       each element of x by a and adding the
X       corresponding element of y
----
Note: corresponds to SAXPY in LINPACK's bla."
X  (map 'vector #'(lambda (u v) (+ (* a u) v)) x y))
X
#|
(setf *x* (vector  1.0  2.0  3.0))
(setf *y* (vector 10.0 20.0 30.0))
(a*x+y pi *x* *y*)
;==> #(13.141592653589793 26.283185307179586 39.424777960769376)
*x*  ;==> #(1.0 2.0 3.0)
*y*  ;==> #(10.0 20.0 30.0)
|#
X
;-------------------------------------------------------------------------------
(defun a*x+y! (a x y &key (result y))
X  "given
X   [1] a (required)
X       ==> a number
X   [2] x (required)
X       ==> a sequence of numbers
X   [3] y (required)
X       ==> a sequence of numbers
X   [4] result (keyword; y)
X       <== a sequence of numbers
returns
X   [1] a sequence of numbers obtained by multiplying
X       each element of x by a and adding the
X       corresponding element of y
----
Note: corresponds to SAXPY in LINPACK's bla."
X  (map-into result #'(lambda (u v) (+ (* a u) v)) x y))
X
#|
(setf *x* (vector  1.0  2.0  3.0))
(setf *y* (vector 10.0 20.0 30.0))
(a*x+y! pi *x* *y*)
;==> #(13.141592653589793 26.283185307179586 39.424777960769376)
*x*  ;==> #(1.0 2.0 3.0)
*y*  ;==> #(13.141592653589793 26.283185307179586 39.424777960769376)
|#
X
;-------------------------------------------------------------------------------
(defun x+y (x y)
X  "given
X   [1] x (required)
X       ==> a sequence of numbers
X   [2] y (required)
X       ==> a sequence of numbers
returns
X   [1] a vector of numbers obtained by adding
X       x and y together element by element"
X  (map 'vector #'+ x y))
X
#|
(setf *x* (vector  1.0  2.0  3.0))
(setf *y* (vector 10.0 20.0 30.0))
(x+y *x* *y*)
;==> #(11.0 22.0 33.0)
*x*  ;==> #(1.0 2.0 3.0)
*y*  ;==> #(10.0 20.0 30.0)
|#
X
;-------------------------------------------------------------------------------
(defun x+y! (x y &key (result x))
X  "given
X   [1] x (required)
X       ==> a sequence of numbers
X   [2] y (required)
X       ==> a sequence of numbers
X   [3] result (keyword; x)
X       <== a sequence of numbers
returns
X   [1] a sequence of numbers obtained by adding
X       x and y together element by element"
X  (map-into result #'+ x y))
X
#|
(setf *x* (vector  1.0  2.0  3.0))
(setf *y* (vector 10.0 20.0 30.0))
(x+y! *x* *y*)
;==> #(11.0 22.0 33.0)
*x*  ;==> #(11.0 22.0 33.0)
*y*  ;==> #(10.0 20.0 30.0)
|#
X
;-------------------------------------------------------------------------------
(defun x+b (x b)
X  "given
X   [1] x (required)
X       ==> a sequence of numbers
X   [2] b (required)
X       ==> a number
returns
X   [1] a vector of numbers obtained by adding
X       b to each element of x"
X  (map 'vector #'(lambda (y) (+ b y)) x))
X
#|
(setf *x* (vector  1.0  2.0  3.0))
(x+b *x* 10.0)
;==> #(11.0 12.0 13.0)
*x*  ;==> #(1.0 2.0 3.0)
|#
X
;-------------------------------------------------------------------------------
(defun x+b! (x b &key (result x))
X  "given
X   [1] x (required)
X       ==> a sequence of numbers
X   [2] b (required)
X       ==> a number
X   [3] result (keyword; x)
X       <== a sequence of numbers
returns
X   [1] a sequence of numbers obtained by adding
X       b to each element of x"
X  (map-into result #'(lambda (y) (+ b y)) x))
X
#|
(setf *x* (vector  1.0  2.0  3.0))
(x+b! *x* 10.0)
;==> #(11.0 12.0 13.0)
*x*  ;==> #(11.0 12.0 13.0)
|#
X
;-------------------------------------------------------------------------------
(defun x-y (x y)
X  "given
X   [1] x (required)
X       ==> a sequence of numbers
X   [2] y (required)
X       ==> a sequence of numbers
returns
X   [1] a vector of numbers obtained by subtracting
X       elements of y from corresponding elements of x"
X  (map 'vector #'- x y))
X
#|
(setf *x* (vector  1.0  2.0  3.0))
(setf *y* (vector 10.0 20.0 30.0))
(x-y *x* *y*)
;==> #(-9.0 -18.0 -27.0)
*x*  ;==> #(1.0 2.0 3.0)
*y*  ;==> #(10.0 20.0 30.0)
|#
X
;-------------------------------------------------------------------------------
(defun x-y! (x y &key (result x))
X  "given
X   [1] x (required)
X       ==> a sequence of numbers
X   [2] y (required)
X       ==> a sequence of numbers
X   [3] result (keyword; x)
X       <== a sequence of numbers
returns
X   [1] a sequence of numbers obtained by subtracting
X       elements of y from corresponding elements of x"
X  (map-into result #'- x y))
X
#|
(setf *x* (vector  1.0  2.0  3.0))
(setf *y* (vector 10.0 20.0 30.0))
(x-y! *x* *y*)
;==> #(-9.0 -18.0 -27.0)
*x*  ;==> #(-9.0 -18.0 -27.0)
*y*  ;==> #(10.0 20.0 30.0)
|#
X
;-------------------------------------------------------------------------------
(defun x*y (x y)
X  "given
X   [1] x (required)
X       ==> a sequence of numbers
X   [2] y (required)
X       ==> a sequence of numbers
returns
X   [1] a vector of numbers obtained by multiplying
X       corresponding elements of x and y"
X  (map 'vector #'* x y))
X
#|
(setf *x* (vector  1.0  2.0  3.0))
(setf *y* (vector 10.0 20.0 30.0))
(x*y *x* *y*)
;==> #(10.0 40.0 90.0)
*x*  ;==> #(1.0 2.0 3.0)
*y*  ;==> #(10.0 20.0 30.0)
|#
X
;-------------------------------------------------------------------------------
(defun x*y! (x y &key (result x))
X  "given
X   [1] x (required)
X       ==> a sequence of numbers
X   [2] y (required)
X       ==> a sequence of numbers
X   [3] result (keyword; x)
X       <== a sequence of numbers
returns
X   [1] a sequence of numbers obtained by multiplying
X       corresponding elements of x and y"
X  (map-into result #'* x y))
X
#|
(setf *x* (vector  1.0  2.0  3.0))
(setf *y* (vector 10.0 20.0 30.0))
(x*y! *x* *y*)
;==> #(10.0 40.0 90.0)
*x*  ;==> #(10.0 40.0 90.0)
*y*  ;==> #(10.0 20.0 30.0)
|#
X
;-------------------------------------------------------------------------------
(defun add-sequences (&rest seqs)
X  "given
X   [1] a arbitrary number of sequences of numbers,
X       all of the same size
returns
X   [1] a new sequence formed by adding the sequences
X       together on an element by element basis"
X  (apply #'map (type-of (car seqs)) #'+ seqs))
X
;;; (add-sequences '(1 2 3) '(4 5 6) '(7 8 9))  ;==>  (12 15 18)
;;; (add-sequences #(1 2 3) #(4 5 6) #(7 8 9))  ;==> #(12 15 18)
X
;-------------------------------------------------------------------------------
(defun multiply-sequences (&rest seqs)
X  "given
X   [1] a arbitrary number of sequences of numbers,
X       all of the same size
returns
X   [1] a new sequence formed by multiplying the sequences
X       together on an element by element basis"
X  (apply #'map (type-of (car seqs)) #'* seqs))
X
;;; (multiply-sequences '(1 2 3) '(4 5 6) '(7 8 9))  ;==>  (28 80 162)
;;; (multiply-sequences #(1 2 3) #(4 5 6) #(7 8 9))  ;==> #(28 80 162)
X
;-------------------------------------------------------------------------------
(defun circular-shift-sequence
X       (a-seq
X        &key
X        (result
X         (make-array (length a-seq))))
X  "given
X   [1] a-seq (required)
X       ==> a sequence
X   [2] result (keyword; vector of same length as a-seq)
X       <== a sequence
returns
X   [1] a `circularly' left-shifted sequence; i.e.,
X       maps   x(0),   x(1), x(2), ..., x(n-2), x(n-1)
X       to     x(n-1), x(0), x(1), ..., x(n-3), x(n-2)
----
Note: result can be the same sequence as result"
X  (let* ((n (length a-seq))
X         (n-1 (1- n))
X         (temp (elt a-seq n-1)))
X    (dotimes (i (1- n))
X      (setf (elt result n-1)
X            (elt a-seq (decf n-1))))
X    (setf (elt result 0) temp)
X    result))
X
#|
(setf *a* (vector 10.0 20.0 30.0 40.0 50.0 60.0 70.0))
(circular-shift-sequence *a* :result *a*)
;==> #(70.0 10.0 20.0 30.0 40.0 50.0 60.0)
(circular-shift-sequence *a*)
;==> #(60.0 70.0 10.0 20.0 30.0 40.0 50.0)
*a*
;==> #(70.0 10.0 20.0 30.0 40.0 50.0 60.0)
|#
X
;-------------------------------------------------------------------------------
;;; The name of this function is a misleading in that it can accept
;;; any type of sequences (not just vectors).
(defun copy-vector
X       (from-vector
X        to-vector
X        &key
X        (start 0)
X        (end (length from-vector)))
X  "given
X   [1] from-vector (required)
X       ==> a vector (actually, works with any sequence)
X   [2] to-vector (required)
X       ==> a vector (actually, works with any sequence)
X   [3] start (keyword; 0)
X       ==> start index of from-vector to be transferred
X   [4] end (keyword; length of from-vector)
X       ==> 1 + end index of from-vector to be transferred
returns
X   [1] to-vector, with elements start to end - 1
X       from-vector copied into corresponding
X       elements 0 to end - start - 1 of to-vector
---
Note: if to-vector is being created on the spot,
X      might consider using CL's copy-seq instead."
X  (replace to-vector from-vector :start2 start :end2 end))
X
#|
(copy-vector #(1 2 3) (make-array 3))
;==> #(1 2 3)
(copy-vector #(1 2 3) (make-array 6 :initial-element 0) :start 1)
;==> #(2 3 0 0 0 0)
(copy-vector #(1 2 3) (make-array 6 :initial-element 0) :end 2)
;==> #(1 2 0 0 0 0)
|#
X
;-------------------------------------------------------------------------------
(defun cumulative-sums
X       (the-seq
X        &key
X        (start 0)
X        (end (length the-seq))
X        (result (make-array (- end start))))
X  "given
X   [1] the-seq (required)
X       ==> a sequence of numbers
X   [2] start (keyword; 0)
X       ==> start index of sequence to be used
X   [3] end (keyword; length of the-seq)
X       ==> 1 + end index of sequence to be used
returns
X   [1] sequence with cumulative sum
X       of elements in specified subsequence"
X  (let ((N (- end start)))
X    (if (plusp N) (setf (elt result 0) (elt the-seq start)))
X    (dotimes (i (1- N) result)
X      (setf (elt result (1+ i)) (+ (elt result i)
X                                   (elt the-seq (incf start)))))))
X
;;; (cumulative-sums #(10 20 30 40 55))  ;==> #(10 30 60 100 155)
;;; (cumulative-sums #(10 20 30 40 55) :start 2)  ;==> #(30 70 125)
X
;-------------------------------------------------------------------------------
(defun difference
X       (sequence-of-numbers
X        &key
X        (start 0)
X        (end (length sequence-of-numbers))
X        (lag 1)
X        (result (make-array (- end start lag))))
X  "given
X   [1] sequence-of-numbers (required)
X       ==> a sequence of numbers
X   [2] start (keyword; 0)
X       ==> start index
X   [3] end (keyword; length of sequence-of-numbers)
X       ==> end index plus one
X   [4] lag (keyword; 1)
X       ==> lag for differencing
X   [5] result (keyword; vector of length (- end start lag))
X       <== results
return
X   [1] a sequence of length (- end start lag)
X       with lag-th order difference of sequence
X       from index start to index end-1
---
Note: See Newton, 1988, page 31."
X  (if (or (not (= start 0))
X          (not (= end (length sequence-of-numbers))))
X    (setf sequence-of-numbers (subseq sequence-of-numbers start end)))
X  (let ((k lag))
X    (dotimes (j (- end start lag) result)
X      (setf (elt result j) (- (elt sequence-of-numbers k)
X                              (elt sequence-of-numbers j)))
X      (incf k))))
X
;;; (difference #(1 2 4 5 7 9 11))  ;==> #(1 2 1 2 2 2)
X
;-------------------------------------------------------------------------------
(defun transform-a-sequence (a-function the-seq)
X  "given
X   [1] a-function (required)
X       ==> a function with one argument
X   [2] the-seq (required)
X       ==> a sequence of numbers
returns
X   [1] a sequence of numbers obtained by applying
X       a-function to each element of the-seq"
X  (map (type-of the-seq) a-function the-seq))
X
;;; (transform-a-sequence #'convert-to-dB '(1 10 100)) ;==>  (0.0 10.0 20.0)
;;; (transform-a-sequence #'convert-to-dB #(1 10 100)) ;==> #(0.0 10.0 20.0)
X
;-------------------------------------------------------------------------------
(defun transform-a-sequence!
X       (a-function
X        the-seq
X        &key
X        (result the-seq))
X  "given
X   [1] a-function (required)
X       ==> a function with one argument
X   [2] the-seq (required)
X       ==> a sequence of numbers
X   [3] result (keyword; the-seq)
X       <== a sequence of numbers
returns
X   [1] a sequence of numbers obtained by applying
X       a-function to each element of the-seq"
X  (map-into result a-function the-seq))
X
#|
(setf *a* #(1 10 100))
(transform-a-sequence! #'convert-to-dB *a*)  ;==> #(0.0 10.0 20.0)
*a*  ;==> #(0.0 10.0 20.0)
|#
X
;-------------------------------------------------------------------------------
(defun linear-interpolation!
X       (a-sequence
X        &key
X        (interpolation-predicate #'zerop))
X  "given
X   [1] a-sequence (required)
X       <=> a sequence of numbers
X   [2] interpolation-predicate (keyword; #'zerop)
X       ==> a unitary predicate which, when true,
X           indicates that a value in a-sequence
X           is to be replaced by a linear interpolation
returns
X   [1] a-sequence, modified with linear interpolates
X       over specified values"
X  (let ((n (length a-sequence)))
X    ;;; check to make sure that the beginning and the end of the sequence
X    ;;; are well-defined and that the sequence is long enough to do something
X    ;;; with
X    (assert (and
X             (> n 2)
X             (not (funcall interpolation-predicate (elt a-sequence 0)))
X             (not (funcall interpolation-predicate (elt a-sequence (1- n))))))
X    (do ((i-start 0)
X         (i-stop (1- n))
X         i-first-bad i-last-bad
X         i-good-before-first-bad
X         i-good-after-last-bad
X         rise run)
X        ((>= i-start i-stop) a-sequence)   ;;; exit test
X      (block this-block
X        (dotimes (i (1+ (- i-stop i-start)) (setf i-start (1+ i-stop)))
X          (when (funcall interpolation-predicate
X                         (elt a-sequence (+ i i-start)))
X            (setf i-first-bad (+ i i-start)
X                  i-good-before-first-bad (1- i-first-bad))
X            (dotimes (j (- i-stop i-first-bad)
X                        (error "algorithm failure!!! --- evacuate building!!!"))
X              (when (not (funcall interpolation-predicate
X                                  (elt a-sequence (+ j i-first-bad 1))))
X                (setf i-last-bad (+ j i-first-bad)
X                      i-good-after-last-bad (1+ i-last-bad)
X                      rise (- (elt a-sequence i-good-after-last-bad)
X                              (elt a-sequence i-good-before-first-bad))
X                      run (- i-good-after-last-bad i-good-before-first-bad))
X                (dotimes (k (- i-good-after-last-bad i-first-bad))
X                  (setf (elt a-sequence (+ i-first-bad k))
X                        (+ (elt a-sequence i-good-before-first-bad)
X                           (* (/ (1+ k) run) rise))))
X                (setf i-start (1+ i-good-after-last-bad))
X                (return-from this-block)))))))))
X
#|
(setf *a* #(1 10 0 100))
(linear-interpolation! *a*)  ;==> #(1 10 55 100)
*a*  ;==> #(1 10 55 100)
(setf *b* #(1 10 nil 100))
(linear-interpolation!
X *b*
X :interpolation-predicate #'null)  ;==> #(1 10 55 100)
*b*  ;==> #(1 10 55 100)
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  dot-product
;;;                 euclidean-norm
;;;                 sum-of-squares
;;;                 sum
;;;                 min-and-max-of-seq
;;;                 min-of-seq
;;;                 max-of-seq
;;;                 binary-search
;;;                 compare-seqs
;;;  all take sequence(s) of numbers as input and return various scalar(s)
;;;  of interest.  For example, the dot-product between two sequences
;;;  of the same length.
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun dot-product (x y)
X  "given
X   [1] x (required)
X       ==> a sequence of numbers
X   [2] y (required)
X       ==> a sequence of numbers
returns
X   [1] the dot product of x and y
----
Note: corresponds to SDOT in LINPACK's bla"
X  (let ((sum 0.0))
X    (dotimes (i (length x) sum)
X      (incf sum (* (elt x i) (elt y i))))))
X
;;; (dot-product #(1 2 3) #(4 5 6))  ;==> 32.0
X
;-------------------------------------------------------------------------------
(defun euclidean-norm (x)
X  "given
X   [1] x (required)
X       ==> a sequence of numbers
returns
X   [1] Euclidean norm of x
----
Note: corresponds to LINPACK's SNRM2,
X      but this is POOR implementation!"
X  (sqrt (dot-product x x)))
X
;;; (euclidean-norm #(4 3))  ;==> 5.0
X
;-------------------------------------------------------------------------------
(defun sum-of-squares
X       (the-seq
X        &key
X        (start 0)
X        (end (length the-seq)))
X  "given
X   [1] the-seq (required)
X       ==> a sequence of numbers
X   [2] start (keyword; 0)
X       ==> start index of sequence to be used
X   [3] end (keyword; length of the-seq)
X       ==> 1 + end index of sequence to be used
returns
X   [1] sum of squares of elements
X       in specified subsequence"
X  #-allegro
X  (reduce #'+ the-seq :start start :end end :key #'(lambda (x) (* x x)))
X  #+allegro
X  (let ((SS 0.0)
X        (j start))
X    (dotimes (i (- end start) SS)
X      (incf SS (expt (elt the-seq j) 2))
X      (incf j)
X      ))
X  )
X
#|
(sum-of-squares #(4 3))                               ;==> 25
(sum-of-squares #(1.0 2.0 3.0))                       ;==> 14.0
(sum-of-squares #(1.0 2.0 3.0 4.0))                   ;==> 30.0
(sum-of-squares #(1.0 2.0 3.0 4.0) :end 3)            ;==> 14.0
(sum-of-squares #(1.0 2.0 3.0 4.0) :start 1 :end 3)   ;==> 13.0
|#
X
;-------------------------------------------------------------------------------
(defun sum
X       (the-seq
X        &key
X        (start 0)
X        (end (length the-seq)))
X  "given
X   [1] the-seq (required)
X       ==> a sequence of numbers
X   [2] start (keyword; 0)
X       ==> start index of sequence to be used
X   [3] end (keyword; length of the-seq)
X       ==> 1 + end index of sequence to be used
returns
X   [1] sum of elements in specified subsequence"
X  (reduce #'+ the-seq :start start :end end))
X
;;; (sum #(4 3))  ;==> 7
;;; (sum #(7 8 3 4 3 10 11))  ;==> 46
;;; (sum #(7 8 3 4 3 10 11) :start 3 :end 5)  ;==> 7
X
;-------------------------------------------------------------------------------
(defun min-and-max-of-seq
X       (the-seq
X        &key
X        (start 0)
X        (end (length the-seq)))
X  "given
X   [1] the-seq (required)
X       ==> a sequence of real-valued numbers
X   [2] start (keyword; 0)
X       ==> start index of sequence to be used
X   [3] end (keyword; length of the-seq)
X       ==> 1 + end index of sequence to be used
returns
X   [1] minimum value in sequence
X   [2] maximum value in sequence"
X  (values  (reduce #'min the-seq :start start :end end)
X           (reduce #'max the-seq :start start :end end)))
X
;;; (min-and-max-of-seq #(7 8 3 4 3 10 11))  ;==> 3 and 11
X
;-------------------------------------------------------------------------------
(defun min-of-seq
X       (the-seq
X        &key
X        (start 0)
X        (end (length the-seq)))
X  "given
X   [1] the-seq (required)
X       ==> a sequence of real-valued numbers
X   [2] start (keyword; 0)
X       ==> start index of sequence to be used
X   [3] end (keyword; length of the-seq)
X       ==> 1 + end index of sequence to be used
returns
X   [1] minimum value in sequence
X   [2] index of minimum value"
X  (let ((the-min (reduce #'min the-seq :start start :end end)))
X    (values the-min (position the-min the-seq :start start :end end))))
X
;;; (min-of-seq #(7 8 3 4 3 10 11))  ;==> 3 and 2
X
;-------------------------------------------------------------------------------
(defun max-of-seq
X       (the-seq
X        &key
X        (start 0)
X        (end (length the-seq)))
X  "given
X   [1] the-seq (required)
X       ==> a sequence of real-valued numbers
X   [2] start (keyword; 0)
X       ==> start index of sequence to be used
X   [3] end (keyword; length of the-seq)
X       ==> 1 + end index of sequence to be used
returns
X   [1] maximum value in sequence
X   [2] index of maximum value"
X  (let ((the-max (reduce #'max the-seq :start start :end end)))
X    (values the-max (position the-max the-seq :start start :end end))))
X
;;; (max-of-seq #(7 8 3 4 3 10 11))  ;==> 11 and 6
X
;-------------------------------------------------------------------------------
(defun binary-search
X       (value
X        ordered-seq
X        &key
X        (lower-test #'<=)
X        (upper-test #'<=))
X  "given
X   [1] value (required)
X       ==> a real-valued number
X   [2] ordered-seq (required)
X       ==> an ordered sequence of real-valued numbers
X   [3] lower-test (keyword; #'<=)
X       ==> a predicate of two arguments
X           that defines the lower test
X   [4] upper-test (keyword; #'<=)
X       ==> a predicate of two arguments
X           that defines the upper test
returns
X   [1] if lower-test and upper-test are set
X       to their default values, this function
X       returns an index i such that v[i] <= value <= v[i+1],
X       where 0 <= i <= n-2 (if there is no such i, nil is returned);
X       if nondefault values are supplied for lower-test and upper-test,
X       then the condition `v[i] <= value <= v[i+1]'
X       gets modified in an obvious way
---
Note: value can be an arbitrary object, and ordered-seq
can be a list of arbitrary objects if the binary predicates
lower-test and upper-test are appropriate set"
X  (binary-search-internal
X   value ordered-seq 0 (1- (length ordered-seq)) lower-test upper-test))
X
;;; (binary-search 8 #(1 2 4 5 7 9 11))  ;==> 4
X
;-------------------------------------------------------------------------------
(defun compare-seqs (one-seq another-seq)
X  "given
X   [1] one-seq (required)
X       ==> any sequence
X   [2] another-seq (required)
X       ==> any sequence (must be same length as one-seq)
returns
X   [1] maximum absolute difference between corresponding
X       elements of the two sequences
X   [2] average absolute difference
---
Note: useful for comparing results that in theory
should be identical"
X  (assert (= (length one-seq) (length another-seq)))
X  (let ((max-abs-diff 0.0)
X        (ave-abs-diff 0.0)
X        (n (length one-seq))
X        abs-diff)
X    (dotimes (i n (values (/ ave-abs-diff n) max-abs-diff))
X      (setf abs-diff (abs (- (elt one-seq i) (elt another-seq i))))
X      (when (> abs-diff max-abs-diff) (setf max-abs-diff abs-diff))
X      (incf ave-abs-diff abs-diff))))
X
#|
(compare-seqs
X #(1 2 3 4 5 6   7 8 9 10)
X #(1 2 3 4 5 6.1 7 8 9 10))
;==>
0.009999999999999964
0.09999999999999964
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  iota
;;;                 sample-from-a-function
;;;  both generate a sequence of numbers.  In the case of iota, the sequence
;;;  of numbers is equally spaced; in the case of sample-from-a-function,
;;;  the sequence is obtained by sampling a function at a specified set
;;;  of points
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun iota
X       (start
X        end
X        &key
X        (type 'vector)
X        (float-p nil)
X        (skip 1)
X        (result (make-sequence type (1+ (/ (- end start) skip)))))
X  "given
X   [1] start (required)
X       ==> an integer
X   [2] end (required)
X       ==> an integer >= start
X   [3] type (keyword; 'vector)
X       ==> type of sequence desired
X   [4] float-p (keyword; nil)
X       ==> if true, sequence consists of floats
X   [5] skip (keyword; 1)
X       ==> desired increment between integers in sequence
X   [6] result (keyword; sequence of length (1+ (- start end)))
X       ==> storage space for result
returns
X   [1] a sequence of length (1+ (/ (- end start) skip))
X       with numbers start, start + skip, ..., end - skip, end
X       (these are stored in whatever result points to)"
X  (if float-p
X    (dotimes (i (1+ (/ (- end start) skip)) result)
X      (setf (elt result i) (float (+ start (* i skip)))))
X    (dotimes (i (1+ (/ (- end start) skip)) result)
X      (setf (elt result i) (+ start (* i skip))))))
X
;;; (iota 7 10 :type 'vector)  ;==> #(7 8 9 10)
X
;-------------------------------------------------------------------------------
(defun sample-from-a-function
X       (a-function
X        &key
X        (seq-of-x-values '())
X        (x0 0.0)
X        (delta-x 1.0)
X        (n (length seq-of-x-values)))
X  "given
X   [1] a-function (required)
X       ==> a function to be sampled from
X   [2] seq-of-x-values (keyword; '())
X       ==> sequence of values at which a-function
X           is to be sampled (default '())
X   [3] x0 (keyword; 0.0)
X       ==> point of first value to be sampled
X   [4] delta-x (keyword; 1.0)
X       ==> increment between points
X   [5] n (keyword; length of seq-of-x-values)
X       ==> number of samples
returns
X   [1] an array with the values of a-function
X       at a specified set of points; and
X   [2] an array with the specified set of points.
---
Note that the specified set of points either is in seq-of-x-values
X -- if it is non-nil -- or is given by
x0, x0 + delta-x, x0 + 2*delta-x, ..., x0 + (n-1)*delta-x"
X  (let ((the-samples (make-array n))
X        (x (if seq-of-x-values
X             (if (arrayp seq-of-x-values)
X               seq-of-x-values
X               (make-array n :initial-contents seq-of-x-values))
X             (make-array n))))
X    (cond
X     (seq-of-x-values
X      (dotimes (i n (values the-samples x))
X        (setf (aref the-samples i)
X              (funcall a-function (elt seq-of-x-values i)))))
X     (t
X      (dotimes (i n (values the-samples x))
X        (setf (aref x i) x0)
X        (setf (aref the-samples i) (funcall a-function x0))
X        (incf x0 delta-x))))))
X
#|
(sample-from-a-function #'cos
X                        :delta-x (/ (* 2 pi) 4)
X                        :n 5)
;==>
#(1.0 6.123031769111886E-17 -1.0 -1.836909530733566E-16 1.0)
#(0.0 1.5707963267948966 3.141592653589793 4.71238898038469 6.283185307179586)
|#
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  The functions  sign
;;;                 sampling-time->Nyquist-frequency
;;;                 Nyquist-frequency->sampling-time
;;;  are, respectively, a Lisp implementation of the Fortran SIGN function
;;;  (useful because it differs from the Common Lisp function signum);
;;;  a simple function that takes the sampling time (called delta t
;;;  in the SAPA book) and returns the Nyquist frequency; and another
;;;  simple function for going from the Nyquist frequency to the
;;;  sampling time
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
(defun sign (a1 a2)
X  "implementation of Fortran SIGN function"
X  (if (minusp a2)
X    (- (abs a1))
X    (abs a1)))
X
;-------------------------------------------------------------------------------
(defun sampling-time->Nyquist-frequency (sampling-time)
X  "given a sampling time, returns the associated Nyquist frequency
---
Note: see page 98 of the SAPA book"
X  (/ (* 2 sampling-time)))
X
;-------------------------------------------------------------------------------
(defun Nyquist-frequency->sampling-time (Nyquist-frequency)
X  "given a Nyquist frequency, returns the associated sampling-time
---
Note: see page 98 of the SAPA book"
X  (/ (* 2 Nyquist-frequency)))
X
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;;  Everything below here consists of internal symbols in the SAPA package
;;;  and should be regarded as "dirty laundry" ...
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
;;; Used by power-of-2
(defun power-of-2+1 (n)
X  (cond
X   ((or (not (numberp n)) (not (plusp n)) (< n 1))
X    nil)
X   ((= n 1) 1)
X   (t
X    (let ((local-answer (power-of-2+1 (/ n 2))))
X      (if (numberp local-answer)
X        (1+ local-answer))))))
X
;-------------------------------------------------------------------------------
;;; NOTE: used by binary-search.  This arrangement was done because the way
;;;       that this function is recursively called makes the lambda list
;;;       bulky as a user interface.  Originally, we had &optional instead of
;;;       &key in binary-search, and there was no binary-search-internal.
;;;       There is probably a slick way of getting around this is Lisp,
;;;       but I (dbp) don't know what it is as of 7/29/91!
;;; NOTE: `end' here is one less than the usual Lisp parameter of that name
(defun binary-search-internal
X       (value
X        ordered-seq
X        start
X        end
X        lower-test
X        upper-test)
X  (cond
X   ((= (- end start) 1)
X    (if (and (funcall lower-test (elt ordered-seq start) value)
X             (funcall upper-test value (elt ordered-seq end)))
X      start))   ;returns nil if termination condition fails
X   (t
X    (let ((middle (truncate (+ start end) 2)))
X      (if (funcall upper-test value (elt ordered-seq middle))
X        (binary-search-internal
X         value ordered-seq start middle lower-test upper-test)
X        (binary-search-internal
X         value ordered-seq middle end lower-test upper-test)
X        )))))
SHAR_EOF
chmod 0644 utilities.lisp ||
echo 'restore of utilities.lisp failed'
Wc_c="`wc -c < 'utilities.lisp'`"
test 40336 -eq "$Wc_c" ||
	echo 'utilities.lisp: original size 40336, current size' "$Wc_c"
fi
exit 0

