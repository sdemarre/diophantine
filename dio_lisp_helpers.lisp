(in-package :maxima)

;; see Solving the generalized Pell equation x^2 - Dy^2 = N, John P. Robertson, July 31, 2004, Pages 4 - 8. http://www.jpr2718.org/pell.pdf
;; computes the simple continued fraction expansion of the quadratic irrational (P+sqrt(D))/Q
;; returns a maxima list of equations with values which are arrays
(defun $dio_pqa (p q d)
  "computes the continued fraction of (p+sqrt(d))/q"
  (if (not (zerop (mod (- d (* p p)) q)))
      ($dio_pqa (* p (abs q)) (* q (abs q)) (* d q q))
      (let ((qd (sqrt (* d 1.0d0)))
	    (p (list p))
	    (q (list q)))
	(let ((a (list 1 0))
	      (b (list 0 1))
	      (g (list (car q) (- (car p))))
	      as
	      xi
	      cxi
	      period-started
	      loop-found
	      integer-partial
	      periodic-partials
	      non-periodic-partials
	      period-initial-p
	      period-initial-q)
	  (do ((i 0 (incf i)))
	      (loop-found)
	    (progn
	      (when (> i 0)
		(push (- (* (car as) (car q)) (car p)) p)
		(push (/ (- d (* (car p) (car p))) (car q)) q)
		(when period-started
		  (setf loop-found (and (= period-initial-p (car p)) (= period-initial-q (car q))))))
	      (push (floor (+ (car p) qd) (car q)) as)
	      (push (+ (* (car as) (car a)) (cadr a)) a)
	      (push (+ (* (car as) (car b)) (cadr b)) b)
	      (push (+ (* (car as) (car g)) (cadr g)) g)
	      (push (/ (+ (car p) qd) (car q)) xi)
	      (push (/ (- (car p) qd) (car q)) cxi)
	      (when (zerop i)
		(push (car as) integer-partial))
	      (when (not period-started)
		(when (and (> (car xi) 1) (< (car cxi) 0) (< -1 (car cxi)))
		  (setf period-started t
			period-initial-p (car p)
			period-initial-q (car q))))
	      (if period-started
		  (push (car as) periodic-partials)
		  (unless (zerop i)
		    (push (car as) non-periodic-partials)))))
	  (macrolet ((mlist (sym l)
		       `(list '(mequal simp) ,sym (cons '(mlist simp) (reverse ,l)))))
	    (list '(mlist simp)
		  (mlist '$a as)
		  (mlist '$ip integer-partial)
		  (mlist '$ap non-periodic-partials)
		  (mlist '$pp (rest periodic-partials))
		  (mlist '$p p)
		  (mlist '$q q)
		  (mlist '$g g)
		  (mlist '$b b)
		  (mlist '$xi xi)
		  (mlist '$cxi cxi)))))))

(defun $dio_find_1 (a b c cf)
  (let ((cf (rest cf)))
    (flet ((p (y z)
	     (let ((r (+ (* a y y) (* b y z) (* c z z))))
	       (format t "p(~a,~a)=~a~%" y z r)
	       r)))
      (let* ((result (list '(mlist simp)))
	     (ynm1 1)
	     (znm1 0)
	     (yn (car cf))
	     (zn 1))
	(loop for cn in (rest cf) do (progn
				       (when (= 1 (p yn zn))
					 (setf result (list '(mlist simp) yn zn))
					 (return result))
				       (format t "cn=~a~%" cn)
				       (psetf yn (+ ynm1 (* cn yn))
					      zn (+ znm1 (* cn zn))
					      ynm1 yn
					      znm1 zn)))
	result))))

(defun $dio_point_ranges (point-data)
  (let ((point-data (rest (mapcar #'rest point-data))))
    (list '(mlist simp)
	  (reduce #'min point-data :key #'first)
	  (reduce #'max point-data :key #'first)
	  (reduce #'min point-data :key #'second)
	  (reduce #'max point-data :key #'second))))

(defun $dio_brute_force (a b c d e f limit &optional pos)
  (flet ((f (x y) (+ (* a x x) (* b x y) (* c y y) (* d x) (* e y) f)))
    (cons '(mlist simp)
          (reverse
           (let (result)
             (loop for x from (if pos 1 (- limit)) to limit do
                  (loop for y from (if pos 1 (- limit)) to limit do
                       (when (zerop (f x y))
                         (push (list '(mlist simp) x y) result))))
             result)))))

(defun $dio_compute_k_power (a b q l)
  "Finds k so that (a+sqrt(q)*b)^k = 1 mod l. Assumes such a k exists. If not,
this function will search forever."
  (let ((ca a)
	(cb b))
    (do* ((k 1 (1+ k)))
	 ((and (= (mod ca l) 1) (= (mod cb l) 0)) k)
      (psetf ca (mod (+ (* a ca) (* b cb q)) l)
	     cb (mod (+ (* a cb) (* b ca)) l)))))

(defun pqa-get (pqa var)
  (let ((var-eq-list (find-if #'(lambda (candidate) (eq (second candidate) var)) (rest pqa))))
    (when var-eq-list (rest (third var-eq-list)))))
(defun $square_divisors (n)
  (cons '(mlist simp) (loop for f from 1 to (isqrt (abs n)) when (zerop (mod n (* f f))) collect f)))
(defun min-pos-pell-values (d n)
  (rest (mcall '$map '$rhs (mcall '$dio_min_pos_pell_solution '$t '$u d n))))
;;  lmm as described in "Solving the generalized Pell equation" 2004 by John P. Robertson, http://www.jpr2718.org/pell.pdf
(defun $dio_lmm (x y d n)
  (let (result)
    (loop for f in (rest ($square_divisors n)) do
	 (let* ((m (/ n (* f f)))
		(absm (abs m))
		(limit (/ absm 2)))
	   (loop for z from (1+ (truncate (- limit))) to (1+ (truncate limit)) do
		(when (= (mod (* z z) absm) (mod d absm))
		  (let ((pqa ($dio_pqa z absm d))
			(found nil))
		    (loop for qi in (rest (pqa-get pqa '$q))
		       and r in (cddr (pqa-get pqa '$g))
		       and s in (cddr (pqa-get pqa '$b))
		       until found
		       do (when (or (= qi 1) (= qi -1))
			    (setf found t)
			    (if (= (- (* r r) (* d s s)) m)
				(push (list '(mlist simp)
					    (list '(mequal simp) x (* f r))
					    (list '(mequal simp) y (* f s)))
				      result)
				(let ((sols (min-pos-pell-values d -1)))
				  (when sols
				    (destructuring-bind (v w) sols
				      (push (list '(mlist simp)
						  (list '(mequal simp) x (* f (+ (* r v) (* s w d))))
						  (list '(mequal simp) y (* f (+ (* r w) (* s v)))))
					    result))))))))))))
    (cons '(mlist simp) result)))



(defun dio_integer_pell_array (px py lx ly d k x-trans-coeffs y-trans-coeffs)
  (destructuring-bind (txx txy txc txl) x-trans-coeffs
    (destructuring-bind (tyx tyy tyc tyl) y-trans-coeffs
      (let ((cx lx)
            (cy ly)
            result)
        (flet ((xt-int-p () (zerop (mod (+ (* txx cx) (* txy cy) txc) txl)))
               (yt-int-p () (zerop (mod (+ (* tyx cy) (* tyy cy) tyc) tyl))))
          (loop for p from 0 to (1- k) do (progn
                                            (push (and (xt-int-p) (yt-int-p)) result)
                                            (psetf cx (+ (* cx px) (* cy py d))
                                                   cy (+ (* cy px) (* cx py))))))))))

(defun check-mod-solution (m a b c d e f)
  (loop for x from 0 below m do
       (let ((z (+ (* a x x) (* d x) f))
	     (u (+ (* b x) e)))
	 (loop for y from 0 below m do
	      (when (zerop (mod (+ z (* y (+ u (* c y)))) m))
		(return-from check-mod-solution)))))
  t)

(defun $dio_check_mod (eq)
  (let ((coeffs (mapcar #'caddr (rest (mcall '$dio_coeffs eq)))))
    (not (or (apply #'check-mod-solution (cons 4 coeffs))
	     (apply #'check-mod-solution (cons 16 coeffs))
	     (apply #'check-mod-solution (cons 25 coeffs))))))

(defun qmul (s u q x y)
  "computes (t+sqrt(q)*u)*(x+sqrt(q)*y)"
  (cons (+ (* s x) (* q u y)) (+ (* s y) (* u x))))

(defmacro qmulf (s u q x y)
  `(psetf ,s (+ (* ,s ,x) (* ,q ,u ,y))
	  ,u (+ (* ,s ,y) (* ,u ,x))))

(defun optimize-stream (gen k)
  (cons gen k))

;; receive i -> add expression i+j*k but check if it can be combined with an existing expression first

