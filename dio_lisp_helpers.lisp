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