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

(defun $dio_point_ranges (point-data)
  (let ((point-data (rest (mapcar #'rest point-data))))
    (list '(mlist simp)
	  (reduce #'min point-data :key #'first)
	  (reduce #'max point-data :key #'first)
	  (reduce #'min point-data :key #'second)
	  (reduce #'max point-data :key #'second))))

;; ;;
;; 	   (+ ,@(unless (zerop a) (list (list '* a 'x 'x)))
;; 	      ,@(unless (zerop b) (list (list '* b 'x 'y)))
;; 	      ,@(unless (zerop c) (list (list '* c 'y 'y)))
;; 	      ,@(unless (zerop d) (list (list '* d 'x)))
;; 	      ,@(unless (zerop e) (list (list '* e 'y)))
;; 	      ,@(unless (zerop f) (list f)))
(defun make-bf-checker-form (a b c d e f)
  `(lambda (limit pos)
     (declare (optimize (speed 3) (safety 0))
	      (type (integer -10000000 10000000) limit))
     (let (result
	   (start (if pos 0 (- limit))))
       (declare (type (integer -10000000 10000000) start))
       (loop for x fixnum from start to limit do
	    (let ((fxy (+ (* ,a (the fixnum (* x x)))
			  (* ,b x start)
			  (* ,c (the fixnum (* start start)))
			  (* ,d x)
			  (* ,e start)
			  ,f)))
	      (declare (type fixnum fxy))
	      (loop for y fixnum from start to limit do
		   (when (zerop fxy)
		     (push (list '(mlist simp) x y) result))
		   (incf fxy (+ ,@(unless (zerop c) (list (list '* (* 2 c) 'y)))
				,@(unless (zerop b) (list (list '* b 'x)))
				,@(unless (zerop e) (list e))
				,@(unless (zerop c) (list c)))))))
       (cons '(mlist simp) result))))
(defun $dio_brute_force (a b c d e f limit &optional pos)
  (with-output-to-string (*error-output*)
    (compile 'bf-fun (make-bf-checker-form a b c d e f)))
  (funcall 'bf-fun limit pos))

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
           (format t " [d=~a, m=~a, range=[~a,~a]]~%" d absm (1+ (truncate (- limit))) (1+ (truncate limit)))
           (when (mcall '$zn_nth_root d 2 absm)
            (loop for z from (1+ (truncate (- limit))) to (1+ (truncate limit)) do
                 (when (= (mod (* z z) absm) (mod d absm))
                   (format t "  lmm found ~a^2=~a (mod ~a)~%" z d absm)
                   (let ((pqa ($dio_pqa z absm d))
                         (found nil))
                     (loop for qi in (rest (pqa-get pqa '$q))
                        and r in (cddr (pqa-get pqa '$g))
                        and s in (cddr (pqa-get pqa '$b))
                        until found
                        do (when (or (= qi 1) (= qi -1))
                             (setf found t)
                             (if (= (- (* r r) (* d s s)) m)
                                 (pushnew (list '(mlist simp)
                                                (list '(mequal simp) x (* f r))
                                                (list '(mequal simp) y (* f s)))
                                          result :test #'equalp)
                                 (let ((sols (min-pos-pell-values d -1)))
                                   (when sols
                                     (destructuring-bind (v w) sols
                                       (pushnew (list '(mlist simp)
                                                      (list '(mequal simp) x (* f (+ (* r v) (* s w d))))
                                                      (list '(mequal simp) y (* f (+ (* r w) (* s v)))))
                                                result :test #'equalp)))))))))))))
    (cons '(mlist simp) (remove-duplicates result))))

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

;; receive i -> add expression i+j*k but check if it can be combined with an existing expression first

(defun all-elements-in-ipowers (k divisor increment ipowers)
  (loop for i from 0 to (1- (/ k divisor)) always (gethash (+ (* i divisor) increment) ipowers)))
(defun unmark-all-elements (k divisor increment ipowers)
  (loop for i from 0 to (1- (/ k divisor)) do (remhash (+ (* i divisor) increment) ipowers)))
(defun $find_min_coverage (integer-powers k k-divisors)
  (let* ((ipowers (make-hash-table :size (length integer-powers)))
	 result)
    (loop for i in (rest integer-powers) do (setf (gethash i ipowers) t))
    (loop for divisor in (rest k-divisors) do
	 (unless (zerop (hash-table-count ipowers))
	   (loop for increment from 0 to (1- divisor) do
		(when (all-elements-in-ipowers k divisor increment ipowers)
		  (unmark-all-elements k divisor increment ipowers)
		  (push (list '(mlist simp) divisor increment) result)
		  (when (zerop (hash-table-count ipowers))
		    (return))))))
    (cons '(mlist simp) (reverse result))))

(defun $find_int_powers_lisp (txx txy txc txd tyx tyy tyc tyd min-pell-x min-pell-y fund-x fund-y q k l sign1 sign2)
  (let ((cx fund-x)
        (cy fund-y)
        result)
    (flet ((transform-x () (/ (+ (* sign1 cx txx) (* sign2 cy txy) txc) txd))
           (transform-y () (/ (+ (* sign1 cx tyx) (* sign2 cy tyy) tyc) tyd)))
     (loop for pow from 0 to (1- k) do
          (let ((actual-x (transform-x))
                (actual-y (transform-y)))
            (when (and (integerp actual-x) (integerp actual-y))
              (push pow result))
            (psetf cx (mod (+ (* q cy min-pell-y) (* cx min-pell-x)) l)
                   cy (mod (+ (* cy min-pell-x) (* cx min-pell-y)) l)))))
    (cons '(mlist simp) (reverse result))))


;; unmarked = length(integer-powers)
;; for d in divisors
;;    for increment from 0 to d-1
;;       if all-elements-integer(k,d,increment, integer-powers)
;;          count-marked = mark-all-elements(k,d,increment, integer-powers)
;;          add-expression-to-result(k,d,increment)
;;          unmarked -= count-marked
;;          if unmarked == 0
;;              return result


(defun check-bf-fun (limit &optional pos)
  (declare (optimize (speed 3) (safety 0))
           (type (integer -100000000 100000000) limit))
  (let (result
	(y-limit (if pos 0 (- limit))))
    (format t "y-limit = ~a~%" y-limit)
    (loop for x fixnum from y-limit to limit
       do (let ((fxy
		 (the fixnum (+ (the fixnum (* -425 (the fixnum (* x x))))
				(* 9 (the fixnum (* y-limit y-limit)))
				(* 42 y-limit)
				(* 5100 x)
				-15226))))
	    (declare (type fixnum fxy))
	    (loop for y fixnum from y-limit to limit
	       do (progn
		    (when (zerop fxy)
		      (push (list '(mlist simp) x y) result))
		   (incf fxy (+ (* 2 9 y) 51))))))
    (cons '(mlist simp) result)))
