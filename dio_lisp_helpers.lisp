(in-package :maxima)
(ql:quickload :iterate)

(defun dio-cf-to-maxima-result (a b result)
  (let (non-repeating-part
	repeating-part
	in-repeating-part)
    (iterate:iter (iterate:for (n na . nb) in result)
		  (unless in-repeating-part
		    (setf in-repeating-part (and (= na a) (= nb b))))
		  (if (not in-repeating-part)
		      (push n non-repeating-part)
		      (push n repeating-part)))
    (list '(mlist simp) (cons '(mlist simp) (reverse non-repeating-part)) (cons '(mlist simp) (reverse repeating-part)))))
(defun dio-next-cf (a b d qd)
  (let ((f (floor (- qd b) (* 2 a))))
    (cons f (cons (+ (* -1 a f f) (* -1 b f) (/ (- d (* b b)) 4 a))
		  (- (* -2 a f) b)))))
(defun $dio_cf_expansion (a b d)
  (let ((dio-cf-hash (make-hash-table :test 'equal))
	(qd (sqrt d))
	result)
    (iterate:iter (iterate:while (not (numberp (gethash (cons a b) dio-cf-hash))))
		  (let ((next-step (dio-next-cf a b d qd)))
		    (setf (gethash (cons a b) dio-cf-hash) (car next-step))
		    (push (cons (car next-step) (cons a b)) result)
		    (setf a (cadr next-step)
			  b (cddr next-step))))
    (dio-cf-to-maxima-result a b (reverse result))))


(defun $dio_point_ranges (point-data)
  (let ((point-data (rest (mapcar #'rest point-data))))
    (list '(mlist simp)
	  (reduce #'min point-data :key #'first)
	  (reduce #'max point-data :key #'first)
	  (reduce #'min point-data :key #'second)
	  (reduce #'max point-data :key #'second))))
