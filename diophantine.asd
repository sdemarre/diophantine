(defsystem diophantine :depends-on ("maxima")
	   :defsystem-depends-on ("maxima-file")
	   :name "diophantine"
	   :version "2016.08.21"
	   :maintainer "Serge De Marre"
	   :author "Serge De Marre"
	   :licence ""
	   :description "diophantine equation solver"
	   :long-description "Maxima program to solve diophantine equations of the form ax^2+bxy+cy^2+dx+ey+f=0 with a,b,c,d,e,f constant integers.
Based on Dario Alpern's solution/code found at https://www.alpertron.com.ar/QUAD.HTM"

	   :components
	   ((:file "dio_lisp_helpers")
	    (:maxima-file "diophantine")
	    (:maxima-file "diophantine-devel")
	    (:maxima-file "diophantine-draw")))
