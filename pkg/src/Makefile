
all:	workall

workall:
	(cd gsl-1.9; ./configure --prefix=${PWD}/gslinstall --exec-prefix=${PWD}/gslinstall; make; make install prefix=${PWD}/gslinstall exec_prefix=${PWD}/gslinstall; cd ../gslinstall; mv include ..; cd lib; cp *.so ../..)
	(cd CLAPACK-3.2.1; mkdir SRC/VARIANTS/LIB; make f2clib; make blaslib; make tmglib; make lapacklib; cp *.so ../; make; cp INCLUDE/*.h ../)
	(cp Makevars MakeSplits; cd MakeSplits; make; cp ttsplits.so ../; make clean; rm Makevars)	
	(cp Makevars Regression; cd Regression;  make; cp shotgun.so ../; make clean; rm Makevars)
	(cp Makevars Cvv; cd Cvv; make; cp cvv.so ../; make clean; rm Makevars)
	(cp Makevars TrainTest; cd TrainTest; make; cp traintest.so ../; make clean; rm Makevars)

clean:
	(rm -f *.so)

cleanall:
	(cd gsl-1.9; make clean)
	(rm -rf gslinstall)
	(cd CLAPACK-3.2.1; make cleanall; rmdir SRC/VARIANTS/LIB)
	(rm -f *.so; rm -f *.h; rm -rf include) 
	(cd Regression; make clean)
	(cd MakeSplits; make clean)
	(cd Cvv; make clean)
	(cd TrainTest; make clean)

