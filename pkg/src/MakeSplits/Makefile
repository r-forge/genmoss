SPLITS_OBJS = ttsplits.so

all:	${SPLITS_OBJS}
	rm -f *.o

model.o: Model.cpp Model.h
	g++ -m64 -g -c -fPIC Model.cpp -o model.o

table.o: newtable.cpp newtable.h
	g++ -m64 -g -c -fPIC newtable.cpp -o table.o

data.o: Data.cpp Data.h
	g++ -m64 -g -c -fPIC Data.cpp -o data.o

ttsplits.so: splits.cpp splits.h data.o model.o table.o
	R CMD SHLIB splits.cpp data.o model.o table.o -o ttsplits.so

clean:
	rm -f *.o
	rm -f ${SPLITS_OBJS}
