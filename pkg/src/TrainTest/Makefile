RJMCMC_OBJS = traintest.so

all:	${RJMCMC_OBJS}
	rm -f *.o

params.o: Params.cpp Params.h
	g++ -m64 -g -c -fPIC Params.cpp -o params.o

table.o: newtable.cpp newtable.h
	g++ -m64 -g -c -fPIC newtable.cpp -o table.o

data.o: Data.cpp Data.h
	g++ -m64 -g -c -fPIC Data.cpp -o data.o

graph.o: newgraph.cpp graph.h
	g++ -m64 -g -c -fPIC newgraph.cpp -o graph.o
		
model.o: model.cpp model.h
	g++ -m64 -g -c -fPIC model.cpp -o model.o
	echo ${R_HOME}
	
traintest.so: grjmcmc.cpp grjmcmc.h graph.o table.o data.o model.o params.o
	${R_HOME}/bin/R CMD SHLIB grjmcmc.cpp graph.o table.o data.o model.o params.o -o traintest.so
		
clean:
		rm -f *.o
		rm -f ${RJMCMC_OBJS}

