CXX = g++
CXXFLAGS = -std=c++11 -Wall -Ofast
PROMG = main_g
PROMM = main_m
OBJSG = Grid.o GVF_BC.o U_BC.o P_BC.o Mixture.o Eigen_Utilities.o\
	THINC.o WENO_1D.o Pressure.o main_g.o
OBJSM = Grid.o GVF_BC.o U_BC.o P_BC.o Mixture.o Eigen_Utilities.o\
	THINC.o WENO_1D.o Pressure.o main_m.o

${PROMG}: ${OBJSG}
	${CXX} -o ${PROMG} ${OBJSG}

${PROMM}: ${OBJSM}
	${CXX} -o ${PROMM} ${OBJSM}

clean:
	rm -f *.o *.ascii ${PROMG} ${PROMM}