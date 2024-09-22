# 文件名: velocity.mk
PETSC_DIR = /share/software/user/open/petsc/3.18.5
PETSC_ARCH =

include ${PETSC_DIR}/lib/petsc/conf/variables
Compiler = gcc
CFLAGS = ${PETSC_CC_INCLUDES} ${CXX_FLAGS} ${CXXFLAGS} ${CPPFLAGS} ${PSOURCECXX} -std=c++11


all: velocity_solver


velocity_solver: frameworkcpp.o settings.o flowstationary.o primary.o
	$(Compiler) frameworkcpp.o settings.o flowstationary.o primary.o -o velocity_solver ${PETSC_LIB} $(CFLAGS)


frameworkcpp.o: frameworkcpp.cpp
	$(Compiler) -c frameworkcpp.cpp $(CFLAGS)

settings.o: settings.cpp
	$(Compiler) -c settings.cpp $(CFLAGS)

flowstationary.o: flowstationary.cpp
	$(Compiler) -c flowstationary.cpp $(CFLAGS)

primary.o: primary.cpp
	$(Compiler) -c primary.cpp $(CFLAGS)


clean:
	rm -f *.o velocity_solver
