PETSC_DIR=/share/software/user/open/petsc/3.18.5
PETSC_ARCH=

include ${PETSC_DIR}/lib/petsc/conf/variables
Compiler=gcc
CFLAGS= ${PETSC_CC_INCLUDES} ${CXX_FLAGS} ${CXXFLAGS} ${CPPFLAGS} ${PSOURCECXX} -std=c++11$

all: primarydiff

primarydiff: BasicDataStructure.o settingdiff.o primarydiff.o main.o 
	$(Compiler) BasicDataStructure.o settingdiff.o primarydiff.o main.o -o primarydiff ${PETSC_LIB} $(CFLAGS)
	
BasicDataStructure.o: BasicDataStructure.cpp
	$(Compiler)	-c BasicDataStructure.cpp $(CFLAGS)
	
settingdiff.o: settingdiff.cpp
	$(Compiler)	-c settingdiff.cpp $(CFLAGS)
	
primarydiff.o: primarydiff.cpp
	$(Compiler)	-c primarydiff.cpp $(CFLAGS)
	
main.o: main.cpp
	$(Compiler)	-c main.cpp $(CFLAGS)

clean:
	rm *.o primarydiff
