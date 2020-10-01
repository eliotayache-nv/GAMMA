INITIAL    = test_GRB
TIMESTEP   = rk3
GEOMETRY   = spherical
HYDRO      = rel_sph
SOLVER     = hllc
DIMENSIONS = 2d
IO         = text2d

OS_NAME := $(shell uname -s | tr A-Z a-z)
HOST_NAME := $(shell hostname | cut -c-6)s

ifeq ($(HOME), /home/t/ehra20)
	CXX     = mpicxx
	CXXFLAGS = -Wall -Wextra -std=c++0x -03 -fopenmp 		#run this line on distant
# 	CXXFLAGS = -Wall -Wextra -std=c++0x -g -fopenmp 		#run this line on distant
	LFLAGS = -fopenmp -lgsl -lgslcblas -lm	#run this line on distant
else
ifeq ($(HOME), /home/rwe-ubuntu)
	CXX     = mpicxx
	CXXFLAGS = -Wall -Wextra -std=c++0x -O3 -fopenmp 		#run this line on distant
# 	CXXFLAGS = -Wall -Wextra -std=c++0x -g -fopenmp 		#run this line on distant
	LFLAGS = -L/usr/local/lib -fopenmp -lgsl -lgslcblas -lm	#run this line on distant
	IFLAGS = -I/usr/local/include -I/usr/include -I/usr/lib/openmpi/include
else
ifeq ($(OS_NAME), linux)
	CXX     = /usr/bin/g++
	CXXFLAGS = -Wall -Wextra -std=c++11 -O3 		#run this line on distant
	LFLAGS = -L/usr/local/lib -fopenmp -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5 -lgsl -lgslcblas -lm	#run this line on distant
else
ifeq ($(OS_NAME), darwin)
	CXX     = /usr/local/bin/mpic++
# 	CXXFLAGS = -Wall -Wextra -g -std=c++11 -O0  	#run this line on local
# 	CXXFLAGS = -Wall -Wextra -Qunused-arguments -std=c++11 -g  	#run this line on local
	CXXFLAGS = -Wall -Wextra -std=c++11 -g -fopenmp	#run this line on local
# 	LFLAGS = -L/usr/local/lib -fopenmp -lhdf5 -lgsl -lm 		#run this line on local
	LFLAGS = -L/usr/local/lib -fopenmp -lgsl -lm	#run this line on local
	IFLAGS = -I/usr/local/include -I/usr/include -I/usr/include/hdf5/serial
endif
endif
endif
endif

IFLAGS = -I/usr/local/include -I/usr/include -I/usr/include/hdf5/serial

BIN    = ./bin
SOURCE = ./src
EXEC   = ./bin/GAMMA
TEST   = ./bin/test

INITIALL   =$(SOURCE)/Initial/$(INITIAL)
TIMESTEPP  =$(SOURCE)/Timestep/$(TIMESTEP)
GEOMETRYY  =$(SOURCE)/Geometry/$(GEOMETRY)
HYDROO     =$(SOURCE)/Hydro/$(HYDRO)
SOLVERR    =$(SOURCE)/Solver/$(SOLVER)
DIMENSIONSS=$(SOURCE)/Dimensions/$(DIMENSIONS)
IOO        =$(SOURCE)/IO/$(IO)

SOURCES = $(wildcard $(SOURCE)/*.cpp) 
OBJECTS = $(SOURCES:.cpp=.o) $(INITIALL).o $(GEOMETRYY).o $(HYDROO).o $(SOLVERR).o $(DIMENSIONSS).o $(IOO).o


all : $(BIN)/GAMMA

test : $(BIN)/test

help :
	@grep -E "^# !" Makefile | sed -e 's/# !/ /g'

os: 
	@echo $(OS_NAME)

clean:
	rm -f $(EXEC) $(OBJECTS)

# ------------------------------------------
# Main Executable
# ------------------------------------------
$(EXEC): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $(EXEC) $(IFLAGS) $(LFLAGS)

$(TEST): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $(EXEC) $(IFLAGS) $(LFLAGS)

# ------------------------------------------
# Temorary files (*.o) (IFLAGS should be added here)
# ------------------------------------------
$(SOURCE)/%.o: $(SOURCE)/%.cpp
	$(CXX) $(CXXFLAGS) $(IFLAGS) -c $< -o $@ 

$(INITIALL).o: $(INITIALL).cpp
	$(CXX) $(CXXFLAGS) $(IFLAGS) -c $< -o $@

$(TIMESTEPP).o: $(TIMESTEPP).cpp
	$(CXX) $(CXXFLAGS) $(IFLAGS) -c $< -o $@

$(GEOMETRYY).o: $(GEOMETRYY).cpp
	$(CXX) $(CXXFLAGS) $(IFLAGS) -c $< -o $@

$(HYDROO).o: $(HYDROO).cpp
	$(CXX) $(CXXFLAGS) $(IFLAGS) -c $< -o $@

$(SOLVERR).o: $(SOLVERR).cpp
	$(CXX) $(CXXFLAGS) $(IFLAGS) -c $< -o $@

$(DIMENSIONSS).o: $(DIMENSIONSS).cpp
	$(CXX) $(CXXFLAGS) $(IFLAGS) -c $< -o $@

$(IOO).o: $(IOO).cpp
	$(CXX) $(CXXFLAGS) $(IFLAGS) -c $< -o $@


