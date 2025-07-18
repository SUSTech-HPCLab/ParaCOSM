########## Compiler #########
# export LD_LIBRARY_PATH=/home/cc/intel/oneapi/tbb/latest/lib:$LD_LIBRARY_PATH
## for gcc with gnu openmp
# OPENMP_FLAGS = -fopenmp
# CC = icpx -fdiagnostics-color=always -Wno-sign-compare -qopenmp -qopenmp-link=static # using library openmp
###########

## for icpx with intel openmp
ONEAPI_ROOT = /home/cc/intel/oneapi
TBB_ROOT = $(ONEAPI_ROOT)/tbb/latest

TBB_LIB_DIR = $(TBB_ROOT)/lib/intel64/gcc4.8
TBB_INC = $(TBB_ROOT)/include

# need to use icpx on intel compiler from oneAPI 
# Do these commands to setup the environment:
# sudo apt install libstdc++-12-dev
# source $(ONEAPI_ROOT)/setvars.sh
CC = icpx -fdiagnostics-color=always -qopenmp -qopenmp-link=static \
	-I$(TBB_INC) -L$(TBB_LIB_DIR) -ltbb -Wl,-rpath=$(TBB_LIB_DIR) 

#############################

FLAGS = -std=c++17 -g -O3 -Wall -I. 
# FLAGS = -std=c++17 -g -pg -fno-inline -Wall -I.
LIBS = -pthread

BUILD = build
OBJ = build/obj

GRAPH = graph
MATCH = matching
UTILS = utils
MODEL = model
DMS = DecisionMakingSystem
BUILD_TOOLS = build/tools
CONFIG = configuration/config.h

all : dir $(BUILD)/csm

dir: $(OBJ)

$(OBJ) :
	mkdir -p $(OBJ)

#################### start ####################

$(BUILD)/csm: $(OBJ)/main.o \
		$(OBJ)/matching.o \
		$(OBJ)/CSMPP.o \
		$(OBJ)/graph.o $(OBJ)/induced_graph.o \
		$(OBJ)/globals.o \
		$(OBJ)/DecisionMakingSystem.o
	$(CC) $(FLAGS) $(OBJ)/main.o \
		$(OBJ)/matching.o \
		$(OBJ)/CSMPP.o \
		$(OBJ)/graph.o $(OBJ)/induced_graph.o \
		$(OBJ)/globals.o \
		$(OBJ)/DecisionMakingSystem.o \
		-o $(BUILD)/csm $(LIBS)

$(OBJ)/main.o: $(MATCH)/main.cpp \
		$(UTILS)/CLI11.hpp \
		$(UTILS)/globals.h $(UTILS)/types.h \
		$(GRAPH)/graph.h \
		$(MATCH)/CSMPP.h
	$(CC) -c $(FLAGS) $(MATCH)/main.cpp -o $(OBJ)/main.o

#################### matching ####################
$(OBJ)/CSMPP.o: $(MATCH)/CSMPP.cpp \
		$(UTILS)/types.h $(UTILS)/utils.h \
		$(UTILS)/globals.h \
		$(GRAPH)/graph.h \
		$(MATCH)/matching.h \
		$(MATCH)/CSMPP.h \
		$(CONFIG)
	$(CC) -c $(FLAGS) $(MATCH)/CSMPP.cpp \
		-o $(OBJ)/CSMPP.o

$(OBJ)/matching.o: $(MATCH)/matching.cpp \
		$(UTILS)/types.h \
		$(GRAPH)/graph.h \
		$(MATCH)/matching.h
	$(CC) -c $(FLAGS) $(MATCH)/matching.cpp \
		-o $(OBJ)/matching.o

#################### graph ####################

$(OBJ)/graph.o: $(GRAPH)/graph.cpp \
		$(UTILS)/types.h $(UTILS)/utils.h \
		$(GRAPH)/graph.h $(DMS)/DecisionMakingSystem.h \
		$(CONFIG)
	$(CC) -c $(FLAGS) $(GRAPH)/graph.cpp \
		-o $(OBJ)/graph.o

$(OBJ)/induced_graph.o: $(GRAPH)/induced_graph.cpp \
		$(UTILS)/types.h \
		$(GRAPH)/induced_graph.h $(GRAPH)/graph.h
	$(CC) -c $(FLAGS) $(GRAPH)/induced_graph.cpp \
		-o $(OBJ)/induced_graph.o

#################### utils ####################

$(OBJ)/globals.o: $(UTILS)/globals.cpp $(UTILS)/globals.h
	$(CC) -c $(FLAGS) $(UTILS)/globals.cpp \
		-o $(OBJ)/globals.o

#################### DMS ####################

$(OBJ)/DecisionMakingSystem.o: $(DMS)/DecisionMakingSystem.cpp $(DMS)/DecisionMakingSystem.h $(UTILS)/types.h $(UTILS)/globals.h 
	$(CC) -c $(FLAGS) $(DMS)/DecisionMakingSystem.cpp \
		-o $(OBJ)/DecisionMakingSystem.o

#################### end ####################

.PHONY: clean

clean: 
	rm -r ${BUILD}
