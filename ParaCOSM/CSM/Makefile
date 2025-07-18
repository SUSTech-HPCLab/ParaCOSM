
########## Compiler #########


## for icpx with intel openmp
ONEAPI_ROOT = /home/cc/intel/oneapi
TBB_ROOT = $(ONEAPI_ROOT)/tbb/latest

TBB_LIB_DIR = $(TBB_ROOT)/lib/intel64/gcc4.8
TBB_INC = $(TBB_ROOT)/include


## for gcc with gnu openmp
# OPENMP_FLAGS = -fopenmp
# CC = g++ -fdiagnostics-color=always -fopenmp -fPIE -I$(TBB_INC) -L$(TBB_LIB_DIR) -ltbb -Wl,-rpath=$(TBB_LIB_DIR)  # using library openmp

# need to use icpx on intel compiler from oneAPI 
# Do these commands to setup the environment:
# sudo apt install libstdc++-12-dev
# source $(ONEAPI_ROOT)/setvars.sh
CC = icpx -fdiagnostics-color=always -qopenmp -qopenmp-link=static -I/home/cc/haibin2/CSM-Benchmark/ContinuousSubgraphMatching/external/taskflow
	

#############################

########## FLAGS ############

FLAGS = -std=c++17 -g -O3  -Wall -I. -I$(TBB_INC) -L$(TBB_LIB_DIR) -ltbb -Wl,-rpath=$(TBB_LIB_DIR) 
FLAGS += -Wno-unused-command-line-argument
# FLAGS += -Wunused-but-set-variable  -Wunused-variable

## To Optimize for speed:
# FLAGS += -pg # using library gprof
# CC = clang++ -O3 -emit-llvm -S

#############################

LIBS = -pthread

BUILD = build
OBJ = build/obj

GRAPH = graph
MATCH = matching
UTILS = utils
BUILD_TOOLS = build/tools

IEDYN = $(MATCH)/Iedyn
GRAPHFLOW = $(MATCH)/GraphFlow
SYMBI = $(MATCH)/SymBi
SJ_TREE = $(MATCH)/SJTree
TURBOFLUX = $(MATCH)/TurboFlux

PARALLEL_SYMBI = $(MATCH)/Parallel_SymBi
PARALLEL_TURBOFLUX = $(MATCH)/Parallel_TurboFlux
PARALLEL_GRAPHFLOW = $(MATCH)/Parallel_GraphFlow
PARALLEL_NEWSP = $(MATCH)/Parallel_NewSP

all : dir $(BUILD)/csm

dir: $(OBJ)

$(OBJ) :
	mkdir -p $(OBJ)

#################### start ####################

$(BUILD)/csm: $(OBJ)/main.o \
		$(OBJ)/matching.o \
		$(OBJ)/sj_tree.o $(OBJ)/graphflow.o \
		$(OBJ)/turboflux.o $(OBJ)/symbi.o \
		$(OBJ)/iedyn.o \
		$(OBJ)/parallel_symbi.o \
		$(OBJ)/parallel_turboflux.o \
		$(OBJ)/parallel_graphflow.o \
		$(OBJ)/graph.o $(OBJ)/induced_graph.o \
		$(OBJ)/globals.o 
# $(OBJ)/parallel_newsp.o \
		
	$(CC) $(FLAGS) $(OBJ)/main.o \
		$(OBJ)/matching.o \
		$(OBJ)/sj_tree.o $(OBJ)/graphflow.o \
		$(OBJ)/turboflux.o $(OBJ)/symbi.o \
		$(OBJ)/iedyn.o $(OBJ)/parallel_symbi.o  $(OBJ)/parallel_turboflux.o $(OBJ)/parallel_graphflow.o \
		$(OBJ)/graph.o $(OBJ)/induced_graph.o \
		$(OBJ)/globals.o \
		-o $(BUILD)/csm $(LIBS)

# $(OBJ)/parallel_newsp.o \

$(OBJ)/main.o: $(MATCH)/main.cpp \
		$(UTILS)/CLI11.hpp \
		$(UTILS)/globals.h $(UTILS)/types.h \
		$(GRAPH)/graph.h \
		$(SJ_TREE)/sj_tree.h $(GRAPHFLOW)/graphflow.h \
		$(TURBOFLUX)/turboflux.h $(SYMBI)/symbi.h \
		$(IEDYN)/iedyn.h \
		$(PARALLEL_TURBOFLUX)/parallel_turboflux.h \
		$(PARALLEL_GRAPHFLOW)/parallel_graphflow.h \
		$(PARALLEL_SYMBI)/parallel.h \

	$(CC) -c $(FLAGS) $(MATCH)/main.cpp -o $(OBJ)/main.o

# $(PARALLEL_NEWSP)/parallel_newsp.h
#################### matching ####################

$(OBJ)/iedyn.o: $(IEDYN)/iedyn.cpp \
		$(UTILS)/types.h $(UTILS)/globals.h \
		$(GRAPH)/graph.h \
		$(MATCH)/matching.h \
		$(IEDYN)/iedyn.h
	$(CC) -c $(FLAGS) $(IEDYN)/iedyn.cpp \
		-o $(OBJ)/iedyn.o

$(OBJ)/symbi.o: $(SYMBI)/symbi.cpp \
		$(UTILS)/types.h $(UTILS)/globals.h \
		$(GRAPH)/graph.h \
		$(MATCH)/matching.h \
		$(SYMBI)/symbi.h
	$(CC) -c $(FLAGS) $(SYMBI)/symbi.cpp \
		-o $(OBJ)/symbi.o

$(OBJ)/turboflux.o: $(TURBOFLUX)/turboflux.cpp \
		$(UTILS)/types.h $(UTILS)/globals.h \
		$(GRAPH)/graph.h \
		$(MATCH)/matching.h \
		$(TURBOFLUX)/turboflux.h
	$(CC) -c $(FLAGS) $(TURBOFLUX)/turboflux.cpp \
		-o $(OBJ)/turboflux.o

$(OBJ)/graphflow.o: $(GRAPHFLOW)/graphflow.cpp \
		$(UTILS)/types.h $(UTILS)/utils.h \
		$(UTILS)/globals.h \
		$(GRAPH)/graph.h \
		$(MATCH)/matching.h \
		$(GRAPHFLOW)/graphflow.h
	$(CC) -c $(FLAGS) $(GRAPHFLOW)/graphflow.cpp \
		-o $(OBJ)/graphflow.o

$(OBJ)/sj_tree.o: $(SJ_TREE)/sj_tree.cpp \
		$(UTILS)/types.h $(UTILS)/globals.h \
		$(GRAPH)/graph.cpp $(GRAPH)/induced_graph.h \
		$(MATCH)/matching.h \
		$(SJ_TREE)/sj_tree.h
	$(CC) -c $(FLAGS) $(SJ_TREE)/sj_tree.cpp \
		-o $(OBJ)/sj_tree.o

$(OBJ)/parallel_symbi.o: $(PARALLEL_SYMBI)/parrallel_symbi.cpp \
		$(UTILS)/types.h $(UTILS)/globals.h \
		$(GRAPH)/graph.cpp $(GRAPH)/induced_graph.h \
		$(MATCH)/matching.h \
		$(PARALLEL_SYMBI)/parallel.h
	$(CC) -c $(FLAGS)  $(PARALLEL_SYMBI)/parrallel_symbi.cpp \
		-o $(OBJ)/parallel_symbi.o
	# $(CC) -c $(FLAGS)  $(PARALLEL_SYMBI)/parrallel_symbi.cpp \
	# 	-S -o $(OBJ)/parallel_symbi.s

$(OBJ)/parallel_turboflux.o: $(PARALLEL_TURBOFLUX)/parallel_turboflux.cpp \
		$(UTILS)/types.h $(UTILS)/globals.h \
		$(GRAPH)/graph.cpp $(GRAPH)/induced_graph.h \
		$(MATCH)/matching.h \
		$(PARALLEL_TURBOFLUX)/parallel_turboflux.h
	$(CC) -c $(FLAGS)  $(PARALLEL_TURBOFLUX)/parallel_turboflux.cpp \
		-o $(OBJ)/parallel_turboflux.o
	# $(CC) -c $(FLAGS)  $(PARALLEL_TURBOFLUX)/parallel_turboflux.cpp \
	# 	-S -o $(OBJ)/parallel_turboflux.s

$(OBJ)/parallel_graphflow.o: $(PARALLEL_GRAPHFLOW)/parallel_graphflow.cpp \
		$(UTILS)/types.h $(UTILS)/globals.h \
		$(GRAPH)/graph.cpp $(GRAPH)/induced_graph.h \
		$(MATCH)/matching.h \
		$(PARALLEL_GRAPHFLOW)/parallel_graphflow.h
	$(CC) -c $(FLAGS)  $(PARALLEL_GRAPHFLOW)/parallel_graphflow.cpp \
		-o $(OBJ)/parallel_graphflow.o

# $(OBJ)/parallel_newsp.o: $(PARALLEL_NEWSP)/parallel_newsp.cpp \
# 		$(UTILS)/types.h $(UTILS)/globals.h \
# 		$(GRAPH)/graph.cpp $(GRAPH)/induced_graph.h \
# 		$(MATCH)/matching.h \
# 		$(PARALLEL_NEWSP)/parallel_newsp.h
# 	$(CC) -c $(FLAGS)  $(PARALLEL_NEWSP)/parallel_newsp.cpp \
# 		-o $(OBJ)/parallel_newsp.o

$(OBJ)/matching.o: $(MATCH)/matching.cpp \
		$(UTILS)/types.h \
		$(GRAPH)/graph.h \
		$(MATCH)/matching.h
	$(CC) -c $(FLAGS) $(MATCH)/matching.cpp \
		-o $(OBJ)/matching.o

#################### graph ####################

$(OBJ)/graph.o: $(GRAPH)/graph.cpp \
		$(UTILS)/types.h $(UTILS)/utils.h \
		$(GRAPH)/graph.h
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

#################### end ####################

.PHONY: clean

clean: 
	rm -r ${BUILD}/obj
	rm ${BUILD}/csm