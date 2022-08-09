ROOTCFLAGS     = $(shell root-config --cflags)
ROOTLIBS       = $(shell root-config --libs)
ROOTGLIBS      = $(shell root-config --glibs) 

INCLUDES       = -I./include 

CXX            = g++
CXXFLAGS       = -fPIC -fno-var-tracking -Wno-deprecated -D_GNU_SOURCE -O2  $(INCLUDES) 
CXXFLAGS      += $(ROOTCFLAGS)

LD             = g++
LDFLAGS        = 

SOFLAGS        = -O -shared  -fPIC #-flat_namespace 
LIBS           = $(ROOTLIBS) 

GLIBS         = $(ROOTGLIBS) -lMinuit -lTreePlayer -lTMVA

SRCS = src/main.cc src/Analyzer_DisplacedMuon.cc
SRCS_CHI2RZ = src/main_Chi2rz.cc src/Analyzer_DisplacedMuon_Chi2rz.cc
SRCS_CHI2BEND = src/main_Chi2bend.cc src/Analyzer_DisplacedMuon_Chi2bend.cc
SRCS_Z0 = src/main_Z0.cc src/Analyzer_DisplacedMuon_Z0.cc
TESTSRCS = src/parallelTest.cc

OBJS =  $(patsubst %.C,%.o,$(SRCS:.cc=.o))
OBJS_CHI2RZ =  $(patsubst %.C,%.o,$(SRCS_CHI2RZ:.cc=.o))
OBJS_CHI2BEND =  $(patsubst %.C,%.o,$(SRCS_CHI2BEND:.cc=.o))
OBJS_Z0 =  $(patsubst %.C,%.o,$(SRCS_Z0:.cc=.o))
TESTOBJS = $(patsubst %.C,%.o,$(TESTSRCS:.cc=.o))

LIB=lib/main.so
LIB_CHI2RZ=lib/main_Chi2rz.so
LIB_CHI2BEND=lib/main_Chi2bend.so
LIB_Z0=lib/main_Z0.so
TEST=lib/parallelTest.so

.SUFFIXES: .cc,.C,.hh,.h

# Rules ====================================
all: $(LIB) RunAll RunTest RunAll_Chi2rz RunAll_Chi2bend RunAll_Z0

lib : $(LIB)

$(LIB): $(OBJS)
	@echo "Creating library $(LIB)"
	mkdir -p lib
	$(LD) $(LDFLAGS) $(GLIBS) $(SOFLAGS) $(OBJS) -o $(LIB) -fopenmp
	@echo "$(LIB) successfully compiled!"

$(LIB_CHI2RZ): $(OBJS_CHI2RZ)
	@echo "Creating library $(LIB_CHI2RZ)"
	mkdir -p lib
	$(LD) $(LDFLAGS) $(GLIBS) $(SOFLAGS) $(OBJS_CHI2RZ) -o $(LIB_CHI2RZ) -fopenmp
	@echo "$(LIB_CHI2RZ) successfully compiled!"

$(LIB_CHI2BEND): $(OBJS_CHI2BEND)
	@echo "Creating library $(LIB_CHI2BEND)"
	mkdir -p lib
	$(LD) $(LDFLAGS) $(GLIBS) $(SOFLAGS) $(OBJS_CHI2BEND) -o $(LIB_CHI2BEND) -fopenmp
	@echo "$(LIB_CHI2BEND) successfully compiled!"

$(LIB_Z0): $(OBJS_Z0)
	@echo "Creating library $(LIB_Z0)"
	mkdir -p lib
	$(LD) $(LDFLAGS) $(GLIBS) $(SOFLAGS) $(OBJS_Z0) -o $(LIB_Z0) -fopenmp
	@echo "$(LIB_Z0) successfully compiled!"

$(TEST): $(TESTOBJS)
	@echo "Creating library $(TEST)"
	mkdir -p lib
	$(LD) $(LDFLAGS) $(GLIBS) $(SOFLAGS) $(TESTOBJS) -o $(TEST) -fopenmp
	@echo "$(TEST) successfully compiled!"

RunAll : src/main.cc $(LIB)	
	mkdir -p bin
	$(CXX) $(CXXFLAGS) -ldl $(LDFLAGS) -o $@ $^ $(GLIBS) -fopenmp -lpthread

RunAll_Chi2rz : src/main_Chi2rz.cc $(LIB_CHI2RZ)	
	mkdir -p bin
	$(CXX) $(CXXFLAGS) -ldl $(LDFLAGS) -o $@ $^ $(GLIBS) -fopenmp -lpthread

RunAll_Chi2bend : src/main_Chi2bend.cc $(LIB_CHI2BEND)	
	mkdir -p bin
	$(CXX) $(CXXFLAGS) -ldl $(LDFLAGS) -o $@ $^ $(GLIBS) -fopenmp -lpthread

RunAll_Z0 : src/main_Z0.cc $(LIB_Z0)	
	mkdir -p bin
	$(CXX) $(CXXFLAGS) -ldl $(LDFLAGS) -o $@ $^ $(GLIBS) -fopenmp -lpthread

RunTest: src/parallelTest.cc $(TEST)
	mkdir -p bin  
	$(CXX) $(CXXFLAGS) -ldl $(LDFLAGS) -o $@ $^ $(GLIBS) -fopenmp -lpthread

clean:
	$(RM) $(OBJS)	
	$(RM) $(OBJS_CHI2RZ)	
	$(RM) $(OBJS_CHI2BEND)	
	$(RM) $(OBJS_Z0)	
	$(RM) $(TESTOBJS)
	$(RM) $(LIB)
	$(RM) $(LIB_CHI2RZ)
	$(RM) $(LIB_CHI2BEND)
	$(RM) $(LIB_Z0)
	$(RM) $(TEST)
	$(RM) RunAll
	$(RM) RunAll_Chi2rz
	$(RM) RunAll_Chi2bend
	$(RM) RunAll_Z0
	$(RM) RunTest

purge:
	$(RM) $(OBJS)
	$(RM) $(OBJS_CHI2RZ)
	$(RM) $(OBJS_CHI2BEND)
	$(RM) $(OBJS_Z0)
	$(RM) $(TESTOBJS)

deps: $(SRCS) $(TESTSRCS) $(SRCS_CHI2RZ) $(SRCS_CHI2BEND) $(SRCS_Z0)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
