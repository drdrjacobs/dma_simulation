# Makefile for tests

# defines library paths, and compiler flags, includes, libs shared with 
# simulation exectables
include ../config.mk

# Directory names
SRCDIR = src
OBJDIR2d_TEST = obj2d_test
OBJDIR3d_TEST = obj3d_test
BINDIR = bin

# Input Names
CPP_SRC_PATHS = $(wildcard $(SRCDIR)/*.cpp)
CPP_FILES = $(CPP_SRC_PATHS:$(SRCDIR)/%=%)

# -----------------------------------------------------------------------------
# Object files
# -----------------------------------------------------------------------------

# C++ Object Files
CPP_OBJ2d = $(addprefix $(OBJDIR2d_TEST)/, $(addsuffix .o, $(CPP_FILES)))
CPP_OBJ3d = $(addprefix $(OBJDIR3d_TEST)/, $(addsuffix .o, $(CPP_FILES)))

MAIN_OBJ = main.cpp.o

# -----------------------------------------------------------------------------
# Make rules
# -----------------------------------------------------------------------------

LINK = $(GPP) $(FLAGS) -o $(BINDIR)/$@ $(INCLUDE) $(SIMULATION_OBJ) $^ $(LIBS)
COMPILE = $(GPP) $(FLAGS) -DDIMENSIONS=$(DIMENSIONS) -c -o $@ $(INCLUDE) $<

# Top level rules

all: simulation_executables tests_2d tests_3d

# Make sure all simulation code is compiled
simulation_executables:
	$(MAKE) -C ../

# Include simulation code in linking of tests
tests_2d: DIMENSIONS = 2
tests_2d: SIMULATION_OBJ = $(filter-out ../$(OBJDIR2d)/$(MAIN_OBJ), \
	$(wildcard ../$(OBJDIR2d)/*.cpp.o))
tests_2d: $(CPP_OBJ2d)
	$(LINK)

tests_3d: DIMENSIONS = 3
tests_3d:  SIMULATION_OBJ = $(filter-out ../$(OBJDIR3d)/$(MAIN_OBJ), \
	$(wildcard ../$(OBJDIR3d)/*.cpp.o))	
tests_3d: $(CPP_OBJ3d)
	$(LINK)

# Compile C++ Source Files
$(CPP_OBJ2d): $(OBJDIR2d_TEST)/%.o : $(SRCDIR)/%
	$(COMPILE)

$(CPP_OBJ3d): $(OBJDIR3d_TEST)/%.o : $(SRCDIR)/%
	$(COMPILE)

# Clean everything including temporary Emacs files
clean:
	rm -f $(BINDIR)/* $(OBJDIR2d_TEST)/*.o $(OBJDIR3d_TEST)/*.o

.PHONY: clean all simulation_executables
