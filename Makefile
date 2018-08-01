# Makefile

# Directory names
SRCDIR = src
OBJDIR2d = obj2d
OBJDIR3d = obj3d
BINDIR = bin

# Input Names
CPP_SRC_PATHS = $(wildcard $(SRCDIR)/*.cpp)
CPP_FILES = $(CPP_SRC_PATHS:$(SRCDIR)/%=%)

# -----------------------------------------------------------------------------

# C++ Compiler and Flags
BOOST_PATH = /central/software/boost/1_66_0
BOOST_INC_PATH = $(BOOST_PATH)/include
BOOST_LIB_PATH = $(BOOST_PATH)/lib

EIGEN_INC_PATH = /home/djacobso/eigen-eigen-b3f3d4950030

GPP = g++
FLAGS = -g -Wall -std=c++11 -pthread -O0
INCLUDE = -I$(BOOST_INC_PATH) -I$(EIGEN_INC_PATH)

# -----------------------------------------------------------------------------
# Object files
# -----------------------------------------------------------------------------

# C++ Object Files
CPP_OBJ2d = $(addprefix $(OBJDIR2d)/, $(addsuffix .o, $(CPP_FILES)))
CPP_OBJ3d = $(addprefix $(OBJDIR3d)/, $(addsuffix .o, $(CPP_FILES)))

# -----------------------------------------------------------------------------
# Make rules
# -----------------------------------------------------------------------------

LINK = $(GPP) $(FLAGS) -o $(BINDIR)/$@ $(INCLUDE) $^
COMPILE = $(GPP) $(FLAGS) -DDIMENSIONS=$(DIMENSIONS) -c -o $@ $(INCLUDE) $<

# Top level rules

all: ga_simulation_2d ga_simulation_3d

ga_simulation_2d: DIMENSIONS = 2
ga_simulation_2d: $(CPP_OBJ2d)
	$(LINK)

ga_simulation_3d: DIMENSIONS = 3
ga_simulation_3d: $(CPP_OBJ3d)
	$(LINK)

# Compile C++ Source Files
$(CPP_OBJ2d): $(OBJDIR2d)/%.o : $(SRCDIR)/%
	$(COMPILE)

$(CPP_OBJ3d): $(OBJDIR3d)/%.o : $(SRCDIR)/%
	$(COMPILE)

# Clean everything including temporary Emacs files
clean:
	rm -f $(BINDIR)/* $(OBJDIR2d)/*.o $(OBJDIR3d)/*.o

.PHONY: clean all
