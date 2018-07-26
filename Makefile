# Makefile

# Directory names
SRCDIR = src
OBJDIR = obj
BINDIR = bin

# Input Names
CPP_SRC_PATHS = $(wildcard $(SRCDIR)/*.cpp)
CPP_FILES = $(CPP_SRC_PATHS:$(SRCDIR)/%=%)

# -----------------------------------------------------------------------------

# C++ Compiler and Flags
BOOST_PATH = /central/software/boost/1_66_0
BOOST_INC_PATH = $(BOOST_PATH)/include
BOOST_LIB_PATH = $(BOOST_PATH)/lib

GPP = g++
FLAGS = -Wall -D_REENTRANT -std=c++11 -pthread -O3
INCLUDE = -I$(CUDA_INC_PATH) -I$(BOOST_INC_PATH)

# -----------------------------------------------------------------------------
# Object files
# -----------------------------------------------------------------------------

# C++ Object Files
CPP_OBJ = $(addprefix $(OBJDIR)/, $(addsuffix .o, $(CPP_FILES)))

# -----------------------------------------------------------------------------
# Make rules
# -----------------------------------------------------------------------------

# Top level rules
all: ga_simulation

ga_simulation: $(CPP_OBJ)
	$(GPP) $(FLAGS) -o $(BINDIR)/$@ $(INCLUDE) $^

# Compile C++ Source Files
$(CPP_OBJ): $(OBJDIR)/%.o : $(SRCDIR)/%
	$(GPP) $(FLAGS) -c -o $@ $(INCLUDE) $<

# Clean everything including temporary Emacs files
clean:
	rm -f $(BINDIR)/* $(OBJDIR)/*.o $(SRCDIR)/*~ *~

.PHONY: clean all
