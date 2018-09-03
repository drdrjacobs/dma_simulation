# Makefile

# defines directory names, library paths, and compiler flags, includes, libs 
# shared with tests
include config.mk

# Input Names
CPP_SRC_PATHS = $(wildcard $(SRCDIR)/*.cpp)
CPP_FILES = $(CPP_SRC_PATHS:$(SRCDIR)/%=%)

# -----------------------------------------------------------------------------
# Object files
# -----------------------------------------------------------------------------

# C++ Object Files
CPP_OBJ2d = $(addprefix $(OBJDIR2d)/, $(addsuffix .o, $(CPP_FILES)))
CPP_OBJ3d = $(addprefix $(OBJDIR3d)/, $(addsuffix .o, $(CPP_FILES)))

# -----------------------------------------------------------------------------
# Make rules
# -----------------------------------------------------------------------------

LINK = $(GPP) $(FLAGS) -o $(BINDIR)/$@ $(INCLUDE) $^ $(LIBS)
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
