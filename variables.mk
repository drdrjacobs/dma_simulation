# Makefile variables for simulation exectables and test executables

# Directory names
SRCDIR = src
OBJDIR2d = obj2d
OBJDIR3d = obj3d
BINDIR = bin

# Libraries
BOOST_PATH = /home/djacobso/boost_1_66_0
BOOST_INC_PATH = $(BOOST_PATH)
BOOST_LIB_PATH = $(BOOST_PATH)/stage/lib

EIGEN_INC_PATH = /home/djacobso/eigen-eigen-b3f3d4950030

# Compiler and flags
GPP = g++
FLAGS = -g -Wall -std=c++17 -pthread -O0
INCLUDE = -I$(BOOST_INC_PATH) -I$(EIGEN_INC_PATH)
LIBS= -L$(BOOST_LIB_PATH) -lboost_serialization -lboost_filesystem \
	-lboost_system
