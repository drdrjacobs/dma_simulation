# making executables for 2d and 3d
# using make targets on same files while defining different macros
# make does not track this, so have to do it manually here
rm obj/simulation.cpp.o
echo "Making 2d"
make ga_simulation_2d
rm obj/simulation.cpp.o
echo -e "\nMaking 3d"
make ga_simulation_3d
rm obj/simulation.cpp.o
