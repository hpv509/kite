PREFIX := build
CC := gcc
CXX := g++
BUILD_TYPE := Release

CPPFLAGS :=   -I/home/henrique/include   -I/usr/local/include/eigen3      -Itools/Src   -ISrc -ISrc/Tools -ISrc/Lattice -ISrc/Hamiltonian -ISrc/Vector -ISrc/Simulation   -DCOMPILE_WAVEPACKET=0

CXXFLAGS := -std=c++17 -O3 -fopenmp
LDFLAGS := -fopenmp -lhdf5_cpp -lhdf5
DESTDIR ?=
