cmake_minimum_required(VERSION 2.6)
project(2drtheta)
enable_language(Fortran)

set (CMAKE_VERBOSE_MAKEFILE OFF)
set (CMAKE_Fortran_FLAGS_RELEASE "-fopenmp -O3")
set (CMAKE_Fortran_FLAGS_DEBUG   "-fopenmp -O2")
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS_DEBUG} -fimplicit-none -Wall -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -Wimplicit-interface -Wunused-parameter -fwhole-file -fcheck=all -std=f2008 -pedantic -fbacktrace")

add_executable(2drtheta.lx main.f90 mydefs.f90 inits.f90 calcs.f90)

add_library(2drthetalib mydefs.f90 inits.f90 calcs.f90)
