cmake_minimum_required(VERSION 3.17)
project(bone-mechanical-characterization Fortran)


enable_language(Fortran)

find_package(LAPACK REQUIRED)

add_executable(main solid_mechanics.f90)
target_link_libraries(main ${LAPACK_LIBRARIES})


