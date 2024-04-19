cmake_minimum_required(VERSION 3.23)
project(proyecto)

set(CMAKE_CXX_STANDARD 14)

add_executable(proyecto main.cpp EKF.cpp EKF_test.cpp include/Matrix.h include/R_z.h src/Matrix.cpp src/R_z.cpp src/R_x.cpp src/R_y.cpp include/R_x.h include/R_y.h src/AccelPointMass.cpp include/AccelPointMass.h)
