cmake_minimum_required(VERSION 3.2)
project(localization)

include_directories(include)
add_definitions(-Wall -O3 -g -DUNIX -DLINUX -D_REENTRANT)

set(SOURCES
    src/gmatrix_linalg.c
    src/gmatrix_sparse.c
    src/gmatrix_statistics.c
    src/gmatrix.c
    src/gps.c
    src/imu.c
    src/kalman.c
    src/localization.c
    src/magnetometer.c
    src/rotation.c
    src/sonar.c
)

add_library(localization STATIC ${SOURCES})
