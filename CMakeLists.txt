cmake_minimum_required(VERSION 3.2)
project(localization)

include_directories(include)
add_definitions(-Wall -O3 -g -DUNIX -DLINUX -D_REENTRANT)

set(SOURCES
    cmex/lib/gmatrix_linalg.c
    cmex/lib/gmatrix_sparse.c
    cmex/lib/gmatrix_statistics.c
    cmex/lib/gmatrix.c
    cmex/lib/gps.c
    cmex/lib/imu.c
    cmex/lib/kalman.c
    cmex/lib/localization.c
    cmex/lib/magnetometer.c
    cmex/lib/rotation.c
    cmex/lib/sonar.c
)

add_library(localization STATIC ${SOURCES})
