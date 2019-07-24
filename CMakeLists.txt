CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(fastreso CXX)

set(CMAKE_MODULE_PATH ${fastreso_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})
# -fopenmp enables the use of multithreading in "#pragma omp parallel for" loops
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -fopenmp -Wall -Wconversion -Wextra -Wmaybe-uninitialized -Wpedantic -std=c++14" )

#------------------------------------------------------------------------
# Turn on GSL
#------------------------------------------------------------------------
find_package(GSL)
include_directories(${GSL_INCLUDE_DIRS})

##------------------------------------------------------------------------
## Add source folder
##------------------------------------------------------------------------

include_directories(./)


# defines target library dist
set(fastreso_source
    TParticle.cpp
    TParticle_THERMINATOR.cpp
    TParticle_AZYHYDRO.cpp
    TFastReso_THERMINATOR.cpp
    TFastReso_AZYHYDRO.cpp
    )

add_library(fastreso ${fastreso_source}) 
target_link_libraries(fastreso ${GSL_LIBRARIES})
add_executable(run_fastreso_therminator.exe run_fastreso_therminator.cpp)
target_link_libraries(run_fastreso_therminator.exe fastreso)
add_executable(run_fastreso_azyhydro.exe run_fastreso_azyhydro.cpp)
target_link_libraries(run_fastreso_azyhydro.exe fastreso)
