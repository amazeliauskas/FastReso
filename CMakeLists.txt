CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(fastreso CXX)

set(CMAKE_MODULE_PATH ${fastreso_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -fopenmp -Wall -std=c++17")

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
    TFastReso_THERMINATOR.cpp
    TFastReso_AZYHYDRO.cpp
    TParticle_THERMINATOR.cpp
    TParticle_AZYHYDRO.cpp
    TKernel.cpp
    )

add_library(fastreso ${fastreso_source}) 
target_link_libraries(fastreso ${GSL_LIBRARIES})

option(FASTRESO_ROOT "Link with ROOT Librariers" OFF)
if (FASTRESO_ROOT)
    message("FASTRESO_ROOT is on") 
   add_definitions("-DEKT_ROOT") 
   find_package(ROOT)
   include_directories(AFTER ${ROOT_INCLUDE_DIR})
   message("FASTRESO_ROOT links with ${ROOT_LIBRARIES}")
endif()



# compile executables
add_executable(run_fastreso_therminator.exe run_fastreso_therminator.cpp)
add_executable(run_fastreso_azyhydro.exe run_fastreso_azyhydro.cpp)
add_executable(run_kernels.exe run_kernels.cpp)

# link libraries
target_link_libraries(run_fastreso_therminator.exe fastreso)
target_link_libraries(run_fastreso_azyhydro.exe fastreso)
target_link_libraries(run_kernels.exe fastreso)


if (FASTRESO_ROOT)
add_executable(run_kernels_root.exe run_kernels_root.cpp)
target_link_libraries(run_kernels_root.exe fastreso)
target_link_libraries(run_kernels_root.exe ${ROOT_LIBRARIES})
endif()




