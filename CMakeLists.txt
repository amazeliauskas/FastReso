CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(fastreso CXX)

set(CMAKE_MODULE_PATH ${fastreso_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -Wall -Wconversion -Wextra -Wmaybe-uninitialized -Wpedantic -std=c++14" )

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
    TFastReso.cpp
    )

add_library(fastreso ${fastreso_source}) 
target_link_libraries(fastreso ${GSL_LIBRARIES})
add_executable(run_fastreso.exe run_fastreso.cpp)
#add_executable(run_temperature.exe run_temperature.cpp)
target_link_libraries(run_fastreso.exe fastreso)
