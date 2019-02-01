cmake_minimum_required(VERSION 3.8)

project(CoRBS)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_ROOT}/Modules/)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/)

find_package(LIBIGL)
if(NOT LIBIGL_FOUND)
  set(FATAL_ERROR "Libigl not found")
endif()

if(UNIX AND NOT APPLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -g -fopenmp")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -g")
endif()

include_directories(${CMAKE_SOURCE_DIR}/include/)

file(GLOB CORBS_SRC src/*.cpp)

# set_target_properties(CoRBS PROPERTIES VERSION ${PROJECT_VERSION})
# set_target_properties(CoRBS PROPERTIES PUBLIC_HEADER include/)

add_executable(${PROJECT_NAME} ${CORBS_SRC})
target_link_libraries(${PROJECT_NAME} igl::core)
