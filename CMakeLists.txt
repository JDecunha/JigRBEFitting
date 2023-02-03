#CMakeLists file for JigRBEFitting
cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project (JigRBE VERSION 1.0.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD "17")


### Find and include ROOT libraries
find_package(ROOT 6.25 CONFIG REQUIRED)
include("${ROOT_USE_FILE}") #include the ROOT CMakelists which defines $ROOT_LIBRARIES

### Find GNU Scientific Library
find_package(GSL REQUIRED)   

### Configure executable
add_executable(JigRBE JigRBEFitting.cc)

#Add our source files
add_subdirectory(src)

#Link my project
target_link_libraries(JigRBE PRIVATE src ${ROOT_LIBRARIES} GSL::gsl)


