#locate our project source files and make a library
file(GLOB_RECURSE sources ${CMAKE_SOURCE_DIR}/src/*.cc ${CMAKE_SOURCE_DIR}/src/survival_fitting/*.cc)
add_library(src ${sources})

#set includes
target_include_directories(src PUBLIC ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/include/survival_fitting ${CMAKE_SOURCE_DIR}/include/survival_plotting ${CMAKE_SOURCE_DIR}/include/proton_spectra) #set the include files

#link ceres to source
target_link_libraries(src Ceres::ceres)


