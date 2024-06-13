cmake_minimum_required(VERSION 3.0)
project(AndorsCMOS CXX)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenCV REQUIRED)
include_directories(include/)
include_directories(${OpenCV_INCLUDE_DIRS})
link_directories(lib/)

add_executable(main.cpp DataToImage/DataToImage.cpp MCR-ALS/linear_algebra.cpp MCR-ALS/MCR-ALS.cpp MCR-ALS/Regressors.cpp MCR-ALS/utils.cpp)
target_link_libraries(acqcam atcore ${OpenCV_LIBRARIES})


