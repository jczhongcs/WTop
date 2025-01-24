cmake_minimum_required(VERSION 3.19)
project(testOutXML)
#set(CMAKE_CXX_COMPILER /opt/homebrew/Cellar/llvm/17.0.6/bin/clang++)


set(CMAKE_CXX_STANDARD 14)
include_directories(${PROJECT_SOURCE_DIR}/include)
add_library(tinyxml
        ${PROJECT_SOURCE_DIR}/tinyxml/tinyxml.cpp
        ${PROJECT_SOURCE_DIR}/tinyxml/tinystr.cpp
        ${PROJECT_SOURCE_DIR}/tinyxml/tinyxmlparser.cpp
        ${PROJECT_SOURCE_DIR}/tinyxml/tinyxmlerror.cpp
        )
add_executable(testOutXML main.cpp)
target_link_libraries(testOutXML tinyxml)

