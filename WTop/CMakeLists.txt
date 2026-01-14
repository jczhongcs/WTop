cmake_minimum_required(VERSION 3.0)
set(CMAKE_BUILD_TYPE Debug)

#
#set(CMAKE_CXX_COMPILER /opt/homebrew/Cellar/llvm/19.1.0/bin/clang++)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -g -fno-omit-frame-pointer -fsanitize-recover=address")

#set(CMAKE_CXX_COMPILER /opt/homebrew/Cellar/llvm/19.1.0/bin/clang++)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address,leak -g -fno-omit-frame-pointer")



project(WTop)

set(CMAKE_CXX_STANDARD 17)




include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()


SET(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

#cwdtw lib
ADD_SUBDIRECTORY(src/wavelib)

include_directories(${PROJECT_SOURCE_DIR}/src/tinyxml)

AUX_SOURCE_DIRECTORY(src/merge DIR_MERGE)
AUX_SOURCE_DIRECTORY(src/mylib DIR_MYLIB)
AUX_SOURCE_DIRECTORY(src/util DIR_UTIL)
AUX_SOURCE_DIRECTORY(src/console DIR_PROCESS)
AUX_SOURCE_DIRECTORY(src/config DIR_CONFIG)
AUX_SOURCE_DIRECTORY(src/ptm DIR_PTM)
AUX_SOURCE_DIRECTORY(src/ms DIR_MS)
AUX_SOURCE_DIRECTORY(src/thread DIR_THREAD)
AUX_SOURCE_DIRECTORY(src/protein DIR_PROTEIN)
AUX_SOURCE_DIRECTORY(src/filter DIR_FILTER)
AUX_SOURCE_DIRECTORY(src/prsm DIR_PRSM)
AUX_SOURCE_DIRECTORY(src/tinyxml DIR_TINYXML)


MESSAGE(STATUS "building protein")
add_library (protein STATIC ${DIR_PROTEIN} )
TARGET_LINK_LIBRARIES(protein)
set_target_properties(protein PROPERTIES OUTPUT_NAME "protein")
set_target_properties(protein PROPERTIES CLEAN_DIRECT_OUTPUT 1)


MESSAGE(STATUS "building filter")
add_library (filter STATIC ${DIR_FILTER} )
TARGET_LINK_LIBRARIES(filter ms protein)
set_target_properties(filter PROPERTIES OUTPUT_NAME "filter")
set_target_properties(filter PROPERTIES CLEAN_DIRECT_OUTPUT 1)



MESSAGE(STATUS "building thread")
add_library (thread STATIC ${DIR_MS} ${DIR_CONFIG} )
TARGET_LINK_LIBRARIES(thread)
set_target_properties(thread PROPERTIES OUTPUT_NAME "thread")
set_target_properties(thread PROPERTIES CLEAN_DIRECT_OUTPUT 1)

MESSAGE(STATUS "building ms")
add_library (ms STATIC ${DIR_MS})
TARGET_LINK_LIBRARIES(ms mylib util protein filter)
set_target_properties(ms PROPERTIES OUTPUT_NAME "ms")
set_target_properties(ms PROPERTIES CLEAN_DIRECT_OUTPUT 1)

MESSAGE(STATUS "building config")
add_library (config STATIC ${DIR_PTM})
TARGET_LINK_LIBRARIES(config)
set_target_properties(config PROPERTIES OUTPUT_NAME "config")
set_target_properties(config PROPERTIES CLEAN_DIRECT_OUTPUT 1)


MESSAGE(STATUS "building ptm")
add_library (ptm STATIC ${DIR_MERGE})
set_target_properties(ptm PROPERTIES OUTPUT_NAME "ptm")
set_target_properties(ptm PROPERTIES CLEAN_DIRECT_OUTPUT 1)

MESSAGE(STATUS "building mylib")
add_library (mylib STATIC ${DIR_MYLIB} ${DIR_MS} ${DIR_PROTEIN})
TARGET_LINK_LIBRARIES(mylib wavelib)
set_target_properties(mylib PROPERTIES OUTPUT_NAME "mylib")
set_target_properties(mylib PROPERTIES CLEAN_DIRECT_OUTPUT 1)


MESSAGE(STATUS "building util")
add_library (util STATIC ${DIR_UTIL} src/util/xmlUtil.cpp src/util/xmlUtil.h)
set_target_properties(util PROPERTIES OUTPUT_NAME "util")
set_target_properties(util PROPERTIES CLEAN_DIRECT_OUTPUT 1)


MESSAGE(STATUS "building tinyxml")
add_library(tinyxml STATIC ${DIR_TINYXML})
set_target_properties(tinyxml PROPERTIES OUTPUT_NAME "tinyxml")
set_target_properties(tinyxml PROPERTIES CLEAN_DIRECT_OUTPUT 1)


MESSAGE(STATUS "building prsm")
add_library (prsm STATIC ${DIR_PRSM} ${DIR_MERGE} ${DIR_MYLIB} src/prsm/prsm_store_ptm.h src/prsm/candi_truncation.cpp src/prsm/candi_truncation.h src/prsm/candidate_prsm.cpp src/prsm/candidate_prsm.hpp)
TARGET_LINK_LIBRARIES(prsm ms protein util tinyxml)
set_target_properties(prsm PROPERTIES OUTPUT_NAME "prsm")
set_target_properties(prsm PROPERTIES CLEAN_DIRECT_OUTPUT 1)

MESSAGE(STATUS "building merge")
add_library (merge STATIC ${DIR_MERGE} ${DIR_MS} ${DIR_FILTER})
TARGET_LINK_LIBRARIES(merge wavelib mylib util prsm)
set_target_properties(merge PROPERTIES OUTPUT_NAME "merge")
set_target_properties(merge PROPERTIES CLEAN_DIRECT_OUTPUT 1)


ADD_EXECUTABLE(wtop ${DIR_PROCESS} ${DIR_CONFIG} ${DIR_THREAD} src/util/xmlUtil.cpp)
TARGET_LINK_LIBRARIES(wtop util mylib merge config -lpthread)

ADD_EXECUTABLE(xmlUtil_exec src/util/xmlUtil.cpp)
TARGET_LINK_LIBRARIES(xmlUtil_exec util tinyxml mylib)

#add_custom_target(run_wtop
#        COMMAND ASAN_OPTIONS=detect_leaks=1:verbosity=1 ${EXECUTABLE_OUTPUT_PATH}/wtop -i r -p r -x bin/config
#        DEPENDS wtop
#        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/bin
#        COMMENT "Running wtop with ASAN leak detection"
#        )