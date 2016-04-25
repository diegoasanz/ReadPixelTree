cmake_minimum_required(VERSION 3.3)
project(ReadPixelTree)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES src/ReadPixelTree.cpp)

add_executable(ReadPixelTree ${SOURCE_FILES})

set(CMAKE_INCLUDE_PATH `root-config --incdir`)
set(COMANDO1 root-config)
EXEC_PROGRAM(${COMANDO1} ARGS --incdir OUTPUT_VARIABLE ROOT_INC_DIR)
include_directories(${ROOT_INC_DIR})

set(CMAKE_LIBRARY_PATH `root-config --libs --glibs --auxlibs`)

set(CMAKE_CXX_FLAGS "-m64 -pthread -stdlib=libc++ -std=c++1y -m64 -pthread -stdlib=libc++ -std=c++1y -m64 -I/usr/local/include -L/usr/local/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lpthread -Wl,-rpath,/usr/local/lib -stdlib=libc++ -lm -ldl -L/usr/local/lib -lGui -lpthread -Wl,-rpath,/usr/local/lib -stdlib=libc++ -lm -ldl")

set(CMAKE_SHARED_LINKER_FLAGS "-m64 -pthread -stdlib=libc++ -std=c++1y -m64 -pthread -stdlib=libc++ -std=c++1y -m64 -I/usr/local/include -L/usr/local/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lpthread -Wl,-rpath,/usr/local/lib -stdlib=libc++ -lm -ldl -L/usr/local/lib -lGui -lpthread -Wl,-rpath,/usr/local/lib -stdlib=libc++ -lm -ldl")
