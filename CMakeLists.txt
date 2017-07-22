cmake_minimum_required(VERSION 3.2)
project(m3vcftools VERSION 1.0.0)

execute_process(COMMAND date OUTPUT_VARIABLE DATE OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND whoami OUTPUT_VARIABLE USER OUTPUT_STRIP_TRAILING_WHITESPACE)

add_definitions(-DVERSION="${PROJECT_VERSION}" -DUSER="${USER}" -DDATE="${DATE}")

find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(ZLIB REQUIRED)
find_library(STATGEN_LIBRARY StatGen)

add_executable(m3vcftools
  src/Concat.cpp
  src/Compress.cpp
  src/Convert.cpp
  src/m3vcfBlock.cpp
  src/m3vcfFile.cpp
  src/m3vcfFileWriter.cpp
  src/m3vcfFileReader.cpp
  src/m3vcfHeader.cpp
  src/m3vcfRecord.cpp
  src/Main.cpp
  src/Tokenize.cpp
  src/Unique.cpp)

target_link_libraries(m3vcftools ${STATGEN_LIBRARY} ${ZLIB_LIBRARIES})

install(TARGETS m3vcftools RUNTIME DESTINATION bin)
