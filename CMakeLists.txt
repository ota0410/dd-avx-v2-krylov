cmake_minimum_required(VERSION 2.8)
project(DD-AVX CXX)
set(serial "2.0.1")
set(version 2.0.1)

set(novec 0)
set(SSE2 0)
set(AVX 0)
set(AVX2 1)

if(novec EQUAL 1)
	set(SSE2 0)
	set(AVX 0)
	set(AVX2 0)
	set(CMAKE_CXX_FLAGS "-no-vec -O3") # コンパイルオプションの指定
elseif(SSE2 EQUAL 1)
	set(novec 0)
	set(AVX 0)
	set(AVX2 0)
	set(CMAKE_CXX_FLAGS "-sse2 -O3") # コンパイルオプションの指定
elseif(AVX EQUAL 1)
	set(novec 0)
	set(SSE2 0)
	set(AVX2 0)
	set(CMAKE_CXX_FLAGS "-mavx -O3") # コンパイルオプションの指定
elseif(AVX2 EQUAL 1)
	set(novec 0)
	set(SSE2 0)
	set(AVX 1)
	set(CMAKE_CXX_FLAGS "-mavx2 -O3") # コンパイルオプションの指定
endif()

configure_file(
	include/DD-AVX_Config.hpp.in
        include/DD-AVX_Config.hpp
	)


find_package(OpenMP REQUIRED)
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")


include_directories(include/)
add_library(ddavx STATIC
      src/DD-AVX.cpp
      src/scalar/DD-AVX_Scalar.cpp
      src/system/DD-AVX_vecIO.cpp
      src/system/DD-AVX_matIO.cpp
      src/vector/DD-AVX_axpy.cpp
      src/vector/DD-AVX_xpay.cpp
      src/vector/DD-AVX_axpyz.cpp
      src/vector/DD-AVX_scale.cpp
      src/vector/DD-AVX_dot.cpp
      src/vector/DD-AVX_nrm2.cpp
      src/matvec/DD-AVX_matvec.cpp
      src/matvec/DD-AVX_crs.cpp
      )
INSTALL(TARGETS ddavx
      LIBRARY DESTINATION lib
      ARCHIVE DESTINATION lib
      )
INSTALL(FILES include/DD-AVX.hpp
      include/DD-AVX_Config.hpp
      include/DD-AVX_MACRO_AVX.hpp
      include/DD-AVX_MACRO_SSE2.hpp
      include/DD-AVX_MACRO_Scalar.hpp
      include/DD-AVX_Config.hpp
      DESTINATION include/
      )

add_executable(sample/main sample/main.cpp)
target_link_libraries(sample/main ddavx)
#add_executable(sample/main2 sample/main2.cpp)
#target_link_libraries(sample/main2 ddavx)
