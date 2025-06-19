cmake_minimum_required (VERSION 2.8)

SET(CMAKE_BUILD_TYPE "Debug")

project (gem_ext)
find_package (OpenCV REQUIRED)
add_executable (gem_ext 
      image.cpp
      img.c
      str.c 
			hash.c
			hash_func.c
      utils.c 
			col.cpp
      main.cpp
			extract_main.cpp
			prep_main.cpp
			part_main.cpp
			merge_main.cpp
			bin_extr_main.cpp
    )
target_link_libraries (gem_ext ${OpenCV_LIBS} tiff raw z m pthread)
