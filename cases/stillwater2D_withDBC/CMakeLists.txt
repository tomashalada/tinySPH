add_definitions(-D_USE_MATH_DEFINES)

include_directories("../../libs/include")
include_directories("../../src_stillWater")
include_directories("../../src_stillWater/general")
include_directories("../../src_stillWater/WCSPH")
include_directories("measuretools")

add_executable(stillwater2D_DBC
    main.cpp
)

find_package(OpenMP REQUIRED)
target_link_libraries(stillwater2D_DBC CompactNSearch OpenMP::OpenMP_CXX)
