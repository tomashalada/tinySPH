add_definitions(-D_USE_MATH_DEFINES)

include_directories("../../libs/include")
include_directories("../../src")
include_directories("../../src/general")
include_directories("../../src/WCSPH_MDBC")
include_directories("measuretools")

add_executable(stillwater2D_mDBC_periodic
    main.cpp
)

find_package(OpenMP REQUIRED)
target_link_libraries(stillwater2D_mDBC_periodic CompactNSearch OpenMP::OpenMP_CXX)
