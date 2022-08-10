add_definitions(-D_USE_MATH_DEFINES)

include_directories("../../libs/include")
include_directories("../../src")
include_directories("../../src/general")
include_directories("../../src/WCSPH_BI")

add_executable(dambreak2D_BI
    main.cpp
)

find_package(OpenMP REQUIRED)
target_link_libraries(dambreak2D_BI CompactNSearch OpenMP::OpenMP_CXX)
