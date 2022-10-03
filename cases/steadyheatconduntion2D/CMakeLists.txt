add_definitions(-D_USE_MATH_DEFINES)

include_directories("../../libs/include")
include_directories("../../src")
include_directories("../../src/general")
include_directories("../../src/heatSIMPLE")

add_executable(heatequation
    main.cpp
)

find_package(OpenMP REQUIRED)
target_link_libraries(heatequation CompactNSearch OpenMP::OpenMP_CXX)
