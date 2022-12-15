
find_package(Eigen3 3.3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

include_directories(../source)

add_executable(rlsfilter rlsfilter.cpp ../source/RLS.h)

target_link_libraries(rlsfilter Eigen3::Eigen)

