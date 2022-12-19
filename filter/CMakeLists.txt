
include_directories(../source)

add_executable(rlsfilter 
   rlsfilter.cpp 
   ../source/RLS.h
)

add_executable(polyrlsfilter 
    polyrlsfilter.cpp 
    ../source/RLS.h 
    ../source/PolyRLS.h
)

