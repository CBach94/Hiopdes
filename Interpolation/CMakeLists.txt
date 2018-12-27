set(Interpolation_Source
  Interpolation.h
  Interpolation.cpp
  LagrangeBasis.h
  LagrangeBasis.cpp)

add_library(Interpolation ${Interpolation_Source})
target_link_libraries(Interpolation ${ARMADILLO_LIBRARIES})
