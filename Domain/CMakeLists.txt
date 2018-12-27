set(Element_Source
  Element.h
  Element.cpp
  Projection.h
  Projection.cpp
  Quadratures.h)

add_library(Element ${Element_Source})
target_link_libraries(Element armadillo)
