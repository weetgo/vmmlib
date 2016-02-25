# Copyright (c) 2011-2016 Daniel Pfeifer <daniel@pfeifer-mail.de>
#                         Stefan.Eilemann@epfl.ch

set(VMMLIB_PUBLIC_HEADERS
  aabb.hpp
  enable_if.hpp
  frustum.hpp
  frustumCuller.hpp
  intersection.hpp
  lowpass_filter.hpp
  math.hpp
  matrix.hpp
  matrix_functors.hpp
  quaternion.hpp
  types.hpp
  util.hpp
  vector.hpp
  visibility.hpp
  vmmlib.hpp
)

common_library(vmmlib)
