# Copyright (c) 2011 Daniel Pfeifer <daniel@pfeifer-mail.de>
#               2011-2014 Stefan.Eilemann@epfl.ch

set(VMMLIB_PUBLIC_HEADERS
  aabb.hpp
  enable_if.hpp
  exception.hpp
  frustum.hpp
  frustum_culler.hpp
  intersection.hpp
  lowpass_filter.hpp
  math.hpp
  matrix.hpp
  matrix_functors.hpp
  qr_decomposition.hpp
  quaternion.hpp
  util.hpp
  vector.hpp
  visibility.hpp
  vmmlib.hpp
  vmmlib_config.hpp
)

common_library(vmmlib)
