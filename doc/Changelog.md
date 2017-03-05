# Changelog {#Changelog}

# git master

# Release 1.12 (09-12-2016)

* [62](https://github.com/Eyescale/vmmlib/pull/62)
  Removed a check in Matrix::computeInverse() that was discarding
  valid matrices

# Release 1.11 (30-06-2016)

* [57](https://github.com/Eyescale/vmmlib/pull/57):
  Fix handling of non-invertible 2x2, 3x3 and 4x4 matrices
* [50](https://github.com/Eyescale/vmmlib/pull/50):
  Improved matrix, quaternion and lowpass filter API
* [49](https://github.com/Eyescale/vmmlib/pull/49):
  Add lookAt support for 4x4 matrices

# Release 1.10 (07-04-2016)

* [49](https://github.com/Eyescale/vmmlib/pull/49):
  Add lookAt support for 4x4 matrices
* Improved API and test coverage

# Release 1.9 (02-11-2015)

* [35](https://github.com/Eyescale/vmmlib/pull/35):
  Fix Vector::rotate, cross and compute_normal, add Matrix-from-iter-ctor
* [26](https://github.com/Eyescale/vmmlib/pull/26):
  Clean up Vector::get_sub_vector()
* Fix AABB default ctor for unsigned values
* [24](https://github.com/Eyescale/vmmlib/pull/24):
  Extended Vector typedefs
* [23](https://github.com/Eyescale/vmmlib/pull/23):
  Implement Vector::product

# Release 1.8 (03-10-2014)

* Sanitization of matrix::get_translation API
* Compilation options and warnings cleanup
* Fix for shadowed member variables
