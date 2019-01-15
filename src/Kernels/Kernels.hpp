// This file wraps all the header files in this directory

#ifndef _KERNELS_
#define _KERNELS_

#include "IsotropicGaussian.hpp"
#include "IsotropicLaplace.hpp"
#include "AnisotropicGaussian.hpp"
#include "ProductLaplace.hpp"
#ifdef USE_MATERN
#include "IsotropicMatern.hpp"
#include "ChordalMatern.hpp"
#endif
#include "InvMultiquadric.hpp"
#include "Chi2.hpp"
#include "Polynomial.hpp"

#endif
