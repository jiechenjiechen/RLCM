#!/bin/bash
#
# Perform unit tests.

#export NumThreads=1
#export Valgrind=

export NumThreads=3
export Valgrind=

#export NumThreads=3
#export Valgrind="valgrind --leak-check=full --show-reachable=yes"

export USE_MATERN=1

export OMP_NUM_THREADS=${NumThreads}

(set -x; $Valgrind Test_Common.ex ${NumThreads})
(set -x; $Valgrind Test_LibSVM_IO.ex ${NumThreads})
(set -x; $Valgrind Test_LibSVM_IO_binary.ex ${NumThreads})
(set -x; $Valgrind Test_METIS_IO.ex ${NumThreads})
(set -x; $Valgrind Test_CSR_IO_binary.ex ${NumThreads})
(set -x; $Valgrind Test_Raw_binary.ex ${NumThreads})
(set -x; $Valgrind Test_spblas.ex ${NumThreads})
(set -x; $Valgrind Test_PCG.ex ${NumThreads} 100)
(set -x; $Valgrind Test_PCG.ex ${NumThreads} 20)
(set -x; $Valgrind Test_GMRES.ex ${NumThreads} 100)
(set -x; $Valgrind Test_GMRES.ex ${NumThreads} 25)
(set -x; $Valgrind Test_Lanczos.ex ${NumThreads})
(set -x; $Valgrind Test_IsotropicGaussian.ex ${NumThreads})
(set -x; $Valgrind Test_IsotropicLaplace.ex ${NumThreads})
if [ "$USE_MATERN" = "1" ]; then
(set -x; $Valgrind Test_IsotropicMatern.ex ${NumThreads})
(set -x; $Valgrind Test_ChordalMatern.ex ${NumThreads})
fi
(set -x; $Valgrind Test_AnisotropicGaussian.ex ${NumThreads})
(set -x; $Valgrind Test_ProductLaplace.ex ${NumThreads})
(set -x; $Valgrind Test_InvMultiquadric.ex ${NumThreads})
(set -x; $Valgrind Test_Chi2.ex ${NumThreads})
(set -x; $Valgrind Test_Polynomial.ex ${NumThreads})
(set -x; $Valgrind Test_SixHumps.ex ${NumThreads})
(set -x; $Valgrind Test_Bird.ex ${NumThreads})
(set -x; $Valgrind Test_NoName.ex ${NumThreads})
(set -x; $Valgrind Test_DVector.ex ${NumThreads})
(set -x; $Valgrind Test_DMatrix.ex ${NumThreads})
(set -x; $Valgrind Test_SMatrix.ex ${NumThreads})
(set -x; $Valgrind Test_SPoint.ex ${NumThreads})
(set -x; $Valgrind Test_DPoint.ex ${NumThreads})
(set -x; $Valgrind Test_DPointArray.ex ${NumThreads})
(set -x; $Valgrind Test_SPointArray.ex ${NumThreads})
(set -x; $Valgrind Test_CMatrix.ex ${NumThreads} RAND)
(set -x; $Valgrind Test_CMatrix.ex ${NumThreads} BBOX)
# NOTE: Nothing is interesting to test about the class Node
(set -x; $Valgrind Test_BMatrix.ex ${NumThreads} RAND)
(set -x; $Valgrind Test_BMatrix.ex ${NumThreads} BBOX)

Test_KRR.sh
Test_KRR_binary.sh

if [ "$USE_MATERN" = "1" ]; then
Test_GP.sh
Test_GP2.sh
fi

