#!/bin/bash
#
# Test the applications under ../app/KRR_binary . If this script is
# called by Test_All.sh, comment out all the exports on NumThreads and
# Valgrind.

#export NumThreads=1
#export Valgrind=

#export NumThreads=3
#export Valgrind=

#export NumThreads=3
#export Valgrind="valgrind --leak-check=full --show-reachable=yes"

export OMP_NUM_THREADS=${NumThreads}

export FileTrain_List=("Toy_regr_binary.train" "Toy_bin_binary.train" "Toy_mult_binary.train")
export FileTest_List=("Toy_regr_binary.test" "Toy_bin_binary.test" "Toy_mult_binary.test")
export NumClasses_List=(1 2 3)

export Dir="../app/KRR_binary"
export KernelType="IsotropicGaussian"
export d=2
export Seed=1
export Rank=6
export Budget=162
export List_sigma=(0.3 1)
export List_lambda=(1 1e-1)
export Par="RAND"
export DiagCorrect=1e-8
export Refinement=0

export N0=${Rank}
export Num_sigma=${#List_sigma[*]}
export Num_lambda=${#List_lambda[*]}

for idx in ${!NumClasses_List[*]}; do
FileTrain=${FileTrain_List[$idx]}
FileTest=${FileTest_List[$idx]}
NumClasses=${NumClasses_List[$idx]}

(set -x; ${Valgrind} ${Dir}/KRR_Standard_${KernelType}_binary.ex ${NumThreads} ${FileTrain} ${FileTest} ${NumClasses} ${d} ${Budget} ${Num_sigma} ${List_sigma[*]} ${Num_lambda} ${List_lambda[*]})

(set -x; ${Valgrind} ${Dir}/KRR_Fourier_${KernelType}_binary.ex ${NumThreads} ${FileTrain} ${FileTest} ${NumClasses} ${d} ${Seed} ${Rank} ${Budget} ${Num_sigma} ${List_sigma[*]} ${Num_lambda} ${List_lambda[*]})

(set -x; ${Valgrind} ${Dir}/KRR_Nystrom_${KernelType}_binary.ex ${NumThreads} ${FileTrain} ${FileTest} ${NumClasses} ${d} ${Seed} ${Rank} ${Budget} ${Num_sigma} ${List_sigma[*]} ${Num_lambda} ${List_lambda[*]})

(set -x; ${Valgrind} ${Dir}/KRR_RLCM_${KernelType}_binary.ex ${NumThreads} ${FileTrain} ${FileTest} ${NumClasses} ${d} ${Seed} ${Rank} ${Budget} ${Num_sigma} ${List_sigma[*]} ${Num_lambda} ${List_lambda[*]} ${Par} ${DiagCorrect} ${Refinement})

(set -x; ${Valgrind} ${Dir}/KRR_BlockDiag_${KernelType}_binary.ex ${NumThreads} ${FileTrain} ${FileTest} ${NumClasses} ${d} ${Seed} ${N0} ${Budget} ${Num_sigma} ${List_sigma[*]} ${Num_lambda} ${List_lambda[*]} ${Par})

done
