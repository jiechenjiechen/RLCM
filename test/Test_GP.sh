#!/bin/bash
#
# Test the applications under ../app/GP : GP_Standard.ex and
# GP_RLCM.ex. If this script is called by Test_All.sh, comment out all
# the exports on NumThreads and Valgrind.

#export NumThreads=1
#export Valgrind=

#export NumThreads=3
#export Valgrind=

#export NumThreads=3
#export Valgrind="valgrind --leak-check=full --show-reachable=yes"
#
# NOTE: Valgrind check is very slow. Reduce the number of ell, nu, and
# tau.

export OMP_NUM_THREADS=${NumThreads}

export Dir="../app/GP"
export d=2
export Dim=(40 50)
export Lower=(-0.8 -1.0)
export Upper=(+0.8 +1.0)
export ell=0.5
export nu=1.5
export tau=-3
export List_ell=(0.1 0.5 1.0)
export List_nu=(0.5 1.5 2.5)
export List_tau=(-5 -3 -1)
#export List_ell=(0.1 1.0)
#export List_nu=(0.5 2.5)
#export List_tau=(-1)
export Seed=1
export PortionTrain=0.5
export IsCheckFiniteDiff=1
export DiffStepSize=(1e-5 1e-5 1e-5)

export Num_ell=${#List_ell[*]}
export Num_nu=${#List_nu[*]}
export Num_tau=${#List_tau[*]}

export r=250
export DiagCorrect=1e-8

export OutputRandomField=1
export RandomFieldFileBasename="Test_GP_Standard_RF"
export OutputLogLik=1
export LogLikFileName="Test_GP_Standard_LogLik.txt"
export ComputeFisher=1
export OutputFisher=1
export FisherFileName="Test_GP_Standard_Fisher.txt"
export OutputKrigedRandomField=1
export KrigedRandomFieldFileBasename="Test_GP_Standard_RF_kriged"
export OutputPredictions=1
export PredictionsFileName="Test_GP_Standard_Pred.txt"

(set -x; ${Valgrind} ${Dir}/GP_Standard.ex ${NumThreads} ${d} ${Dim[*]} ${Lower[*]} ${Upper[*]} ${ell} ${nu} ${tau} ${Num_ell} ${List_ell[*]} ${Num_nu} ${List_nu[*]} ${Num_tau} ${List_tau[*]} ${Seed} ${PortionTrain} ${IsCheckFiniteDiff} ${DiffStepSize[*]} ${OutputRandomField} ${RandomFieldFileBasename} ${OutputLogLik} ${LogLikFileName} ${ComputeFisher} ${OutputFisher} ${FisherFileName} ${OutputKrigedRandomField} ${KrigedRandomFieldFileBasename} ${OutputPredictions} ${PredictionsFileName})

export OutputRandomField=1
export RandomFieldFileBasename="Test_GP_RLCM_RF"
export OutputLogLik=1
export LogLikFileName="Test_GP_RLCM_LogLik.txt"
export ComputeFisher=1
export OutputFisher=1
export FisherFileName="Test_GP_RLCM_Fisher.txt"
export OutputKrigedRandomField=1
export KrigedRandomFieldFileBasename="Test_GP_RLCM_RF_kriged"
export OutputPredictions=1
export PredictionsFileName="Test_GP_RLCM_Pred.txt"

(set -x; ${Valgrind} ${Dir}/GP_RLCM.ex ${NumThreads} ${d} ${Dim[*]} ${Lower[*]} ${Upper[*]} ${ell} ${nu} ${tau} ${Num_ell} ${List_ell[*]} ${Num_nu} ${List_nu[*]} ${Num_tau} ${List_tau[*]} ${r} ${DiagCorrect} ${Seed} ${PortionTrain} ${IsCheckFiniteDiff} ${DiffStepSize[*]} ${OutputRandomField} ${RandomFieldFileBasename} ${OutputLogLik} ${LogLikFileName} ${ComputeFisher} ${OutputFisher} ${FisherFileName} ${OutputKrigedRandomField} ${KrigedRandomFieldFileBasename} ${OutputPredictions} ${PredictionsFileName})
