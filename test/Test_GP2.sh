#!/bin/bash
#
# Test the applications under ../app/GP : GP_Standard_TestFunction.ex
# and GP_RLCM_TestFunction.ex. If this script is called by
# Test_All.sh, comment out all the exports on NumThreads and Valgrind.

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
export Dim=(40 40)
export Lower=(0 0)
export Upper=(1 1)
export noise=1e-2
export List_sigma1=(0.05 0.1 0.15)
export List_sigma2=(0.05 0.1 0.15)
export List_tau=(-5 -3 -1)
#export List_sigma1=(0.1)
#export List_sigma2=(0.1)
#export List_tau=(-4)
export Seed=1
export PortionTrain=0.5
export IsCheckFiniteDiff=1
export DiffStepSize=(1e-5 1e-5 1e-5)

export Num_sigma1=${#List_sigma1[*]}
export Num_sigma2=${#List_sigma2[*]}
export Num_tau=${#List_tau[*]}

export r=125
export DiagCorrect=1e-8

export OutputRandomField=1
export RandomFieldFileBasename="Test_GP2_Standard_RF"
export OutputLogLik=1
export LogLikFileName="Test_GP2_Standard_LogLik.txt"
export ComputeFisher=1
export OutputFisher=1
export FisherFileName="Test_GP2_Standard_Fisher.txt"
export OutputKrigedRandomField=1
export KrigedRandomFieldFileBasename="Test_GP2_Standard_RF_kriged"
export OutputPredictions=1
export PredictionsFileName="Test_GP2_Standard_Pred.txt"

(set -x; ${Valgrind} ${Dir}/GP_Standard_TestFunction.ex ${NumThreads} ${Dim[*]} ${Lower[*]} ${Upper[*]} ${noise} ${Num_sigma1} ${List_sigma1[*]} ${Num_sigma2} ${List_sigma2[*]} ${Num_tau} ${List_tau[*]} ${Seed} ${PortionTrain} ${IsCheckFiniteDiff} ${DiffStepSize[*]} ${OutputRandomField} ${RandomFieldFileBasename} ${OutputLogLik} ${LogLikFileName} ${ComputeFisher} ${OutputFisher} ${FisherFileName} ${OutputKrigedRandomField} ${KrigedRandomFieldFileBasename} ${OutputPredictions} ${PredictionsFileName})

export OutputRandomField=1
export RandomFieldFileBasename="Test_GP2_RLCM_RF"
export OutputLogLik=1
export LogLikFileName="Test_GP2_RLCM_LogLik.txt"
export ComputeFisher=1
export OutputFisher=1
export FisherFileName="Test_GP2_RLCM_Fisher.txt"
export OutputKrigedRandomField=1
export KrigedRandomFieldFileBasename="Test_GP2_RLCM_RF_kriged"
export OutputPredictions=1
export PredictionsFileName="Test_GP2_RLCM_Pred.txt"

(set -x; ${Valgrind} ${Dir}/GP_RLCM_TestFunction.ex ${NumThreads} ${Dim[*]} ${Lower[*]} ${Upper[*]} ${noise} ${Num_sigma1} ${List_sigma1[*]} ${Num_sigma2} ${List_sigma2[*]} ${Num_tau} ${List_tau[*]} ${r} ${DiagCorrect} ${Seed} ${PortionTrain} ${IsCheckFiniteDiff} ${DiffStepSize[*]} ${OutputRandomField} ${RandomFieldFileBasename} ${OutputLogLik} ${LogLikFileName} ${ComputeFisher} ${OutputFisher} ${FisherFileName} ${OutputKrigedRandomField} ${KrigedRandomFieldFileBasename} ${OutputPredictions} ${PredictionsFileName})

