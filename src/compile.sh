#!/bin/bash


#Source folders
zeroD="./0D"
integrator="./integrator"
adapters="./include/adapters"
chem_config="./."

chemgen="$HOME/abhijeet/05_chemgen/draft_1/src"

sundials_include="$HOME/abhijeet/10_CVODES/sundials/include"
sundials_build_include="$HOME/abhijeet/10_CVODES/sundials_build_dir/include"
sundials_build_src_sundials="$HOME/abhijeet/10_CVODES/sundials_build_dir/src/sundials"
sundials_build_src_cvodes="$HOME/abhijeet/10_CVODES/sundials_build_dir/src/cvodes"
sundials_build_src_sunlinsol="$HOME/abhijeet/10_CVODES/sundials_build_dir/src/sunlinsol/dense"
sundials_build_src_sunmatrix="$HOME/abhijeet/10_CVODES/sundials_build_dir/src/sunmatrix/dense"
sundials_build_src_nvector_serial="$HOME/abhijeet/10_CVODES/sundials_build_dir/src/nvector/serial"

#Source files
main="main.cpp"

zeroD_reactor="$zeroD/IdealGasConstPressureAdiabaticReactor.cpp"
integrator_CVODESSerialIntegrator="$integrator/CVODESSerialIntegrator.cpp"


#Executable
exec_name="reactor_test"

#Compile command
g++ -I $zeroD -I $integrator -I $adapters -I $chem_config  \
                                            \
    -I $chemgen                             \
                                            \
    -I $sundials_build_include              \
    -I $sundials_include                    \
                                            \
    -Wfatal-errors                          \
    -o $exec_name                           \
                                            \
    -L "$sundials_build_src_sundials"         \
    -L "$sundials_build_src_cvodes"         \
    -L "$sundials_build_src_nvector_serial" \
    -L "$sundials_build_src_sunlinsol"      \
    -L "$sundials_build_src_sunmatrix"      \
    -lsundials_core                         \
    -lsundials_cvodes                       \
    -lsundials_sunmatrixdense               \
    -lsundials_nvecserial                   \
    -lsundials_sunlinsoldense               \
                                            \
    $main                                   \
    $zeroD_reactor                          \
    $integrator_CVODESSerialIntegrator      \
    

