add_executable (OpenSMOKE_GasKinetics tests/OpenSMOKE_GasKinetics.cpp)
target_link_libraries (OpenSMOKE_GasKinetics PRIVATE opensmoke)

add_executable (OpenSMOKE_ODESolver tests/OpenSMOKE_ODESolver.cpp)
target_link_libraries (OpenSMOKE_ODESolver PRIVATE opensmoke)

add_executable (OpenSMOKE_SolidKinetics tests/OpenSMOKE_SolidKinetics.cpp)
target_link_libraries (OpenSMOKE_SolidKinetics PRIVATE opensmoke)

add_executable (OpenSMOKE_SurfaceKinetics tests/OpenSMOKE_SurfaceKinetics.cpp)
target_link_libraries (OpenSMOKE_SurfaceKinetics PRIVATE opensmoke)

add_test (NAME OpenSMOKE_GasKinetics
          COMMAND bash -c "
            ./OpenSMOKE_GasKinetics 2> OpenSMOKE_GasKinetics.log
            diff OpenSMOKE_GasKinetics.log ../tests/OpenSMOKE_GasKinetics.ref")

add_test (NAME OpenSMOKE_ODESolver
          COMMAND bash -c "
            ./OpenSMOKE_ODESolver 2> OpenSMOKE_ODESolver.log
            diff OpenSMOKE_ODESolver.log ../tests/OpenSMOKE_ODESolver.ref")
        
add_test (NAME OpenSMOKE_SolidKinetics
          COMMAND bash -c "
            ./OpenSMOKE_SolidKinetics 2> OpenSMOKE_SolidKinetics.log
            diff OpenSMOKE_SolidKinetics.log ../tests/OpenSMOKE_SolidKinetics.ref")

add_test (NAME OpenSMOKE_SurfaceKinetics
          COMMAND bash -c "
            ./OpenSMOKE_SurfaceKinetics 2> OpenSMOKE_SurfaceKinetics.log
            diff OpenSMOKE_SurfaceKinetics.log ../tests/OpenSMOKE_SurfaceKinetics.ref")

add_executable (OpenSMOKE_DRG tests/OpenSMOKE_DRG.cpp)
target_link_libraries (OpenSMOKE_DRG PRIVATE opensmoke)

add_test (NAME OpenSMOKE_DRG
          COMMAND bash -c "
            ./OpenSMOKE_DRG 2> OpenSMOKE_DRG.log
            diff OpenSMOKE_DRG.log ../tests/OpenSMOKE_DRG.ref")

add_test (NAME OpenSMOKE_PythonInterface
          COMMAND bash -c "python3 ../tests/OpenSMOKE_PythonInterface.py")

add_test (NAME OpenSMOKE_PythonInterfaceSolid
          COMMAND bash -c "python3 ../tests/OpenSMOKE_PythonInterfaceSolid.py")


add_executable (OpenSMOKE_GasJacobian tests/OpenSMOKE_GasJacobian.cpp)
target_link_libraries (OpenSMOKE_GasJacobian PRIVATE opensmoke)

add_test (NAME OpenSMOKE_GasJacobian
          COMMAND bash -c "
            ./OpenSMOKE_GasJacobian 2> OpenSMOKE_GasJacobian.log
            diff OpenSMOKE_GasJacobian.log ../tests/OpenSMOKE_GasJacobian.ref")

add_executable (OpenSMOKE_ODESolverJacobian tests/OpenSMOKE_ODESolverJacobian.cpp)
target_link_libraries (OpenSMOKE_ODESolverJacobian PRIVATE opensmoke)

add_test (NAME OpenSMOKE_ODESolverJacobian
          COMMAND bash -c "
            ./OpenSMOKE_ODESolverJacobian 2> OpenSMOKE_ODESolverJacobian.log
            diff OpenSMOKE_ODESolverJacobian.log ../tests/OpenSMOKE_ODESolverJacobian.ref")
