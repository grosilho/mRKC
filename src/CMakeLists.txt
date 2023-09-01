
include_directories(BEFORE ../include)

add_library(MainHeader INTERFACE)
target_link_libraries(MainHeader INTERFACE Eigen3::Eigen)
#target_include_directories(MainHeader INTERFACE)

add_subdirectory(lean_vtk)
add_subdirectory(ionic_models)

add_library(chebyshev_methods ChebyshevMethods.cpp)
target_include_directories(chebyshev_methods PUBLIC MainHeader)
target_link_libraries(chebyshev_methods PUBLIC MainHeader)

add_library(legendre_methods LegendreMethods.cpp)
target_include_directories(legendre_methods PUBLIC MainHeader)
target_link_libraries(legendre_methods PUBLIC MainHeader)

# ODEs
add_library(ode Ode.cpp)
target_link_libraries(ode PUBLIC MainHeader)
target_include_directories(ode PUBLIC MainHeader)

add_library(ode_problems OdeProblems.cpp)
target_include_directories(ode_problems PUBLIC ../include/lean_vtk ../include/ionic_models)
target_link_libraries(ode_problems PUBLIC ode lean_vtk myokit_ionic_model)

add_library(multirate_ode_problems MultirateOdeProblems.cpp)
target_link_libraries(multirate_ode_problems PUBLIC ode_problems)

# Others
add_library(parameters Parameters.cpp)
target_include_directories(parameters PUBLIC ../include/lean_vtk ../include/ionic_models ../include/getpot)
target_link_libraries(parameters PUBLIC                             
                            multirate_ode_rungekutta_integrator
                            multirate_ode_problems
                            stabilized_ode_rungekutta_integrator
                            classical_ode_rungekutta_integrator
                            )


# ODE Integrators
add_library(time_integrator TimeIntegrator.cpp)
target_link_libraries(time_integrator PUBLIC ode parameters)

add_library(error_controller ErrorController.cpp)
target_link_libraries(error_controller PUBLIC parameters)

add_library(ode_rungekutta_integrator OdeRungeKuttaIntegrator.cpp)
target_link_libraries(ode_rungekutta_integrator PUBLIC time_integrator error_controller)

add_library(multirate_ode_rungekutta_integrator MultirateOdeRungeKuttaIntegrator.cpp)
target_link_libraries(multirate_ode_rungekutta_integrator PUBLIC ode_rungekutta_integrator)

add_library(stabilized_ode_rungekutta_integrator StabilizedOdeRungeKuttaIntegrators.cpp)
target_link_libraries(stabilized_ode_rungekutta_integrator PUBLIC ode_rungekutta_integrator multirate_ode_rungekutta_integrator chebyshev_methods legendre_methods)

add_library(classical_ode_rungekutta_integrator ClassicalOdeRungeKuttaIntegrators.cpp)
target_link_libraries(classical_ode_rungekutta_integrator PUBLIC ode_rungekutta_integrator)

