include_directories(BEFORE ../../include/ionic_models)

add_library(utils INTERFACE)
target_link_libraries(utils INTERFACE Eigen3::Eigen)

# Ionic models
add_library(ionic_model IonicModel.cpp)
target_link_libraries(ionic_model PUBLIC Eigen3::Eigen)

add_library(myokit_ionic_model MyoKitIonicModel.cpp)
target_link_libraries(myokit_ionic_model PUBLIC ionic_model)
