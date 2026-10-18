cmake_minimum_required(VERSION 3.20)
project(ThetaOptimality VERSION 0.1.0 LANGUAGES CXX)
include(GNUInstallDirs)
include(CTest)
add_library(theta_opt INTERFACE)
add_library(ThetaOptimality::theta_opt ALIAS theta_opt)
target_compile_features(theta_opt INTERFACE cxx_std_20)
target_include_directories(theta_opt INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
add_executable(theta_decay examples/decay.cpp)
target_link_libraries(theta_decay PRIVATE theta_opt)
if(BUILD_TESTING)
    add_executable(theta_validation tests/validation.cpp)
    target_link_libraries(theta_validation PRIVATE theta_opt)
    add_test(NAME discrete_optimality_validation COMMAND theta_validation)
endif()
install(TARGETS theta_opt EXPORT ThetaOptimalityTargets)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT ThetaOptimalityTargets NAMESPACE ThetaOptimality::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ThetaOptimality)
include(CMakePackageConfigHelpers)
configure_package_config_file(cmake/ThetaOptimalityConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/ThetaOptimalityConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ThetaOptimality)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/ThetaOptimalityConfigVersion.cmake
    VERSION ${PROJECT_VERSION} COMPATIBILITY SameMajorVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/ThetaOptimalityConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/ThetaOptimalityConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ThetaOptimality)
