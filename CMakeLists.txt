cmake_minimum_required(VERSION 3.20)
project(ThetaOptimality VERSION 0.6.2 LANGUAGES CXX)
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
add_executable(decay_qualification examples/decay_qualification.cpp)
target_link_libraries(decay_qualification PRIVATE theta_opt)
add_executable(two_temperature_qualification examples/two_temperature_qualification.cpp)
target_link_libraries(two_temperature_qualification PRIVATE theta_opt)
add_executable(oscillator examples/oscillator.cpp)
target_link_libraries(oscillator PRIVATE theta_opt)
# The optional mixture example links GPL-3.0-or-later ThermoGPU code.
option(THETA_BUILD_GAS_EXAMPLE "Build the ThermoGPU composition example" ON)
set(THERMOGPU_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/third_party/ThermoGPU" CACHE PATH "ThermoGPU source tree")
if(THETA_BUILD_GAS_EXAMPLE)
    add_executable(peng_robinson_inverse examples/peng_robinson_inverse.cpp
        ${THERMOGPU_ROOT}/src/cpu/pr_scalar.cpp)
    target_include_directories(peng_robinson_inverse PRIVATE ${THERMOGPU_ROOT}/include)
    target_link_libraries(peng_robinson_inverse PRIVATE theta_opt)
    if(BUILD_TESTING)
        add_test(NAME peng_robinson_composition_qualification COMMAND peng_robinson_inverse ${CMAKE_CURRENT_BINARY_DIR}/gas_results)
    endif()
endif()
if(BUILD_TESTING)
    add_executable(theta_validation tests/validation.cpp)
    target_link_libraries(theta_validation PRIVATE theta_opt)
    add_test(NAME discrete_optimality_validation COMMAND theta_validation)
    add_test(NAME oscillator_hidden_state_validation COMMAND oscillator ${CMAKE_CURRENT_BINARY_DIR}/oscillator_results)
    add_test(NAME decay_state_qualification COMMAND decay_qualification ${CMAKE_CURRENT_BINARY_DIR}/qualification_results)
    add_test(NAME two_temperature_state_qualification COMMAND two_temperature_qualification ${CMAKE_CURRENT_BINARY_DIR}/thermal_results)
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
