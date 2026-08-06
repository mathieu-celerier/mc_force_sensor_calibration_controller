add_fsm_state_simple(PressureCheck)
add_fsm_state_simple(ChooseTransition)
add_fsm_state_simple(ShowForces)
add_fsm_state_simple(InitialPosture)
add_fsm_state_simple(CalibrationMotion)
add_fsm_state_simple(CalibrationMotionLogging)
add_fsm_state_simple(RunCalibrationScript)
add_fsm_state_simple(CheckResults)

add_fsm_data_directory(data)

target_include_directories(RunCalibrationScript
                           PRIVATE "${PROJECT_BINARY_DIR}/src")

# These states use boost::filesystem directly (temp_directory_path, exists, create_directory,
# copy_file, ...) which add_fsm_state_simple does not link by default
target_link_libraries(CalibrationMotionLogging PRIVATE Boost::filesystem)
target_link_libraries(RunCalibrationScript PRIVATE Boost::filesystem)
target_link_libraries(CheckResults PRIVATE Boost::filesystem)
