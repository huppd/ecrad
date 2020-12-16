set(driver_report_c_SRC Report.cpp)

PREPEND(driver_report_c_SRC ${CMAKE_CURRENT_SOURCE_DIR} ${driver_report_c_SRC})

target_sources(driver_report_c PUBLIC ${driver_report_c_SRC})

set(driver_report_fortran_wrapper_SRC ecrad_driver_report.f90)

PREPEND(driver_report_fortran_wrapper_SRC ${CMAKE_CURRENT_SOURCE_DIR} ${driver_report_fortran_wrapper_SRC})

target_sources(driver_report_fortran_wrapper PUBLIC ${driver_report_fortran_wrapper_SRC})

