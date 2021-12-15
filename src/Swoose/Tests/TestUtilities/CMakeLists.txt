# Build a shared library module from MockModule and the model
add_library(MockModule SHARED
            ${CMAKE_CURRENT_SOURCE_DIR}/MockModule.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/MockModule.h
            ${CMAKE_CURRENT_SOURCE_DIR}/MockQmCalculator.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/MockQmCalculator.h
            )
set_target_properties(MockModule PROPERTIES
                      POSITION_INDEPENDENT_CODE ON
                      OUTPUT_NAME mockmodule
                      )
target_link_libraries(MockModule PUBLIC Scine::UtilsOS Scine::Swoose)
