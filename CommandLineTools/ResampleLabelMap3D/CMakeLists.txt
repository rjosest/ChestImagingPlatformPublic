cmake_minimum_required(VERSION 2.8.0)

project(ResampleLabelMap3D )

set(MODULE_NAME ResampleLabelMap3D)


find_package(SlicerExecutionModel REQUIRED)
include(${SlicerExecutionModel_USE_FILE})

 FIND_PACKAGE( ITK REQUIRED )

 FIND_PACKAGE( CIP REQUIRED )

set(MODULE_INCLUDE_DIRECTORIES
  ${CIP_INCLUDE_DIRECTORIES}
  )

set(MODULE_TARGET_LIBRARIES
  ${ITK_LIBRARIES}
  ${CIP_LIBRARIES}
  )

set(MODULE_SRCS
  )

SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  INCLUDE_DIRECTORIES ${MODULE_INCLUDE_DIRECTORIES}
  ADDITIONAL_SRCS ${MODULE_SRCS}
  RUNTIME_OUTPUT_DIRECTORY        "${CMAKE_BINARY_DIR}/bin"
  EXECUTABLE_ONLY
  )
