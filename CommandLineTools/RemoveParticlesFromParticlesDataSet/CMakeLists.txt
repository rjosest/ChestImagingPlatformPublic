CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT( RemoveParticlesFromParticlesDataSet )

SET ( MODULE_NAME RemoveParticlesFromParticlesDataSet )

FIND_PACKAGE ( SlicerExecutionModel REQUIRED )
INCLUDE ( ${SlicerExecutionModel_USE_FILE} )

FIND_PACKAGE( CIP REQUIRED )

SET ( MODULE_TARGET_LIBRARIES
  ${VTK_LIBRARIES}
  )

SET ( MODULE_SRCS RemoveParticlesFromParticlesDataSet.cxx )

SEMMacroBuildCLI (
  NAME ${MODULE_NAME}
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  INCLUDE_DIRECTORIES ${MODULE_INCLUDE_DIRECTORIES}
  ADDITIONAL_SRCS ${MODULE_SRCS}
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
  EXECUTABLE_ONLY
  )




