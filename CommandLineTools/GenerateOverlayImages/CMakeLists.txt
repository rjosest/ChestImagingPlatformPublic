CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT( GenerateOverlayImages )

SET ( MODULE_NAME GenerateOverlayImages )

FIND_PACKAGE ( SlicerExecutionModel REQUIRED )
INCLUDE ( ${SlicerExecutionModel_USE_FILE} )

FIND_PACKAGE( CIP REQUIRED )

SET ( MODULE_INCLUDE_DIRECTORIES
    ${CMAKE_SOURCE_DIR}/Common
  )

SET ( MODULE_TARGET_LIBRARIES
  ${ITK_LIBRARIES}
  CIPCommon
  )

SET ( MODULE_SRCS GenerateOverlayImages.cxx )

SEMMacroBuildCLI (
  NAME ${MODULE_NAME}
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  INCLUDE_DIRECTORIES ${MODULE_INCLUDE_DIRECTORIES}
  ADDITIONAL_SRCS ${MODULE_SRCS}
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
  EXECUTABLE_ONLY
  )
