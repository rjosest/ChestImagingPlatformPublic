CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT( ConvertDicom )

FIND_PACKAGE ( ITK )
IF ( ITK_FOUND )
  INCLUDE(${ITK_USE_FILE})
ELSE ( ITK_FOUND )
  MESSAGE ( FATAL_ERROR "Cannot build without ITK" )
ENDIF ( ITK_FOUND )

INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/Common )
INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/Utilities )

SET ( ConvertDicom_SOURCE ConvertDicom.cxx )

ADD_EXECUTABLE ( ConvertDicom ${ConvertDicom_SOURCE} )
TARGET_LINK_LIBRARIES ( ConvertDicom ${ITK_LIBRARIES}  )

IF ( CIP_BUILD_TESTING )
 SUBDIRS ( Testing )
ENDIF( CIP_BUILD_TESTING ) 
