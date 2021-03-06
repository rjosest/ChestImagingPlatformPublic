CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT( EditVesselParticles )

SET( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}" )

FIND_PACKAGE( VTK )
IF ( VTK_FOUND )  
   INCLUDE ( ${VTK_USE_FILE} )    
ELSE ( VTK_FOUND )
  MESSAGE ( FATAL_ERROR "Cannot build without VTK" )
ENDIF( VTK_FOUND )

FIND_PACKAGE ( ITK )
IF ( ITK_FOUND )
  INCLUDE(${ITK_USE_FILE})
ELSE ( ITK_FOUND )
  MESSAGE ( FATAL_ERROR "Cannot build without ITK" )
ENDIF ( ITK_FOUND )

INCLUDE_DIRECTORIES(
${CMAKE_SOURCE_DIR}/Utilities
${CMAKE_SOURCE_DIR}/Common
)

SET ( VTK_LIBS
vtkCommon 
vtkImaging 
vtkFiltering 
vtkGraphics 
vtkRendering 
vtkWidgets 
vtkHybrid
)

SET ( EditVesselParticles_SOURCE EditVesselParticles.cxx )

ADD_EXECUTABLE ( EditVesselParticles ${EditVesselParticles_SOURCE})
TARGET_LINK_LIBRARIES ( EditVesselParticles ${ITK_LIBRARIES} ${VTK_LIBS} CIPUtilities CIPCommon )
