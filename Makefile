# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canoncical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /projects/lmi/people/jross/Downloads/cmake/cmake-2.8.6-Linux-i386/bin/cmake

# The command to remove a file.
RM = /projects/lmi/people/jross/Downloads/cmake/cmake-2.8.6-Linux-i386/bin/cmake -E remove -f

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /projects/lmi/people/jross/Downloads/cmake/cmake-2.8.6-Linux-i386/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /projects/lmi/people/jross/Downloads/ChestImagingPlatformPrivate

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /projects/lmi/people/jross/Downloads/ChestImagingPlatformPrivate

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/projects/lmi/people/jross/Downloads/cmake/cmake-2.8.6-Linux-i386/bin/ccmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/projects/lmi/people/jross/Downloads/cmake/cmake-2.8.6-Linux-i386/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /projects/lmi/people/jross/Downloads/ChestImagingPlatformPrivate/CMakeFiles /projects/lmi/people/jross/Downloads/ChestImagingPlatformPrivate/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /projects/lmi/people/jross/Downloads/ChestImagingPlatformPrivate/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean
.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named CIPUtilities

# Build rule for target.
CIPUtilities: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 CIPUtilities
.PHONY : CIPUtilities

# fast build rule for target.
CIPUtilities/fast:
	$(MAKE) -f Utilities/CMakeFiles/CIPUtilities.dir/build.make Utilities/CMakeFiles/CIPUtilities.dir/build
.PHONY : CIPUtilities/fast

#=============================================================================
# Target rules for targets named CIPCommon

# Build rule for target.
CIPCommon: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 CIPCommon
.PHONY : CIPCommon

# fast build rule for target.
CIPCommon/fast:
	$(MAKE) -f Common/CMakeFiles/CIPCommon.dir/build.make Common/CMakeFiles/CIPCommon.dir/build
.PHONY : CIPCommon/fast

#=============================================================================
# Target rules for targets named CIPIO

# Build rule for target.
CIPIO: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 CIPIO
.PHONY : CIPIO

# fast build rule for target.
CIPIO/fast:
	$(MAKE) -f IO/CMakeFiles/CIPIO.dir/build.make IO/CMakeFiles/CIPIO.dir/build
.PHONY : CIPIO/fast

#=============================================================================
# Target rules for targets named ConvertDicom

# Build rule for target.
ConvertDicom: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 ConvertDicom
.PHONY : ConvertDicom

# fast build rule for target.
ConvertDicom/fast:
	$(MAKE) -f CommandLineTools/ConvertDicom/CMakeFiles/ConvertDicom.dir/build.make CommandLineTools/ConvertDicom/CMakeFiles/ConvertDicom.dir/build
.PHONY : ConvertDicom/fast

#=============================================================================
# Target rules for targets named MergeParticleDataSets

# Build rule for target.
MergeParticleDataSets: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 MergeParticleDataSets
.PHONY : MergeParticleDataSets

# fast build rule for target.
MergeParticleDataSets/fast:
	$(MAKE) -f CommandLineTools/MergeParticleDataSets/CMakeFiles/MergeParticleDataSets.dir/build.make CommandLineTools/MergeParticleDataSets/CMakeFiles/MergeParticleDataSets.dir/build
.PHONY : MergeParticleDataSets/fast

#=============================================================================
# Target rules for targets named RemoveParticlesFromParticlesDataSet

# Build rule for target.
RemoveParticlesFromParticlesDataSet: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 RemoveParticlesFromParticlesDataSet
.PHONY : RemoveParticlesFromParticlesDataSet

# fast build rule for target.
RemoveParticlesFromParticlesDataSet/fast:
	$(MAKE) -f CommandLineTools/RemoveParticlesFromParticlesDataSet/CMakeFiles/RemoveParticlesFromParticlesDataSet.dir/build.make CommandLineTools/RemoveParticlesFromParticlesDataSet/CMakeFiles/RemoveParticlesFromParticlesDataSet.dir/build
.PHONY : RemoveParticlesFromParticlesDataSet/fast

#=============================================================================
# Target rules for targets named ExtractParticlesFromChestRegionChestType

# Build rule for target.
ExtractParticlesFromChestRegionChestType: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 ExtractParticlesFromChestRegionChestType
.PHONY : ExtractParticlesFromChestRegionChestType

# fast build rule for target.
ExtractParticlesFromChestRegionChestType/fast:
	$(MAKE) -f CommandLineTools/ExtractParticlesFromChestRegionChestType/CMakeFiles/ExtractParticlesFromChestRegionChestType.dir/build.make CommandLineTools/ExtractParticlesFromChestRegionChestType/CMakeFiles/ExtractParticlesFromChestRegionChestType.dir/build
.PHONY : ExtractParticlesFromChestRegionChestType/fast

#=============================================================================
# Target rules for targets named GenerateOtsuLungCast

# Build rule for target.
GenerateOtsuLungCast: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 GenerateOtsuLungCast
.PHONY : GenerateOtsuLungCast

# fast build rule for target.
GenerateOtsuLungCast/fast:
	$(MAKE) -f CommandLineTools/GenerateOtsuLungCast/CMakeFiles/GenerateOtsuLungCast.dir/build.make CommandLineTools/GenerateOtsuLungCast/CMakeFiles/GenerateOtsuLungCast.dir/build
.PHONY : GenerateOtsuLungCast/fast

#=============================================================================
# Target rules for targets named GenerateAtlasConvexHull

# Build rule for target.
GenerateAtlasConvexHull: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 GenerateAtlasConvexHull
.PHONY : GenerateAtlasConvexHull

# fast build rule for target.
GenerateAtlasConvexHull/fast:
	$(MAKE) -f CommandLineTools/GenerateAtlasConvexHull/CMakeFiles/GenerateAtlasConvexHull.dir/build.make CommandLineTools/GenerateAtlasConvexHull/CMakeFiles/GenerateAtlasConvexHull.dir/build
.PHONY : GenerateAtlasConvexHull/fast

#=============================================================================
# Target rules for targets named RegisterLungAtlas

# Build rule for target.
RegisterLungAtlas: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 RegisterLungAtlas
.PHONY : RegisterLungAtlas

# fast build rule for target.
RegisterLungAtlas/fast:
	$(MAKE) -f CommandLineTools/RegisterLungAtlas/CMakeFiles/RegisterLungAtlas.dir/build.make CommandLineTools/RegisterLungAtlas/CMakeFiles/RegisterLungAtlas.dir/build
.PHONY : RegisterLungAtlas/fast

#=============================================================================
# Target rules for targets named ResampleLabelMap

# Build rule for target.
ResampleLabelMap: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 ResampleLabelMap
.PHONY : ResampleLabelMap

# fast build rule for target.
ResampleLabelMap/fast:
	$(MAKE) -f CommandLineTools/ResampleLabelMap/CMakeFiles/ResampleLabelMap.dir/build.make CommandLineTools/ResampleLabelMap/CMakeFiles/ResampleLabelMap.dir/build
.PHONY : ResampleLabelMap/fast

#=============================================================================
# Target rules for targets named QualityControl

# Build rule for target.
QualityControl: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 QualityControl
.PHONY : QualityControl

# fast build rule for target.
QualityControl/fast:
	$(MAKE) -f CommandLineTools/QualityControl/CMakeFiles/QualityControl.dir/build.make CommandLineTools/QualityControl/CMakeFiles/QualityControl.dir/build
.PHONY : QualityControl/fast

#=============================================================================
# Target rules for targets named RemoveChestTypeFromLabelMapUsingParticles

# Build rule for target.
RemoveChestTypeFromLabelMapUsingParticles: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 RemoveChestTypeFromLabelMapUsingParticles
.PHONY : RemoveChestTypeFromLabelMapUsingParticles

# fast build rule for target.
RemoveChestTypeFromLabelMapUsingParticles/fast:
	$(MAKE) -f CommandLineTools/RemoveChestTypeFromLabelMapUsingParticles/CMakeFiles/RemoveChestTypeFromLabelMapUsingParticles.dir/build.make CommandLineTools/RemoveChestTypeFromLabelMapUsingParticles/CMakeFiles/RemoveChestTypeFromLabelMapUsingParticles.dir/build
.PHONY : RemoveChestTypeFromLabelMapUsingParticles/fast

#=============================================================================
# Target rules for targets named SegmentLungLobes

# Build rule for target.
SegmentLungLobes: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 SegmentLungLobes
.PHONY : SegmentLungLobes

# fast build rule for target.
SegmentLungLobes/fast:
	$(MAKE) -f CommandLineTools/SegmentLungLobes/CMakeFiles/SegmentLungLobes.dir/build.make CommandLineTools/SegmentLungLobes/CMakeFiles/SegmentLungLobes.dir/build
.PHONY : SegmentLungLobes/fast

#=============================================================================
# Target rules for targets named ReadDicomWriteTags

# Build rule for target.
ReadDicomWriteTags: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 ReadDicomWriteTags
.PHONY : ReadDicomWriteTags

# fast build rule for target.
ReadDicomWriteTags/fast:
	$(MAKE) -f CommandLineTools/ReadDicomWriteTags/CMakeFiles/ReadDicomWriteTags.dir/build.make CommandLineTools/ReadDicomWriteTags/CMakeFiles/ReadDicomWriteTags.dir/build
.PHONY : ReadDicomWriteTags/fast

#=============================================================================
# Target rules for targets named SplitLeftLungRightLung

# Build rule for target.
SplitLeftLungRightLung: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 SplitLeftLungRightLung
.PHONY : SplitLeftLungRightLung

# fast build rule for target.
SplitLeftLungRightLung/fast:
	$(MAKE) -f CommandLineTools/SplitLeftLungRightLung/CMakeFiles/SplitLeftLungRightLung.dir/build.make CommandLineTools/SplitLeftLungRightLung/CMakeFiles/SplitLeftLungRightLung.dir/build
.PHONY : SplitLeftLungRightLung/fast

#=============================================================================
# Target rules for targets named ExtractChestLabelMap

# Build rule for target.
ExtractChestLabelMap: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 ExtractChestLabelMap
.PHONY : ExtractChestLabelMap

# fast build rule for target.
ExtractChestLabelMap/fast:
	$(MAKE) -f CommandLineTools/ExtractChestLabelMap/CMakeFiles/ExtractChestLabelMap.dir/build.make CommandLineTools/ExtractChestLabelMap/CMakeFiles/ExtractChestLabelMap.dir/build
.PHONY : ExtractChestLabelMap/fast

#=============================================================================
# Target rules for targets named LabelParticlesByChestRegionChestType

# Build rule for target.
LabelParticlesByChestRegionChestType: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 LabelParticlesByChestRegionChestType
.PHONY : LabelParticlesByChestRegionChestType

# fast build rule for target.
LabelParticlesByChestRegionChestType/fast:
	$(MAKE) -f CommandLineTools/LabelParticlesByChestRegionChestType/CMakeFiles/LabelParticlesByChestRegionChestType.dir/build.make CommandLineTools/LabelParticlesByChestRegionChestType/CMakeFiles/LabelParticlesByChestRegionChestType.dir/build
.PHONY : LabelParticlesByChestRegionChestType/fast

#=============================================================================
# Target rules for targets named GenerateStenciledLabelMapFromParticles

# Build rule for target.
GenerateStenciledLabelMapFromParticles: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 GenerateStenciledLabelMapFromParticles
.PHONY : GenerateStenciledLabelMapFromParticles

# fast build rule for target.
GenerateStenciledLabelMapFromParticles/fast:
	$(MAKE) -f CommandLineTools/GenerateStenciledLabelMapFromParticles/CMakeFiles/GenerateStenciledLabelMapFromParticles.dir/build.make CommandLineTools/GenerateStenciledLabelMapFromParticles/CMakeFiles/GenerateStenciledLabelMapFromParticles.dir/build
.PHONY : GenerateStenciledLabelMapFromParticles/fast

#=============================================================================
# Target rules for targets named GenerateDistanceMapFromLabelMap

# Build rule for target.
GenerateDistanceMapFromLabelMap: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 GenerateDistanceMapFromLabelMap
.PHONY : GenerateDistanceMapFromLabelMap

# fast build rule for target.
GenerateDistanceMapFromLabelMap/fast:
	$(MAKE) -f CommandLineTools/GenerateDistanceMapFromLabelMap/CMakeFiles/GenerateDistanceMapFromLabelMap.dir/build.make CommandLineTools/GenerateDistanceMapFromLabelMap/CMakeFiles/GenerateDistanceMapFromLabelMap.dir/build
.PHONY : GenerateDistanceMapFromLabelMap/fast

#=============================================================================
# Target rules for targets named FilterVesselParticleData

# Build rule for target.
FilterVesselParticleData: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 FilterVesselParticleData
.PHONY : FilterVesselParticleData

# fast build rule for target.
FilterVesselParticleData/fast:
	$(MAKE) -f CommandLineTools/FilterVesselParticleData/CMakeFiles/FilterVesselParticleData.dir/build.make CommandLineTools/FilterVesselParticleData/CMakeFiles/FilterVesselParticleData.dir/build
.PHONY : FilterVesselParticleData/fast

#=============================================================================
# Target rules for targets named FilterAirwayParticleData

# Build rule for target.
FilterAirwayParticleData: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 FilterAirwayParticleData
.PHONY : FilterAirwayParticleData

# fast build rule for target.
FilterAirwayParticleData/fast:
	$(MAKE) -f CommandLineTools/FilterAirwayParticleData/CMakeFiles/FilterAirwayParticleData.dir/build.make CommandLineTools/FilterAirwayParticleData/CMakeFiles/FilterAirwayParticleData.dir/build
.PHONY : FilterAirwayParticleData/fast

#=============================================================================
# Target rules for targets named ScourParticleData

# Build rule for target.
ScourParticleData: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 ScourParticleData
.PHONY : ScourParticleData

# fast build rule for target.
ScourParticleData/fast:
	$(MAKE) -f CommandLineTools/ScourParticleData/CMakeFiles/ScourParticleData.dir/build.make CommandLineTools/ScourParticleData/CMakeFiles/ScourParticleData.dir/build
.PHONY : ScourParticleData/fast

#=============================================================================
# Target rules for targets named ReadNRRDsWriteVTK

# Build rule for target.
ReadNRRDsWriteVTK: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 ReadNRRDsWriteVTK
.PHONY : ReadNRRDsWriteVTK

# fast build rule for target.
ReadNRRDsWriteVTK/fast:
	$(MAKE) -f CommandLineTools/ReadNRRDsWriteVTK/CMakeFiles/ReadNRRDsWriteVTK.dir/build.make CommandLineTools/ReadNRRDsWriteVTK/CMakeFiles/ReadNRRDsWriteVTK.dir/build
.PHONY : ReadNRRDsWriteVTK/fast

#=============================================================================
# Target rules for targets named ComputeCrossSectionalArea

# Build rule for target.
ComputeCrossSectionalArea: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 ComputeCrossSectionalArea
.PHONY : ComputeCrossSectionalArea

# fast build rule for target.
ComputeCrossSectionalArea/fast:
	$(MAKE) -f CommandLineTools/ComputeCrossSectionalArea/CMakeFiles/ComputeCrossSectionalArea.dir/build.make CommandLineTools/ComputeCrossSectionalArea/CMakeFiles/ComputeCrossSectionalArea.dir/build
.PHONY : ComputeCrossSectionalArea/fast

#=============================================================================
# Target rules for targets named ComputeIntensityStatistics

# Build rule for target.
ComputeIntensityStatistics: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 ComputeIntensityStatistics
.PHONY : ComputeIntensityStatistics

# fast build rule for target.
ComputeIntensityStatistics/fast:
	$(MAKE) -f CommandLineTools/ComputeIntensityStatistics/CMakeFiles/ComputeIntensityStatistics.dir/build.make CommandLineTools/ComputeIntensityStatistics/CMakeFiles/ComputeIntensityStatistics.dir/build
.PHONY : ComputeIntensityStatistics/fast

#=============================================================================
# Target rules for targets named FilterFissureParticleData

# Build rule for target.
FilterFissureParticleData: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 FilterFissureParticleData
.PHONY : FilterFissureParticleData

# fast build rule for target.
FilterFissureParticleData/fast:
	$(MAKE) -f CommandLineTools/FilterFissureParticleData/CMakeFiles/FilterFissureParticleData.dir/build.make CommandLineTools/FilterFissureParticleData/CMakeFiles/FilterFissureParticleData.dir/build
.PHONY : FilterFissureParticleData/fast

#=============================================================================
# Target rules for targets named GenerateFissureShapeModels

# Build rule for target.
GenerateFissureShapeModels: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 GenerateFissureShapeModels
.PHONY : GenerateFissureShapeModels

# fast build rule for target.
GenerateFissureShapeModels/fast:
	$(MAKE) -f CommandLineTools/GenerateFissureShapeModels/CMakeFiles/GenerateFissureShapeModels.dir/build.make CommandLineTools/GenerateFissureShapeModels/CMakeFiles/GenerateFissureShapeModels.dir/build
.PHONY : GenerateFissureShapeModels/fast

#=============================================================================
# Target rules for targets named FitFissureSurfaceModelToParticleData

# Build rule for target.
FitFissureSurfaceModelToParticleData: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 FitFissureSurfaceModelToParticleData
.PHONY : FitFissureSurfaceModelToParticleData

# fast build rule for target.
FitFissureSurfaceModelToParticleData/fast:
	$(MAKE) -f CommandLineTools/FitFissureSurfaceModelToParticleData/CMakeFiles/FitFissureSurfaceModelToParticleData.dir/build.make CommandLineTools/FitFissureSurfaceModelToParticleData/CMakeFiles/FitFissureSurfaceModelToParticleData.dir/build
.PHONY : FitFissureSurfaceModelToParticleData/fast

#=============================================================================
# Target rules for targets named FitFissureSurfaceModelsToParticleData

# Build rule for target.
FitFissureSurfaceModelsToParticleData: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 FitFissureSurfaceModelsToParticleData
.PHONY : FitFissureSurfaceModelsToParticleData

# fast build rule for target.
FitFissureSurfaceModelsToParticleData/fast:
	$(MAKE) -f CommandLineTools/FitFissureSurfaceModelsToParticleData/CMakeFiles/FitFissureSurfaceModelsToParticleData.dir/build.make CommandLineTools/FitFissureSurfaceModelsToParticleData/CMakeFiles/FitFissureSurfaceModelsToParticleData.dir/build
.PHONY : FitFissureSurfaceModelsToParticleData/fast

#=============================================================================
# Target rules for targets named ClassifyFissureParticles

# Build rule for target.
ClassifyFissureParticles: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 ClassifyFissureParticles
.PHONY : ClassifyFissureParticles

# fast build rule for target.
ClassifyFissureParticles/fast:
	$(MAKE) -f CommandLineTools/ClassifyFissureParticles/CMakeFiles/ClassifyFissureParticles.dir/build.make CommandLineTools/ClassifyFissureParticles/CMakeFiles/ClassifyFissureParticles.dir/build
.PHONY : ClassifyFissureParticles/fast

#=============================================================================
# Target rules for targets named CropLung

# Build rule for target.
CropLung: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 CropLung
.PHONY : CropLung

# fast build rule for target.
CropLung/fast:
	$(MAKE) -f CommandLineTools/CropLung/CMakeFiles/CropLung.dir/build.make CommandLineTools/CropLung/CMakeFiles/CropLung.dir/build
.PHONY : CropLung/fast

#=============================================================================
# Target rules for targets named GenerateModel

# Build rule for target.
GenerateModel: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 GenerateModel
.PHONY : GenerateModel

# fast build rule for target.
GenerateModel/fast:
	$(MAKE) -f CommandLineTools/GenerateModel/CMakeFiles/GenerateModel.dir/build.make CommandLineTools/GenerateModel/CMakeFiles/GenerateModel.dir/build
.PHONY : GenerateModel/fast

#=============================================================================
# Target rules for targets named GenerateImageSubVolumes

# Build rule for target.
GenerateImageSubVolumes: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 GenerateImageSubVolumes
.PHONY : GenerateImageSubVolumes

# fast build rule for target.
GenerateImageSubVolumes/fast:
	$(MAKE) -f CommandLineTools/GenerateImageSubVolumes/CMakeFiles/GenerateImageSubVolumes.dir/build.make CommandLineTools/GenerateImageSubVolumes/CMakeFiles/GenerateImageSubVolumes.dir/build
.PHONY : GenerateImageSubVolumes/fast

#=============================================================================
# Target rules for targets named GenerateMedianFilteredImage

# Build rule for target.
GenerateMedianFilteredImage: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 GenerateMedianFilteredImage
.PHONY : GenerateMedianFilteredImage

# fast build rule for target.
GenerateMedianFilteredImage/fast:
	$(MAKE) -f CommandLineTools/GenerateMedianFilteredImage/CMakeFiles/GenerateMedianFilteredImage.dir/build.make CommandLineTools/GenerateMedianFilteredImage/CMakeFiles/GenerateMedianFilteredImage.dir/build
.PHONY : GenerateMedianFilteredImage/fast

#=============================================================================
# Target rules for targets named GenerateStatisticsForAirwayGenerationLabeling

# Build rule for target.
GenerateStatisticsForAirwayGenerationLabeling: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 GenerateStatisticsForAirwayGenerationLabeling
.PHONY : GenerateStatisticsForAirwayGenerationLabeling

# fast build rule for target.
GenerateStatisticsForAirwayGenerationLabeling/fast:
	$(MAKE) -f CommandLineTools/GenerateStatisticsForAirwayGenerationLabeling/CMakeFiles/GenerateStatisticsForAirwayGenerationLabeling.dir/build.make CommandLineTools/GenerateStatisticsForAirwayGenerationLabeling/CMakeFiles/GenerateStatisticsForAirwayGenerationLabeling.dir/build
.PHONY : GenerateStatisticsForAirwayGenerationLabeling/fast

#=============================================================================
# Target rules for targets named EvaluateLungLobeSegmentationResults

# Build rule for target.
EvaluateLungLobeSegmentationResults: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 EvaluateLungLobeSegmentationResults
.PHONY : EvaluateLungLobeSegmentationResults

# fast build rule for target.
EvaluateLungLobeSegmentationResults/fast:
	$(MAKE) -f CommandLineTools/EvaluateLungLobeSegmentationResults/CMakeFiles/EvaluateLungLobeSegmentationResults.dir/build.make CommandLineTools/EvaluateLungLobeSegmentationResults/CMakeFiles/EvaluateLungLobeSegmentationResults.dir/build
.PHONY : EvaluateLungLobeSegmentationResults/fast

#=============================================================================
# Target rules for targets named ConvertLabelMapValueToChestRegionChestType

# Build rule for target.
ConvertLabelMapValueToChestRegionChestType: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 ConvertLabelMapValueToChestRegionChestType
.PHONY : ConvertLabelMapValueToChestRegionChestType

# fast build rule for target.
ConvertLabelMapValueToChestRegionChestType/fast:
	$(MAKE) -f CommandLineTools/ConvertLabelMapValueToChestRegionChestType/CMakeFiles/ConvertLabelMapValueToChestRegionChestType.dir/build.make CommandLineTools/ConvertLabelMapValueToChestRegionChestType/CMakeFiles/ConvertLabelMapValueToChestRegionChestType.dir/build
.PHONY : ConvertLabelMapValueToChestRegionChestType/fast

#=============================================================================
# Target rules for targets named GenerateRegionHistogramsAndParenchymaPhenotypes

# Build rule for target.
GenerateRegionHistogramsAndParenchymaPhenotypes: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 GenerateRegionHistogramsAndParenchymaPhenotypes
.PHONY : GenerateRegionHistogramsAndParenchymaPhenotypes

# fast build rule for target.
GenerateRegionHistogramsAndParenchymaPhenotypes/fast:
	$(MAKE) -f CommandLineTools/GenerateRegionHistogramsAndParenchymaPhenotypes/CMakeFiles/GenerateRegionHistogramsAndParenchymaPhenotypes.dir/build.make CommandLineTools/GenerateRegionHistogramsAndParenchymaPhenotypes/CMakeFiles/GenerateRegionHistogramsAndParenchymaPhenotypes.dir/build
.PHONY : GenerateRegionHistogramsAndParenchymaPhenotypes/fast

#=============================================================================
# Target rules for targets named ComputeDistanceMap

# Build rule for target.
ComputeDistanceMap: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 ComputeDistanceMap
.PHONY : ComputeDistanceMap

# fast build rule for target.
ComputeDistanceMap/fast:
	$(MAKE) -f CommandLineTools/ComputeDistanceMap/CMakeFiles/ComputeDistanceMap.dir/build.make CommandLineTools/ComputeDistanceMap/CMakeFiles/ComputeDistanceMap.dir/build
.PHONY : ComputeDistanceMap/fast

#=============================================================================
# Target rules for targets named TransferFieldDataToFromPointData

# Build rule for target.
TransferFieldDataToFromPointData: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 TransferFieldDataToFromPointData
.PHONY : TransferFieldDataToFromPointData

# fast build rule for target.
TransferFieldDataToFromPointData/fast:
	$(MAKE) -f CommandLineTools/TransferFieldDataToFromPointData/CMakeFiles/TransferFieldDataToFromPointData.dir/build.make CommandLineTools/TransferFieldDataToFromPointData/CMakeFiles/TransferFieldDataToFromPointData.dir/build
.PHONY : TransferFieldDataToFromPointData/fast

#=============================================================================
# Target rules for targets named PerturbParticles

# Build rule for target.
PerturbParticles: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 PerturbParticles
.PHONY : PerturbParticles

# fast build rule for target.
PerturbParticles/fast:
	$(MAKE) -f CommandLineTools/PerturbParticles/CMakeFiles/PerturbParticles.dir/build.make CommandLineTools/PerturbParticles/CMakeFiles/PerturbParticles.dir/build
.PHONY : PerturbParticles/fast

#=============================================================================
# Target rules for targets named ReadWriteRegionAndTypePoints

# Build rule for target.
ReadWriteRegionAndTypePoints: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 ReadWriteRegionAndTypePoints
.PHONY : ReadWriteRegionAndTypePoints

# fast build rule for target.
ReadWriteRegionAndTypePoints/fast:
	$(MAKE) -f CommandLineTools/ReadWriteRegionAndTypePoints/CMakeFiles/ReadWriteRegionAndTypePoints.dir/build.make CommandLineTools/ReadWriteRegionAndTypePoints/CMakeFiles/ReadWriteRegionAndTypePoints.dir/build
.PHONY : ReadWriteRegionAndTypePoints/fast

#=============================================================================
# Target rules for targets named PerformMorphological

# Build rule for target.
PerformMorphological: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 PerformMorphological
.PHONY : PerformMorphological

# fast build rule for target.
PerformMorphological/fast:
	$(MAKE) -f CommandLineTools/PerformMorphological/CMakeFiles/PerformMorphological.dir/build.make CommandLineTools/PerformMorphological/CMakeFiles/PerformMorphological.dir/build
.PHONY : PerformMorphological/fast

#=============================================================================
# Target rules for targets named GenerateOverlayImages

# Build rule for target.
GenerateOverlayImages: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 GenerateOverlayImages
.PHONY : GenerateOverlayImages

# fast build rule for target.
GenerateOverlayImages/fast:
	$(MAKE) -f CommandLineTools/GenerateOverlayImages/CMakeFiles/GenerateOverlayImages.dir/build.make CommandLineTools/GenerateOverlayImages/CMakeFiles/GenerateOverlayImages.dir/build
.PHONY : GenerateOverlayImages/fast

#=============================================================================
# Target rules for targets named ReadParticlesWriteConnectedParticles

# Build rule for target.
ReadParticlesWriteConnectedParticles: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 ReadParticlesWriteConnectedParticles
.PHONY : ReadParticlesWriteConnectedParticles

# fast build rule for target.
ReadParticlesWriteConnectedParticles/fast:
	$(MAKE) -f CommandLineTools/ReadParticlesWriteConnectedParticles/CMakeFiles/ReadParticlesWriteConnectedParticles.dir/build.make CommandLineTools/ReadParticlesWriteConnectedParticles/CMakeFiles/ReadParticlesWriteConnectedParticles.dir/build
.PHONY : ReadParticlesWriteConnectedParticles/fast

#=============================================================================
# Target rules for targets named ReadWriteImageData

# Build rule for target.
ReadWriteImageData: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 ReadWriteImageData
.PHONY : ReadWriteImageData

# fast build rule for target.
ReadWriteImageData/fast:
	$(MAKE) -f CommandLineTools/ReadWriteImageData/CMakeFiles/ReadWriteImageData.dir/build.make CommandLineTools/ReadWriteImageData/CMakeFiles/ReadWriteImageData.dir/build
.PHONY : ReadWriteImageData/fast

#=============================================================================
# Target rules for targets named ACILAssistant

# Build rule for target.
ACILAssistant: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 ACILAssistant
.PHONY : ACILAssistant

# fast build rule for target.
ACILAssistant/fast:
	$(MAKE) -f InteractiveTools/ACILAssistant/CMakeFiles/ACILAssistant.dir/build.make InteractiveTools/ACILAssistant/CMakeFiles/ACILAssistant.dir/build
.PHONY : ACILAssistant/fast

#=============================================================================
# Target rules for targets named ViewChestData

# Build rule for target.
ViewChestData: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 ViewChestData
.PHONY : ViewChestData

# fast build rule for target.
ViewChestData/fast:
	$(MAKE) -f InteractiveTools/ViewChestData/CMakeFiles/ViewChestData.dir/build.make InteractiveTools/ViewChestData/CMakeFiles/ViewChestData.dir/build
.PHONY : ViewChestData/fast

#=============================================================================
# Target rules for targets named EditVesselParticles

# Build rule for target.
EditVesselParticles: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 EditVesselParticles
.PHONY : EditVesselParticles

# fast build rule for target.
EditVesselParticles/fast:
	$(MAKE) -f InteractiveTools/EditVesselParticles/CMakeFiles/EditVesselParticles.dir/build.make InteractiveTools/EditVesselParticles/CMakeFiles/EditVesselParticles.dir/build
.PHONY : EditVesselParticles/fast

#=============================================================================
# Target rules for targets named EditAirwayParticles

# Build rule for target.
EditAirwayParticles: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 EditAirwayParticles
.PHONY : EditAirwayParticles

# fast build rule for target.
EditAirwayParticles/fast:
	$(MAKE) -f InteractiveTools/EditAirwayParticles/CMakeFiles/EditAirwayParticles.dir/build.make InteractiveTools/EditAirwayParticles/CMakeFiles/EditAirwayParticles.dir/build
.PHONY : EditAirwayParticles/fast

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... rebuild_cache"
	@echo "... CIPUtilities"
	@echo "... CIPCommon"
	@echo "... CIPIO"
	@echo "... ConvertDicom"
	@echo "... MergeParticleDataSets"
	@echo "... RemoveParticlesFromParticlesDataSet"
	@echo "... ExtractParticlesFromChestRegionChestType"
	@echo "... GenerateOtsuLungCast"
	@echo "... GenerateAtlasConvexHull"
	@echo "... RegisterLungAtlas"
	@echo "... ResampleLabelMap"
	@echo "... QualityControl"
	@echo "... RemoveChestTypeFromLabelMapUsingParticles"
	@echo "... SegmentLungLobes"
	@echo "... ReadDicomWriteTags"
	@echo "... SplitLeftLungRightLung"
	@echo "... ExtractChestLabelMap"
	@echo "... LabelParticlesByChestRegionChestType"
	@echo "... GenerateStenciledLabelMapFromParticles"
	@echo "... GenerateDistanceMapFromLabelMap"
	@echo "... FilterVesselParticleData"
	@echo "... FilterAirwayParticleData"
	@echo "... ScourParticleData"
	@echo "... ReadNRRDsWriteVTK"
	@echo "... ComputeCrossSectionalArea"
	@echo "... ComputeIntensityStatistics"
	@echo "... FilterFissureParticleData"
	@echo "... GenerateFissureShapeModels"
	@echo "... FitFissureSurfaceModelToParticleData"
	@echo "... FitFissureSurfaceModelsToParticleData"
	@echo "... ClassifyFissureParticles"
	@echo "... CropLung"
	@echo "... GenerateModel"
	@echo "... GenerateImageSubVolumes"
	@echo "... GenerateMedianFilteredImage"
	@echo "... GenerateStatisticsForAirwayGenerationLabeling"
	@echo "... EvaluateLungLobeSegmentationResults"
	@echo "... ConvertLabelMapValueToChestRegionChestType"
	@echo "... GenerateRegionHistogramsAndParenchymaPhenotypes"
	@echo "... ComputeDistanceMap"
	@echo "... TransferFieldDataToFromPointData"
	@echo "... PerturbParticles"
	@echo "... ReadWriteRegionAndTypePoints"
	@echo "... PerformMorphological"
	@echo "... GenerateOverlayImages"
	@echo "... ReadParticlesWriteConnectedParticles"
	@echo "... ReadWriteImageData"
	@echo "... ACILAssistant"
	@echo "... ViewChestData"
	@echo "... EditVesselParticles"
	@echo "... EditAirwayParticles"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system
