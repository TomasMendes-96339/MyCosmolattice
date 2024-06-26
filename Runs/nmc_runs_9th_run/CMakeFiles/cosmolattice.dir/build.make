# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/tomas_mendes/projects/cosmolattice

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/tomas_mendes/projects/cosmolattice/nmc_runs_9th_run

# Include any dependencies generated for this target.
include CMakeFiles/cosmolattice.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/cosmolattice.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/cosmolattice.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/cosmolattice.dir/flags.make

CMakeFiles/cosmolattice.dir/src/cosmolattice.cpp.o: CMakeFiles/cosmolattice.dir/flags.make
CMakeFiles/cosmolattice.dir/src/cosmolattice.cpp.o: ../src/cosmolattice.cpp
CMakeFiles/cosmolattice.dir/src/cosmolattice.cpp.o: CMakeFiles/cosmolattice.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/tomas_mendes/projects/cosmolattice/nmc_runs_9th_run/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/cosmolattice.dir/src/cosmolattice.cpp.o"
	/usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/cosmolattice.dir/src/cosmolattice.cpp.o -MF CMakeFiles/cosmolattice.dir/src/cosmolattice.cpp.o.d -o CMakeFiles/cosmolattice.dir/src/cosmolattice.cpp.o -c /home/tomas_mendes/projects/cosmolattice/src/cosmolattice.cpp

CMakeFiles/cosmolattice.dir/src/cosmolattice.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cosmolattice.dir/src/cosmolattice.cpp.i"
	/usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/tomas_mendes/projects/cosmolattice/src/cosmolattice.cpp > CMakeFiles/cosmolattice.dir/src/cosmolattice.cpp.i

CMakeFiles/cosmolattice.dir/src/cosmolattice.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cosmolattice.dir/src/cosmolattice.cpp.s"
	/usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/tomas_mendes/projects/cosmolattice/src/cosmolattice.cpp -o CMakeFiles/cosmolattice.dir/src/cosmolattice.cpp.s

# Object files for target cosmolattice
cosmolattice_OBJECTS = \
"CMakeFiles/cosmolattice.dir/src/cosmolattice.cpp.o"

# External object files for target cosmolattice
cosmolattice_EXTERNAL_OBJECTS =

GB_wMass: CMakeFiles/cosmolattice.dir/src/cosmolattice.cpp.o
GB_wMass: CMakeFiles/cosmolattice.dir/build.make
GB_wMass: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_cxx.so
GB_wMass: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi.so
GB_wMass: ../dependencies/MyLibs/lib/libpfft.a
GB_wMass: ../dependencies/MyLibs/lib/libfftw3_threads.a
GB_wMass: ../dependencies/MyLibs/lib/libfftw3_mpi.a
GB_wMass: ../dependencies/MyLibs/lib/libfftw3.a
GB_wMass: CMakeFiles/cosmolattice.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/tomas_mendes/projects/cosmolattice/nmc_runs_9th_run/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable GB_wMass"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/cosmolattice.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/cosmolattice.dir/build: GB_wMass
.PHONY : CMakeFiles/cosmolattice.dir/build

CMakeFiles/cosmolattice.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/cosmolattice.dir/cmake_clean.cmake
.PHONY : CMakeFiles/cosmolattice.dir/clean

CMakeFiles/cosmolattice.dir/depend:
	cd /home/tomas_mendes/projects/cosmolattice/nmc_runs_9th_run && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/tomas_mendes/projects/cosmolattice /home/tomas_mendes/projects/cosmolattice /home/tomas_mendes/projects/cosmolattice/nmc_runs_9th_run /home/tomas_mendes/projects/cosmolattice/nmc_runs_9th_run /home/tomas_mendes/projects/cosmolattice/nmc_runs_9th_run/CMakeFiles/cosmolattice.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/cosmolattice.dir/depend

