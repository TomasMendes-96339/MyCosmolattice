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
CMAKE_BINARY_DIR = /home/tomas_mendes/projects/cosmolattice/nmc_runs_12th_run

# Utility rule file for clean-cmake.

# Include any custom commands dependencies for this target.
include CMakeFiles/clean-cmake.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/clean-cmake.dir/progress.make

CMakeFiles/clean-cmake:
	/usr/bin/cmake -P /home/tomas_mendes/projects/cosmolattice/src/cmake/clean-all.cmake

clean-cmake: CMakeFiles/clean-cmake
clean-cmake: CMakeFiles/clean-cmake.dir/build.make
.PHONY : clean-cmake

# Rule to build all files generated by this target.
CMakeFiles/clean-cmake.dir/build: clean-cmake
.PHONY : CMakeFiles/clean-cmake.dir/build

CMakeFiles/clean-cmake.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/clean-cmake.dir/cmake_clean.cmake
.PHONY : CMakeFiles/clean-cmake.dir/clean

CMakeFiles/clean-cmake.dir/depend:
	cd /home/tomas_mendes/projects/cosmolattice/nmc_runs_12th_run && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/tomas_mendes/projects/cosmolattice /home/tomas_mendes/projects/cosmolattice /home/tomas_mendes/projects/cosmolattice/nmc_runs_12th_run /home/tomas_mendes/projects/cosmolattice/nmc_runs_12th_run /home/tomas_mendes/projects/cosmolattice/nmc_runs_12th_run/CMakeFiles/clean-cmake.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/clean-cmake.dir/depend

