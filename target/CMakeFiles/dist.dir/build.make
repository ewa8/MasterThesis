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
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.22.3/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.22.3/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/ewarusiecka/Documents/University/DTU/Special_course/MSN/target

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/ewarusiecka/Documents/University/DTU/Special_course/MSN/target

# Utility rule file for dist.

# Include any custom commands dependencies for this target.
include CMakeFiles/dist.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/dist.dir/progress.make

CMakeFiles/dist:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/Users/ewarusiecka/Documents/University/DTU/Special_course/MSN/target/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Creating a source distribution from nestml_7572c325c24847f4afc90a8b52d50e81_module..."
	/usr/bin/make package_source

dist: CMakeFiles/dist
dist: CMakeFiles/dist.dir/build.make
.PHONY : dist

# Rule to build all files generated by this target.
CMakeFiles/dist.dir/build: dist
.PHONY : CMakeFiles/dist.dir/build

CMakeFiles/dist.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/dist.dir/cmake_clean.cmake
.PHONY : CMakeFiles/dist.dir/clean

CMakeFiles/dist.dir/depend:
	cd /Users/ewarusiecka/Documents/University/DTU/Special_course/MSN/target && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/ewarusiecka/Documents/University/DTU/Special_course/MSN/target /Users/ewarusiecka/Documents/University/DTU/Special_course/MSN/target /Users/ewarusiecka/Documents/University/DTU/Special_course/MSN/target /Users/ewarusiecka/Documents/University/DTU/Special_course/MSN/target /Users/ewarusiecka/Documents/University/DTU/Special_course/MSN/target/CMakeFiles/dist.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/dist.dir/depend

