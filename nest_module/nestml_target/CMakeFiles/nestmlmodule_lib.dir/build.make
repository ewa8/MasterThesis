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
CMAKE_SOURCE_DIR = /Users/ewarusiecka/Documents/University/DTU/Special_course/MSN/nest_module/nestml_target

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/ewarusiecka/Documents/University/DTU/Special_course/MSN/nest_module/nestml_target

# Include any dependencies generated for this target.
include CMakeFiles/nestmlmodule_lib.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/nestmlmodule_lib.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/nestmlmodule_lib.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/nestmlmodule_lib.dir/flags.make

CMakeFiles/nestmlmodule_lib.dir/nestmlmodule.o: CMakeFiles/nestmlmodule_lib.dir/flags.make
CMakeFiles/nestmlmodule_lib.dir/nestmlmodule.o: nestmlmodule.cpp
CMakeFiles/nestmlmodule_lib.dir/nestmlmodule.o: CMakeFiles/nestmlmodule_lib.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/ewarusiecka/Documents/University/DTU/Special_course/MSN/nest_module/nestml_target/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/nestmlmodule_lib.dir/nestmlmodule.o"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/nestmlmodule_lib.dir/nestmlmodule.o -MF CMakeFiles/nestmlmodule_lib.dir/nestmlmodule.o.d -o CMakeFiles/nestmlmodule_lib.dir/nestmlmodule.o -c /Users/ewarusiecka/Documents/University/DTU/Special_course/MSN/nest_module/nestml_target/nestmlmodule.cpp

CMakeFiles/nestmlmodule_lib.dir/nestmlmodule.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/nestmlmodule_lib.dir/nestmlmodule.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/ewarusiecka/Documents/University/DTU/Special_course/MSN/nest_module/nestml_target/nestmlmodule.cpp > CMakeFiles/nestmlmodule_lib.dir/nestmlmodule.i

CMakeFiles/nestmlmodule_lib.dir/nestmlmodule.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/nestmlmodule_lib.dir/nestmlmodule.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/ewarusiecka/Documents/University/DTU/Special_course/MSN/nest_module/nestml_target/nestmlmodule.cpp -o CMakeFiles/nestmlmodule_lib.dir/nestmlmodule.s

CMakeFiles/nestmlmodule_lib.dir/msn.o: CMakeFiles/nestmlmodule_lib.dir/flags.make
CMakeFiles/nestmlmodule_lib.dir/msn.o: msn.cpp
CMakeFiles/nestmlmodule_lib.dir/msn.o: CMakeFiles/nestmlmodule_lib.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/ewarusiecka/Documents/University/DTU/Special_course/MSN/nest_module/nestml_target/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/nestmlmodule_lib.dir/msn.o"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/nestmlmodule_lib.dir/msn.o -MF CMakeFiles/nestmlmodule_lib.dir/msn.o.d -o CMakeFiles/nestmlmodule_lib.dir/msn.o -c /Users/ewarusiecka/Documents/University/DTU/Special_course/MSN/nest_module/nestml_target/msn.cpp

CMakeFiles/nestmlmodule_lib.dir/msn.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/nestmlmodule_lib.dir/msn.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/ewarusiecka/Documents/University/DTU/Special_course/MSN/nest_module/nestml_target/msn.cpp > CMakeFiles/nestmlmodule_lib.dir/msn.i

CMakeFiles/nestmlmodule_lib.dir/msn.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/nestmlmodule_lib.dir/msn.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/ewarusiecka/Documents/University/DTU/Special_course/MSN/nest_module/nestml_target/msn.cpp -o CMakeFiles/nestmlmodule_lib.dir/msn.s

# Object files for target nestmlmodule_lib
nestmlmodule_lib_OBJECTS = \
"CMakeFiles/nestmlmodule_lib.dir/nestmlmodule.o" \
"CMakeFiles/nestmlmodule_lib.dir/msn.o"

# External object files for target nestmlmodule_lib
nestmlmodule_lib_EXTERNAL_OBJECTS =

libnestmlmodule.dylib: CMakeFiles/nestmlmodule_lib.dir/nestmlmodule.o
libnestmlmodule.dylib: CMakeFiles/nestmlmodule_lib.dir/msn.o
libnestmlmodule.dylib: CMakeFiles/nestmlmodule_lib.dir/build.make
libnestmlmodule.dylib: CMakeFiles/nestmlmodule_lib.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/ewarusiecka/Documents/University/DTU/Special_course/MSN/nest_module/nestml_target/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX shared library libnestmlmodule.dylib"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/nestmlmodule_lib.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/nestmlmodule_lib.dir/build: libnestmlmodule.dylib
.PHONY : CMakeFiles/nestmlmodule_lib.dir/build

CMakeFiles/nestmlmodule_lib.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/nestmlmodule_lib.dir/cmake_clean.cmake
.PHONY : CMakeFiles/nestmlmodule_lib.dir/clean

CMakeFiles/nestmlmodule_lib.dir/depend:
	cd /Users/ewarusiecka/Documents/University/DTU/Special_course/MSN/nest_module/nestml_target && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/ewarusiecka/Documents/University/DTU/Special_course/MSN/nest_module/nestml_target /Users/ewarusiecka/Documents/University/DTU/Special_course/MSN/nest_module/nestml_target /Users/ewarusiecka/Documents/University/DTU/Special_course/MSN/nest_module/nestml_target /Users/ewarusiecka/Documents/University/DTU/Special_course/MSN/nest_module/nestml_target /Users/ewarusiecka/Documents/University/DTU/Special_course/MSN/nest_module/nestml_target/CMakeFiles/nestmlmodule_lib.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/nestmlmodule_lib.dir/depend

