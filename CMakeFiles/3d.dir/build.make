# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/qiheng/phd/obstacles

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/qiheng/phd/obstacles

# Include any dependencies generated for this target.
include CMakeFiles/3d.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/3d.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/3d.dir/flags.make

CMakeFiles/3d.dir/3dregions.cpp.o: CMakeFiles/3d.dir/flags.make
CMakeFiles/3d.dir/3dregions.cpp.o: 3dregions.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qiheng/phd/obstacles/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/3d.dir/3dregions.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/3d.dir/3dregions.cpp.o -c /home/qiheng/phd/obstacles/3dregions.cpp

CMakeFiles/3d.dir/3dregions.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/3d.dir/3dregions.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/qiheng/phd/obstacles/3dregions.cpp > CMakeFiles/3d.dir/3dregions.cpp.i

CMakeFiles/3d.dir/3dregions.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/3d.dir/3dregions.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/qiheng/phd/obstacles/3dregions.cpp -o CMakeFiles/3d.dir/3dregions.cpp.s

# Object files for target 3d
3d_OBJECTS = \
"CMakeFiles/3d.dir/3dregions.cpp.o"

# External object files for target 3d
3d_EXTERNAL_OBJECTS =

3d: CMakeFiles/3d.dir/3dregions.cpp.o
3d: CMakeFiles/3d.dir/build.make
3d: CMakeFiles/3d.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/qiheng/phd/obstacles/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable 3d"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/3d.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/3d.dir/build: 3d

.PHONY : CMakeFiles/3d.dir/build

CMakeFiles/3d.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/3d.dir/cmake_clean.cmake
.PHONY : CMakeFiles/3d.dir/clean

CMakeFiles/3d.dir/depend:
	cd /home/qiheng/phd/obstacles && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/qiheng/phd/obstacles /home/qiheng/phd/obstacles /home/qiheng/phd/obstacles /home/qiheng/phd/obstacles /home/qiheng/phd/obstacles/CMakeFiles/3d.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/3d.dir/depend

