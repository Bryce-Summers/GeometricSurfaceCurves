# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.2

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
CMAKE_SOURCE_DIR = /home/bryce/Documents/research

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/bryce/Documents/research/build

# Include any dependencies generated for this target.
include src/CMakeFiles/meshedit.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/meshedit.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/meshedit.dir/flags.make

src/CMakeFiles/meshedit.dir/UnionFind/UF_ADT.cpp.o: src/CMakeFiles/meshedit.dir/flags.make
src/CMakeFiles/meshedit.dir/UnionFind/UF_ADT.cpp.o: ../src/UnionFind/UF_ADT.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/bryce/Documents/research/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/CMakeFiles/meshedit.dir/UnionFind/UF_ADT.cpp.o"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/meshedit.dir/UnionFind/UF_ADT.cpp.o -c /home/bryce/Documents/research/src/UnionFind/UF_ADT.cpp

src/CMakeFiles/meshedit.dir/UnionFind/UF_ADT.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/meshedit.dir/UnionFind/UF_ADT.cpp.i"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/bryce/Documents/research/src/UnionFind/UF_ADT.cpp > CMakeFiles/meshedit.dir/UnionFind/UF_ADT.cpp.i

src/CMakeFiles/meshedit.dir/UnionFind/UF_ADT.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/meshedit.dir/UnionFind/UF_ADT.cpp.s"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/bryce/Documents/research/src/UnionFind/UF_ADT.cpp -o CMakeFiles/meshedit.dir/UnionFind/UF_ADT.cpp.s

src/CMakeFiles/meshedit.dir/UnionFind/UF_ADT.cpp.o.requires:
.PHONY : src/CMakeFiles/meshedit.dir/UnionFind/UF_ADT.cpp.o.requires

src/CMakeFiles/meshedit.dir/UnionFind/UF_ADT.cpp.o.provides: src/CMakeFiles/meshedit.dir/UnionFind/UF_ADT.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/meshedit.dir/build.make src/CMakeFiles/meshedit.dir/UnionFind/UF_ADT.cpp.o.provides.build
.PHONY : src/CMakeFiles/meshedit.dir/UnionFind/UF_ADT.cpp.o.provides

src/CMakeFiles/meshedit.dir/UnionFind/UF_ADT.cpp.o.provides.build: src/CMakeFiles/meshedit.dir/UnionFind/UF_ADT.cpp.o

src/CMakeFiles/meshedit.dir/UnionFind/UF_Serial.cpp.o: src/CMakeFiles/meshedit.dir/flags.make
src/CMakeFiles/meshedit.dir/UnionFind/UF_Serial.cpp.o: ../src/UnionFind/UF_Serial.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/bryce/Documents/research/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/CMakeFiles/meshedit.dir/UnionFind/UF_Serial.cpp.o"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/meshedit.dir/UnionFind/UF_Serial.cpp.o -c /home/bryce/Documents/research/src/UnionFind/UF_Serial.cpp

src/CMakeFiles/meshedit.dir/UnionFind/UF_Serial.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/meshedit.dir/UnionFind/UF_Serial.cpp.i"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/bryce/Documents/research/src/UnionFind/UF_Serial.cpp > CMakeFiles/meshedit.dir/UnionFind/UF_Serial.cpp.i

src/CMakeFiles/meshedit.dir/UnionFind/UF_Serial.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/meshedit.dir/UnionFind/UF_Serial.cpp.s"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/bryce/Documents/research/src/UnionFind/UF_Serial.cpp -o CMakeFiles/meshedit.dir/UnionFind/UF_Serial.cpp.s

src/CMakeFiles/meshedit.dir/UnionFind/UF_Serial.cpp.o.requires:
.PHONY : src/CMakeFiles/meshedit.dir/UnionFind/UF_Serial.cpp.o.requires

src/CMakeFiles/meshedit.dir/UnionFind/UF_Serial.cpp.o.provides: src/CMakeFiles/meshedit.dir/UnionFind/UF_Serial.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/meshedit.dir/build.make src/CMakeFiles/meshedit.dir/UnionFind/UF_Serial.cpp.o.provides.build
.PHONY : src/CMakeFiles/meshedit.dir/UnionFind/UF_Serial.cpp.o.provides

src/CMakeFiles/meshedit.dir/UnionFind/UF_Serial.cpp.o.provides.build: src/CMakeFiles/meshedit.dir/UnionFind/UF_Serial.cpp.o

src/CMakeFiles/meshedit.dir/scene.cpp.o: src/CMakeFiles/meshedit.dir/flags.make
src/CMakeFiles/meshedit.dir/scene.cpp.o: ../src/scene.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/bryce/Documents/research/build/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/CMakeFiles/meshedit.dir/scene.cpp.o"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/meshedit.dir/scene.cpp.o -c /home/bryce/Documents/research/src/scene.cpp

src/CMakeFiles/meshedit.dir/scene.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/meshedit.dir/scene.cpp.i"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/bryce/Documents/research/src/scene.cpp > CMakeFiles/meshedit.dir/scene.cpp.i

src/CMakeFiles/meshedit.dir/scene.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/meshedit.dir/scene.cpp.s"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/bryce/Documents/research/src/scene.cpp -o CMakeFiles/meshedit.dir/scene.cpp.s

src/CMakeFiles/meshedit.dir/scene.cpp.o.requires:
.PHONY : src/CMakeFiles/meshedit.dir/scene.cpp.o.requires

src/CMakeFiles/meshedit.dir/scene.cpp.o.provides: src/CMakeFiles/meshedit.dir/scene.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/meshedit.dir/build.make src/CMakeFiles/meshedit.dir/scene.cpp.o.provides.build
.PHONY : src/CMakeFiles/meshedit.dir/scene.cpp.o.provides

src/CMakeFiles/meshedit.dir/scene.cpp.o.provides.build: src/CMakeFiles/meshedit.dir/scene.cpp.o

src/CMakeFiles/meshedit.dir/camera.cpp.o: src/CMakeFiles/meshedit.dir/flags.make
src/CMakeFiles/meshedit.dir/camera.cpp.o: ../src/camera.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/bryce/Documents/research/build/CMakeFiles $(CMAKE_PROGRESS_4)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/CMakeFiles/meshedit.dir/camera.cpp.o"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/meshedit.dir/camera.cpp.o -c /home/bryce/Documents/research/src/camera.cpp

src/CMakeFiles/meshedit.dir/camera.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/meshedit.dir/camera.cpp.i"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/bryce/Documents/research/src/camera.cpp > CMakeFiles/meshedit.dir/camera.cpp.i

src/CMakeFiles/meshedit.dir/camera.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/meshedit.dir/camera.cpp.s"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/bryce/Documents/research/src/camera.cpp -o CMakeFiles/meshedit.dir/camera.cpp.s

src/CMakeFiles/meshedit.dir/camera.cpp.o.requires:
.PHONY : src/CMakeFiles/meshedit.dir/camera.cpp.o.requires

src/CMakeFiles/meshedit.dir/camera.cpp.o.provides: src/CMakeFiles/meshedit.dir/camera.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/meshedit.dir/build.make src/CMakeFiles/meshedit.dir/camera.cpp.o.provides.build
.PHONY : src/CMakeFiles/meshedit.dir/camera.cpp.o.provides

src/CMakeFiles/meshedit.dir/camera.cpp.o.provides.build: src/CMakeFiles/meshedit.dir/camera.cpp.o

src/CMakeFiles/meshedit.dir/light.cpp.o: src/CMakeFiles/meshedit.dir/flags.make
src/CMakeFiles/meshedit.dir/light.cpp.o: ../src/light.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/bryce/Documents/research/build/CMakeFiles $(CMAKE_PROGRESS_5)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/CMakeFiles/meshedit.dir/light.cpp.o"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/meshedit.dir/light.cpp.o -c /home/bryce/Documents/research/src/light.cpp

src/CMakeFiles/meshedit.dir/light.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/meshedit.dir/light.cpp.i"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/bryce/Documents/research/src/light.cpp > CMakeFiles/meshedit.dir/light.cpp.i

src/CMakeFiles/meshedit.dir/light.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/meshedit.dir/light.cpp.s"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/bryce/Documents/research/src/light.cpp -o CMakeFiles/meshedit.dir/light.cpp.s

src/CMakeFiles/meshedit.dir/light.cpp.o.requires:
.PHONY : src/CMakeFiles/meshedit.dir/light.cpp.o.requires

src/CMakeFiles/meshedit.dir/light.cpp.o.provides: src/CMakeFiles/meshedit.dir/light.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/meshedit.dir/build.make src/CMakeFiles/meshedit.dir/light.cpp.o.provides.build
.PHONY : src/CMakeFiles/meshedit.dir/light.cpp.o.provides

src/CMakeFiles/meshedit.dir/light.cpp.o.provides.build: src/CMakeFiles/meshedit.dir/light.cpp.o

src/CMakeFiles/meshedit.dir/mesh.cpp.o: src/CMakeFiles/meshedit.dir/flags.make
src/CMakeFiles/meshedit.dir/mesh.cpp.o: ../src/mesh.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/bryce/Documents/research/build/CMakeFiles $(CMAKE_PROGRESS_6)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/CMakeFiles/meshedit.dir/mesh.cpp.o"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/meshedit.dir/mesh.cpp.o -c /home/bryce/Documents/research/src/mesh.cpp

src/CMakeFiles/meshedit.dir/mesh.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/meshedit.dir/mesh.cpp.i"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/bryce/Documents/research/src/mesh.cpp > CMakeFiles/meshedit.dir/mesh.cpp.i

src/CMakeFiles/meshedit.dir/mesh.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/meshedit.dir/mesh.cpp.s"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/bryce/Documents/research/src/mesh.cpp -o CMakeFiles/meshedit.dir/mesh.cpp.s

src/CMakeFiles/meshedit.dir/mesh.cpp.o.requires:
.PHONY : src/CMakeFiles/meshedit.dir/mesh.cpp.o.requires

src/CMakeFiles/meshedit.dir/mesh.cpp.o.provides: src/CMakeFiles/meshedit.dir/mesh.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/meshedit.dir/build.make src/CMakeFiles/meshedit.dir/mesh.cpp.o.provides.build
.PHONY : src/CMakeFiles/meshedit.dir/mesh.cpp.o.provides

src/CMakeFiles/meshedit.dir/mesh.cpp.o.provides.build: src/CMakeFiles/meshedit.dir/mesh.cpp.o

src/CMakeFiles/meshedit.dir/material.cpp.o: src/CMakeFiles/meshedit.dir/flags.make
src/CMakeFiles/meshedit.dir/material.cpp.o: ../src/material.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/bryce/Documents/research/build/CMakeFiles $(CMAKE_PROGRESS_7)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/CMakeFiles/meshedit.dir/material.cpp.o"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/meshedit.dir/material.cpp.o -c /home/bryce/Documents/research/src/material.cpp

src/CMakeFiles/meshedit.dir/material.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/meshedit.dir/material.cpp.i"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/bryce/Documents/research/src/material.cpp > CMakeFiles/meshedit.dir/material.cpp.i

src/CMakeFiles/meshedit.dir/material.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/meshedit.dir/material.cpp.s"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/bryce/Documents/research/src/material.cpp -o CMakeFiles/meshedit.dir/material.cpp.s

src/CMakeFiles/meshedit.dir/material.cpp.o.requires:
.PHONY : src/CMakeFiles/meshedit.dir/material.cpp.o.requires

src/CMakeFiles/meshedit.dir/material.cpp.o.provides: src/CMakeFiles/meshedit.dir/material.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/meshedit.dir/build.make src/CMakeFiles/meshedit.dir/material.cpp.o.provides.build
.PHONY : src/CMakeFiles/meshedit.dir/material.cpp.o.provides

src/CMakeFiles/meshedit.dir/material.cpp.o.provides.build: src/CMakeFiles/meshedit.dir/material.cpp.o

src/CMakeFiles/meshedit.dir/texture.cpp.o: src/CMakeFiles/meshedit.dir/flags.make
src/CMakeFiles/meshedit.dir/texture.cpp.o: ../src/texture.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/bryce/Documents/research/build/CMakeFiles $(CMAKE_PROGRESS_8)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/CMakeFiles/meshedit.dir/texture.cpp.o"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/meshedit.dir/texture.cpp.o -c /home/bryce/Documents/research/src/texture.cpp

src/CMakeFiles/meshedit.dir/texture.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/meshedit.dir/texture.cpp.i"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/bryce/Documents/research/src/texture.cpp > CMakeFiles/meshedit.dir/texture.cpp.i

src/CMakeFiles/meshedit.dir/texture.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/meshedit.dir/texture.cpp.s"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/bryce/Documents/research/src/texture.cpp -o CMakeFiles/meshedit.dir/texture.cpp.s

src/CMakeFiles/meshedit.dir/texture.cpp.o.requires:
.PHONY : src/CMakeFiles/meshedit.dir/texture.cpp.o.requires

src/CMakeFiles/meshedit.dir/texture.cpp.o.provides: src/CMakeFiles/meshedit.dir/texture.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/meshedit.dir/build.make src/CMakeFiles/meshedit.dir/texture.cpp.o.provides.build
.PHONY : src/CMakeFiles/meshedit.dir/texture.cpp.o.provides

src/CMakeFiles/meshedit.dir/texture.cpp.o.provides.build: src/CMakeFiles/meshedit.dir/texture.cpp.o

src/CMakeFiles/meshedit.dir/collada.cpp.o: src/CMakeFiles/meshedit.dir/flags.make
src/CMakeFiles/meshedit.dir/collada.cpp.o: ../src/collada.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/bryce/Documents/research/build/CMakeFiles $(CMAKE_PROGRESS_9)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/CMakeFiles/meshedit.dir/collada.cpp.o"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/meshedit.dir/collada.cpp.o -c /home/bryce/Documents/research/src/collada.cpp

src/CMakeFiles/meshedit.dir/collada.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/meshedit.dir/collada.cpp.i"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/bryce/Documents/research/src/collada.cpp > CMakeFiles/meshedit.dir/collada.cpp.i

src/CMakeFiles/meshedit.dir/collada.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/meshedit.dir/collada.cpp.s"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/bryce/Documents/research/src/collada.cpp -o CMakeFiles/meshedit.dir/collada.cpp.s

src/CMakeFiles/meshedit.dir/collada.cpp.o.requires:
.PHONY : src/CMakeFiles/meshedit.dir/collada.cpp.o.requires

src/CMakeFiles/meshedit.dir/collada.cpp.o.provides: src/CMakeFiles/meshedit.dir/collada.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/meshedit.dir/build.make src/CMakeFiles/meshedit.dir/collada.cpp.o.provides.build
.PHONY : src/CMakeFiles/meshedit.dir/collada.cpp.o.provides

src/CMakeFiles/meshedit.dir/collada.cpp.o.provides.build: src/CMakeFiles/meshedit.dir/collada.cpp.o

src/CMakeFiles/meshedit.dir/halfEdgeMesh.cpp.o: src/CMakeFiles/meshedit.dir/flags.make
src/CMakeFiles/meshedit.dir/halfEdgeMesh.cpp.o: ../src/halfEdgeMesh.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/bryce/Documents/research/build/CMakeFiles $(CMAKE_PROGRESS_10)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/CMakeFiles/meshedit.dir/halfEdgeMesh.cpp.o"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/meshedit.dir/halfEdgeMesh.cpp.o -c /home/bryce/Documents/research/src/halfEdgeMesh.cpp

src/CMakeFiles/meshedit.dir/halfEdgeMesh.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/meshedit.dir/halfEdgeMesh.cpp.i"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/bryce/Documents/research/src/halfEdgeMesh.cpp > CMakeFiles/meshedit.dir/halfEdgeMesh.cpp.i

src/CMakeFiles/meshedit.dir/halfEdgeMesh.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/meshedit.dir/halfEdgeMesh.cpp.s"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/bryce/Documents/research/src/halfEdgeMesh.cpp -o CMakeFiles/meshedit.dir/halfEdgeMesh.cpp.s

src/CMakeFiles/meshedit.dir/halfEdgeMesh.cpp.o.requires:
.PHONY : src/CMakeFiles/meshedit.dir/halfEdgeMesh.cpp.o.requires

src/CMakeFiles/meshedit.dir/halfEdgeMesh.cpp.o.provides: src/CMakeFiles/meshedit.dir/halfEdgeMesh.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/meshedit.dir/build.make src/CMakeFiles/meshedit.dir/halfEdgeMesh.cpp.o.provides.build
.PHONY : src/CMakeFiles/meshedit.dir/halfEdgeMesh.cpp.o.provides

src/CMakeFiles/meshedit.dir/halfEdgeMesh.cpp.o.provides.build: src/CMakeFiles/meshedit.dir/halfEdgeMesh.cpp.o

src/CMakeFiles/meshedit.dir/meshResampler.cpp.o: src/CMakeFiles/meshedit.dir/flags.make
src/CMakeFiles/meshedit.dir/meshResampler.cpp.o: ../src/meshResampler.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/bryce/Documents/research/build/CMakeFiles $(CMAKE_PROGRESS_11)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/CMakeFiles/meshedit.dir/meshResampler.cpp.o"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/meshedit.dir/meshResampler.cpp.o -c /home/bryce/Documents/research/src/meshResampler.cpp

src/CMakeFiles/meshedit.dir/meshResampler.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/meshedit.dir/meshResampler.cpp.i"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/bryce/Documents/research/src/meshResampler.cpp > CMakeFiles/meshedit.dir/meshResampler.cpp.i

src/CMakeFiles/meshedit.dir/meshResampler.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/meshedit.dir/meshResampler.cpp.s"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/bryce/Documents/research/src/meshResampler.cpp -o CMakeFiles/meshedit.dir/meshResampler.cpp.s

src/CMakeFiles/meshedit.dir/meshResampler.cpp.o.requires:
.PHONY : src/CMakeFiles/meshedit.dir/meshResampler.cpp.o.requires

src/CMakeFiles/meshedit.dir/meshResampler.cpp.o.provides: src/CMakeFiles/meshedit.dir/meshResampler.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/meshedit.dir/build.make src/CMakeFiles/meshedit.dir/meshResampler.cpp.o.provides.build
.PHONY : src/CMakeFiles/meshedit.dir/meshResampler.cpp.o.provides

src/CMakeFiles/meshedit.dir/meshResampler.cpp.o.provides.build: src/CMakeFiles/meshedit.dir/meshResampler.cpp.o

src/CMakeFiles/meshedit.dir/patchDrawer.cpp.o: src/CMakeFiles/meshedit.dir/flags.make
src/CMakeFiles/meshedit.dir/patchDrawer.cpp.o: ../src/patchDrawer.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/bryce/Documents/research/build/CMakeFiles $(CMAKE_PROGRESS_12)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/CMakeFiles/meshedit.dir/patchDrawer.cpp.o"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/meshedit.dir/patchDrawer.cpp.o -c /home/bryce/Documents/research/src/patchDrawer.cpp

src/CMakeFiles/meshedit.dir/patchDrawer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/meshedit.dir/patchDrawer.cpp.i"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/bryce/Documents/research/src/patchDrawer.cpp > CMakeFiles/meshedit.dir/patchDrawer.cpp.i

src/CMakeFiles/meshedit.dir/patchDrawer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/meshedit.dir/patchDrawer.cpp.s"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/bryce/Documents/research/src/patchDrawer.cpp -o CMakeFiles/meshedit.dir/patchDrawer.cpp.s

src/CMakeFiles/meshedit.dir/patchDrawer.cpp.o.requires:
.PHONY : src/CMakeFiles/meshedit.dir/patchDrawer.cpp.o.requires

src/CMakeFiles/meshedit.dir/patchDrawer.cpp.o.provides: src/CMakeFiles/meshedit.dir/patchDrawer.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/meshedit.dir/build.make src/CMakeFiles/meshedit.dir/patchDrawer.cpp.o.provides.build
.PHONY : src/CMakeFiles/meshedit.dir/patchDrawer.cpp.o.provides

src/CMakeFiles/meshedit.dir/patchDrawer.cpp.o.provides.build: src/CMakeFiles/meshedit.dir/patchDrawer.cpp.o

src/CMakeFiles/meshedit.dir/curveTracer.cpp.o: src/CMakeFiles/meshedit.dir/flags.make
src/CMakeFiles/meshedit.dir/curveTracer.cpp.o: ../src/curveTracer.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/bryce/Documents/research/build/CMakeFiles $(CMAKE_PROGRESS_13)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/CMakeFiles/meshedit.dir/curveTracer.cpp.o"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/meshedit.dir/curveTracer.cpp.o -c /home/bryce/Documents/research/src/curveTracer.cpp

src/CMakeFiles/meshedit.dir/curveTracer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/meshedit.dir/curveTracer.cpp.i"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/bryce/Documents/research/src/curveTracer.cpp > CMakeFiles/meshedit.dir/curveTracer.cpp.i

src/CMakeFiles/meshedit.dir/curveTracer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/meshedit.dir/curveTracer.cpp.s"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/bryce/Documents/research/src/curveTracer.cpp -o CMakeFiles/meshedit.dir/curveTracer.cpp.s

src/CMakeFiles/meshedit.dir/curveTracer.cpp.o.requires:
.PHONY : src/CMakeFiles/meshedit.dir/curveTracer.cpp.o.requires

src/CMakeFiles/meshedit.dir/curveTracer.cpp.o.provides: src/CMakeFiles/meshedit.dir/curveTracer.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/meshedit.dir/build.make src/CMakeFiles/meshedit.dir/curveTracer.cpp.o.provides.build
.PHONY : src/CMakeFiles/meshedit.dir/curveTracer.cpp.o.provides

src/CMakeFiles/meshedit.dir/curveTracer.cpp.o.provides.build: src/CMakeFiles/meshedit.dir/curveTracer.cpp.o

src/CMakeFiles/meshedit.dir/meshEdit.cpp.o: src/CMakeFiles/meshedit.dir/flags.make
src/CMakeFiles/meshedit.dir/meshEdit.cpp.o: ../src/meshEdit.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/bryce/Documents/research/build/CMakeFiles $(CMAKE_PROGRESS_14)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/CMakeFiles/meshedit.dir/meshEdit.cpp.o"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/meshedit.dir/meshEdit.cpp.o -c /home/bryce/Documents/research/src/meshEdit.cpp

src/CMakeFiles/meshedit.dir/meshEdit.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/meshedit.dir/meshEdit.cpp.i"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/bryce/Documents/research/src/meshEdit.cpp > CMakeFiles/meshedit.dir/meshEdit.cpp.i

src/CMakeFiles/meshedit.dir/meshEdit.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/meshedit.dir/meshEdit.cpp.s"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/bryce/Documents/research/src/meshEdit.cpp -o CMakeFiles/meshedit.dir/meshEdit.cpp.s

src/CMakeFiles/meshedit.dir/meshEdit.cpp.o.requires:
.PHONY : src/CMakeFiles/meshedit.dir/meshEdit.cpp.o.requires

src/CMakeFiles/meshedit.dir/meshEdit.cpp.o.provides: src/CMakeFiles/meshedit.dir/meshEdit.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/meshedit.dir/build.make src/CMakeFiles/meshedit.dir/meshEdit.cpp.o.provides.build
.PHONY : src/CMakeFiles/meshedit.dir/meshEdit.cpp.o.provides

src/CMakeFiles/meshedit.dir/meshEdit.cpp.o.provides.build: src/CMakeFiles/meshedit.dir/meshEdit.cpp.o

src/CMakeFiles/meshedit.dir/png.cpp.o: src/CMakeFiles/meshedit.dir/flags.make
src/CMakeFiles/meshedit.dir/png.cpp.o: ../src/png.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/bryce/Documents/research/build/CMakeFiles $(CMAKE_PROGRESS_15)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/CMakeFiles/meshedit.dir/png.cpp.o"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/meshedit.dir/png.cpp.o -c /home/bryce/Documents/research/src/png.cpp

src/CMakeFiles/meshedit.dir/png.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/meshedit.dir/png.cpp.i"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/bryce/Documents/research/src/png.cpp > CMakeFiles/meshedit.dir/png.cpp.i

src/CMakeFiles/meshedit.dir/png.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/meshedit.dir/png.cpp.s"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/bryce/Documents/research/src/png.cpp -o CMakeFiles/meshedit.dir/png.cpp.s

src/CMakeFiles/meshedit.dir/png.cpp.o.requires:
.PHONY : src/CMakeFiles/meshedit.dir/png.cpp.o.requires

src/CMakeFiles/meshedit.dir/png.cpp.o.provides: src/CMakeFiles/meshedit.dir/png.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/meshedit.dir/build.make src/CMakeFiles/meshedit.dir/png.cpp.o.provides.build
.PHONY : src/CMakeFiles/meshedit.dir/png.cpp.o.provides

src/CMakeFiles/meshedit.dir/png.cpp.o.provides.build: src/CMakeFiles/meshedit.dir/png.cpp.o

src/CMakeFiles/meshedit.dir/main.cpp.o: src/CMakeFiles/meshedit.dir/flags.make
src/CMakeFiles/meshedit.dir/main.cpp.o: ../src/main.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/bryce/Documents/research/build/CMakeFiles $(CMAKE_PROGRESS_16)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/CMakeFiles/meshedit.dir/main.cpp.o"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/meshedit.dir/main.cpp.o -c /home/bryce/Documents/research/src/main.cpp

src/CMakeFiles/meshedit.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/meshedit.dir/main.cpp.i"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/bryce/Documents/research/src/main.cpp > CMakeFiles/meshedit.dir/main.cpp.i

src/CMakeFiles/meshedit.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/meshedit.dir/main.cpp.s"
	cd /home/bryce/Documents/research/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/bryce/Documents/research/src/main.cpp -o CMakeFiles/meshedit.dir/main.cpp.s

src/CMakeFiles/meshedit.dir/main.cpp.o.requires:
.PHONY : src/CMakeFiles/meshedit.dir/main.cpp.o.requires

src/CMakeFiles/meshedit.dir/main.cpp.o.provides: src/CMakeFiles/meshedit.dir/main.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/meshedit.dir/build.make src/CMakeFiles/meshedit.dir/main.cpp.o.provides.build
.PHONY : src/CMakeFiles/meshedit.dir/main.cpp.o.provides

src/CMakeFiles/meshedit.dir/main.cpp.o.provides.build: src/CMakeFiles/meshedit.dir/main.cpp.o

# Object files for target meshedit
meshedit_OBJECTS = \
"CMakeFiles/meshedit.dir/UnionFind/UF_ADT.cpp.o" \
"CMakeFiles/meshedit.dir/UnionFind/UF_Serial.cpp.o" \
"CMakeFiles/meshedit.dir/scene.cpp.o" \
"CMakeFiles/meshedit.dir/camera.cpp.o" \
"CMakeFiles/meshedit.dir/light.cpp.o" \
"CMakeFiles/meshedit.dir/mesh.cpp.o" \
"CMakeFiles/meshedit.dir/material.cpp.o" \
"CMakeFiles/meshedit.dir/texture.cpp.o" \
"CMakeFiles/meshedit.dir/collada.cpp.o" \
"CMakeFiles/meshedit.dir/halfEdgeMesh.cpp.o" \
"CMakeFiles/meshedit.dir/meshResampler.cpp.o" \
"CMakeFiles/meshedit.dir/patchDrawer.cpp.o" \
"CMakeFiles/meshedit.dir/curveTracer.cpp.o" \
"CMakeFiles/meshedit.dir/meshEdit.cpp.o" \
"CMakeFiles/meshedit.dir/png.cpp.o" \
"CMakeFiles/meshedit.dir/main.cpp.o"

# External object files for target meshedit
meshedit_EXTERNAL_OBJECTS =

meshedit: src/CMakeFiles/meshedit.dir/UnionFind/UF_ADT.cpp.o
meshedit: src/CMakeFiles/meshedit.dir/UnionFind/UF_Serial.cpp.o
meshedit: src/CMakeFiles/meshedit.dir/scene.cpp.o
meshedit: src/CMakeFiles/meshedit.dir/camera.cpp.o
meshedit: src/CMakeFiles/meshedit.dir/light.cpp.o
meshedit: src/CMakeFiles/meshedit.dir/mesh.cpp.o
meshedit: src/CMakeFiles/meshedit.dir/material.cpp.o
meshedit: src/CMakeFiles/meshedit.dir/texture.cpp.o
meshedit: src/CMakeFiles/meshedit.dir/collada.cpp.o
meshedit: src/CMakeFiles/meshedit.dir/halfEdgeMesh.cpp.o
meshedit: src/CMakeFiles/meshedit.dir/meshResampler.cpp.o
meshedit: src/CMakeFiles/meshedit.dir/patchDrawer.cpp.o
meshedit: src/CMakeFiles/meshedit.dir/curveTracer.cpp.o
meshedit: src/CMakeFiles/meshedit.dir/meshEdit.cpp.o
meshedit: src/CMakeFiles/meshedit.dir/png.cpp.o
meshedit: src/CMakeFiles/meshedit.dir/main.cpp.o
meshedit: src/CMakeFiles/meshedit.dir/build.make
meshedit: CMU462/src/libCMU462.a
meshedit: CMU462/deps/glew/libglew.a
meshedit: CMU462/deps/glfw/src/libglfw3.a
meshedit: /usr/lib/x86_64-linux-gnu/libX11.so
meshedit: /usr/lib/x86_64-linux-gnu/libXrandr.so
meshedit: /usr/lib/x86_64-linux-gnu/libXinerama.so
meshedit: /usr/lib/x86_64-linux-gnu/libXi.so
meshedit: /usr/lib/x86_64-linux-gnu/libXxf86vm.so
meshedit: /usr/lib/x86_64-linux-gnu/librt.so
meshedit: /usr/lib/x86_64-linux-gnu/libm.so
meshedit: /usr/lib/x86_64-linux-gnu/libXcursor.so
meshedit: /usr/lib/x86_64-linux-gnu/libGL.so
meshedit: /usr/lib/x86_64-linux-gnu/libGLU.so
meshedit: /usr/lib/x86_64-linux-gnu/libGL.so
meshedit: /usr/lib/x86_64-linux-gnu/libfreetype.so
meshedit: /usr/lib/x86_64-linux-gnu/libGLU.so
meshedit: /usr/lib/x86_64-linux-gnu/libfreetype.so
meshedit: src/CMakeFiles/meshedit.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../meshedit"
	cd /home/bryce/Documents/research/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/meshedit.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/meshedit.dir/build: meshedit
.PHONY : src/CMakeFiles/meshedit.dir/build

src/CMakeFiles/meshedit.dir/requires: src/CMakeFiles/meshedit.dir/UnionFind/UF_ADT.cpp.o.requires
src/CMakeFiles/meshedit.dir/requires: src/CMakeFiles/meshedit.dir/UnionFind/UF_Serial.cpp.o.requires
src/CMakeFiles/meshedit.dir/requires: src/CMakeFiles/meshedit.dir/scene.cpp.o.requires
src/CMakeFiles/meshedit.dir/requires: src/CMakeFiles/meshedit.dir/camera.cpp.o.requires
src/CMakeFiles/meshedit.dir/requires: src/CMakeFiles/meshedit.dir/light.cpp.o.requires
src/CMakeFiles/meshedit.dir/requires: src/CMakeFiles/meshedit.dir/mesh.cpp.o.requires
src/CMakeFiles/meshedit.dir/requires: src/CMakeFiles/meshedit.dir/material.cpp.o.requires
src/CMakeFiles/meshedit.dir/requires: src/CMakeFiles/meshedit.dir/texture.cpp.o.requires
src/CMakeFiles/meshedit.dir/requires: src/CMakeFiles/meshedit.dir/collada.cpp.o.requires
src/CMakeFiles/meshedit.dir/requires: src/CMakeFiles/meshedit.dir/halfEdgeMesh.cpp.o.requires
src/CMakeFiles/meshedit.dir/requires: src/CMakeFiles/meshedit.dir/meshResampler.cpp.o.requires
src/CMakeFiles/meshedit.dir/requires: src/CMakeFiles/meshedit.dir/patchDrawer.cpp.o.requires
src/CMakeFiles/meshedit.dir/requires: src/CMakeFiles/meshedit.dir/curveTracer.cpp.o.requires
src/CMakeFiles/meshedit.dir/requires: src/CMakeFiles/meshedit.dir/meshEdit.cpp.o.requires
src/CMakeFiles/meshedit.dir/requires: src/CMakeFiles/meshedit.dir/png.cpp.o.requires
src/CMakeFiles/meshedit.dir/requires: src/CMakeFiles/meshedit.dir/main.cpp.o.requires
.PHONY : src/CMakeFiles/meshedit.dir/requires

src/CMakeFiles/meshedit.dir/clean:
	cd /home/bryce/Documents/research/build/src && $(CMAKE_COMMAND) -P CMakeFiles/meshedit.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/meshedit.dir/clean

src/CMakeFiles/meshedit.dir/depend:
	cd /home/bryce/Documents/research/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/bryce/Documents/research /home/bryce/Documents/research/src /home/bryce/Documents/research/build /home/bryce/Documents/research/build/src /home/bryce/Documents/research/build/src/CMakeFiles/meshedit.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/meshedit.dir/depend
