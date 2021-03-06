cmake_minimum_required(VERSION 2.8)

# Collada viewer source
set(COLLADA_VIEWER_SOURCE
    UnionFind/UF_ADT.cpp
    UnionFind/UF_Serial.cpp
    scene.cpp
    camera.cpp
    light.cpp
    mesh.cpp
    material.cpp
    texture.cpp
    collada.cpp
    style.cpp
    #critical_points.cpp
    svg_exporter.cpp
    halfEdgeMesh.cpp
    PatchFunctions.cpp
    pointCurve.cpp
    meshResampler.cpp
    bezierPatch.cpp
    patchDrawer.cpp
    curveTracer.cpp
    meshEdit.cpp
    png.cpp
    testing.cpp
    main.cpp
)

# Collada viewer header
set(COLLADA_VIEWER_HEADER
    UnionFind/UF_ADT.h
    UnionFind/UF_Serial.h
    scene.h
    camera.h
    light.h
    mesh.h
    material.h
    texture.h
    collada.h
    style.h
    svg_exporter.h
    critical_points.h
    PatchFunctions.h
    halfEdgeMesh.h
    pointCurve.h
    bezierPatch.h
    patchDrawer.h
    curveTracer.h
    meshResampler.h
    meshEdit.h
    png.h
    testing.h
)


#-------------------------------------------------------------------------------
# Set include directories
#-------------------------------------------------------------------------------
include_directories(
  ${GLEW_INCLUDE_DIRS}
  ${GLFW_INCLUDE_DIRS}
  ${CMU462_INCLUDE_DIRS}
  ${FREETYPE_INCLUDE_DIRS}
)

#-------------------------------------------------------------------------------
# Set link directories
#-------------------------------------------------------------------------------
link_directories(
  ${GLEW_LIBRARY_DIRS}
  ${GLFW_LIBRARY_DIRS}
  ${CMU462_LIBRARY_DIRS}
  ${FREETYPE_LIBRARY_DIRS}
)

#-------------------------------------------------------------------------------
# Add executable
#-------------------------------------------------------------------------------
add_executable( meshedit
    ${COLLADA_VIEWER_SOURCE}
    ${COLLADA_VIEWER_HEADER}
)

target_link_libraries( meshedit
    CMU462 ${CMU462_LIBRARIES}
    glew ${GLEW_LIBRARIES}
    glfw ${GLFW_LIBRARIES}
    ${OPENGL_LIBRARIES}
    ${FREETYPE_LIBRARIES}
)

#-------------------------------------------------------------------------------
# Platform-specific configurations for target
#-------------------------------------------------------------------------------
if(APPLE)
  set_property( TARGET meshedit APPEND_STRING PROPERTY COMPILE_FLAGS
                "-Wno-deprecated-declarations -Wno-c++11-extensions")
endif(APPLE)


# Put executable in build directory root
set(EXECUTABLE_OUTPUT_PATH ..)

# Install to project root
install(TARGETS meshedit DESTINATION ${ColladaViewer_SOURCE_DIR})
