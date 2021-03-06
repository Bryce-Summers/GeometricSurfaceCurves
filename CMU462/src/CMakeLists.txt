cmake_minimum_required(VERSION 2.8)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# CMU462 library source files
set(CMU462_SOURCE
    polynomial.cpp
    polynomialVector3D.cpp
    vector2D.cpp
    vector3D.cpp
    vector4D.cpp
    matrix3x3.cpp
    matrix4x4.cpp
    quaternion.cpp
    complex.cpp
    color.cpp
    osdtext.cpp
    osdfont.c
    viewer.cpp
    base64.cpp
    tinyxml2.cpp
)

# CMU462 library header files
set(CMU462_HEADER
    CMU462.h
    polynomial.h
    polynomialVector3D.h
    vector2D.h
    vector3D.h
    vector4D.h
    matrix3x3.h
    matrix4x4.h
    quaternion.h
    complex.h
    color.h
    osdtext.h
    viewer.h
    base64.h
    tinyxml2.h
    renderer.h
)

#-------------------------------------------------------------------------------
# Building static library (always)
#-------------------------------------------------------------------------------
add_library(CMU462 STATIC ${CMU462_SOURCE} ${CMU462_HEADER})

target_link_libraries(
  CMU462
  ${GLEW_LIBRARIES}
  ${GLFW_LIBRARIES}
  ${OPENGL_LIBRARIES}
  ${FREETYPE_LIBRARIES}
)

#-------------------------------------------------------------------------------
# Building shared library
#-------------------------------------------------------------------------------
if(CMU462_BUILD_SHARED)
  add_library(CMU462_SHARED SHARED ${CMU462_SOURCE} ${CMU462_HEADER})
  target_link_libraries(
    CMU462_SHARED
    ${GLEW_LIBRARIES}
    ${GLFW_LIBRARIES}
    ${OPENGL_LIBRARIES}
    ${FREETYPE_LIBRARIES}
  )
endif()

#-------------------------------------------------------------------------------
# Platform-specific target configurations
#-------------------------------------------------------------------------------

###################
# Building on OSX #
###################
if (APPLE)

  # Clang options
  target_compile_options(CMU462 PRIVATE -Wno-constant-conversion)

  # Framework dependencies
  target_link_libraries(CMU462
    ${COCOA_LIBRARIES}
    ${IOKIT_LIBRARIES}
    ${COREVIDEO_LIBRARIES})

  # Output name
  if(CMU462_BUILD_DEBUG)
    set_target_properties(CMU462 PROPERTIES OUTPUT_NAME CMU462_osx_d)
    if(CMU462_BUILD_SHARED)
      set_target_properties(CMU462_SHARED PROPERTIES OUTPUT_NAME CMU462_osx_d)
    endif()
  else(CMU462_BUILD_DEBUG)
    set_target_properties(CMU462 PROPERTIES OUTPUT_NAME CMU462_osx)
    if(CMU462_BUILD_SHARED)
      set_target_properties(CMU462_SHARED PROPERTIES OUTPUT_NAME CMU462_osx)
    endif()
  endif(CMU462_BUILD_DEBUG)
  
endif(APPLE)

##################
# Build on Linux #
##################
set(LINUX UNIX AND NOT APPLE)
if(LINUX)
  if(CMU462_BUILD_DEBUG)
    set_target_properties(CMU462 PROPERTIES OUTPUT_NAME CMU462_d)
    if(CMU462_BUILD_SHARED)
      set_target_properties(CMU462_SHARED PROPERTIES OUTPUT_NAME CMU462_d)
    endif()
  else(CMU462_BUILD_DEBUG)
    set_target_properties(CMU462 PROPERTIES OUTPUT_NAME CMU462)
    if(CMU462_BUILD_SHARED)
      set_target_properties(CMU462_SHARED PROPERTIES OUTPUT_NAME CMU462)
    endif()
  endif(CMU462_BUILD_DEBUG)
endif(LINUX)

####################
# Build on Windows #
####################
if(WIN32)
  
  if(MSVC)
    if(CMU462_BUILD_DEBUG)
      set_target_properties(CMU462 PROPERTIES OUTPUT_NAME CMU462_d)
      if(CMU462_BUILD_SHARED)
        set_target_properties(CMU462_SHARED PROPERTIES OUTPUT_NAME CMU462_d)
      endif()
    else(CMU462_BUILD_DEBUG)
      set_target_properties(CMU462 PROPERTIES OUTPUT_NAME CMU462)
      if(CMU462_BUILD_SHARED)
        set_target_properties(CMU462_SHARED PROPERTIES OUTPUT_NAME CMU462)
      endif()
    endif(CMU462_BUILD_DEBUG)
  endif(MSVC)

  if(MINGW)
    if(CMU462_BUILD_DEBUG)
      set_target_properties(CMU462 PROPERTIES OUTPUT_NAME CMU462_d)
      if(CMU462_BUILD_SHARED)
        set_target_properties(CMU462_SHARED PROPERTIES OUTPUT_NAME CMU462_d)
      endif()
    else(CMU462_BUILD_DEBUG)
      set_target_properties(CMU462 PROPERTIES OUTPUT_NAME CMU462)
      if(CMU462_BUILD_SHARED)
        set_target_properties(CMU462_SHARED PROPERTIES OUTPUT_NAME CMU462)
      endif()
    endif(CMU462_BUILD_DEBUG)
  endif(MINGW)

endif(WIN32)

#-------------------------------------------------------------------------------
# Installation configurations
#-------------------------------------------------------------------------------
install(TARGETS CMU462 DESTINATION lib)
if(CMU462_BUILD_SHARED)
  install(TARGETS CMU462_SHARED DESTINATION lib)
endif()

# install headers
install(FILES ${CMU462_HEADER} DESTINATION include/CMU462)