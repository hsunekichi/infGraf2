cmake_minimum_required (VERSION 3.2)
project(nori)

add_subdirectory(ext ext_build)

#set(CMAKE_CXX_COMPILER "nvc++")

include_directories(
  # Nori include files
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  # tinyformat string formatting library
  ${TFM_INCLUDE_DIR}
  # Eigen linear algebra library
  SYSTEM ${EIGEN_INCLUDE_DIR}
  # OpenEXR high dynamic range bitmap library
  SYSTEM ${OPENEXR_INCLUDE_DIRS}
  # Intel Thread Building Blocks
  SYSTEM ${TBB_INCLUDE_DIR}
  # Pseudorandom number generator
  ${PCG32_INCLUDE_DIR}
  # PugiXML parser
  ${PUGIXML_INCLUDE_DIR}
  # Helper functions for statistical hypothesis tests
  ${HYPOTHESIS_INCLUDE_DIR}
  # GLFW library for OpenGL context creation
  SYSTEM ${GLFW_INCLUDE_DIR}
  # GLEW library for accessing OpenGL functions
  SYSTEM ${GLEW_INCLUDE_DIR}
  # NanoVG drawing library
  SYSTEM ${NANOVG_INCLUDE_DIR}
  # NanoGUI user interface library
  SYSTEM ${NANOGUI_INCLUDE_DIR}
  SYSTEM ${NANOGUI_EXTRA_INCS}
  # Portable filesystem API
  SYSTEM ${FILESYSTEM_INCLUDE_DIR}
  # STB Image Write
  SYSTEM ${STB_IMAGE_WRITE_INCLUDE_DIR}
)

# The following lines build the main executable. If you add a source
# code file to Nori, be sure to include it in this list.
add_executable(nori

  # Header files
  include/nori/common.h
  include/nori/accel.h
  include/nori/bbox.h
  include/nori/bitmap.h
  include/nori/block.h
  include/nori/bsdf.h
  include/nori/camera.h
  include/nori/color.h
  include/nori/dpdf.h
  include/nori/frame.h
  include/nori/gui.h
  include/nori/integrator.h
  include/nori/emitter.h
  include/nori/mesh.h
  include/nori/object.h
  include/nori/parser.h
  include/nori/proplist.h
  include/nori/math.h
  include/nori/ray.h
  include/nori/rfilter.h
  include/nori/sampler.h
  include/nori/scene.h
  include/nori/texture.h
  include/nori/timer.h
  include/nori/transform.h
  include/nori/vector.h
  include/nori/warp.h
  include/nori/pathtracing.h
  include/nori/kdtree.h

  # Source code files
  src/common.cpp
  src/accel.cpp
  src/area.cpp
  src/bitmap.cpp
  src/block.cpp
  src/chi2test.cpp
  src/dielectric.cpp
  src/diffuse.cpp
  src/environment.cpp  
  src/gui.cpp
  src/independent.cpp
  src/parser.cpp
  src/main.cpp
  src/mesh.cpp
  src/math.cpp
  src/microfacet.cpp
  src/mirror.cpp
  src/obj.cpp
  src/object.cpp
  src/perspective.cpp
  src/proplist.cpp
  src/rfilter.cpp
  src/scene.cpp
  src/texture.cpp
  src/ttest.cpp
  src/warp.cpp
  src/normals.cpp
  src/simpleIntegrator.cpp
  src/AOintegrator.cpp
  src/whitted.cpp
  src/GGX.cpp
  src/area.cpp
  src/orenNayar.cpp
  src/stratified.cpp
  src/pointLight.cpp
  src/path_ems.cpp
  src/path_mats.cpp
  src/path_mis.cpp
  src/beckmann.cpp
  src/lightProbeMipMap.cpp
  src/pathtracing.cpp
  src/directWhitted.cpp
  src/depthIntegrator.cpp
  src/dipole.cpp
  src/subsurfacePbrt.cpp
  src/normalizedDiffusion.cpp
  src/indirectWhitted.cpp
  src/agregate.cpp
  src/skin.cpp
)

add_definitions(${NANOGUI_EXTRA_DEFS})

# The following lines build the warping test application
add_executable(warptest
  include/nori/warp.h
  include/nori/math.h
  src/warp.cpp
  src/warptest.cpp
  src/microfacet.cpp
  src/object.cpp
  src/proplist.cpp
  src/common.cpp
  src/math.cpp
)

if (WIN32)
  target_link_libraries(nori tbb_static pugixml IlmImf nanogui ${NANOGUI_EXTRA_LIBS} zlibstatic)
else()
  target_link_libraries(nori tbb_static pugixml IlmImf nanogui ${NANOGUI_EXTRA_LIBS})
endif()

target_link_libraries(warptest tbb_static nanogui ${NANOGUI_EXTRA_LIBS})

set_property(TARGET nori PROPERTY CXX_STANDARD 17)


if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endif()

# Force colored output for the ninja generator
if (CMAKE_GENERATOR STREQUAL "Ninja")
  if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fcolor-diagnostics")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fcolor-diagnostics")
  elseif (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fdiagnostics-color=always")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-color=always")
  endif()
endif()

# vim: set et ts=2 sw=2 ft=cmake nospell:
