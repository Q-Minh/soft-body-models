cmake_minimum_required(VERSION 3.14)
project(sbm VERSION 0.0.1 LANGUAGES CXX)

set_property(GLOBAL PROPERTY USE_FOLDER ON)

include(FetchContent)

set(LIBIGL_USE_STATIC_LIBRARY     OFF CACHE STRING   "Use libigl as static library" )
set(LIBIGL_WITH_CGAL              ON  CACHE INTERNAL "Use CGAL"                     )
set(LIBIGL_WITH_COMISO            OFF CACHE INTERNAL "Use CoMiso"                   )
set(LIBIGL_WITH_CORK              OFF CACHE INTERNAL "Use Cork"                     )
set(LIBIGL_WITH_EMBREE            OFF CACHE INTERNAL "Use Embree"                   )
set(LIBIGL_WITH_MATLAB            OFF CACHE INTERNAL "Use Matlab"                   )
set(LIBIGL_WITH_MOSEK             OFF CACHE INTERNAL "Use MOSEK"                    )
set(LIBIGL_WITH_PNG               OFF CACHE INTERNAL "Use PNG"                      )
set(LIBIGL_WITH_TETGEN            ON  CACHE INTERNAL "Use Tetgen"                   )
set(LIBIGL_WITH_TRIANGLE          OFF CACHE INTERNAL "Use Triangle"                 )
set(LIBIGL_WITH_PREDICATES        OFF CACHE INTERNAL "Use exact predicates"         )
set(LIBIGL_WITH_XML               OFF CACHE INTERNAL "Use XML"                      )
set(LIBIGL_WITHOUT_COPYLEFT       OFF CACHE INTERNAL "Disable Copyleft libraries"   )
set(LIBIGL_WITH_OPENGL            OFF CACHE INTERNAL "download libigl with opengl"  )
set(LIBIGL_WITH_OPENGL_GLFW       OFF CACHE INTERNAL "download libigl with glfw"    )
set(LIBIGL_WITH_OPENGL_GLFW_IMGUI OFF CACHE INTERNAL "download libigl with imgui"   )
FetchContent_Declare(
  _libigl
  GIT_REPOSITORY https://github.com/libigl/libigl
  GIT_TAG        v2.2.0
)
FetchContent_MakeAvailable(_libigl)

function (create_program _exec_name _folder _sources)
    add_executable(${_exec_name})
    set_target_properties(${_exec_name} PROPERTIES FOLDER ${_exec_name})
    target_compile_features(${_exec_name} PRIVATE cxx_std_17)

    target_sources(${_exec_name}
    PRIVATE
        # source files
        ${_sources}
    )

    target_link_libraries(${_exec_name} 
    PRIVATE 
        igl::core 
        igl::cgal 
        igl::tetgen
    )    
endfunction()

create_program("surface-from-tet" "surface-from-tet" "surface_from_tet.cpp")
create_program("check-mesh" "check-mesh" "check_mesh.cpp")
create_program(
    "convert-triangle-mesh-format" 
    "convert-triangle-mesh-format" 
    "convert_triangle_mesh_format.cpp"
)