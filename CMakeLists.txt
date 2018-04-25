###################################################################
# Package: SuperRest
###################################################################

atlas_subdir( SuperRest )

set( extra_deps )
set( extra_libs )

atlas_depends_on_subdirs(
    PUBLIC
    Superflow
    SusyNtuple
    ${extra_deps}
)

find_package( ROOT COMPONENTS Gpad Graf Graf3d Core Tree MathCore Hist RIO )

atlas_add_library( SuperRestLib
    SuperRest/*.h Root/*.cxx
    PUBLIC_HEADERS SuperRest
    PRIVATE_INCLUDE_DIRS ${ROOT_INCLUDE_DIRS} $ENV{TestArea}/../RestFrames/include
    LINK_LIBRARIES SuperflowLib SusyNtupleLib ${ROOT_LIBRARIES}
)

set( extra_libs )

# define a function to call the add_executable macro on our executables
function(SuperRestExec filename)
    set(execname)
    get_filename_component(execname ${filename} NAME_WE)
    atlas_add_executable( ${execname} "util/${execname}.cxx"
        INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
        LINK_LIBRARIES ${ROOT_LIBRARIES} SuperRestLib ${extra_libs}
    )
endfunction(SuperRestExec)

# build everything in util/ as an executable
file(GLOB files "util/*.cxx")
foreach(file ${files})
    SuperRestExec(${file})
endforeach()

atlas_install_data( data/* )
