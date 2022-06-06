cmake_minimum_required(VERSION 3.10)

enable_testing()

get_filename_component(TEMPLATE_DPU_PROJECT_RELATIVE_PATH ${CMAKE_CURRENT_LIST_DIR}/.cmakeDpuProject ABSOLUTE)


function(add_dpu_project name HOST_SRC_LIST
                              HOST_INCLUDE_DIRS
                              DPU_SRC_LIST
                              DPU_INCLUDE_DIRS
                              COMMON_INCLUDE_DIRS
                              NR_TASKLETS
                              NR_DPUS
                              DPU_STACK_SIZE_BYTE
                              HOSTDPU_XFER_BUFFER_SIZE_LOG2_BYTE)

set(SOLVED_HOST_SRC_LIST "")
set(SOLVED_HOST_INCLUDE_DIRS "")
set(SOLVED_DPU_SRC_LIST "")
set(SOLVED_DPU_INCLUDE_DIRS "")
set(SOLVED_COMMON_INCLUDE_DIRS "")

foreach (SRC ${HOST_SRC_LIST})
     get_filename_component(SOLVED_SRC ${SRC} ABSOLUTE)
     list(APPEND SOLVED_HOST_SRC_LIST ${SOLVED_SRC})
endforeach()
foreach (SRC ${DPU_SRC_LIST})
     get_filename_component(SOLVED_SRC ${SRC} ABSOLUTE)
     list(APPEND SOLVED_DPU_SRC_LIST ${SOLVED_SRC})
endforeach()
foreach (SRC ${DPU_INCLUDE_DIRS})
     get_filename_component(SOLVED_SRC ${SRC} ABSOLUTE)
     list(APPEND SOLVED_DPU_INCLUDE_DIRS ${SOLVED_SRC})
endforeach()
foreach (SRC ${HOST_INCLUDE_DIRS})
     get_filename_component(SOLVED_SRC ${SRC} ABSOLUTE)
     list(APPEND SOLVED_HOST_INCLUDE_DIRS ${SOLVED_SRC})
endforeach()
foreach (SRC ${COMMON_INCLUDE_DIRS})
     get_filename_component(SOLVED_SRC ${SRC} ABSOLUTE)
     list(APPEND SOLVED_COMMON_INCLUDE_DIRS ${SOLVED_SRC})
endforeach()

set(DPU_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/${name}-build-dpu)
set(DPU_BINARY_NAME ${name}-dpu)
if (NOT DEFINED UPMEM_HOME)
        if ( "$ENV{UPMEM_HOME}" STREQUAL "")
                set(UPMEM_HOME "/usr")
        else ()
                set(UPMEM_HOME $ENV{UPMEM_HOME})
        endif ()
endif ()
include(${UPMEM_HOME}/share/upmem/cmake/include/host/DpuHost.cmake)
link_directories("${DPU_HOST_LINK_DIRECTORIES}")

include(ExternalProject)

string(REPLACE ";" "|" SOLVED_DPU_SRC_LIST    "${SOLVED_DPU_SRC_LIST}")
string(REPLACE ";" "|" SOLVED_DPU_INCLUDE_DIRS    "${SOLVED_DPU_INCLUDE_DIRS}")
string(REPLACE ";" "|" SOLVED_COMMON_INCLUDE_DIRS    "${SOLVED_COMMON_INCLUDE_DIRS}")

# pass ARGN to dpu project
set(DPU_ARGN ${ARGN})
ExternalProject_Add(
        ${DPU_BINARY_NAME}
	BINARY_DIR ${DPU_BINARY_DIR}
	SOURCE_DIR ${TEMPLATE_DPU_PROJECT_RELATIVE_PATH}
        CMAKE_ARGS -DNAME=${DPU_BINARY_NAME}
                   -DDPU_SRC_LIST=${SOLVED_DPU_SRC_LIST}
                   -DDPU_INCLUDE_DIRS=${SOLVED_DPU_INCLUDE_DIRS}
                   -DCOMMON_INCLUDE_DIRS=${SOLVED_COMMON_INCLUDE_DIRS}
                   -DCMAKE_TOOLCHAIN_FILE=${UPMEM_HOME}/share/upmem/cmake/dpu.cmake
                   -DHOSTDPU_XFER_BUFFER_SIZE_LOG2_BYTE=${HOSTDPU_XFER_BUFFER_SIZE_LOG2_BYTE}
                   -DUPMEM_HOME=${UPMEM_HOME}
                   -DNR_TASKLETS=${NR_TASKLETS}
                   -DDPU_STACK_SIZE_BYTE=${DPU_STACK_SIZE_BYTE}
		    -DDPU_ARGN=${DPU_ARGN}
        BUILD_ALWAYS TRUE
        INSTALL_COMMAND ""
)

set(${name}-path  "${DPU_BINARY_DIR}/${DPU_BINARY_NAME}" PARENT_SCOPE)

if(SOLVED_HOST_SRC_LIST)
        add_executable(${name}-host ${SOLVED_HOST_SRC_LIST})
        target_compile_definitions(${name}-host PUBLIC -DHOSTDPU_XFER_BUFFER_SIZE_LOG2_BYTE=${HOSTDPU_XFER_BUFFER_SIZE_LOG2_BYTE})
        target_compile_definitions(${name}-host PUBLIC -DMULTITHREAD) 
        target_compile_definitions(${name}-host PUBLIC -DNR_TASKLETS=${NR_TASKLETS}) 
        target_compile_definitions(${name}-host PUBLIC -DNR_DPUS=${NR_DPUS}) 
        target_compile_definitions(${name}-host PUBLIC -DDPU_BINARY="${DPU_BINARY_DIR}/${DPU_BINARY_NAME}")
        target_include_directories(${name}-host PUBLIC ${DPU_HOST_INCLUDE_DIRECTORIES}
                                                       ${SOLVED_HOST_INCLUDE_DIRS}
                                                       ${SOLVED_COMMON_INCLUDE_DIRS})
        target_link_libraries(${name}-host ${DPU_HOST_LIBRARIES})
        add_dependencies(${name}-host ${name}-dpu)
endif()

endfunction()
