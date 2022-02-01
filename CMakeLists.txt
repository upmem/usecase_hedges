cmake_minimum_required(VERSION 3.10)
project(hedges C CXX)
include(dpu_recipes/dpu.cmake)

set(NR_DPUS 1)
set(NR_TASKLETS 2)
set(HOSTDPU_XFER_BUFFER_SIZE_LOG2_BYTE 24)
set(DPU_STACK_SIZE_BYTE 8000)
set(ENCODER_STACK_MAXSEQ_DPU_NBYTES 300)
set(MAX_DECODED_VBITS 400)
set(BUDDY_ALLOCATOR_SIZE_BYTE 8000)
set(ENCODER_NPACKET_PER_DPU 2)
set(ENABLE_TEST_DPU_HEAP  0)
set(ENABLE_HEAP_HOST_DEBUGGING 0)
set(ENABLE_DPU_PRINT 0)

add_dpu_project(
	dpuencoder
	""
	""
	"hedges/dpu_hedges_encoder.c;"
	"."
	"hedges/"
	${NR_TASKLETS}
	${NR_DPUS}
	${DPU_STACK_SIZE_BYTE}
	${HOSTDPU_XFER_BUFFER_SIZE_LOG2_BYTE}
	"-DENABLE_DPU_PRINT=${ENABLE_DPU_PRINT} \
	 -DENABLE_TEST_DPU_HEAP=${ENABLE_TEST_DPU_HEAP} \
	 -DENABLE_HEAP_HOST_DEBUGGING=${ENABLE_HEAP_HOST_DEBUGGING} \
	 -DENCODER_STACK_MAXSEQ_DPU_NBYTES=${ENCODER_STACK_MAXSEQ_DPU_NBYTES} \
	 -DBUDDY_ALLOCATOR_SIZE_BYTE=${BUDDY_ALLOCATOR_SIZE_BYTE}"
)

add_dpu_project(
 	dpudecoder
 	""
 	""
 	"hedges/dpu_hedges_decoder.c;"
 	"dpu_recipes/heap;.;hedges/"
 	""
 	${NR_TASKLETS}
 	${NR_DPUS}
 	${DPU_STACK_SIZE_BYTE}
 	${HOSTDPU_XFER_BUFFER_SIZE_LOG2_BYTE}
 	"-DENABLE_DPU_PRINT=${ENABLE_DPU_PRINT} \
	 -DENABLE_TEST_DPU_HEAP=${ENABLE_TEST_DPU_HEAP} \
	 -DENABLE_HEAP_HOST_DEBUGGING=${ENABLE_HEAP_HOST_DEBUGGING} \
	 -DMAX_DECODED_VBITS=${MAX_DECODED_VBITS} \
	 -DENCODER_STACK_MAXSEQ_DPU_NBYTES=${ENCODER_STACK_MAXSEQ_DPU_NBYTES}\
     -DBUDDY_ALLOCATOR_SIZE_BYTE=${BUDDY_ALLOCATOR_SIZE_BYTE}"
)

set(CMAKE_CXX_FLAGS "-O3   -g3")
if (NOT DEFINED UPMEM_HOME)
        if ( "$ENV{UPMEM_HOME}" STREQUAL "")
                set(UPMEM_HOME "/usr")
        else ()
                set(UPMEM_HOME $ENV{UPMEM_HOME})
        endif ()
endif ()
include(${UPMEM_HOME}/share/upmem/cmake/include/host/DpuHost.cmake)
get_filename_component(HEDGES ${CMAKE_CURRENT_SOURCE_DIR}/hedges ABSOLUTE)
get_filename_component(SCHIFRA ${CMAKE_CURRENT_SOURCE_DIR}/hedges/schifra ABSOLUTE)
get_filename_component(NUMPY_1_13_INCLUDE /usr/local/lib/python2.7/dist-packages/numpy-1.13.3-py2.7-linux-x86_64.egg/numpy/core/include ABSOLUTE)
get_filename_component(PYTHON_27_INCLUDE /usr/local/include/python2.7/ ABSOLUTE)
add_library(NRpyDNAcode SHARED hedges/NRpyDNAcode.cpp)
target_include_directories(NRpyDNAcode PUBLIC  hedges
									           hedges/common
									           "${DPU_HOST_INCLUDE_DIRECTORIES}"
									           .        
									           ${PYTHON_27_INCLUDE}
									           ${NUMPY_1_13_INCLUDE}
)
target_link_libraries(NRpyDNAcode ${DPU_HOST_LIBRARIES})
target_compile_definitions( NRpyDNAcode PUBLIC  -DMAX_DECODED_VBITS=${MAX_DECODED_VBITS}                                                                      )             
target_compile_definitions( NRpyDNAcode PUBLIC  -DENABLE_DPU_PRINT=${ENABLE_DPU_PRINT}) 
target_compile_definitions( NRpyDNAcode PUBLIC  -DENABLE_TEST_DPU_HEAP=${ENABLE_TEST_DPU_HEAP})
target_compile_definitions( NRpyDNAcode PUBLIC  -DENABLE_HEAP_HOST_DEBUGGING=${ENABLE_HEAP_HOST_DEBUGGING} )
target_compile_definitions( NRpyDNAcode PUBLIC  -DENCODER_NPACKET_PER_DPU=${ENCODER_NPACKET_PER_DPU}  )
target_compile_definitions( NRpyDNAcode PUBLIC  -DNR_DPUS=${NR_DPUS} )
target_compile_definitions( NRpyDNAcode PUBLIC  -DNR_TASKLETS=${NR_TASKLETS})
target_compile_definitions( NRpyDNAcode PUBLIC  -DDPU_DECODE_BINARY="${dpudecoder-path}")
target_compile_definitions( NRpyDNAcode PUBLIC  -DDPU_ENCODE_BINARY="${dpuencoder-path}")     
set_target_properties(NRpyDNAcode PROPERTIES PREFIX "")

add_library(NRpyRS SHARED ../hedges/NRpyRS.cpp)
target_include_directories(NRpyRS PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
 										  ${PYTHON_27_INCLUDE} 
 										  ${NUMPY_1_13_INCLUDE}
)
set_target_properties(NRpyRS PROPERTIES PREFIX "")