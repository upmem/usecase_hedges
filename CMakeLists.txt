cmake_minimum_required(VERSION 3.10)
project(hedges C CXX)
include(dpu.cmake)

# number of DPU (WARNING : NRpyDNACode support only one DPU in the current state
set(NR_DPUS 1)
# number of TASKLETS for each DPU
set(NR_TASKLETS 20)
# STATIC DPU BUFFER SIZE for HOST/DPU Tensor xfer
set(HOSTDPU_XFER_BUFFER_SIZE_LOG2_BYTE 19)
# stack size for each DPU
set(DPU_STACK_SIZE_BYTE 2000)
# limit for sequence of decoded vbits
set(MAX_DECODED_VBITS 400)
# DPU buddy allocator space in byte
set(BUDDY_ALLOCATOR_SIZE_BYTE 8000)
# enable (compile time) diff of CPU model
# shoveled HEAP and DPU shoveled HEAP
set(ENABLE_HEAP_HOST_DEBUGGING 0)

# enable (compile time) printf() from dpus
set(ENABLE_DPU_PRINT 0)

# enable (compile time) DPU perf measurement
set(MESURE_PERF 1)

# enable (compile time) DPU DMA engine bandwidth measurement
set(MESURE_BW 1)

# number of decimal bits for hypothesis score quantification
set(DECODER_QUANT_FRAC_BITS 16)

# HEAP_MAX_ITEM gives the maximum number of hypothesis
set(HEAP_MAX_ITEM 50000)

# bits buffer, and uncoded bit buffer. Theses two buffer
# are declared in TASKLET stack, so we must limit its size
# NOTE that the effective size of both buffer are computed as
# uint64_t codetext_size = vbitlen(8 /*8 bit per packed bit*/ * packed_message_size) + RPRIMER;
# but we must use constants for C array SIZE, so we fix the limit
# at the smallest possible value right now.
# WARNING : 300 works only for coderate = 0.5
set(ENCODER_STACK_MAXSEQ_DPU_NBYTES 300)

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
	"-DMESURE_BW=${MESURE_BW} \
	 -DMESURE_PERF=${MESURE_PERF} \
	 -DENABLE_DPU_PRINT=${ENABLE_DPU_PRINT} \
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
 	"heap;.;hedges/"
 	""
 	${NR_TASKLETS}
 	${NR_DPUS}
 	${DPU_STACK_SIZE_BYTE}
 	${HOSTDPU_XFER_BUFFER_SIZE_LOG2_BYTE}
	"-DMESURE_BW=${MESURE_BW} \
	 -DMESURE_PERF=${MESURE_PERF} \
 	 -DENABLE_DPU_PRINT=${ENABLE_DPU_PRINT} \
	 -DENABLE_TEST_DPU_HEAP=${ENABLE_TEST_DPU_HEAP} \
	 -DENABLE_HEAP_HOST_DEBUGGING=${ENABLE_HEAP_HOST_DEBUGGING} \
	 -DMAX_DECODED_VBITS=${MAX_DECODED_VBITS} \
	 -DENCODER_STACK_MAXSEQ_DPU_NBYTES=${ENCODER_STACK_MAXSEQ_DPU_NBYTES} \
	 -DBUDDY_ALLOCATOR_SIZE_BYTE=${BUDDY_ALLOCATOR_SIZE_BYTE}\
	 -DHEAP_MAX_ITEM=${HEAP_MAX_ITEM}\
	 -DDECODER_QUANT_FRAC_BITS=${DECODER_QUANT_FRAC_BITS}"
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
get_filename_component(NUMPY_1_13_INCLUDE  $ENV{HOME}/.local/lib/python2.7/site-packages/numpy/core/include  ABSOLUTE)
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

target_compile_definitions(NRpyDNAcode PUBLIC -DMAX_DECODED_VBITS=${MAX_DECODED_VBITS})             
target_compile_definitions(NRpyDNAcode PUBLIC -DHEAP_MAX_ITEM=${HEAP_MAX_ITEM})             
target_compile_definitions(NRpyDNAcode PUBLIC -DENABLE_DPU_PRINT=${ENABLE_DPU_PRINT}) 
target_compile_definitions(NRpyDNAcode PUBLIC -DENABLE_TEST_DPU_HEAP=${ENABLE_TEST_DPU_HEAP})
target_compile_definitions(NRpyDNAcode PUBLIC -DENABLE_HEAP_HOST_DEBUGGING=${ENABLE_HEAP_HOST_DEBUGGING})
target_compile_definitions(NRpyDNAcode PUBLIC -DNR_DPUS=${NR_DPUS})
target_compile_definitions(NRpyDNAcode PUBLIC -DNR_TASKLETS=${NR_TASKLETS})
target_compile_definitions(NRpyDNAcode PUBLIC -DDPU_DECODE_BINARY="${dpudecoder-path}")
target_compile_definitions(NRpyDNAcode PUBLIC -DDPU_ENCODE_BINARY="${dpuencoder-path}")     
target_compile_definitions(NRpyDNAcode PUBLIC -DMESURE_PERF=${MESURE_PERF})     
target_compile_definitions(NRpyDNAcode PUBLIC -DMESURE_BW=${MESURE_BW})     
target_compile_definitions(NRpyDNAcode PUBLIC -DDECODER_QUANT_FRAC_BITS=${DECODER_QUANT_FRAC_BITS})

add_library(NRpyRS SHARED ../hedges/NRpyRS.cpp)
target_include_directories(NRpyRS PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
                                         ${PYTHON_27_INCLUDE} 
                                         ${NUMPY_1_13_INCLUDE}
)
set_target_properties(NRpyDNAcode PROPERTIES PREFIX "")
set_target_properties(NRpyRS PROPERTIES PREFIX "")

