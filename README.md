# usecase_dnastorage
DNA storage error-correcting code pipeline (inner HEDGES, outer RS/LDPC) running on UPMEM PIMM DPU.

# docker fetch

For convenience, we provide a packaged docker image available on ghrc.io.
This image is compatible with UPMEM PIMM cloud server environement.
No extra packages are required to build and run the project on DPU.

```
# login to ghcr.io
docker login -u [GITHUB LOGIN] -p [GITHUB GHCR TOKEN] ghcr.io
# docker pull HEDGES image
make build_docker
```

# usage

## build project
```
make clean && make
```

## run on DPU HW
```
make testDpuProgrammHw
```

## run on DPU simulator
```
make testDpuProgrammSimulator
```

## run with global DPU profiling mode

This mode is used to get detailed information about dpu_programm
performances, for the current run.

in test_programm.py
```
test_dpu_decoder = True
test_dpu_statistics = False
```

## run with  statistics profiling mode

See (https://sdk.upmem.com/2021.2.0/260_Profiling.html) for more infotmations about statistics profiling mode.
NOTE : For statics mode, DPU binaries needs to be compiled with '-pg' flag for statistics profiling.

in dpu_recipes/.cmakeDpuProject/CMakeLists.txt
, add -pg to CMAKE_C_FLAGS
```
set(CMAKE_C_FLAGS "-O3 -pg -fstack-size-section -fshort-enums   -DNR_TASKLETS=${NR_TASKLETS} -DSTACK_SIZE_DEFAULT=${DPU_STACK_SIZE_BYTE} ${DPU_ARGN}")
```
then rebuild project
```
make clean && make
```

and set test_dpu_statistics=True in test_programm.py
```
test_dpu_decoder = True
test_dpu_statistics = True
```

Then run programm on HW with dpu-statistics profiler
```
make testStatisticsProfiler
```
