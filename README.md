# usecase_dnastorage
DNA storage error-correcting code pipeline (inner HEDGES, outer RS/LDPC) running on UPMEM PIMM DPU.

It is based on (2004821117, https://www.pnas.org/doi/full/10.1073/pnas.2004821117), the original HEDGES Research
Article. (Williamm H.Press, John A.Hawkins and all, Texas University, June 6, 2020).

```
     ____________________________________________________________
    |   ___________      _______________      _______________    |
    |  | synthetic |    | outer encoder |    | inner encoder |   |
    |  | data      |--->| (R.S)         |--->| (HEDGES)      |   |
    |  |           |    |     [CPU]     |    |   [CPU/DPU]   |   |
    |  |___________|    |_______________|    |_______________|   |
    |                                               |            |
    |                                          [A.C.G.T...]      |
    |                                      (encoded DNA strands) |
    |   ___________      _______________      ______|________    |
    |  | synthetic |    | outer decoder |    | inner decoder |   |
    |  | data      |<---| (R.S)         |<---| (HEDGES)      |   |
    |  |           |    |     [CPU]     |    |    [CPU/DPU]  |   |
    |  |___________|    |_______________|    |_______________|   |
    |____________________________________________________________|
```

The current implementation support only one DPU as it consists of a POC.


# performances summary
(measured with 1 DPU with some estimated multi DPU perfs)


| CR    | time (sec/DPU)  | decoding throughput (seq/sec/DPU) | DPU pipeline efficiency (%) |
|-------|:----------------|:-----------------------------:|:-------------------------------:|
| 0.5   |  36.12          |      14.12                    |               69                |
| 0.33  |  10.0           |      50.51                    |               62                |
| 0.25  |  5.31           |      95.92                    |               69                |

# docker fetch

For convenience, we provide a packaged docker image available on ghrc.io.
This image is compatible with UPMEM PIMM cloud server environement.
No extra packages are required to build and run the project on DPU.

login to ghcr.io
```
docker login -u [GITHUB LOGIN] -p [GITHUB GHCR TOKEN] ghcr.io
```
docker pull HEDGES image
```
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

This mode is used to get detailed information about dpu program
performances, for the current run.

in test_programm.py
```
test_dpu_decoder = True
test_dpu_statistics = False
```

## run with statistics profiling mode

See (https://sdk.upmem.com/2021.2.0/260_Profiling.html) for more information about statistics profiling mode.
NOTE : DPU binaries needs to be compiled with '-pg' flag for statistics profiling.

in dpu_recipes/.cmakeDpuProject/CMakeLists.txt, add -pg to CMAKE_C_FLAGS
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
