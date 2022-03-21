# usecase_dnastorage
DNA storage error-correcting code pipeline (inner HEDGES, outer RS/LDPC) running on UPMEM PIMM DPU.

It is based on (2004821117, https://www.pnas.org/doi/full/10.1073/pnas.2004821117), the original HEDGES Research
Article. (Williamm H.Press, John A.Hawkins and all, Texas University, June 6, 2020).

The codebase relative to 2004821117 is (https://github.com/whpress/hedges).

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
(measured with 1 DPU with some estimated PIMM server with 2560 DPUs perfs)


| CR    |  srate    |    drate       |  irate        | strands/DPU run      | time (sec/DPU)  | decoding throughput (seq/sec/DPU) | (estimated) decoding throughput (seq/sec/2560 DPUs) | DPU pipeline efficiency (%) |
|-------|:----------|:---------------|:--------------|:---------------------|:----------------|:---------------------------------:|:---------------------------------:|:-------------------------------:|
| 0.5   |  0.03105  |   0.00945      |     0.00585   |  510                 |  29.72          |      17                           |      43,925                       |               69                |
| 0.33  |  0.03105  |   0.00945      |     0.00585   |  510                 |  8.31           |      61.4                         |     157,092                       |               62                |
| 0.25  |  0.03105  |   0.00945      |     0.00585   |  510                 |  4.37           |      116.6                        |     298,569                       |               69                |

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
