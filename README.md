## usecase_dnastorage
DNA storage error-correcting code pipeline (inner -> HEDGES, outer -> RS (Reed-Solomon)) running on UPMEM PIMM DPU.

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

## performances summary
(measured with 1 DPU with some estimated PIMM server with 2560 DPUs perfs)


| CR    |  srate    |    drate       |  irate        | strands/DPU run      | time (sec/DPU)  | decoding throughput (seq/sec/DPU) | (estimated) decoding throughput (seq/sec/2560 DPUs) | DPU pipeline efficiency (%) |
|-------|:----------|:---------------|:--------------|:---------------------|:----------------|:---------------------------------:|:---------------------------------:|:-------------------------------:|
| 0.5   |  0.03105  |   0.00945      |     0.00585   |  510                 |  29.72          |      17                           |      43,925                       |               69                |
| 0.33  |  0.03105  |   0.00945      |     0.00585   |  510                 |  8.31           |      61.4                         |     157,092                       |               62                |
| 0.25  |  0.03105  |   0.00945      |     0.00585   |  510                 |  4.37           |      116.6                        |     298,569                       |               69                |

## usage

### install dependencies (numpy1.13 + python2.7)

HEDGES project is built on nr3python library which depends
on numpy1.13 and python2.7.

```
./install.sh
```

### build project
```
make clean && make
```

### run on DPU HW
in hedges_pipeline.py, make sure that dpu is enabled for decoding step
```
test_dpu_decoder = True
```
```
make test
```
