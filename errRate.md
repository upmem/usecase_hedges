# twist parameters
# (srate, drate, irate) = array([0.02, 0.03, 0.01])
# totstrandlen = 200  # total length of DNA strand


## pure CPU Error Rate validation


# scenario 0

[HEDGES TEST] for each packet, these statistics are shown in two groups:
[HEDGES TEST] 1.1 HEDGES decode failures, 1.2 HEDGES bytes thus declared as erasures
[HEDGES TEST] 1.3 R-S total errors detected in packet, 1.4 max errors detected in a single decode
[HEDGES TEST] 2.1 R-S reported as initially-uncorrected-but-recoverable total, 2.2 same, but max in single decode
[HEDGES TEST] 2.3 R-S total error codes; if zero, then R-S corrected all errors
[HEDGES TEST] 2.4 Actual number of byte errors when compared to known plaintext input
[HEDGES TEST] CR index 3 : 0.5
[HEDGES TEST] real inner CR for hedges (inner) level  0.38
[HEDGES TEST] test_dpu_encoder 0
[HEDGES TEST] test_dpu_decoder 0
[HEDGES TEST] test_dpu_statistics 0
[HEDGES TEST] npackets 400
[HEDGES TEST] hlimit 90000
[HEDGES TEST] dpu_profiling 0
[MAX HEAP HOST] 89839
[DPU HEAP_MAX_ITEM] 93000
[WARINIG][DPU] hlimit too big for DPU heap limit HEAP_MAX_ITEM
 

[END][CPU] (baddecodes 10298, erasures 36923, tot_detect 93454, max_detect 8034) ( totbytes 1338000 , tot_uncorect  531, max_uncorect  467, toterrcodes 1326, badbytes 27123)
[END][CPU] DNA Perr 0.06 , byteErrorRate 2.03e-02


# scenario 1 (500K Hlimit)

[HEDGES TEST] for each packet, these statistics are shown in two groups:
[HEDGES TEST] 1.1 HEDGES decode failures, 1.2 HEDGES bytes thus declared as erasures
[HEDGES TEST] 1.3 R-S total errors detected in packet, 1.4 max errors detected in a single decode
[HEDGES TEST] 2.1 R-S reported as initially-uncorrected-but-recoverable total, 2.2 same, but max in single decode
[HEDGES TEST] 2.3 R-S total error codes; if zero, then R-S corrected all errors
[HEDGES TEST] 2.4 Actual number of byte errors when compared to known plaintext input
[HEDGES TEST] CR index 3 : 0.5
[HEDGES TEST] real inner CR for hedges (inner) level  0.38
[HEDGES TEST] test_dpu_encoder 0
[HEDGES TEST] test_dpu_decoder 0
[HEDGES TEST] test_dpu_statistics 0
[HEDGES TEST] npackets 400
[HEDGES TEST] hlimit 500000
[HEDGES TEST] dpu_profiling 0
[MAX HEAP HOST] 499942
[DPU HEAP_MAX_ITEM] 93000


[END][CPU] (baddecodes 4208, erasures 12532, tot_detect 57477, max_detect 5168) ( totbytes 1338000 , tot_uncorect    6, max_uncorect    6, toterrcodes   17, badbytes  283)
[END][CPU] DNA Perr 0.06 , byteErrorRate 2.12e-04

