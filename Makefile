.PHONY: clean build

PWD = $(shell pwd)

# docker login -u dgerinmem -p ghp_... ghcr.io
# docker pull ghcr.io/upmem/usecase_dnastorage/hedges/debian10:sw-20-release
# DOCKER_IMG = ghcr.io/upmem/usecase_dnastorage/hedges/debian10:sw-20-release
DOCKER_IMG = ghcr.io/upmem/usecase_dnastorage/hedges/ubuntu18:sw-20-release
# FOR PERF
# sudo mount -o remount,mode=755 /sys/kernel/debug
# sudo mount -o remount,mode=755 /sys/kernel/debug/tracing
# sudo chmod a+w /sys/kernel/debug/tracing/uprobe_events
# sudo chmod a+w /sys/kernel/debug/tracing/kprobe_events
# sudo bash -c "echo -1 > /proc/sys/kernel/perf_event_paranoid"
# perf install at doccker boot (necessary for perf with debian based image)
# testAppFuncProfilerWithPerfInstall:
#   docker run --rm  -i --privileged --device=/dev/dpu_rank0    --security-opt seccomp=default.json   -v ${PWD}:/work ${DOCKER_IMG} bash -c "git clone --depth 1 --branch \"v$(uname -r | sed 's/\([0-9]\+\.[0-9]\+\).*/\1/g')\" https://git.kernel.org/pub/scm/linux/kernel/git/torvalds/linux.git && cd linux/tools/perf && make install && cp ./perf /usr/bin/perf  && cd /work && PYTHONPATH=/work/build  dpu-profiling functions -A -- python -u test_program.py"



build_docker:
	docker pull ${DOCKER_IMG}

clean:
	docker run --rm -v ${PWD}:/work ${DOCKER_IMG} bash -c "rm -rf /work/build/"
	rm -rf dpu/build
build: ${SRC}
	docker run --rm  -v ${PWD}:/work ${DOCKER_IMG} bash -c "mkdir -p /work/build && cd /work/build && cmake -GNinja /work/ && ninja"

testDpuProgrammSimulator:
	docker run --rm  -v ${PWD}:/work ${DOCKER_IMG} bash -c "cd /work/ && UPMEM_PROFILE_BASE=\"backend=simulator\"  PYTHONPATH=/work/build python -u test_program.py"

testDpuProgrammSimulatorStackAnalyzer:
	docker run --rm  -v ${PWD}:/work ${DOCKER_IMG} dpu_stack_analyzer work/build/dpu/hedges_encode_dpu
	docker run --rm  -v ${PWD}:/work ${DOCKER_IMG} dpu_stack_analyzer work/build/dpu/hedges_decode_dpu

testDpuProgrammSimulatorTrace:
	docker run --rm  -v ${PWD}:/work ${DOCKER_IMG} bash -c "cd /work/ && UPMEM_TRACE_DIR=/work  UPMEM_PROFILE_BASE=\"backend=simulator\"  PYTHONPATH=/work/build python -u test_program.py"

testDpuProgrammHw:
	docker run --rm  -i --privileged --device=/dev/dpu_rank0   -v ${PWD}:/work ${DOCKER_IMG} bash -c "cd /work &&   python -u test_program.py"

testAppMemProfiler:
	docker run --rm  -i --privileged --device=/dev/dpu_rank0    --security-opt seccomp=default.json   -v ${PWD}:/work ${DOCKER_IMG} bash -c "cd /work && PYTHONPATH=/work/build dpu-profiling memory-transfer -- python -u test_program.py"

testAppFuncProfiler:
	docker run --rm  -i --privileged --device=/dev/dpu_rank0    --security-opt seccomp=default.json   -v ${PWD}:/work ${DOCKER_IMG} bash -c "  dpu-profiling functions -A -- python -u test_program.py"

testSectionProfiler:
	dpu-profiling dpu-sections -- dpu-lldb --batch --one-line run -- ./build/dpu/hedges_dpu

testStatisticsProfiler:
	dpu-profiling dpu-statistics -- dpu-lldb --batch --one-line run -- ./build/dpu/hedges_dpu

testAppChromeProfiler:
	docker run --rm  -i --privileged --device=/dev/dpu_rank0   --security-opt seccomp=default.json  -v ${PWD}:/work ${DOCKER_IMG} bash -c "cd /work && UPMEM_PROFILE_BASE=\"backend=simulator\" PYTHONPATH=/work/build dpu-profiling functions -o chrometf.json -A  -- python -u test_program.py"

print_module_help_files:
	docker run --rm  -v ${PWD}:/work ${DOCKER_IMG} bash -c "cd /work/ && PYTHONPATH=/work/build/hedges python -u print_module_help_files.py"

.DEFAULT_GOAL := build
