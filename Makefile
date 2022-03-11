.PHONY: clean build

PWD = $(shell pwd)

DOCKER_IMG = ghcr.io/upmem/usecase_dnastorage/hedges/ubuntu18:sw-20-release

build_docker:
	docker pull ${DOCKER_IMG}

clean:
	docker run --rm -v ${PWD}:/work ${DOCKER_IMG} bash -c "rm -rf /work/build/"
	rm -rf dpu/build
build: ${SRC}
	docker run --rm  -v ${PWD}:/work ${DOCKER_IMG} bash -c "mkdir -p /work/build && cd /work/build && cmake -GNinja /work/ && ninja"

testDpuProgrammSimulator:
	docker run --rm  -v ${PWD}:/work ${DOCKER_IMG} bash -c "cd /work/ && UPMEM_PROFILE_BASE=\"backend=simulator\"  PYTHONPATH=/work/build python -u hedges_pipeline.py"

testDpuProgrammSimulatorStackAnalyzer:
	docker run --rm  -v ${PWD}:/work ${DOCKER_IMG} dpu_stack_analyzer /work/build/dpuencoder-build-dpu/dpuencoder-dpu
	docker run --rm  -v ${PWD}:/work ${DOCKER_IMG} dpu_stack_analyzer /work/build/dpudecoder-build-dpu/dpudecoder-dpu

testDpuProgrammSimulatorTrace:
	docker run --rm  -v ${PWD}:/work ${DOCKER_IMG} bash -c "cd /work/ && UPMEM_TRACE_DIR=/work  UPMEM_PROFILE_BASE=\"backend=simulator\"  PYTHONPATH=/work/build python -u hedges_pipeline.py && dputrace -i trace-0000-00 -t 1 -no-color -no-tree trace-0000-00 > trace.out"

testDpuProgrammHw:
	docker run --rm  -i --privileged --device=/dev/dpu_rank0   -v ${PWD}:/work ${DOCKER_IMG} bash -c "cd /work && PYTHONPATH=/work/build   python -u hedges_pipeline.py"


testDpuProgrammHwPostlldb:
	docker run --rm  -i --privileged --device=/dev/dpu_rank0   -v ${PWD}:/work ${DOCKER_IMG} bash -c "cd /work && PYTHONPATH=/work/build   python -u hedges_pipeline.py && dpu-lldb-attach-dpu 0.0.0 ./build/dpudecoder-build-dpu/dpudecoder-dpu"


testAppMemProfiler:
	docker run --rm  -i --device=/dev/dpu_rank0    --security-opt seccomp=dpu_recipes/docker/docker_with_perf.json --mount type=bind,source=/sys/kernel/tracing,target=/sys/kernel/tracing -v ${PWD}:/work ${DOCKER_IMG} bash -c "cd /work && PYTHONPATH=/work/build dpu-profiling memory-transfer -- python -u hedges_pipeline.py"

testAppFuncProfiler:
	branch=4.19
	docker run --rm -it --security-opt seccomp=dpu_recipes/docker/docker_with_perf.json  --mount type=bind,source=/sys/kernel/tracing,target=/sys/kernel/tracing  -v  ${PWD}:/work  ghcr.io/upmem/usecase_dnastorage/hedges/ubuntu18:sw-20-release bash -c "apt install libelf-dev -y && python3 -m pip install --upgrade pip && pip3 install paramiko && rm -rf /linux && cd / && git clone --depth 1 --branch v4.19 https://git.kernel.org/pub/scm/linux/kernel/git/torvalds/linux.git && cd /linux/tools/perf && make install && cp ./perf /usr/bin/perf && cd /work && PYTHONPATH=/work/build/host:/usr/lib/python3/dist-packages:/usr/local/lib/python3.6/dist-packages/ dpu-profiling functions -A -- python -u hedges_pipeline.py"

testSectionProfiler:
	dpu-profiling dpu-sections -- dpu-lldb --batch --one-line run -- ./build/dpudecoder-build-dpu/dpudecoder-dpu

testStatisticsProfiler:
	docker run --rm  -i --privileged --device=/dev/dpu_rank0   -v ${PWD}:/work ${DOCKER_IMG} bash -c "cd /work && PYTHONPATH=/work/build dpu-profiling dpu-statistics --    python -u hedges_pipeline.py"

testAppChromeProfiler:
	docker run --rm  -i --privileged --device=/dev/dpu_rank0   --security-opt seccomp=dpu_recipes/docker/docker_with_perf.json  -v ${PWD}:/work ${DOCKER_IMG} bash -c "cd /work && UPMEM_PROFILE_BASE=\"backend=simulator\" PYTHONPATH=/work/build dpu-profiling functions -o chrometf.json -A  -- python -u hedges_pipeline.py"

print_module_help_files:
	docker run --rm  -v ${PWD}:/work ${DOCKER_IMG} bash -c "cd /work/ && PYTHONPATH=/work/build/hedges python -u print_module_help_files.py"

.DEFAULT_GOAL := build
