.PHONY: clean build

PWD = $(shell pwd)

DOCKER_IMG = ghcr.io/upmem/usecase_dnastorage/hedges/ubuntu18:sw-20-release

build_docker:
	docker pull ${DOCKER_IMG}

clean:
	rm -rf ${PWD}/build
build: ${SRC}
	cmake -B build -GNinja ${PWD} && cd build && ninja

testDpuProgrammSimulator:
	PYTHONPATH=build UPMEM_PROFILE_BASE=\"backend=simulator\"   python -u hedges_pipeline.py

testDpuProgrammSimulatorStackAnalyzer:
	dpu_stack_analyzer build/dpuencoder-build-dpu/dpuencoder-dpu
	dpu_stack_analyzer build/dpudecoder-build-dpu/dpudecoder-dpu

testDpuProgrammSimulatorTrace:
	UPMEM_PROFILE_BASE=\"backend=simulator\"  PYTHONPATH=build python -u hedges_pipeline.py && dputrace -i trace-0000-00 -t 1 -no-color -no-tree trace-0000-00 > trace.out

testDpuProgrammHw:
	PYTHONPATH=build   python -u hedges_pipeline.py

testDpuProgrammHwPostlldb:
	PYTHONPATH=build   python -u hedges_pipeline.py && dpu-lldb-attach-dpu 0.0.0 ./build/dpudecoder-build-dpu/dpudecoder-dpu


testAppMemProfiler:
	PYTHONPATH=build dpu-profiling memory-transfer -- python -u hedges_pipeline.py

testAppFuncProfiler:
	PYTHONPATH=build dpu-profiling functions -A -- python -u hedges_pipeline.py

testSectionProfiler:
	dpu-profiling dpu-sections -- dpu-lldb --batch --one-line run -- ./build/dpudecoder-build-dpu/dpudecoder-dpu

testStatisticsProfiler:
	PYTHONPATH=/work/build dpu-profiling dpu-statistics --    python -u hedges_pipeline.py

testAppChromeProfiler:
	UPMEM_PROFILE_BASE=\"backend=simulator\" PYTHONPATH=/work/build dpu-profiling functions -o chrometf.json -A  -- python -u hedges_pipeline.py

.DEFAULT_GOAL := build
