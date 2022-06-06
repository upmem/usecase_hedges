.PHONY: clean build

PWD = $(shell pwd)

DOCKER_IMG = ghcr.io/upmem/usecase_dnastorage/hedges/ubuntu18:sw-20-release

build_docker:
	docker pull ${DOCKER_IMG}

clean:
	rm -rf ${PWD}/build
build: ${SRC}
	cmake -B build -GNinja ${PWD} && cd build && ninja

test:
	PYTHONPATH=build   python -u hedges_pipeline.py

profiling:
	UPMEM_PROFILE_BASE=\"backend=simulator\" PYTHONPATH=/work/build dpu-profiling functions -o chrometf.json -A  -- python -u hedges_pipeline.py

.DEFAULT_GOAL := build
