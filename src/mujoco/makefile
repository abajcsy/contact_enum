MUJOCO_PATH := /home/abajcsy/Documents/mujoco/mjpro150
COMMON=-O2 -isystem ${MUJOCO_PATH}/include -L ${MUJOCO_PATH}/bin -std=c++11 -mavx

demo: mujoco
	g++ $(COMMON) demo.cpp -lmujoco150nogl -lmujoco150 -lGL -lglew $(MUJOCO_PATH)/bin/libglfw.so.3 -o demo

.PHONY: all mojoco
mujoco:
	g++ $(COMMON) ${MUJOCO_PATH}/sample/compile.cpp -lmujoco150nogl -o $(MUJOCO_PATH)/bin/compile
	g++ $(COMMON) ${MUJOCO_PATH}/sample/derivative.cpp -lmujoco150nogl -fopenmp -o $(MUJOCO_PATH)/bin/derivative
	g++ $(COMMON) ${MUJOCO_PATH}/sample/simulate.cpp -lmujoco150 -lGL -lglew $(MUJOCO_PATH)/bin/libglfw.so.3 -o $(MUJOCO_PATH)/bin/simulate
	g++ $(COMMON) ${MUJOCO_PATH}/sample/record.cpp -lmujoco150 -lGL -lglew $(MUJOCO_PATH)/bin/libglfw.so.3 -o $(MUJOCO_PATH)/bin/record
	g++ $(COMMON) ${MUJOCO_PATH}/sample/basic.cpp -lmujoco150 -lGL -lglew $(MUJOCO_PATH)/bin/libglfw.so.3 -o $(MUJOCO_PATH)/bin/basic

all: mujoco demo

.DEFAULT_TARGET: all