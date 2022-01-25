CUDA_DIR = /usr/local/cuda-11.5/
sm = 75

build:
	mpicxx -g -fopenmp -c ./letter.cpp -o ./letter.o
	mpicxx -g -fopenmp -c ./sequence.cpp -o ./sequence.o
	mpicxx -g -fopenmp -c ./program.cpp -o ./program.o
	nvcc -g --std=c++17 -I./inc -c ./cuda_func.cu -o ./cuda_func.o
	mpicxx -g -I./inc -fopenmp -o ./program ./letter.o ./sequence.o ./program.o ./cuda_func.o $(CUDA_DIR)lib64/libcudart_static.a -ldl -lrt
clean:
	rm -f ./*.o
run:
	mpiexec -n 16 ./program
