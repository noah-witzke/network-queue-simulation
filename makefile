main: main.cpp
	g++ -std=c++17 -o main.out main.cpp
	./main.out
	python3 figures.py