rootPath = .
CFLAG = -std=c++11 -I$(rootPath)/Mesh-release
LIBS =  -L$(rootPath)/Mesh-release -lZMesh 
a.exe : MainCable.o removeComments.o main.cpp
	g++ -g -o a.exe $(CFLAG) main.cpp MainCable.o removeComments.o $(LIBS) 
MainCable.o: MainCable.cpp 
	g++  -g $(CFLAG) -c MainCable.cpp 
removeComments.o: removeComments.cpp
	g++  -c removeComments.cpp
