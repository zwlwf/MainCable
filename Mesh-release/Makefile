CFLAG = -g -std=c++11 

MESHOBJECTS_small = feMesh.o ZElem.o ZNode.o ZEntity.o \
		 ZPatch.o ZLine.o littlePiece.o  ZBdNode.o alinalg.o

all: libZMesh.dll libZMesh.a

libZMesh.dll: $(MESHOBJECTS_small)
	g++ -shared -fPIC -o libZMesh.dll $(MESHOBJECTS_small)

libZMesh.a: $(MESHOBJECTS_small)
	ar -rcf libZMesh.a $(MESHOBJECTS_small)
	ranlib libZMesh.a

$(MESHOBJECTS_small): %.o : %.cpp
	g++ $(CFLAG) $< -c -o $@
