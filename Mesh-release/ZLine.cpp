#include "ZLink.h"
#include "feMesh.h"

namespace ZMesh{

ZLine3Q::ZLine3Q( int n, int node1 , int node2, int node3) : ZElem(n, 10)
{ 
	int inds[3] = {node1, node2,node3};
	if (mesh->getNodeList().size() == 0) 
	{
		std::cout<<" create node list at first! \n";
		return ;
	}
	setNodes( inds, 3, mesh->getNodeList() );
	getLength();
}

ZLine2::ZLine2( int n, int node1 , int node2) : ZElem(n, 9)
{ 
	int inds[2] = {node1, node2};
	if (mesh->getNodeList().size() == 0) 
	{
		std::cout<<" create node list at first! \n";
		return ;
	}

	setNodes( inds, 2, mesh->getNodeList() );
	double dl[3];
	dl[0] = _nodes[1]->x() - _nodes[0]->x();
	dl[1] = _nodes[1]->y() - _nodes[0]->y();
	dl[2] = _nodes[1]->z() - _nodes[0]->z();
	_L = sqrt(dl[1]*dl[1]+dl[0]*dl[0]+ dl[2]*dl[2]);
}
} // end namespace ZMesh
