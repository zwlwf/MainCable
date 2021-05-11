#include "ZLink.h"
#include "feMesh.h"

namespace ZMesh{

ZLink3D::ZLink3D( int n, int node1 , int node2) : ZElem(n, 3)
{ 
	int inds[2] = {node1, node2};
	if (mesh->getNodeList().size() == 0) 
	{
		std::cout<<" create node list at first! \n";
		return ;
	}
	setNodes( inds, 2, mesh->getNodeList() );
}

ZLink2D::ZLink2D( int n, int node1 , int node2) : ZElem(n, 1)
{ 
	int inds[2] = {node1, node2};
	if (mesh->getNodeList().size() == 0) 
	{
		std::cout<<" create node list at first! \n";
		return ;
	}

	setNodes( inds, 2, mesh->getNodeList() );
}

int ZLink2D::getNDofPerNode() const
{
	return 2;
}

int ZLink3D::getNDofPerNode() const
{
	return 3;
}

std::vector<int> ZLink2D::getDofInd() const
{
	int a[2]={0,1};
	return std::vector<int>(a,a+2);
}

std::vector<int> ZLink3D::getDofInd() const
{
	int a[3]={0,1,2};
	return std::vector<int>(a,a+3);
}

amat::Matrix<double> ZLink2D::elemKMatrix()
{
	double K,L;
	double dl[2];
	dl[0] = _nodes[1]->x() - _nodes[0]->x();
	dl[1] = _nodes[1]->y() - _nodes[0]->y();
	L = sqrt(dl[1]*dl[1]+dl[0]*dl[0]);
	double E, Area;
	E = _mat->E; 
	Area = _real->at(0);
	K = E*Area/L; // Mat have no num, and it start from 1
	double epsilon0 = _real->at(1);
	double A[2];
	A[0] = dl[0]/L; // cos(theta)
	A[1] = dl[1]/L; // sin(theta)
	amat::Matrix<double> B ;
	B.setup(4,4);
	// calculate B = A^T KA
	for(int i1=0; i1<2; i1++)
		for(int j1=0; j1<2; j1++)
			for(int i2=0; i2<2; i2++)
				for(int j2=0; j2<2; j2++)
					B(i1*2+i2,j1*2+j2) = K*(A[i2]*A[j2] + epsilon0*(( i2==j2? 1 : 0) - A[i2]*A[j2] ))*
						( i1==j1 ? 1 : -1 );

	return B;
}

amat::Matrix<double> ZLink3D::elemKMatrix()
{
	double K,L;
	double dl[3];
	dl[0] = _nodes[1]->x() - _nodes[0]->x();
	dl[1] = _nodes[1]->y() - _nodes[0]->y();
	dl[2] = _nodes[1]->z() - _nodes[0]->z();
	L = sqrt(dl[1]*dl[1]+dl[0]*dl[0]+ dl[2]*dl[2]);
	double E, Area;
	E = _mat->E; 
	Area = _real->at(0);
	double epsilon0 = _real->at(1);
	K = E*Area/L; // Mat have no num, and it start from 1
	double A[3];
	for(int i=0; i<3; i++)
		A[i] = dl[i]/L;
	amat::Matrix<double> B ;
	B.setup(6,6);
	for(int i1=0; i1<2; i1++)
		for(int j1=0; j1<2; j1++)
			for(int i2=0; i2<3; i2++)
				for(int j2=0; j2<3; j2++)
					B(i1*3+i2,j1*3+j2) = K*(A[i2]*A[j2] + epsilon0*(( i2==j2? 1 : 0) - A[i2]*A[j2] ))*
						( i1==j1 ? 1 : -1 );

	return B;
}
} // end namespace ZMesh
