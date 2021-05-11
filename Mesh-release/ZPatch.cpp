#include "ZPatch.h"
#include "feMesh.h"
#include <cmath>

using namespace ZMesh;
	
ZPatch::ZPatch( feMesh *m, amat::Matrix<int> Nodes) 
{
	mesh = m;
	nFace = Nodes.rows();
	dofs.setup(nFace,Nodes.cols());
	areas.setup(nFace,1);
	std::vector<ZNode*> faceTmp(Nodes.cols());
	for(int i=0; i<nFace; i++)
	{
		for (int j=0; j<Nodes.cols(); j++)
		{
			faceTmp[j] = mesh->getNodeByTag( Nodes(i,j));
			// nodeDof has to be created
			dofs(i,j) = mesh->getNodeDofByTag(Nodes(i,j));
		}
		faces.push_back(faceTmp); // has to compile with -std=c++11
		if ( Nodes.cols()>1)
		{
			double dx = faceTmp[1]->x() - faceTmp[0]->x(); 
			double dy = faceTmp[1]->y() - faceTmp[0]->y(); 
			areas(i) = sqrt(dx*dx + dy*dy);
		}
		else
			areas(i) = 0;
	}
}
