#include "feMesh.h"

namespace ZMesh{

std::map<int,int> feMesh::gmsh2FEM_type = {
	{1,2},   // 2-node line ---> ZBeam2D
    {2,6},   // 3-node triangle ---> ZTri2D1V
    {3,8},   // 4-node quadrangle ---> ZQuad4Q
    {4,7},   // 4-node tetrahedron ---> ZTet
    {5,0},   // 8-node hexahedron ---> ?
    {6,0},   // 6-node prim ---> ?
    {7,0},   // 5-node pyramid ---> ?
    {8,10},   // 3-node second order line  ---> ZLine3Q
    {9,0},   // 6-node second order triangle  ---> ?
    {10,0},   // 9-node second order quadrangle  ---> ?
    {11,0},   // 10-node second order tetrahedron  ---> ?
    {15,13},  // 1-node point  ---> ZPoint
    {16,11}   // 8-node second order quadrangle ---> ZQuad8Q
	};   

std::map<int,int> feMesh::FEM2gmsh_type = {
	{FEM_ZLink2D , GMSH_2_NODE_LINE},
	{FEM_ZBeam2D , GMSH_2_NODE_LINE},
	{FEM_ZLink3D , GMSH_2_NODE_LINE},
	{FEM_ZTri2D , GMSH_3_NODE_TRI},
	{FEM_ZBeam3D , GMSH_2_NODE_LINE},
	{FEM_ZTri2D1V , GMSH_3_NODE_TRI},
	{FEM_ZTet , GMSH_4_NODE_TET},
	{FEM_ZQuad4Q , GMSH_4_NODE_QUAD},
	{FEM_ZLine2 ,  GMSH_2_NODE_LINE},
	{FEM_ZLine3Q , GMSH_3_NODE_LINE},
	{FEM_ZQuad8Q , GMSH_8_NODE_QUAD},
	{FEM_ZTri3D , GMSH_3_NODE_TRI},
	{FEM_ZPoint , GMSH_1_NODE_POINT}
	};   

// more see vtk-file-format
std::map<int,int> feMesh::FEM2vtk_type = {
	{FEM_ZLink2D ,3}, // VTK_LINE
	{FEM_ZBeam2D ,3}, // VTK_LINE
	{FEM_ZLink3D , 3}, // VTK_LINE
	{FEM_ZTri2D , 5}, // VTK_TRIANGLE(=5)
	{FEM_ZBeam3D , 3}, // VTK_LINE
	{FEM_ZTri2D1V , 5}, //VTK_TRIANGLE(=5)
	{FEM_ZTet , 10}, //VTK_TETRA (=10)
	{FEM_ZQuad4Q , 9}, //VTK_QUAD (=9)
	{FEM_ZLine2 , 3},
	{FEM_ZLine3Q , 4}, // VTK_POLY_LINE (=4)
	{FEM_ZQuad8Q , 7}, // VTK_POLYGON (=7)
	{FEM_ZTri3D , 5},
	{FEM_ZPoint , 1}  // VTK_VERTEX (=1)
	};   

void feMesh::readRealData( std::string mfile )
{
	std::ifstream is( mfile.c_str() );
	if ( !is )
	{
		std::cout<<" File not opened! \n";
		return;
	}
	std::string line;
	while ( getline(is, line) )
	{
		Real *tmp;
		tmp = new Real[1];
		tmp->set( line );
		realList.insert( tmp );
	}
}

void feMesh::nlist() 
{
	std::cout<<"List of Node"<<std::endl;
	for( std::set<ZNode *, ZEntityLessThen>::iterator it =_nlist.begin();
			it != _nlist.end();
			it++ )
	{
		std::cout<<**it<<std::endl;
	}
}

amat::Matrix<double> feMesh::Move( double dx, double dy, double dz=0.0 )
{
	amat::Matrix<double> deltaDisp(_nlist.size(), 3);
	int i=0;
	for( std::set<ZNode *, ZEntityLessThen>::iterator it =_nlist.begin();
			it != _nlist.end();
			it++ )
	{
		ZNode * np = *it;
		deltaDisp(i,0) = dx;
		deltaDisp(i,1) = dy;
		deltaDisp(i,2) = dz;
		i++;
	}
	return deltaDisp;
}

amat::Matrix<double> feMesh::Rotate( amat::Matrix<double> center, double theta )
{
	amat::Matrix<double> deltaDisp(_nlist.size(), 3);
	int i=0;
	for( std::set<ZNode *, ZEntityLessThen>::iterator it =_nlist.begin();
			it != _nlist.end();
			it++ )
	{
		ZNode * np = *it;
		double rx = np->x() - center(0);
		double ry = np->y() - center(1);
		deltaDisp(i,0) = rx*cos(theta) - ry*sin(theta) - rx;
		deltaDisp(i,1) = rx*sin(theta) + ry*cos(theta) - ry;
		deltaDisp(i,2) = 0.0;
		i++;
	}
	return deltaDisp;
}

void feMesh::elist()
{
	std::cout<<"List of Element"<<std::endl;
	for( std::set<ZElem *, ZEntityLessThen>::iterator it =_elist.begin();
			it != _elist.end();
			it++ )
	{
		std::cout<<**it<<std::endl;
	}
}
void feMesh::listVolume()
{
	typedef std::set<ZElem *, ZEntityLessThen> elemSet;
	std::cout<<"List of element Volume"<<std::endl;
	for( elemSet::iterator it = _elist.begin();
			it != _elist.end();
			it++)
	{
		std::cout<<(*it)->getVolume()<<std::endl;
	}
}

void feMesh::rlist()
{
	for ( std::set<Real*, ZEntityLessThen>::iterator it = realList.begin();
			it != realList.end();
			it++ )
		(*it)->output();
}

void feMesh::bdlist()
{
	std::cout<<"\n\n Boundary Node List \n";
	for (int i=0; i<boundaryNodes.size(); i++)
	{
		boundaryNodes[i]->output();
	}
}

void feMesh::patchList()
{
	std::cout<<" Patch List \n";
	for (
			str2ZPatchMap::const_iterator itm = patches.begin();
			itm != patches.end();
			itm++
		)
	{
		std::cout<< itm->first<<std::endl;
	}
}

void feMesh::bdfacelist()
{
	std::cout<<" list of boundaryfaces \n";
	for (int i=0; i<boundaryFaces.size(); i++)
	{
		std::cout<<i+1<<"#BF ";
		boundaryFaces[i]->disp();
		std::cout<<std::endl;
	}
}

/* create the dof number for each node:
 * look throught _elist, and determinate the number of dof for each node
 */
void feMesh::createNodeDof()
{
	int nInd=0;
	for(std::set<ZNode*,ZEntityLessThen>::iterator itn = this->nBegin();
			itn != this->nEnd();
			itn++)
	{
		(*itn) -> setDof(nInd++);
	}
}
void feMesh::createNodeDof2()
{
	// if created node dofs, clear the data arrays
	if ( nodeDofLocal.size() )
		nodeDofLocal.clear();
	if ( nodeDof.size() )
		nodeDof.clear();
	if ( elemNodeDofLocal.size() )
		elemNodeDofLocal.clear();

	for( std::set<ZElem*,ZEntityLessThen>::iterator it = _elist.begin();
			it != _elist.end();
			it++ )
	{
		ZElem * ep = (*it);
		for (int i=0; i<ep->getNNode(); i++)
		{
			int nodeTag = (*ep)[i]->getNum();
			for( int idof : ep->getDofInd() )
				nodeDofLocal[nodeTag].insert(idof);
		}
	}

	int ind=0; // start index of Dof
	bool toDispDOfList = false;
	int nInd=0;
	for( std::set<ZNode*,ZEntityLessThen >::const_iterator it = _nlist.begin();
			it != _nlist.end();
			it++ )
	{
		int nodeTag = (*it)->getNum();
		nDofPerNode = nodeDofLocal[nodeTag].size();
		for(int i=0; i<nDofPerNode; i++)
			nodeDof[nodeTag].push_back(ind+i);
		(*it)->setDof(nInd++);
		ind += nDofPerNode;
	}
	totalDof = ind;

	// the same node for different element has different dofs, store it in array elemNodeDofLocal
	for( ZElem* ep : *this )
	{
		elemNodeDofLocal[ep->getNum()] = amat::Matrix<int> ( ep->getNNode(), ep->getNDofPerNode() );
		std::vector<int> dofInElemTmp = ep->getDofInd();
		for(int inode = 0; inode < ep->getNNode(); inode++)
		{
			ZNode *np = (*ep)[inode];
			int localDofInElem = -1;
			for( int idof : nodeDofLocal[np->getNum()] )
			{
				// find idof in ep->getDofInd()
				localDofInElem++;
				for(int iii=0; iii<dofInElemTmp.size(); iii++)
					if( dofInElemTmp[iii] == idof )
					{
						elemNodeDofLocal[ep->getNum()](inode, iii) = localDofInElem;
						break;
					}
			}
		}
	}
}

std::list<ZElem *> createMap( feMesh *m, feMesh *base )
{
	std::list<ZElem*> elist;
	int iCount=0;
	for( std::set<ZNode*, ZEntityLessThen>::const_iterator itn = m->nBegin();	
			  itn != m-> nEnd();
			  itn ++ )
	{
		iCount++;
		if (iCount%10000 == 0 ) std::cout<<iCount/10000<<" w nodes done" <<std::endl;
		ZElem * closestE = *(base->eBegin()); // elem that is closest to current node
		double minDist=1.0e10;
		double disTmp;

		std::set<ZElem*, ZEntityLessThen>::const_iterator ite; 
		for(     ite = base->eBegin();
				 ite != base->eEnd();
				 ite ++ )
		{
			disTmp = (*ite)->distance( *itn );
			if ( disTmp <= 0.0 )
			{
				elist.push_back( *ite );
				break;
			}
			else
			{
				if( disTmp < minDist ) 
				{
					minDist = disTmp;
					closestE = *ite;
				}
			}
		}
		if ( ite == base->eEnd() )
		{
			std::cout<<"Warning: Node "<<(*itn)->getNum()<<" not in one of the base element !"<<std::endl;
			elist.push_back(closestE);
		}
	}
	return elist;
}

} // end namespace
