#ifndef _FEMESH_H
#define _FEMESH_H
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <list>
#include <stdlib.h>
#include <omp.h>
#include <set>
#include <map>
#include <unordered_map>
#include "amatrix.hpp" // for Matrix class in amat namespace
#include "aconstants.h"  // for amc_real like..., and have included in amatrix already.
#include "ZElem.h"
#include "ZBdNode.h"
#include "ZPatch.h"
#include "ZPatch2.h"
#include "ZNode.h"
#include "ZLine.h" 
#include "myMap.h"

#ifndef MESHDATA_DOUBLE_PRECISION
#define MESHDATA_DOUBLE_PRECISION 14
#endif

/// \cond
namespace ZMesh {
/// \endcond

class fem;
class fec;
class Runtime;
class nodeField;
class elemField;
class Field;

enum gmsh_ElemType
{
GMSH_2_NODE_LINE=1,
GMSH_3_NODE_TRI=2,
GMSH_4_NODE_QUAD=3,
GMSH_4_NODE_TET=4,
GMSH_8_NODE_HEX=5,
GMSH_6_NODE_PRIM=6,
GMSH_5_NODE_PYRAMID=7,
GMSH_3_NODE_LINE=8,
GMSH_6_NODE_TRI=9,
GMSH_9_NODE_QUAD=10,
GMSH_10_NODE_TET=11,
GMSH_1_NODE_POINT=15,
GMSH_8_NODE_QUAD=16 };

enum FEM_ElemType
{
FEM_NONE=0,
FEM_ZLink2D = 1,
FEM_ZBeam2D = 2,
FEM_ZLink3D = 3,
FEM_ZTri2D = 4,
FEM_ZBeam3D = 5,
FEM_ZTri2D1V = 6,
FEM_ZTet = 7,
FEM_ZQuad4Q = 8,
FEM_ZLine2 = 9,
FEM_ZLine3Q = 10,
FEM_ZQuad8Q = 11,
FEM_ZTri3D = 12,
FEM_ZPoint = 13,
FEM_ZBeamNA2D = 14,
FEM_ZBeamNA3D = 15 
};

/** A general mesh class mainly for beam mass and link*/
class feMesh
{
	friend fem;
	friend fec;
	friend Runtime;
	friend nodeField;
	friend elemField;
	friend Field;
protected:
	int nNodes;		//< Number of nodes
	int nElems;		//< Number of elements
	int nDim;		//< Dimension of the mesh
	int nLoadNodes; //< Number of load nodes
	int nBdNodes;   //< Number of boundary nodes with fixed displacement, Force is not a really boundary
	int totalDof;
	std::set<ZElem*, ZEntityLessThen> _elist;		//< number of nodes to an element
	std::set<ZNode*, ZEntityLessThen> _nlist;    	//< Specifies coordinates of each node
	std::set<Mat*, ZEntityLessThen> matList;
	std::set<Real*, ZEntityLessThen> realList;

	/// Boundary nodes 
	std::vector<ZBdNode* > boundaryNodes; 
	std::vector<ZBdFace2D* > boundaryFaces; 

	// Dof of each node, this is just for multi-element type problem
	std::unordered_map<int ,std::vector<int> > nodeDof;
	// local dof index for node, mainly used for multi-type element problem,
	// and can be used for static assemble function
	// elemNodeDofLocal[elemNum] (i,j) : dof ind in nodeDof[i node's Tag] of j dof of i node of elem elemNum, 
	std::unordered_map<int ,amat::Matrix<int> > elemNodeDofLocal;
	std::unordered_map<int,std::set<int> > nodeDofLocal;
	str2ZPatchMap patches; // store the nodes/faces for the  
	str2ZPatch2Map patches2; // store the boundary patch element 
	// following variable is for single element type problem
	int nNodePerElem;
	int nDofPerNode;

public:

	static std::map<int, int> gmsh2FEM_type ;
	static std::map<int, int> FEM2gmsh_type ; // seldom use
	static std::map<int, int> FEM2vtk_type ; // usually used for output to vtk format.

	feMesh() 
	{
		std::cout<<"set the dimension to 2 by default\n";
		nDim = 2;
	} // empty constructor

	feMesh( int dim) { nDim = dim; } // initial with dimension


	int lastElemNum()
	{
		std::set<ZElem*,ZEntityLessThen>::iterator ite = eEnd();
		ite--;
		return (*ite)->getNum();
	}

	int lastNodeNum()
	{
		std::set<ZNode*,ZEntityLessThen>::iterator itn = nEnd();
		itn--;
		return (*itn)->getNum();
	}
	
	void createDefaultMat() // create default Mat with num 1
	{
		Mat *matFactory;
		matFactory = new Mat(1); 
		matFactory->setE( 1.0 );// matList create temprary here
		matFactory->setNu( 0.3 );// set nu to air
		matList.insert( matFactory );
	}

	void createDefaultReal() // create default Real with num 1
	{
		Real *realFactory;
		realFactory = new Real[1]; 
		realFactory->set("r 1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 ");// realList create temprary here
		realList.insert( realFactory );
	}

    void insertMat( Mat * mat )
    {   matList.insert( mat );  }

    void insertReal( Real * real )
    {   realList.insert( real );  }
	
	void insertElem( ZElem* ep )
	{
		std::pair<std::set<ZElem*,ZEntityLessThen>::iterator, bool> result = _elist.insert( ep );
		if ( result.second )
			nElems++;
		else
			std::cout<<"Insert elements failed!"<<std::endl;
	}

	void insertNode( ZNode* np )
	{
		std::pair<std::set<ZNode*,ZEntityLessThen>::iterator, bool> result = _nlist.insert( np );
		if ( result.second )
			nNodes++;
		else
			std::cout<<"Insert Nodes failed!"<<std::endl;
	}

	int getNNode()
	{ return _nlist.size(); }

	int getNElem()
	{ return _elist.size(); }

	int getDim()
	{ return nDim; }

	std::set<ZNode*, ZEntityLessThen>::const_iterator nBegin()
	{
		return _nlist.begin();
	}

	std::set<ZNode*, ZEntityLessThen>::const_iterator nEnd()
	{
		return _nlist.end();
	}

	std::set<ZElem*, ZEntityLessThen>::const_iterator eBegin()
	{
		return _elist.begin();
	}

	std::set<ZElem*, ZEntityLessThen>::const_iterator eEnd()
	{
		return _elist.end();
	}

    std::set<ZElem*, ZEntityLessThen>::const_iterator begin() 
    {                                                         
        return _elist.begin();                                
    }                                                         
                                                              
    std::set<ZElem*, ZEntityLessThen>::const_iterator end()   
    {                                                         
        return _elist.end();                                  
    }    

	str2ZPatchMap & getPatches() 
	{
		return patches;
	}

	str2ZPatch2Map & getPatches2() 
	{
		return patches2;
	}

//! feMesh(const feMesh& other) // copy constructor

//!	feMesh& operator=(const feMesh& other)

	~feMesh() { }

	ZElem * getElemByTag( int ind ) 
	{
		ZEntity tmp(ind);
		std::set<ZElem*, ZEntityLessThen>::const_iterator it = _elist.find((ZElem*)&tmp);
		if(it != _elist.end())
			return *it;
		else
			return 0;
	}

	ZNode * getNodeByTag( int ind ) 
	{
		ZEntity tmp(ind);
		std::set<ZNode*, ZEntityLessThen>::const_iterator it = _nlist.find((ZNode*)&tmp);
		if(it != _nlist.end())
			return *it;
		else
			return NULL;
	}

	int getNodeDofByTag( int tag)
	{
		return nodeDof[tag][0];
	}

	// return the index in dofs of nTag node of UVWorRotXYZ 
	int getNodeDofLocalByTag( int nTag, int UVWorRotXYZ )
	{
		std::set<int>::iterator it = nodeDofLocal[nTag].find(UVWorRotXYZ);
		assert( it != nodeDofLocal[nTag].end() );
		return std::distance( nodeDofLocal[nTag].begin(), it );
	}

	Real * getRealByTag( int ind ) 
	{
		ZEntity tmp(ind);
		std::set<Real*, ZEntityLessThen>::const_iterator it = realList.find((Real*)&tmp);
		if(it != realList.end())
			return *it;
		else
			return 0;
	}

	Mat * getMatByTag( int ind ) 
	{
		ZEntity tmp(ind);
		std::set<Mat*, ZEntityLessThen>::const_iterator it = matList.find((Mat*)&tmp);
		if(it != matList.end())
			return *it;
		else
			return 0;
	}

	void createMesh(const char* _coordFile, const char* const _elem_NodeListFile, int Dim); // default is triangular, 2D, node has only x,y coordinate

	/// Reads mesh from Gmsh 2 format file for elem of Tri2D type
	void readGmsh(std::string mfile);

	void readGmsh2(std::string mfile);

	void readGmsh3(std::string mfile);

	void readRealData( std::string mfile);

	void readTest(std::string mfile, int );

	void readBoundaryNode( std::string mfile );
	
	void readBoundaryFace( std::string mfile );

	void readPatches ( std::string mfile);

	void readMesh(std::string pathname);

	void write2Foam(std::string pathName);

	void write2Tikz(std::string pathName);

	// create the dof number for each node, this should be called as soon as
	// _elist and _nlist is created
	void createNodeDof();

	// for different element type
	void createNodeDof2();
	
	// like the nlist command in ansys
	void nlist();

	// like the elist command in ansys
	void elist();

	// like the rlist command in ANSYS
	void rlist();

	// list the boundary 
	void bdlist();

	void bdfacelist();

	void patchList();

	amat::Matrix<double> Move( double dx, double dy, double dz );

	amat::Matrix<double> Rotate( amat::Matrix<double> center, double theta );

	const std::set<ZNode*, ZEntityLessThen> & getNodeList() const
	{ return _nlist; }

	void setNodeList( const std::set<ZNode*, ZEntityLessThen> &nlist )
	{ 
		nNodes = nlist.size();
		_nlist = nlist ; 
	}

	void setElemList( const std::set<ZElem*, ZEntityLessThen> &elist )
	{ 
		nElems = elist.size();
		_elist = elist ; 
	}

	int getTotalDof() { return totalDof; }

	const std::set<ZElem*, ZEntityLessThen> & getElemList() const
	{ return _elist; }

	const std::set<Real*, ZEntityLessThen> & getRealList() const
	{ return realList; }

	void setRealList( std::set<Real*, ZEntityLessThen> &rlist )
	{ 
		realList = rlist ; 
	}

	const std::set<Mat*, ZEntityLessThen> & getMatList() const
	{ return matList; }

	void setMatList( std::set<Mat*, ZEntityLessThen> &mlist )
	{ 
		matList = mlist ; 
	}

	void setBdList( std::vector<ZBdNode*> &bdlist )
	{ 
		boundaryNodes = bdlist ; 
	}

	void addBdNode( ZBdNode* bd)
	{
		boundaryNodes.push_back(bd);
	}

	std::set< ZNode*, ZEntityLessThen> getPatchNodes(std::string pname )
	{
		std::set<ZNode*, ZEntityLessThen> nset;
		ZPatch2 &zp = patches2[ pname ];
		for( std::list< ZElem*>::iterator ite = zp.begin();
				ite != zp.end();
				ite++ )
		{
			ZElem* ep = *ite;
			for(int i=0; i<ep->getNNode(); i++)
				nset.insert( (*ep)[i]);
		}
		return nset;
	}

	void listVolume();

    void write2vtk(const std::string pathname);
};

std::list<ZElem *> createMap( feMesh *m, feMesh *base );

} // end namespace ZMesh
#endif
