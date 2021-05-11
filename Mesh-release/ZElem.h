#ifndef _ZELEM_H
#define _ZELEM_H
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <assert.h>
#include <vector>
#include <list>
#include <stdlib.h>
#include "amatrix.hpp" // for Matrix class in amat namespace
#include "aconstants.h"  // for amc_real like..., and have included in amatrix already.
#include "ZEntity.h"
#include "ZNode.h"

#ifndef MESHDATA_DOUBLE_PRECISION
#define MESHDATA_DOUBLE_PRECISION 14
#endif

/// \cond
namespace ZMesh {
/// \endcond
 
	enum MassMode {
	LUMPM,
	CONSISTENT
	};

class feMesh;

class GaussP
{
	int _n;
	std::vector<double> _W; // weight
	std::vector<double> _C; // coordinate
	public:
	GaussP(int n);

	double W(int i)
	{ return _W[i]; }

	double C(int i)
	{ return _C[i]; }

	int getN()
	{ return _n; }
};

class ZElem : public ZEntity
{
	int _type; // set type to private, can not changed by its sons ;XD
	// ZLink2D -> 1, ZBeam2D -> 2, ZLink3D ->3, ZTri2D ->4, ZBeam3D -> 5, ZTri2D1V ->6 , ZTet -> 7 , ZQuad4Q ->8, ZLine2 ->9, ZLine3Q->10, ZQuad8Q->11, ZTri3D->12, ZPoint->13
	protected:
	Real* _real;
	Mat* _mat; // material
	int _nNode;
	std::vector< ZNode*> _nodes;
	void setType(int tt) { _type=tt; }
	public:
	static feMesh *mesh; // globel mesh of model 
	static int defaultMat;
	static int defaultReal;
	ZElem( int n);
	ZElem( int n , int type) ;

	Real * getReal() { return _real; }
	Mat * getMat() { return _mat; }

	void setNodes( int inds[], 
				   int nNodePerElem,
				   const std::set<ZNode*, ZEntityLessThen> & nList);

	ZNode * operator[] ( int n) ;

	ZNode * getNode( int n);

	int getType() { return _type; }
	int getNNode() { return _nNode; }

	std::vector< ZNode*>::const_iterator begin()
	{
		return _nodes.begin();
	}

	std::vector< ZNode*>::const_iterator end()
	{
		return _nodes.end();
	}

	int getNNode() const { return _nNode; } 

	void setRealByTag ( const std::set<Real*, ZEntityLessThen> &rlist, int n)
	{
		_real = getByTag( rlist, n ); 
	}

	void setMatByTag (const std::set<Mat*, ZEntityLessThen> &mlist, int n ) ;

	std::vector<int> getAllNodeDofs();

	friend std::ostream & operator<<( std::ostream & os, ZElem & e);
	virtual void DisplayTypeName() const = 0;
	virtual std::string typeName() const = 0;
	virtual std::vector<double> xiI() const = 0; // return the local coordinate xi of each nodes
	virtual std::vector<double> etaI() const = 0; // return the local coordinate eta of each nodes
	virtual amat::Matrix<double> xi_eta_zeta() const = 0; // return the local coordinate (xi,eta,zeta) of each nodes
	virtual int getNDofPerNode() const = 0;
	virtual double getVolume() const=0;
	virtual amat::Matrix<double> GetJacobi(double xi, double eta) const = 0;

	virtual std::vector<double> vCenter() const=0;
	virtual std::vector<int> getDofInd() const=0; // return dof ind (like 0 1 3 for beam2D )
	virtual amat::Matrix<double> elemKMatrix() = 0;
	virtual amat::Matrix<double> elemMMatrix( MassMode mmode= LUMPM ) = 0;
	virtual amat::Matrix<double> basis( double x, double y=0.0, double z=0.0) const = 0; // shape function of element
	// shape function's dN_i/dx_j at point (x,y,z)
	virtual amat::Matrix<double> dNi( double x, double y=0.0, double z=0.0) const = 0; 
	// shape function \int Ni Nj d\Omega at the region of element
	virtual amat::Matrix<double> int_NiNj() const = 0; 
	// shape function \int Ni d\Omega at the region of element
	virtual amat::Matrix<double> int_Ni() const = 0; 
	// shape function \int \nabla Ni \cdot \nabla Nj d\Omega at the region of element
	virtual amat::Matrix<double> int_dNidNj() const = 0; 

	virtual bool isIn( ZNode * np) =0; 
	virtual void calcVolume() = 0; 

	// distance of node np to this element
	virtual double distance( ZNode * np) =0; 
	virtual amat::Matrix<double> elementForce( amat::Matrix<double> u )=0;
	virtual amat::Matrix<double> elementGlobalForce( amat::Matrix<double> u )=0;

	virtual void get6DofDisp( double xi, double* elemU, //input
							double d6[], double pos[], //ouput
						 double xi_y, double xi_z )=0;  // optional

	virtual void applyNodeForce( double xi, double f_in[], double *f_out) = 0;

};

ZElem* getByTag(std::set<ZElem*, ZEntityLessThen> & list, int n);

// distant between nodes np1 and np2
double dist( ZNode *np1, ZNode *np2 );

// distant between node np1 and line ( np2, np3 )
double dist( ZNode *np1, ZNode *np2 , ZNode *np3);

// distant between node np1 and plane surface ( np2, np3, np4 )
double dist( ZNode *np1, ZNode *np2 , ZNode *np3, ZNode *np4);

double det( double d[3][3] );

} // end namespace ZMesh
#endif
