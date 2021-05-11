#ifndef _ZPATCH_H
#define _ZPATCH_H
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <list>
#include <stdlib.h>
#include <omp.h>
#include <set>
#include <unordered_map>
#include "amatrix.hpp" // for Matrix class in amat namespace
#include "aconstants.h"  // for amc_real like..., and have included in amatrix already.
#include "ZNode.h"


#ifndef MESHDATA_DOUBLE_PRECISION
#define MESHDATA_DOUBLE_PRECISION 14
#endif

/// \cond
namespace ZMesh {
/// \endcond
class feMesh;
class ZPatch
{
	friend feMesh;
	std::list<std::vector<ZNode*> > faces; //in fact this is not important at all
	int nFace; // number of faces consist of
	feMesh *mesh; // which mesh this patch belongs to
	amat::Matrix<int> dofs;
	amat::Matrix<double> areas; // in 2D case, it is the length of this boundary face
	public:
	ZPatch( feMesh *m)
	{
		mesh = m;
	}

	ZPatch() {} // empty constructor is needed is other constructor is defined
	ZPatch( const ZPatch & zp)
	{
		faces = zp.faces;
		nFace = zp.nFace;
		mesh = zp.mesh;
		dofs = zp.dofs;
		areas = zp.areas;
	}
	
	ZPatch( feMesh *m, amat::Matrix<int> Nodes);

	amat::Matrix<double> & getAreas() {	return areas;  }

	amat::Matrix<int> & getDofs() { return dofs; }

	void setMesh( feMesh *m)
	{
		mesh = m;
	}

};

} // end namespace ZMesh
#endif
