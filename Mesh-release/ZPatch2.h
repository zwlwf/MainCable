#ifndef _ZPATCH2_H
#define _ZPATCH2_H
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
#include "ZElem.h"


#ifndef MESHDATA_DOUBLE_PRECISION
#define MESHDATA_DOUBLE_PRECISION 14
#endif

/// \cond
namespace ZMesh {
/// \endcond
class feMesh;
class ZPatch2
{
	friend feMesh;
	std::list< ZElem* > faces; 
	int nFace; // number of faces consist of
	feMesh *mesh; // which mesh this patch belongs to
	public:
	ZPatch2( feMesh *m)
	{
		mesh = m;
		nFace=0;
	}

	ZPatch2() // empty constructor is needed is other constructor is defined
	{
		nFace=0;
	} 

	ZPatch2( const ZPatch2 & zp)
	{
		faces = zp.faces;
		nFace = zp.nFace;
		mesh = zp.mesh;
	}

	std::list< ZElem *>::iterator begin()
	{
		return faces.begin();
	}
	
	std::list< ZElem *>::iterator end()
	{
		return faces.end();
	}

	void add( ZElem * ep)
	{
		faces.push_back( ep );
		nFace++;
	}

	int getNFaces()
	{  return nFace; }

	void setMesh( feMesh *m)
	{
		mesh = m;
	}

};

} // end namespace ZMesh
#endif
