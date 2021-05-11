#ifndef _ZNODE_H
#define _ZNODE_H

#include "ZEntity.h"
#include <iostream>
#include <set>

namespace ZMesh{

class ZNode : public ZEntity
{
	double _x,_y,_z,_dof;
	public:
	ZNode(): ZEntity(-1) // initial a empty node with tag -1
	{}

	ZNode( int n, double x=0.0, double y=0.0, double z=0.0 )
		: ZEntity(n)
	{
		_x = x;
		_y = y;
		_z = z;
		_dof = n-1; // default dof, should be edit when read in mesh in feMesh.
	}

	ZNode( const ZNode & n1) : ZEntity( n1._num )
	{
		_x = n1._x;
		_y = n1._y;
		_z = n1._z;
		_dof = n1._dof;
	}
	
	ZNode operator=( const ZNode & n1)
	{
		_num = n1._num;
		_x = n1._x;
		_y = n1._y;
		_z = n1._z;
		_dof = n1._dof;
	}
		
	void setCoord( double x,double y, double z=0.0)
	{
		_x = x;
		_y = y;
		_z = z;
	}


	double& x() { return _x; }
	double& y() { return _y; }
	double& z() { return _z; }
	friend std::ostream & operator<<( std::ostream & os, ZNode & n);
	int dof() { return _dof; }
	int setDof(int dd) {  _dof = dd; }

};

ZNode* getByTag( std::set<ZNode*, ZEntityLessThen> & list, int n);
} // end namespace ZMesh
#endif
