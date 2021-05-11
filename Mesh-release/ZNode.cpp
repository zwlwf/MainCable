#include "ZNode.h"

namespace ZMesh{

std::ostream & operator<<( std::ostream & os, ZNode &n)
{ 
	os <<n._num << "#N ["<< n._x <<", "<<n._y<<", "<<n._z <<"]"; 
	return os;
}

ZNode* getByTag( std::set<ZNode*, ZEntityLessThen> & list, int n)
{
	ZEntity tmp(n);
	std::set<ZNode*, ZEntityLessThen>::const_iterator it = list.find((ZNode *)&tmp);
	if(it != list.end())
		return *it;
	else
	{
		return NULL;
	}
}

}

