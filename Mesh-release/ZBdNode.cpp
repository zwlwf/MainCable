#include "ZBdNode.h"
namespace ZMesh 
{

ZBdNode* getByTag( const std::set<ZBdNode*, ZEntityLessThen> & list, int n)
{
	ZEntity tmp(n);
	std::set<ZBdNode*, ZEntityLessThen>::iterator it = list.find( (ZBdNode*) &tmp);
	if ( it != list.end() )
	{
		return *it;
	}
	else
		return NULL;
}

}
