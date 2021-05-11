#include "ZEntity.h"

namespace ZMesh {

Real* getByTag(const std::set<Real*, ZEntityLessThen> & list, int n)
{
	if ( ! &list )
	{
		std::cout<<" Real list have not been filled!\n";
		return NULL;
	}
	ZEntity tmp(n);

	std::set<Real*, ZEntityLessThen>::iterator it = list.find((Real *)&tmp);
	if(it != list.end())
		return *it;
	else
	{
		std::cout<<" Real list have not been filled!\n";
		return 0;
	}
}

Mat* getByTag(const std::set<Mat*, ZEntityLessThen> & list, int n)
{
	if ( ! &list )
	{
		std::cout<<" Real list have not been filled!\n";
		return NULL;
	}
	ZEntity tmp(n);
	std::set<Mat*, ZEntityLessThen>::iterator it = list.find((Mat *)&tmp);
	if(it != list.end())
		return *it;
	else
	{
		std::cout<<" Mat list have not been filled!\n";
		return 0;
	}
}

} // end namespace ZMesh
