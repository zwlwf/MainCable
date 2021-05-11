#ifndef _ZENTITY_H
#define _ZENTITY_H
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <list>
#include <stdlib.h>
#include <set>
namespace ZMesh{
// defines the very basic class of ZEntity

class ZEntityLessThen;
class ZEntity
{
	friend ZEntityLessThen;
	protected:
	int _num; // number of ZEntity
	public:
	ZEntity( int n )
	{
		_num = n;
	}
			
	int getNum() const
	{
		return _num;
	}
};

class ZEntityLessThen
{
	public:
	ZEntityLessThen(){}
	bool operator() (ZEntity * const z1, ZEntity * const z2) const
	{
		return z1->getNum() < z2->getNum();
	}
};

class Mat : public ZEntity
{
	public:
	double E; //Yang's module
	double nu; // Possion ratio
	double alpha; // temperature epsilon ratio
	double density; // density of material
	Mat(int n) :ZEntity(n)
	{ 
		// set default material to steel
		E=2.1e11; 
		nu = 0.3;
		density = 7850;
		alpha=0.001;
	}
	void setE( double E1) { E = E1; }
	void setNu( double nu1) { nu = nu1; }
	void setAlpha( double alpha1) { alpha = alpha1; }
	void setDensity( double density1) { density = density1; }

	~Mat() {}
};

class Real: public ZEntity // as I store real in a set, For real data may too much
{
	std::vector<double> data ; 
	public:
	Real() :ZEntity(-1) { }
	Real( int n ) : ZEntity( n)
	{}

	void set( std::string command ) 
	{
		std::istringstream is(command);
		double r;
		char c;
		is >> c;
		if (c != 'r')
		{
			std::cout<<" command must start with r (delimeter is blank)!"<<std::endl;
			return ;
		}
		is>> r;
		_num = (int) r;
		while ( is>> r )
		{
			data.push_back( r);
		}
	}

	void output()
	{
		std::cout<<_num<<"#R [";
		for (int i=0; i<data.size(); i++)
			std::cout<<data[i]<<" ";
		std::cout<<"]\n";
	}

	double& at( int ind)
	{
		return data[ind];
	}

	void append(double ri)
	{
		data.push_back(ri);
	}

	int size()
	{
		return data.size();
	}

	~Real() {}
};

Real* getByTag(const std::set<Real*, ZEntityLessThen> & list, int n) ;
Mat* getByTag(const std::set<Mat*, ZEntityLessThen> & list, int n);
} // end namespace ZMesh
#endif
