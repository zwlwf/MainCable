#ifndef _ZLINE_H
#define _ZLINE_H
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <list>
#include <stdlib.h>
#include "ZElem.h"
#include "feMesh.h"

namespace ZMesh{

class ZLine2: public ZElem
{
	double _L;
	public:
	ZLine2( int n) : ZElem(n, 9)
	{ }
	ZLine2( int n, int node1 , int node2);
	virtual void DisplayTypeName() const
	{
		std::cout<<"Element Type: 2 node linear line element"<< std::endl;
	}

	virtual std::string typeName() const
	{
		return std::string("ZLine2");
	}

	virtual int getNDofPerNode() const
	{
		return 1;
	}

	virtual std::vector<int> getDofInd() const
	{
		int a[1]={0};
		return std::vector<int>(a,a+1);
	}

	virtual amat::Matrix<double> basis( double x, double y=0.0, double z=0.0) const {} // shape function of element
	// shape function's dN_i/dx_j at point (x,y,z)
	virtual amat::Matrix<double> dNi( double x, double y=0.0, double z=0.0) const {}
	virtual double getVolume() const 
	{ return _L; }
	// shape function \int Ni Nj d\Omega at the region of element
	virtual amat::Matrix<double> int_NiNj() const {}
	// shape function \int Ni d\Omega at the region of element
	virtual amat::Matrix<double> int_Ni() const 
	{
		amat::Matrix<double> intN(2,1);
		intN(0,0) = 0.5*_L;
		intN(1,0) = 0.5*_L;
		return intN;
	}
	// shape function \int \nabla Ni \cdot \nabla Nj d\Omega at the region of element
	virtual amat::Matrix<double> int_dNidNj() const {}

	virtual amat::Matrix<double> GetJacobi(double xi, double eta) const {}

	virtual std::vector<double> xiI() const
	{
		std::vector<double> xis(2,-1);
		xis[1] = 1;
		return xis;
	}

	virtual std::vector<double> etaI() const {}

	virtual amat::Matrix<double> xi_eta_zeta() const 
	{	std::cout<<"xi_eta_zeta function not for this element"<<std::endl;
	}

	virtual double distance( ZNode *np)
	{
		bool hasImplemented = false;
		assert( hasImplemented );
		return 0.0;
	}

	virtual bool isIn(ZNode * np)
	{
		std::cout<<"Not implemented for element type "<<getType() <<std::endl;
		return false;
	}

	virtual std::vector<double> vCenter() const
	{
		std::vector<double> C(3,0.0);
		for (int i=0; i<getNNode(); i++)
		{
			C[0] += _nodes[i]->x();
			C[1] += _nodes[i]->y();
			C[2] += _nodes[i]->z();
		}
		for(int i=0; i<3; i++)
			C[i] /= getNNode();
		return C;
	}

	virtual amat::Matrix<double> elemKMatrix()
	{ }



	#include "commonCode.h"
};

class ZLine3Q: public ZElem
{
	double _L;
	void getLength()
	{
		double xiI[3] = {-0.774596669241483,0,0.774596669241483};
		double W[3] = {0.5555555555555,0.8888888888888,0.5555555555555};
		int nGaussP=3;
		_L = 0.0;
		for(int i=0; i<nGaussP; i++)
		{
			double dx=0.0, dy=0.0;
			dx += _nodes[0]->x() * (xiI[i] - 0.5);
			dx += _nodes[1]->x() * (-2.0) *xiI[i];
			dx += _nodes[2]->x() * (xiI[i] + 0.5);
			dy += _nodes[0]->y() * (xiI[i] - 0.5);
			dy += _nodes[1]->y() * (-2.0) *xiI[i];
			dy += _nodes[2]->y() * (xiI[i] + 0.5);
			_L += std::sqrt(dx*dx + dy*dy)*W[i];
		}
	}

	public:
	ZLine3Q( int n) : ZElem(n, 10)
	{ }
	ZLine3Q( int n, int node1 , int node2, int node3);

	virtual void DisplayTypeName() const
	{
		std::cout<<"Element Type: 3-node quadratic line element"<< std::endl;
	}

	virtual std::string typeName() const
	{
		return std::string("ZLine3Q");
	}

	virtual int getNDofPerNode() const
	{
		return 1;
	}

	virtual std::vector<int> getDofInd() const
	{
		int a[1]={0};
		return std::vector<int>(a,a+1);
	}

	virtual amat::Matrix<double> basis( double x, double y=0.0, double z=0.0) const {} // shape function of element
	// shape function's dN_i/dx_j at point (x,y,z)
	virtual amat::Matrix<double> dNi( double x, double y=0.0, double z=0.0) const {}
	// shape function \int Ni Nj d\Omega at the region of element
	virtual amat::Matrix<double> int_NiNj() const {}
	// shape function \int Ni d\Omega at the region of element
	virtual amat::Matrix<double> int_Ni() const
	{
		amat::Matrix<double> intN(3,1);
		intN(0,0) = 1.0/6*_L;
		intN(1,0) = 4.0/6*_L;
		intN(2,0) = 1.0/6*_L;
		return intN;
	}
	// shape function \int \nabla Ni \cdot \nabla Nj d\Omega at the region of element
	virtual double getVolume() const 
	{ return _L; }
	virtual amat::Matrix<double> int_dNidNj() const {}

	virtual amat::Matrix<double> GetJacobi(double xi, double eta) const {}

	virtual std::vector<double> xiI() const
	{
		std::vector<double> xis(2,-1);
		xis[1] = 0.0;
		xis[2] = 1.0;
		return xis;
	}
	virtual std::vector<double> etaI() const {}

	virtual amat::Matrix<double> xi_eta_zeta() const 
	{	std::cout<<"xi_eta_zeta function not for this element"<<std::endl;
	}

	virtual double distance( ZNode *np)
	{
		bool hasImplemented = false;
		assert( hasImplemented );
		return 0.0;
	}

	virtual bool isIn(ZNode * np)
	{
		std::cout<<"Not implemented for element type "<<getType() <<std::endl;
		return false;
	}

	virtual std::vector<double> vCenter() const
	{
		std::vector<double> C(3,0.0);
		for (int i=0; i<getNNode(); i++)
		{
			C[0] += _nodes[i]->x();
			C[1] += _nodes[i]->y();
			C[2] += _nodes[i]->z();
		}
		for(int i=0; i<3; i++)
			C[i] /= getNNode();
		return C;
	}

	virtual amat::Matrix<double> elemKMatrix()
	{
		// do nothing
	}



	#include "commonCode.h"
};

} // end namespace ZMesh
#endif
