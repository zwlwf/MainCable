#ifndef _ZLINK_H
#define _ZLINK_H
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

class ZLink2D: public ZElem
{
	public:
	ZLink2D( int n) : ZElem(n, 1)
	{ }
	ZLink2D( int n, int node1 , int node2);
	virtual void DisplayTypeName() const
	{
		std::cout<<"Element Type: 2D link"<< std::endl;
	}

	virtual std::string typeName() const
	{
		return std::string("ZLink2D");
	}

	virtual int getNDofPerNode() const;

	virtual std::vector<int> getDofInd() const;

	virtual amat::Matrix<double> basis( double x, double y=0.0, double z=0.0) const {} // shape function of element
	// shape function's dN_i/dx_j at point (x,y,z)
	virtual amat::Matrix<double> dNi( double x, double y=0.0, double z=0.0) const {}
	virtual double getVolume() const {}
	// shape function \int Ni Nj d\Omega at the region of element
	virtual amat::Matrix<double> int_NiNj() const {}
	// shape function \int Ni d\Omega at the region of element
	virtual amat::Matrix<double> int_Ni() const {}
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

	//! if the output need a 6*6 or 5*5 ( for node shared with a beam )
	//! you need to reshape the output Matrix by just add zeros line and column.
	virtual amat::Matrix<double> elemKMatrix();
#define elemMMatrix_fun

	virtual amat::Matrix<double> elemMMatrix( MassMode mmode )
	{
		// lumped mass
		amat::Matrix<double> M(4,4);
		double L,c,s;
		double dl[2];
		dl[0] = _nodes[1]->x() - _nodes[0]->x();
		dl[1] = _nodes[1]->y() - _nodes[0]->y();
		L = std::sqrt(dl[1]*dl[1]+dl[0]*dl[0]);
		double rho, Area;
		rho = _mat->density; 
		Area = _real->at(0);
		double Me=rho*Area*L/2;
		double Kdata[16] = {
			Me, 0 ,0, 0,
			0, Me, 0, 0,
			0, 0, Me, 0,
			0, 0, 0, Me };
		M.setdata(Kdata, 4*4);
		return M;
	}

	#include "commonCode.h"
};
#undef elemMMatrix_fun

class ZLink3D: public ZElem
{
	public:
	ZLink3D( int n) : ZElem(n, 3)
	{ }
	ZLink3D( int n, int node1 , int node2);

	virtual void DisplayTypeName() const
	{
		std::cout<<"Element Type: 3D link"<< std::endl;
	}

	virtual std::string typeName() const
	{
		return std::string("ZLink3D");
	}

	virtual int getNDofPerNode() const;

	virtual std::vector<int> getDofInd() const;

	virtual amat::Matrix<double> basis( double x, double y=0.0, double z=0.0) const {} // shape function of element
	// shape function's dN_i/dx_j at point (x,y,z)
	virtual amat::Matrix<double> dNi( double x, double y=0.0, double z=0.0) const {}
	// shape function \int Ni Nj d\Omega at the region of element
	virtual amat::Matrix<double> int_NiNj() const {}
	// shape function \int Ni d\Omega at the region of element
	virtual amat::Matrix<double> int_Ni() const {}
	// shape function \int \nabla Ni \cdot \nabla Nj d\Omega at the region of element
	virtual double getVolume() const {}
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

	//! if the output need a 6*6 or 5*5 ( for node shared with a beam )
	//! you need to reshape the output Matrix by just add zeros line and column.
	virtual amat::Matrix<double> elemKMatrix();
#define elemMMatrix_fun

	virtual amat::Matrix<double> elemMMatrix( MassMode mmode )
	{
		// lumped mass
		amat::Matrix<double> M(6,6);
		double L,c,s;
		double dl[2];
		dl[0] = _nodes[1]->x() - _nodes[0]->x();
		dl[1] = _nodes[1]->y() - _nodes[0]->y();
		L = std::sqrt(dl[1]*dl[1]+dl[0]*dl[0]);
		double rho, Area;
		rho = _mat->density; 
		Area = _real->at(0);
		double Me=rho*Area*L/2;
		double Kdata[36] = {
			Me, 0 ,0, 0, 0, 0,
			0, Me, 0, 0, 0, 0,
			0, 0, Me, 0, 0, 0,
			0, 0, 0, Me, 0 ,0,
			0, 0, 0, 0, Me, 0,
			0, 0, 0, 0, 0 , Me };
		M.setdata(Kdata, 6*6);
		return M;
	}
	
	#include "commonCode.h"
};
#undef elemMMatrix_fun

} // end namespace ZMesh
#endif
