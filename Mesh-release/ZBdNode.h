#ifndef _ZBDNODE_H
#define _ZBDNODE_H

#include "ZNode.h"
#include <iostream>
#include <vector>
#include <cmath>

namespace ZMesh{

class ZBdNode : public ZNode
{
	/* in a more general way in FEM boundary node may be both loadBd and displaceBd,
	 * and the dofDisp and dofForce both exist.
	 * and loadBd means boundary that contribute to rhs item, 
	 * displaceBd means boundary that contribute to fixed node value.
	 * So I do not create two more derived class for loadBd and displaceBd
	 */
	public:
	enum BdType
	{
		loadBd,
		displaceBd
	};
	private:
	std::vector<double> dofDisp; 
	std::vector<double> dofForce;
	std::vector<BdType> _type; // 0 for loadBd, 1 for displaceBd
	std::string labelDisp[6] = { "ux", "uy", "uz", 
								  "rotx","roty", "rotz" };
	std::string labelLoad[6] = { "fx", "fy", "fz", 
								  "Mx","My", "Mz" };
	public:
	ZBdNode( int n )
		: ZNode(n)
	{ 
		dofDisp = std::vector<double> (6, 0.0) ;
		dofForce = std::vector<double> (6, 0.0) ;
		_type = std::vector<BdType> (6, loadBd) ;
	}
	
	ZBdNode( const ZNode* nd )
		: ZNode(*nd)
	{ 
		dofDisp = std::vector<double> (6, 0.0) ;
		dofForce = std::vector<double> (6, 0.0) ;
		_type = std::vector<BdType> (6, loadBd) ;
	}
	
	void setType(BdType type, int i) // set type of dof 
	{
		_type[i] = type;
	}

	BdType getType(int i)
	{
		return _type[i];
	}

	double  ux() const { return dofDisp[0]; }
	double  uy() const { return dofDisp[1]; }
	double  uz() const { return dofDisp[2]; }
	double  rotx() const { return dofDisp[3]; }
	double  roty() const { return dofDisp[4]; }
	double  rotz() const { return dofDisp[5]; }

	double & ux()  { return dofDisp[0]; }
	double & uy()  { return dofDisp[1]; }
	double & uz()  { return dofDisp[2]; }
	double & rotx()  { return dofDisp[3]; }
	double & roty()  { return dofDisp[4]; }
	double & rotz()  { return dofDisp[5]; }
	// return pointer , so that their element can be set directly.
	double & fx() { return dofForce[0]; }
	double & fy() { return dofForce[1]; }
	double & fz() { return dofForce[2]; }
	double & Mx() { return dofForce[3]; }
	double & My() { return dofForce[4]; }
	double & Mz() { return dofForce[5]; }

	double getDofDisp( int ind )
	{ return dofDisp[ind]; }
	double getDofForce( int ind )
	{ return dofForce[ind]; }

	// set the dofDisp and dofForce directly
	void setDofDisp( std::vector<double> d) 
	{
		dofDisp = d; 
	}

	void setDofDisp(  int ind , double d) 
	{
		dofDisp[ind] = d; 
		_type[ind] = displaceBd;
	}

	void setDofForce( std::vector<double > f) 
	{
		dofForce = f; 
	}

	void allFixed() // set all dof of this node to 0
	{
		for ( int i=0; i<6; i++)
		{
			dofDisp[i] = 0;
			_type[i] = displaceBd;
		}
	}

	void output()
	{
		std::cout<<getNum()<<"#BdNode :";
		for (int i=0; i<6; i++)
		{
			if ( _type[i] == displaceBd ) 
				std::cout<<labelDisp[i]<<"="<< dofDisp[i]<<" " ;
			else if ( abs( dofForce[i] ) > 1.0e-6 )
				std::cout<<labelLoad[i]<<"="<< dofForce[i]<<" " ;
		}
		std::cout<<std::endl;
	}
};

// boundary face class, mainly for Neumman condition, so please difine Dirichlet 
// condition on ZBdnode directly.
class ZBdFace2D 
{
	double L;
	double _grad; // const gradient in this boundary face
	ZNode *_n1;
	ZNode *_n2;
	public:
	ZBdFace2D(){} // empty constructor
//	ZBdFace2D(int n1, int n2)
//	{
//	}

	ZBdFace2D(ZNode *n1, ZNode *n2)
	{
		_n1 = n1;
		_n2 = n2;
		L = std::sqrt(
			std::pow((_n1->x() - _n2->x()),2) + 
			std::pow((_n1->y() - _n2->y()),2) + 
			std::pow((_n1->z() - _n2->z()),2));
	}

	ZBdFace2D(ZNode *n1, ZNode *n2, double grad) : ZBdFace2D( n1, n2)
	{
		_grad = grad;
	}

	void setGrad(double grad)
	{
		_grad = grad;
	}

	void disp()
	{
		std::cout<<_n1->getNum()<<" "<<_n2->getNum()<< " " <<_grad <<" "<<L;
	}

	ZNode* getNode1()
	{
		return _n1;
	}

	ZNode* getNode2()
	{
		return _n2;
	}

	double getLength()
	{	return L; }

	double getGrad()
	{	return _grad; }
};

// ZBdFace3D is a real face, it almost 3 node triangle or 4 node polygon
// you can do a switch the calculate
class ZBdFace3D;

ZBdNode* getByTag( const std::set<ZBdNode*, ZEntityLessThen> & list, int n);

} // end namespace ZMesh
#endif
