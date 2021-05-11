#ifndef _MAINCABLE_H
#define _MAINCABLE_H
#include <string>
#include "feMesh.h"

class Hanger
{
	public :
		int num; // # of the node connected to in main cable
		double deck_y; // y coordinate at deck
		double deck_z; // z coordinate at deck
		double Py; // force at y direction
		Hanger( int _num, double _deck_y, double _deck_z, double _Py )
		{
			num = _num;
			deck_y = _deck_y;
			deck_z = _deck_z;
			Py = _Py;
		}

		Hanger( const Hanger & h )
		{
			num = h.num;
			deck_y = h.deck_y;
			deck_z = h.deck_z;
			Py = h.Py;
		}

};

class MainCable
{
	public: 
		MainCable()
		{}

		void initialize(std::string configName);

		void checkInfo();

		void guessH0();

		void guessH0( double H0 );

		void solveY(); // solve when H0 is known, solve a non-linear eqation about y ,  2D main cable

		void solveY(double H); // solve when H0=H is known, solve a non-linear eqation about y, 2D main cable

		void solveYZ(); // solve when H0 is known, solve a non-linear eqation about y and z, that is 3D main cable

		void solveYZ(double H); // solve when H0=H is known, solve a non-linear eqation about y and z

		void raiseDeck(double); 

		void solveH0(); // solve when H0 is unknown, iterate to get a suitable H0 for a given slag f.

		double middleSpanY()
		{
			int middlePointDof = getByTag(this->nodes, middlePointId)->dof();
			return my(middlePointDof,0);
		}

		double middleSpanZ()
		{
			int middlePointDof = getByTag(this->nodes, middlePointId)->dof();
			return mz(middlePointDof,0);
		}

		amat::Matrix<double>& getSolutionY();
		amat::Matrix<double>& getSolutionZ();
		void printStatus(std::string fname);
		void printAPDL(std::string fname);
		double getH0() 
		{
			return H0;
		}

	private :
		int dim=3;
		double q; // weight per length of main cable
		double q_hanger; // weight per length of hanger
		double H0; // axis force at x dir
		int middlePointId;
		double middlePointY;

		std::set<ZMesh::ZNode*,ZMesh::ZEntityLessThen> nodes;
		std::set<ZMesh::ZElem*,ZMesh::ZEntityLessThen> elems;
		std::vector<Hanger> hangers;
		std::vector<int> fixedPoints;
		//nodeField u(2);
		amat::Matrix<double> my;
		amat::Matrix<double> mz;
};

std::string removeComments(std::string prgm) ;

template <class T> 
std::istringstream & operator>>(std::istringstream & is,std::vector<T> &v)
{
	char c;
	T d;
	is>>c;
	std::string dataStr ;
	std::string line;
	std::getline(is, dataStr, '}');
	std::istringstream dataStrStream( dataStr );
	while( dataStrStream>>d )
	{
		v.push_back(d);
	}
	std::getline(is, line, ';'); // read to end of this (key, value)
	return is;
}

#endif
