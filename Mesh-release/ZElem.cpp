#include "ZElem.h"
#include "feMesh.h"

namespace ZMesh {
/// \endcond
int ZElem::defaultMat = 1;
int ZElem::defaultReal = 1;
feMesh *ZElem::mesh;

GaussP::GaussP(int n)
{
	_n = n;
	switch (n)
	{
		case 1:
			_W.push_back(2.0);
			_C.push_back(0.0);
			break;
		case 2:
			_W.push_back(1.0); 
			_W.push_back(1.0);
			_C.push_back(-0.577350269189626); 
			_C.push_back(0.577350269189626);
			break;
		case 3:
			_W.push_back(0.5555555555555); 
			_W.push_back(0.8888888888888);
			_W.push_back(0.5555555555555);
			_C.push_back(-0.774596669241483); 
			_C.push_back(0.0); 
			_C.push_back(0.774596669241483);
			break;
		case 4:
			_W.push_back(0.347854845137454);
			_W.push_back(0.652145154862546);
			_W.push_back(0.652145154862546);
			_W.push_back(0.347854845137454);
			_C.push_back(-0.861136311594053);
			_C.push_back(-0.339981043584856);
			_C.push_back(0.339981043584856);
			_C.push_back(0.861136311594053);
			break;
		default :
			std::cout<<"Error: Add higher order Gauss";
	}
}

ZNode * ZElem::operator[] ( int n)
{
	if ( n>=_nodes.size() )
	{
		std::cout<< " The index of local node is too large!\n";
		return NULL;
	}
	return _nodes[n];
}

ZNode * ZElem::getNode( int n)
{
	if ( n>=_nodes.size() )
	{
		std::cout<< " The index of local node is too large!\n";
		return NULL;
	}
	return _nodes[n];
}

void ZElem::setMatByTag (const std::set<Mat*, ZEntityLessThen> &mlist, int n ) 
{
	_mat = getByTag( mlist ,n) ;
}

ZElem::ZElem( int n) : ZEntity(n)
{
	setRealByTag(mesh->getRealList(), defaultReal);
	setMatByTag(mesh->getMatList(), defaultMat);
}

ZElem::ZElem( int n , int type) : ZEntity(n)
{	
	_type = type ;
	setRealByTag(mesh->getRealList(), defaultReal);
	setMatByTag(mesh->getMatList(), defaultMat);
}

void ZElem::setNodes( int inds[], 
			   int nNodePerElem,
			   const std::set<ZNode*, ZEntityLessThen> & nList)
{
	_nNode = nNodePerElem;
	if ( _nodes.size() != 0)
		_nodes.clear();
	for (int i=0; i<nNodePerElem; i++)
	{
		ZNode tmp(inds[i]);
		std::set<ZNode*, ZEntityLessThen>::const_iterator it = nList.find(&tmp);
		if(it != nList.end())
		{
			_nodes.push_back( *it );
		}
		else
		{
			std::cout<<" Node Number Not Found ! \n";
			return ;
		}
	}
}
	
// re-edit on 2019/03/24
std::vector<int> ZElem::getAllNodeDofs()
{
	std::vector<int> dofs;
	for( ZNode* np : *this )
		dofs.push_back( np->dof() );
	return dofs;
}


std::ostream & operator<<( std::ostream & os, ZElem &e)
{ 
	os <<e._num << "#E nodes: [";
	for (int i=0; i<e._nodes.size(); i++)
		os<<e._nodes[i]->getNum()<<" ";
	os<< "],"; 
	os  << " elemType : " << e._type <<","
		<< " realTag : "<< e._real->getNum()<<","
		<< " matTag : "<<e._mat->getNum() <<"\n";
	return os;
}

ZElem* getByTag( std::set<ZElem*, ZEntityLessThen> & list, int n)
{
	ZEntity tmp(n);
	std::set<ZElem*, ZEntityLessThen>::const_iterator it = list.find((ZElem *)&tmp);
	if(it != list.end())
		return *it;
	else
	{
		std::cout<<" Element list have not been filled!\n";
		return 0;
	}
}
	
double dist( ZNode *np1, ZNode *np2 )
{
	return sqrt( (np1->x() - np2->x())*(np1->x() - np2->x()) +
				 (np1->y() - np2->y())*(np1->y() - np2->y()) +
				 (np1->z() - np2->z())*(np1->z() - np2->z()) );
}

// the distance of node np1 to the line of (np2,np3)
double dist( ZNode *np1, ZNode *np2 , ZNode *np3)
{
	double dr[3];
	double lambda;
	lambda = 
		( // (r3 - r2).(r1-r2)
		 ( np3->x() - np2->x())*(np1->x() - np2->x()) + 
		 ( np3->y() - np2->y())*(np1->y() - np2->y()) + 
		 ( np3->z() - np2->z())*(np1->z() - np2->z())
		)/
		( // (r3 - r2).(r3-r2)
		 ( np3->x() - np2->x())*(np3->x() - np2->x()) + 
		 ( np3->y() - np2->y())*(np3->y() - np2->y()) + 
		 ( np3->z() - np2->z())*(np3->z() - np2->z())
		);
	dr[0] = np1->x() - np2->x() - lambda*( np3->x() - np2->x());
	dr[1] = np1->y() - np2->y() - lambda*( np3->y() - np2->y());
	dr[2] = np1->z() - np2->z() - lambda*( np3->z() - np2->z());
	return sqrt(  dr[0]*dr[0] + dr[1]*dr[1] +  dr[2]*dr[2]);
}


// the distance of node np1 to the surface of (np2,np3,np4)
double dist( ZNode *np1, ZNode *np2 , ZNode *np3, ZNode *np4)
{
	double r1[3] = {np3->x() - np2->x(),
					np3->y() - np2->y(),
					np3->z() - np2->z() };
	double r2[3] = {np4->x() - np2->x(),
					np4->y() - np2->y(),
					np4->z() - np2->z() };
	double n[3] = { r1[1]*r2[2] - r2[1]*r1[2],
					r1[2]*r2[0] - r1[0]*r2[2],
					r1[0]*r2[1] - r2[0]*r1[1] };
	double norm_n = sqrt(  n[0]*n[0] + n[1]*n[1] +  n[2]*n[2]);
	for(int i=0; i<3; i++) n[i] /= norm_n;
	return abs(
			(np1->x()-np2->x())*n[0] +
			(np1->y()-np2->y())*n[1] +
			(np1->z()-np2->z())*n[2] );
}


	double det( double A[3][3] ) 
	{
		double d =
	 A[0][0]*A[1][1]*A[2][2] - 
				A[0][0]*A[1][2]*A[2][1] -
				A[0][1]*A[1][0]*A[2][2] +
				A[0][1]*A[1][2]*A[2][0] +
				A[0][2]*A[1][0]*A[2][1] -
				A[0][2]*A[1][1]*A[2][0];
		return d;
	}


} // end namespace ZMesh
