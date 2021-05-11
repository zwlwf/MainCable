#include "MainCable.h"
#include "alinalg.hpp"
#include "vtkFile.h"

using namespace ZMesh;
using namespace std;

void MainCable::initialize(std::string configName)
{ 
	std::ostringstream configStr;
	std::ifstream is(configName.c_str());
	configStr<<is.rdbuf();
	is.close();
	std::string cleanString = removeComments(configStr.str()); 
	std::istringstream cleanStringStream(cleanString);
	std::vector<double> nodeData;
	std::vector<double> hangerData;
	std::vector<double> middlePointData;
	std::vector<int> elemData;
	std::string itemName, line;
	char c;
	while ( cleanStringStream>> itemName )
	{
		std::cout<<"read a item named "<<itemName<<std::endl;
		if ( itemName == string("nodes") )
		{
			cleanStringStream>>nodeData;
		} else if ( itemName == string("elems") )
		{
			cleanStringStream>>elemData;
		} else if ( itemName == string("hanger") )
		{
			cleanStringStream>>hangerData;
		} else if ( itemName == string("q") )
		{
			cleanStringStream>>this->q;
			std::getline(cleanStringStream, line, ';');
		} else if ( itemName == string("fixedPoints") )
		{
			cleanStringStream>>fixedPoints;
		} else if ( itemName == string("q_hanger") )
		{
			cleanStringStream>>this->q_hanger;
			std::getline(cleanStringStream, line, ';');
		} else if ( itemName == string("middlePoint") )
		{
			cleanStringStream>>middlePointData;
		} else if ( itemName == string("dimension") )
		{
			cleanStringStream>>this->dim;
			std::getline(cleanStringStream, line, ';');
		}
		else
			std::getline(cleanStringStream, line, ';');
	}

	int nodeCol = 4;
	int nNode = nodeData.size()/nodeCol; // only given x is meaningful, y z as initial guess
	int tag;
   	double x,y,z;
	for(int i=0; i<nNode; i++)
	{
		tag = round(nodeData[i*nodeCol+0]);
		x = nodeData[i*nodeCol+1];
		y = nodeData[i*nodeCol+2];
		z = nodeData[i*nodeCol+3];
		this->nodes.insert( new ZNode( tag, x, y, z ) );
	}

	my.setup(nNode,1);
	mz.setup(nNode,1);

	int idof=0;
	for( ZNode* np : nodes)
	{
		my(idof,0) = np->y();
		mz(idof,0) = np->z();
		np->setDof(idof++);
	}

	int elemCol = 3;
	int nElem = elemData.size()/elemCol;

	feMesh tmp(2);
	tmp.createDefaultMat();
	tmp.createDefaultReal();
	tmp.setNodeList( nodes);
	ZElem::mesh = &tmp;

	tag = 1;
	for(int i=0; i<nElem; i++)
	{
		this->elems.insert(new ZLine2( elemData[i*elemCol], elemData[i*elemCol+1], elemData[i*elemCol+2]) );
	}

	int hangerCol = 4;
	int nHinger = hangerData.size()/hangerCol;
	hangers.reserve(nHinger);
	for(int i=0; i<nHinger; i++)
	{
		hangers.push_back( Hanger( round(hangerData[i*hangerCol]),
					hangerData[i*hangerCol+2], 
					hangerData[i*hangerCol+3], 
					hangerData[i*hangerCol+1] ) );
	}

	if ( middlePointData.size() )
	{
		this->middlePointId = round( middlePointData[0] );
		this->middlePointY = middlePointData[1] ;
	}

}

void MainCable::checkInfo()
{
	std::cout<<" ------------ check the information of main cable -----------  "<<std::endl;
	std::cout<<"There are "<<std::endl
			 <<this->nodes.size()<<" nodes ,"<<std::endl
			 <<this->elems.size()<<" elements,"<<std::endl
  			 <<"q = " << q<< "N/m" << std::endl
			 <<"H0( guessed ) = "<<H0<<" N"<<std::endl;
}

void MainCable::printStatus(std::string fname)
{
	vtkFile os( fname.c_str() );

	// output to format of vtk
	for( ZNode* np : this->nodes )
	{
		np->y() = my(np->dof(),0);
		np->z() = mz(np->dof(),0);
	}

	os.printHeader();
	os.printNodes( this->nodes);
	os.printElems( this->elems);
}

void MainCable::printAPDL(std::string fname)
{
	ofstream os( fname.c_str() );
	if ( !os ) 
	{
		std::cout<<"file "<<fname<<" can not be opened!"<<std::endl;
		return ;
	}
	// output to format of vtk
	for( ZNode* np : this->nodes )
	{
		os<<"n,"<<np->getNum()<<","<< np->x() <<","<< my(np->dof(),0) <<","<< mz(np->dof(),0)<< std::endl;
	}

	double Area = this->q/7850/9.8;
	double ds, dx;

	for( ZElem* ep : elems )
	{
		// r,1,Area,epsilon
		ds = sqrt( 
		  pow( ep->getNode(1)->x() - ep->getNode(0)->x(),2) + 
		  pow( my(ep->getNode(1)->dof(),0)-my(ep->getNode(0)->dof(),0), 2) + 
		  pow( mz(ep->getNode(1)->dof(),0)-mz(ep->getNode(0)->dof(),0), 2) );
		dx = fabs(ep->getNode(1)->x() - ep->getNode(0)->x());
	    	
		os<<"r,"<<ep->getNum()<<","<< Area <<","<< this->H0*ds/dx/1.95e11/Area<<std::endl;
	}
	
	
	os.close();
}

void MainCable::guessH0( double H0)
{
	this->H0 = H0;
}

void MainCable::guessH0( )
{
	// add on 2019/02/19
	// 确定垂度点两侧IP点
	double sumF, L, f;
	sumF=0;
	multimap<double, int> x2fixedPoints;
	double xmid = getByTag(this->nodes, middlePointId)->x();
	for(int ip : this->fixedPoints )
	{
		x2fixedPoints.insert( {fabs(getByTag(this->nodes, ip)->x()-xmid), ip} );
	}
	multimap<double,int>::iterator it;
	it = x2fixedPoints.begin();
	ZMesh::ZNode* n1 = getByTag(this->nodes, it->second);
	it++;
	ZMesh::ZNode* n2 = getByTag(this->nodes, it->second);
	L = fabs( n1->x() - n2->x() );
	f = ( n1->y() + n2->y() )/2 - middlePointY;
	sumF += this->q*2*sqrt(L*L/4+f*f);
	for( Hanger hh : hangers )
	{
		double xh = getByTag(this->nodes, hh.num)->x();
		if( ( xh - n1->x() )*(xh - n2->x()) < 0 )
			sumF += hh.Py;
	}

	this->H0 = sumF*L/8.0/f;
}

void MainCable::solveY(double HH0)
{
	int itr = 0;
	int N = nodes.size();
	amat::Matrix<double> K(N, N), rhs(N,1);
	bool converged = false;
	double dL;
	double dphi1,dphi2;
	amat::Matrix<double> du(N,1);
	amat::Matrix<double> u(N,1);
	// set initial value of y for node coordinate
	for ( ZNode* np : nodes )
	{
		u(np->dof()) = my(np->dof(),0);
	}

#ifdef DEBUG
	cout<<"H = "<<HH0<<endl;
#endif
	double omega = 1.0;
	while ( !converged )
	{
		++itr;
		K.zeros();
		rhs.zeros();
		double J_i;
		for( ZElem * ep : elems )
		{
			ZNode* n0 = ep->getNode(0);
			ZNode* n1 = ep->getNode(1);
			int dof1 = n0->dof();
			int dof2 = n1->dof();
			dL = fabs(ep->getNode(1)->x() - ep->getNode(0)->x());
		    J_i= dL/2;
			dphi1 =-1/dL;
			dphi2 = 1/dL; 
			double lambda=0.0;
			double lambda2 = dphi1*my(dof1) + dphi2*my(dof2);
			// Kyy_ij
			K(dof1,dof1) += HH0*dphi1*dphi1*dL + this->q*lambda2/sqrt(1+lambda*lambda+lambda2*lambda2)*J_i*dphi1; 
			K(dof1,dof2) += HH0*dphi1*dphi2*dL + this->q*lambda2/sqrt(1+lambda*lambda+lambda2*lambda2)*J_i*dphi2; 
			K(dof2,dof1) += HH0*dphi2*dphi1*dL + this->q*lambda2/sqrt(1+lambda*lambda+lambda2*lambda2)*J_i*dphi1; 
			K(dof2,dof2) += HH0*dphi2*dphi2*dL + this->q*lambda2/sqrt(1+lambda*lambda+lambda2*lambda2)*J_i*dphi2; 

			rhs(dof1) -= (HH0*lambda2*dL*dphi1+this->q*sqrt(1+lambda*lambda+lambda2*lambda2)*J_i ); 
			rhs(dof2) -= (HH0*lambda2*dL*dphi2+this->q*sqrt(1+lambda*lambda+lambda2*lambda2)*J_i);
		}

		for( Hanger hh : hangers )
		{
			ZNode* np = getByTag( this->nodes, hh.num );
			int dof1 = np->dof();
			rhs(dof1) -= ( hh.Py + this->q_hanger*sqrt((my(dof1)-hh.deck_y)*(my(dof1)-hh.deck_y) + (mz(dof1)-hh.deck_z)*(mz(dof1)-hh.deck_z))); // 考虑吊杆
		}

		// introduce the boundary , for fixed points, du=0 
		for(int id : fixedPoints )
		{
			ZNode* np = getByTag( this->nodes, id );
			int dof1 = np->dof();
			rhs(np->dof()) = 0; 
			for ( int i=0; i<N; i++)
			{
				K(dof1,i) = 0;
			}
			K(dof1,dof1) = 1.0;
		}

		gausselim(K,rhs,du);
		double res = du.l2norm();
		cout<<"iteration = "<<itr<<", res = "<<res<<endl;
#ifdef DEBUG
		char outFile[1024];
		int id;
		cout<<"current itr = "<<itr<<", input [id] of out>>";
		cin>>id;
		sprintf(outFile,"out_%d.vtk", id);
		this->printStatus(outFile);
#endif
		if ( res < 1.0e-6)
		{
			converged = true;
			break;
		}
		u = u + du*omega;
		for(int i=0; i<N; i++)
		{
			my(i,0) = u(i,0);
		}
	}
}

void MainCable::raiseDeck( double htmp)
{
	for ( Hanger & hh : this->hangers )
		hh.deck_y += htmp;
}

void MainCable::solveY( )
{
	this->solveY( this->H0);
}

void MainCable::solveYZ(double HH0)
{
	int itr = 0;
	int N = nodes.size();
	amat::Matrix<double> K(2*N, 2*N), rhs(2*N,1);
	bool converged = false;
	double dL;
	double dphi1, dphi2;
	amat::Matrix<double> du(2*N,1);
	amat::Matrix<double> u(2*N,1);
	// set initial value of y for node coordinate
	for ( ZNode* np : nodes )
	{
		u(2*np->dof()) = my(np->dof(),0);
		u(2*np->dof()+1) = mz(np->dof(),0);
	}

#ifdef DEBUG
	cout<<"H = "<<HH0<<endl;
#endif
	double omega = 1.0;
	while ( !converged )
	{
		++itr;
		K.zeros();
		rhs.zeros();
		double J_i;
		for( ZElem * ep : elems )
		{
			ZNode* n0 = ep->getNode(0);
			ZNode* n1 = ep->getNode(1);
			int dof1 = n0->dof();
			int dof2 = n1->dof();
			dL = fabs(ep->getNode(1)->x() - ep->getNode(0)->x());
		    J_i= dL/2;
			dphi1 =-1/dL;
			dphi2 = 1/dL; 
			double lambda = dphi1*mz(dof1) + dphi2*mz(dof2);
			double lambda2 = dphi1*my(dof1) + dphi2*my(dof2);
			// Kyy_ij
			K(2*dof1,2*dof1) += HH0*dphi1*dphi1*dL + this->q*lambda2/sqrt(1+lambda*lambda+lambda2*lambda2)*J_i*dphi1; 
			K(2*dof1,2*dof2) += HH0*dphi1*dphi2*dL + this->q*lambda2/sqrt(1+lambda*lambda+lambda2*lambda2)*J_i*dphi2; 
			K(2*dof2,2*dof1) += HH0*dphi2*dphi1*dL + this->q*lambda2/sqrt(1+lambda*lambda+lambda2*lambda2)*J_i*dphi1; 
			K(2*dof2,2*dof2) += HH0*dphi2*dphi2*dL + this->q*lambda2/sqrt(1+lambda*lambda+lambda2*lambda2)*J_i*dphi2; 

			// Kyz_ij
			K(2*dof1,2*dof1+1) += this->q*lambda/sqrt(1+lambda*lambda+lambda2*lambda2)*J_i*dphi1; 
			K(2*dof1,2*dof2+1) += this->q*lambda/sqrt(1+lambda*lambda+lambda2*lambda2)*J_i*dphi2; 
			K(2*dof2,2*dof1+1) += this->q*lambda/sqrt(1+lambda*lambda+lambda2*lambda2)*J_i*dphi1; 
			K(2*dof2,2*dof2+1) += this->q*lambda/sqrt(1+lambda*lambda+lambda2*lambda2)*J_i*dphi2; 
			rhs(2*dof1) -= (HH0*lambda2*dL*dphi1+this->q*sqrt(1+lambda*lambda+lambda2*lambda2)*J_i ); 
			//rhs(2*i) -= (HH0*lambda2*dL*dphi1+this->q*sqrt(1+lambda*lambda+lambda2*lambda2)*J_i + Fy(i)/2 );
			rhs(2*dof2) -= (HH0*lambda2*dL*dphi2+this->q*sqrt(1+lambda*lambda+lambda2*lambda2)*J_i);


			// Kzz_ij
			K(2*dof1+1,2*dof1+1) += HH0*dphi1*dphi1*dL;
			K(2*dof1+1,2*dof2+1) += HH0*dphi1*dphi2*dL;
			K(2*dof2+1,2*dof1+1) += HH0*dphi2*dphi1*dL;
			K(2*dof2+1,2*dof2+1) += HH0*dphi2*dphi2*dL;
			rhs(2*dof1+1) -= (HH0*lambda*dL*dphi1);
			rhs(2*dof2+1) -= (HH0*lambda*dL*dphi2);

		}

		for( Hanger hh : hangers )
		{
			ZNode* np = getByTag( this->nodes, hh.num );
			int dof1 = np->dof();
			rhs(2*dof1) -= ( hh.Py + this->q_hanger*sqrt((my(dof1)-hh.deck_y)*(my(dof1)-hh.deck_y) + (mz(dof1)-hh.deck_z)*(mz(dof1)-hh.deck_z))); // 考虑吊杆
			rhs(2*dof1+1) -= ( hh.Py*(mz(dof1)-hh.deck_z)/(my(dof1)-hh.deck_y) ); // 考虑吊杆
			// Kzy_ij
			K(2*dof1+1,2*dof1) -= hh.Py*(mz(dof1) - hh.deck_z)/(my(dof1)-hh.deck_y)/(my(dof1)-hh.deck_y) ;
			// Kzz_ij
			K(2*dof1+1, 2*dof1+1) += hh.Py/(my(dof1) - hh.deck_y);

			K(2*dof1, 2*dof1) += this->q_hanger*(my(dof1)-hh.deck_y)/sqrt((my(dof1)-hh.deck_y)*(my(dof1)-hh.deck_y) + (mz(dof1)-hh.deck_z)*(mz(dof1)-hh.deck_z));
			K(2*dof1, 2*dof1+1) += this->q_hanger*(mz(dof1)-hh.deck_z)/sqrt((my(dof1)-hh.deck_y)*(my(dof1)-hh.deck_y) + (mz(dof1)-hh.deck_z)*(mz(dof1)-hh.deck_z));
		}

		// introduce the boundary , for fixed points, du=0 
		for(int id : fixedPoints )
		{
			ZNode* np = getByTag( this->nodes, id );
			int dof1 = np->dof();
			rhs(2*np->dof()) = 0; 
			rhs(2*np->dof()+1) = 0; 
			for ( int i=0; i<2*N; i++)
			{
				K(2*dof1,i) = 0;
				K(2*dof1+1,i) = 0;
			}
			K(2*dof1,2*dof1) = 1.0;
			K(2*dof1+1,2*dof1+1) = 1.0;
		}

		gausselim(K,rhs,du);
		double res = du.l2norm();
		cout<<"iteration = "<<itr<<", res = "<<res<<endl;
#ifdef DEBUG
		char outFile[1024];
		int id;
		cout<<"current itr = "<<itr<<", input [id] of out>>";
		cin>>id;
		sprintf(outFile,"out_%d.vtk", id);
		this->printStatus(outFile);
#endif
		if ( res < 1.0e-6)
		{
			converged = true;
			break;
		}
		u = u + du*omega;
		for(int i=0; i<N; i++)
		{
			my(i,0) = u(2*i);
			mz(i,0) = u(2*i+1);
		}
	}
	 
}

void MainCable::solveYZ()
{
	this->solveYZ( this->H0 );
}

void MainCable::solveH0()
{
	int middlePointDof = getByTag(this->nodes, middlePointId)->dof();

	double f0, f1, f2;
	double H1, H2;
	H1 = this->H0*1.1;

	if ( this->dim == 2 )
		this->solveY();
	else
		this->solveYZ();

	f0 = my(middlePointDof) - middlePointY;
#ifdef DEBUG
	cout<<"G = "<<f0<<endl;
#endif

	if ( this->dim == 2 )
		this->solveY(H1);
	else
		this->solveYZ(H1);

	f1 = my(middlePointDof) - middlePointY;
#ifdef DEBUG
	cout<<"G = "<<f1<<endl;
#endif

	H2 = (H0*f1 - H1*f0)/(f1 - f0);

	if ( this->dim == 2 )
		this->solveY(H2);
	else
		this->solveYZ( H2 );
	f2 = my(middlePointDof) - middlePointY;
#ifdef DEBUG
	cout<<"G = "<<f2<<endl;
#endif

	while( abs(f2)>1.0e-4 )
	{
		H0 = H1;
		H1 = H2;
		f0 = f1;
		f1 = f2;
		H2 = (H0*f1 - H1*f0)/(f1 - f0);
		if ( this->dim == 2 )
			this->solveY(H2);
		else
			this->solveYZ(H2);
		f2 = my(middlePointDof) - middlePointY;
#ifdef DEBUG
	cout<<"G = "<<f2<<endl;
#endif
	}
	this->H0 = H2;

	std::cout<<"H0 = "<<H0<<std::endl;
}

amat::Matrix<double> & MainCable::getSolutionZ()
{ return mz; }

amat::Matrix<double> & MainCable::getSolutionY()
{ return my; }
