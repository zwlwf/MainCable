#include "feMesh.h"

namespace ZMesh{

void feMesh::readTest(std::string mfile, int elemType)
{
	double E; // Yang's Module
	// here there are only one material
	Mat *matFactory;
	std::ostream & WRITE_IO = std::cout;
	std::ifstream READ_IN( mfile.c_str() );
	if ( !READ_IN ) 
	{
		std::cout<<" file fo "<<mfile <<" can not open !"<<std::endl;
		return ;
	}
	READ_IN >> nNodes >> nElems >> nLoadNodes >> nBdNodes >> E;
	matFactory = new Mat(1); // matList create temprary here
	matFactory->setE( E );
	matList.insert( matFactory );
	WRITE_IO << "\n\n Internal Data \n\n" 
			 << " Number of Nodes         :" << std::setw(5) << nNodes << "\n" 
			 << " Number of Elements      :" << std::setw(5)<< nElems << "\n"
			 << " Number of Loaded Nodes  :" << std::setw(5) << nLoadNodes << "\n"
			 << " Number of Support Nodes :" << std::setw(5) << nBdNodes << "\n" 
			 << " Modulus of Elasticity :" << std::setw(15) << std::setprecision(0) << E << "\n\n";
	for ( int i=0; i< nNodes; i++)
	{
		int NUM; double X,Y,Z;
		X = Y = Z = 0.0;
		READ_IN >> NUM >> X >> Y;
		if ( nDim == 3 ) READ_IN >> Z;
		_nlist.insert( new ZNode(NUM,X,Y,Z) ); // the address located for _nlist may not continuous
	}
	int nodeLocal[20];// max node num = 20
	nNodePerElem = 2; // number of node per element, for both link and beam , 2D and 3D case
	ZElem *p;
	for ( int i=0; i<nElems; i++)
	{
		//int elemType = 1; // element type, this should be input in a more general format.
		int NUM;
		READ_IN >> NUM;
		for ( int j=0; j<nNodePerElem; j++)
			READ_IN >> nodeLocal[j];
		double W[2];
		switch ( elemType)
		{
			case 1: // 2D link
				{
					ZLink2D * tmp;
					tmp = new ZLink2D(NUM);
					p = tmp;
					READ_IN >> W[0] ; // read Area 
					break;
				}
			case 2:
				{
					ZBeam2D * tmp;
					tmp = new ZBeam2D(NUM);
					p = tmp;
					READ_IN >> W[0] >>W[1]; // read Area and I of beam 
					break;
				}
			default :
				std::cout<<" Such elemType not implemented! "<<std::endl;
		}
		p->setNodes(nodeLocal, nNodePerElem, _nlist);
		p->setRealByTag(realList,1);
		p->setMatByTag(matList,1);
		_elist.insert(p);
	}
	// read boundary

	for (int i=0; i<nLoadNodes; i++)
	{
		int NUM;
		double W[3];
		switch ( elemType )
		{
			case 1:
				{
					READ_IN >> NUM >> W[0] >> W[1];
					ZBdNode & tmp = *( new ZBdNode(NUM) ); 
					tmp.fx() = W[0];
					tmp.fy() = W[1];
					boundaryNodes.push_back(&tmp);
					break;
				}
			case 2:
				{
					READ_IN >> NUM >> W[0] >> W[1] >> W[2];
					ZBdNode & tmp = *( new ZBdNode(NUM) ); 
					tmp.fx() = W[0];
					tmp.fy() = W[1];
					tmp.Mx() = W[2];
					boundaryNodes.push_back(&tmp);
					break;
				}
			default :
				std::cout<<"Such Element not implemented !\n";
				break;
		}
	}	 

	for (int i=0; i<nBdNodes; i++)
	{
		int NUM;
		double W[3];
		int flag[3];
		switch ( elemType )
		{
			case 1:
				{
					READ_IN >> NUM >> flag[0] >> flag[1] >> W[0] >> W[1];
					ZBdNode & tmp = *( new ZBdNode(NUM) ); 
					for (int j=0; j<2; j++)
					{
						if ( flag[j] == 0 ) tmp.setDofDisp(j, W[j]);
					}
					boundaryNodes.push_back(&tmp);
					break;
				}
			case 2:
				{
					READ_IN >> NUM >> flag[0] >> flag[1] >> flag[2] >> W[0] >> W[1] >>W[2];
					ZBdNode & tmp = *( new ZBdNode(NUM) ); 
					for (int j=0; j<2; j++)
					{
						if ( flag[j] == 0 ) tmp.setDofDisp(j, W[j]);
					}
					if ( flag[2] == 0 ) tmp.setDofDisp(3, W[2]); // rotx
					boundaryNodes.push_back(&tmp);
					break;
				}
			default :
				std::cout<<"Such Element not implemented !\n";
				break;
		}
	}
	createNodeDof(); // create nodeDof for mesh
	READ_IN.close();
}

// only accept Tri2D type, old one, no physical entity
void feMesh::readGmsh(std::string mfile)
{
	nDim = 2;
    std::cout << "feMesh: readGmsh(): Reading mesh file...\n";
    int dum; double dummy; std::string dums; char ch;
	int nPoint; // may not be nNodes

    std::ifstream infile(mfile.c_str());
	if ( !infile ) 
	{
		std::cout<<" file fo "<<mfile <<" can not open !"<<std::endl;
		return ;
	}
    for(int i = 0; i < 4; i++)		//skip 4 lines
        do
            ch = infile.get();
        while(ch != '\n');
    infile >> nPoint; // point may be too many, so not write to nNodes immediently
    // read coords of points
	std::set<ZNode*, ZEntityLessThen> _nlistTmp;    
    for(int i = 0; i < nPoint; i++)
    {
		double X,Y,Z;
        if ( (infile >> dum).fail() ) break; // read not a int ,dum is the node number
         infile >> X >>Y >>Z;
		_nlistTmp.insert( new ZNode(dum,X,Y,Z) );
    }
	if ( infile.fail() )
	 {
		infile.clear();
	}
    infile >> dums;
    infile >> dums;
    int nelm, elmtype, nBoundaryTags, ntags, nElemTag;
	char c;
    /// elmtype is the standard element type in the Gmsh 2 mesh format - of either faces or elements
    nElemTag = 0; 
    infile >> nelm;
    int elms[20];
    //std::cout << "feMesh: readGmsh2(): Total number of elms is " << nelm << std::endl;

	ZElem *p;
    for(int i = 0; i < nelm; i++)
    {
		bool isElemOk = true; // is all node of this elem in the _nlistTmp
        infile >> dum;
        infile >> elmtype;
		if ( elmtype == 2 ) // input only Tri2D element
		{
			nNodePerElem = 3;
			infile >> ntags;
			for(int j = 0; j < ntags; j++)
				infile >> elms[j+nNodePerElem];		// get tags
			for(int j = 0; j < nNodePerElem; j++)
			{
				infile >> elms[j];			// get node numbers
				// check whether all node is in the _nlistTmp input
				if ( getByTag( _nlistTmp, elms[j]) == NULL )
				{
					isElemOk = false;
				}
			}
			if ( isElemOk )
			{
				ZTri2D *tmp;
				tmp = new ZTri2D(dum);
				p = tmp;
				p->setNodes(elms, nNodePerElem, _nlistTmp);
				p->setRealByTag(realList,1);
				p->setMatByTag(matList,1);
				_elist.insert(p);
			}
		}
		else
		{
			while ( (infile.get(c)).good() )
			{
				if ( c == '\n' ) break; // drop this line
			}
		}

	}
	// insert just used node to _nlist.
	for ( 
			std::set<ZElem*>::const_iterator ite = _elist.begin();
			ite != _elist.end();
			++ite
		)
	{
		ZElem* p1 = *ite;
		for ( int i=0; i< p1->getNNode(); i++)
			_nlist.insert( (*p1)[i] );
	}

	infile.close();
	nNodes = _nlist.size();
	nElems = _elist.size();
	std::cout<<" Number of Nodes = "<<nNodes<<" read,\n"
			 <<" Number of Elements = " <<nElems<<std::endl;

} 

void feMesh::readGmsh2(std::string mfile)
{
	std::cout << "readGmsh2(): Reading mesh file...\n";
	std::ifstream is(mfile.c_str());
	if ( !is )
	{
		std::cout<<"Error: Gmsh file not exist!"<<std::endl;
		exit(-1);
	}
	std::string str1,str2;
	nDim=0;
	int id,idim;
	int nPhy,nEle;
	int* phyDims;
	std::string* phyNames;
	while ( is>>str1)
	{
		if (str1 == "$PhysicalNames")
		{
			is>>nPhy;
			phyDims = new int[nPhy+1];
			phyNames = new std::string[nPhy+1];
			for (int i=0; i<nPhy; i++)
			{
				is>> idim >> id >> str2; 
				std::cout<<" dimension = "<<idim <<", id =" <<id<<", name = "<<str2.substr(1,str2.size()-2)<<std::endl;
				phyDims[id] = idim; 
				phyNames[id] = str2.substr(1,str2.size()-2); 
				if (idim>nDim) nDim = idim;
			}
			break;
		}
	}

	if ( is.eof() )
	{
		std::cout<<"Error! ---> Physicals not exist !"<<std::endl;
		return ;
	}
	std::unordered_map<int , ZPatch2 > zpatches;
	// read node
	std::set<ZNode*, ZEntityLessThen > nlist;
	std::set<ZElem*, ZEntityLessThen > elist;
	while ( is>> str1 )
	{
		if (str1 == "$Nodes")
		{
			is>>nNodes;
			double d1,d2,d3;
			for (int i=0; i<nNodes; i++)
			{
				is>> id ;
				is>> d1 >> d2>>d3;
				nlist.insert( new ZNode(id, d1,d2,d3) );
			}
			std::cout<< nNodes << " Nodes input "<<std::endl;
			break;
		}
	}
    setNodeList(nlist); // node must be filled at first

	if ( is.eof() )
	{
		std::cout<<"Error! ---> Nodes not exist !"<<std::endl;
		return ;
	}
	// read elems
	nElems = 0;
	while ( is>> str1 )
	{
		if (str1 == "$Elements")
		{
			is>>nEle;
			int elmtype,nTags,PhyId,geoId,otherId;
			int d1,d2,d3,d4,d5,d6,d7,d8;
			for (int i=0; i<nEle; i++)
			{
				is>> id ;

				is >> elmtype >> nTags >> PhyId;
				elmtype = gmsh2FEM_type[elmtype];
				is >> geoId; // read out geometry elementary id
				
				for ( int j=0; j<nTags-2; j++)
				{
					is>>otherId;
				}
				if ( phyDims[PhyId] < nDim ) // boundary patch
				{
					switch ( elmtype )
					{
					   case(13): // 1-node point
					   		is>>d1;
							zpatches[PhyId].add( new ZPoint( id, d1) );
							break;
					   case(9): // linear edge
							is>>d1>>d2;
							zpatches[PhyId].add( new ZLine2( id, d1, d2) );
							break;
						case(10): // quadratic edge
							is>>d1>>d2>>d3;
							zpatches[PhyId].add( new ZLine3Q( id, d1, d3,d2) );
							break;
						case(FEM_ZTri2D):
						case (FEM_ZTri2D1V):
							std::cout<<"warning : ZTri3D instead here ";
						case(12): // linear triangles
							is>>d1>>d2>>d3;
							zpatches[PhyId].add( new ZTri3D( id, d1, d2,d3) );
							break;
						case(8):	// linear quads
							is>>d1>>d2>>d3>>d4;
							zpatches[PhyId].add( new ZQuad4Q( id, d1, d2,d3,d4 ) );
							break;
						/*
						case(9):	// quadratic triangles
							std::cout<<"Warning:  quadratic triangles not implemented for patches\n";
							getline(is, str2);
							break;
						case(16):	// quadratic quad (8 nodes)
							std::cout<<"Warning:  quadratic quad (8 nodes) not implemented for patches\n";
							getline(is, str2);
							break;
						case(10):	// quadratic quad (9 nodes)
							std::cout<<"Warning: quadratic quad (9 nodes) not implemented for patches\n";
							getline(is, str2);
							break;
						*/
						default:
							std::cout <<"Warning: readGmsh2(): Element type not recognized for patches" << std::endl;
							getline(is, str2);
					} 
				}
				else  // inner elements
				{
					switch ( elmtype )
					{
					    case(1): 
							is>>d1>>d2;
							elist.insert( new ZLink2D(id, d1,d2) );
					 	    nElems++;
							break;
					    case(2): 
							is>>d1>>d2;
							elist.insert( new ZBeam2D(id, d1,d2) );
					 	    nElems++;
							break;
					    case(3): 
							is>>d1>>d2;
							elist.insert( new ZLink3D(id, d1,d2) );
					 	    nElems++;
							break;
						case(4): 
							is>>d1>>d2>>d3;
							elist.insert( new ZTri2D(id, d1,d2,d3) );
							nElems++;
							break;
					    case(5): 
							is>>d1>>d2;
							elist.insert( new ZBeam3D(id, d1,d2) );
					 	    nElems++;
							break;
						case(6): 
							is>>d1>>d2>>d3;
							elist.insert( new ZTri2D1V(id, d1,d2,d3) );
							nElems++;
							break;
						case(7):
							is>>d1>>d2>>d3>>d4;
							elist.insert ( new ZTet(id, d1,d2,d3,d4) );
							nElems++;
							break;
						case(8):
							is>>d1>>d2>>d3>>d4;
							elist.insert ( new ZQuad4Q(id, d1,d2,d3,d4) );
							nElems++;
							break;
					    case(9): 
							is>>d1>>d2;
							elist.insert( new ZLine2(id, d1,d2) );
					 	    nElems++;
							break;
					    case(10): 
							is>>d1>>d2>>d3;
							elist.insert( new ZLine3Q( id, d1, d3,d2) );
					 	    nElems++;
							break;
						case(11):
							is>>d1>>d2>>d3>>d4>>d5>>d6>>d7>>d8;
							elist.insert ( new ZQuad8Q(id, d1,d2,d3,d4,d5,d6,d7,d8) );
							nElems++;
							break;
						case(12):
							is>>d1>>d2>>d3;
							elist.insert( new ZTri3D(id, d1,d2,d3) );
							nElems++;
							break;
						case(13):
							is>>d1;
							elist.insert( new ZPoint(id, d1) );
							nElems++;
							break;
						default:
							std::cout << " Warning: readGmsh2(): Element type not recognized." << std::endl;
						   getline(is,str2);
						}
				}
			}
			std::cout<< nElems <<" Elements input"<<std::endl; 
			break;
		}
	}


	if ( is.eof() )
	{
		std::cout<<"Error! ---> Elements not exist !"<<std::endl;
		return ;
	}
	is.close();

	for(int i=1; i<=nPhy; i++)
	{
		if( phyDims[i] < nDim )
		{
			zpatches[i].setMesh(this);
			patches2.insert( {phyNames[i], zpatches[i]} );
			std::cout<<zpatches[i].getNFaces()<<" faces in patch "<< phyNames[i]<<" inserted!"<<std::endl;
		}
	}

	setElemList(elist);
	createNodeDof2();
}

// mainly for 2D triangle, which was created by amc, only read in node and mesh
void feMesh::readMesh(std::string pathname)
{
	// read in nodes
	std::set<ZNode*, ZEntityLessThen> nlist;
	std::set<ZNode*, ZEntityLessThen> NewNlist;
	std::set<ZElem*, ZEntityLessThen> elist;
	double d1,d2,d3;
	int i1,i2,i3,i4;
	int tag=0;
	std::ifstream is( (pathname+"/nodes").c_str() );
	if ( !is ) 
	{
		std::cout<<" file can not open !"<<std::endl;
		return ;
	}
	/* 2D point also has z coordinate in data file.
	if ( nDim == 2 )
	{
		while ( (is>>d1>>d2).good() )
		{
			nlist.insert( new ZNode(++tag, d1,d2 ) );
		}
	}
	else if ( nDim == 3)
	{
	*/
		while ( (is>>d1>>d2>>d3).good() )
		{
			nlist.insert( new ZNode(++tag, d1,d2,d3) );
		}
	//}

	tag=0;
	is.close();
	setNodeList(nlist); // node must be filled at first
	is.open( (pathname+"/elems").c_str() );
	if ( !is ) 
	{
		std::cout<<" file can not open !"<<std::endl;
		return ;
	}
	while ( (is>>i1>>i2>>i3).good() )
	{
		// note: 2D triangle nodes
		if ( nDim == 2 )
			elist.insert( new ZTri2D1V(++tag, i1,i2,i3 ) );
		else
		{
			is>>i4;
			elist.insert( new ZTet(++tag, i1,i2,i3,i4 ) );
		}
	}
	tag=0;
	is.close();

	// insert just used node to _nlist.
	for ( 
			std::set<ZElem*>::const_iterator ite = elist.begin();
			ite != elist.end();
			++ite
		)
	{
		ZElem* p1 = *ite;
		for ( int i=0; i< p1->getNNode(); i++)
			NewNlist.insert( (*p1)[i] );
	}
	setNodeList( NewNlist );
	setElemList(elist);
	createNodeDof();
}

// main for ZTri2D1V element
/* the file format for boundary note is 
 * node_1, val_1
 * ...
 * node_i, val_i
 */
void feMesh::readBoundaryNode( std::string mfile )
{
	std::ifstream is;
	std::vector<ZBdNode*> bnds;
	int i1; 
	double d1;
	is.open( mfile.c_str() );
	if ( !is ) 
	{
		std::cout<<"Boundary Node file can not be openned !"<<std::endl;
		return ;
	}
	while ( (is>>i1>>d1).good() )
	{
		ZBdNode * tmp = new ZBdNode(i1); 
		tmp->setDofDisp(0,d1); 
		bnds.push_back( tmp);
	}
	is.close();
	setBdList(bnds);
}

void feMesh::readBoundaryFace( std::string mfile )
{
	std::ifstream is( mfile.c_str() );
	if ( !is ) 
	{
		std::cout<< "File of boundary face can not be opened! "<<std::endl;
		return ;
	}
	int n1, n2;
	double grad;
	ZNode *np1, *np2;
	ZBdFace2D *nf;
	while ( (is>>n1>>n2>> grad).good() )
	{
		np1 = getNodeByTag(n1);
		if ( np1 == NULL ) 
		{
			std::cout<<n1<<" node not found in the list! Fill the nodelist at first!"<<std::endl;
			return ;
		}
		np2 = getNodeByTag(n2);
		if ( np2 == NULL ) 
		{
			std::cout<<n2<<" node not found in the list! Fill the nodelist at first!"<<std::endl;
			return ;
		}
		nf = new ZBdFace2D(np1, np2, grad);
		boundaryFaces.push_back( nf );
	}
	is.close();
}

void feMesh::readPatches ( std::string mfile)
{
	std::ifstream is( mfile.c_str() );
	if ( !is ) 
	{
		std::cout<< "File of boundary patches can not be opened! "<<std::endl;
		return ;
	}
	std::string str;
	char c;
	int rows, cols;
	amat::Matrix<int> nodeTmp;
	while ( true )
	{
		is>>str; // read in name of patch
		is>>rows>>cols; // read in node of faces list
		while ( is.get(c) )
			if (c== '{') break; 
		nodeTmp.setup(rows,cols);
		for (int i=0; i<rows; i++)
			for (int j=0; j<cols; j++)
				is>>nodeTmp(i,j);
		while ( is.get(c) )
			if (c== '}') break; 
		
		if ( is.eof() ) break;
		ZPatch(this, nodeTmp); 
		std::cout<<str<<" inserted "<<rows << " faces!\n";
		patches.insert( {str, ZPatch(this, nodeTmp)} );
	}
	
	is.close();
}
} //end namespace
