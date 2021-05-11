#include "feMesh.h"
// just read in the gmsh and store elements into ZPatch2 by its physical name, so you can deal with the each physical name

namespace ZMesh{

void feMesh::readGmsh3(std::string mfile)
{
	std::cout << "readGmsh3(): Reading mesh file...\n";
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
	int* phyIds;
	std::map<int,std::string> phyId2Name;
	std::map<std::string,int> phyName2Id;
	while ( is>>str1)
	{
		if (str1 == "$PhysicalNames")
		{
			is>>nPhy;
			phyDims = new int[nPhy];
			phyIds = new int[nPhy];
			for (int i=0; i<nPhy; i++)
			{
				is>> idim >> id >> str2; 
				std::cout<<" dimension = "<<idim <<", id =" <<id<<", name = "<<str2.substr(1,str2.size()-2)<<std::endl;
				phyDims[id] = idim; 
				phyIds[i] = id; 
				phyId2Name[id] = str2.substr(1,str2.size()-2); 
				phyName2Id[str2.substr(1,str2.size()-2)] = id; 
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
				//if ( phyDims[PhyId] < nDim ) // boundary patch
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
						case(12): // linear triangles
							is>>d1>>d2>>d3;
							zpatches[PhyId].add( new ZTri3D( id, d1, d2,d3) );
							break;
						case(8):	// linear quads
							is>>d1>>d2>>d3>>d4;
							zpatches[PhyId].add( new ZQuad4Q( id, d1, d2,d3,d4 ) );
							break;
						case(7):
							is>>d1>>d2>>d3>>d4;
							zpatches[PhyId].add( new ZTet(id, d1,d2,d3,d4) );
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

	for(int i=0; i<nPhy; i++)
	{
		//if( phyDims[i] < nDim )
		{
			zpatches[phyIds[i]].setMesh(this);
			patches2.insert( {phyId2Name[ phyIds[i] ], zpatches[phyIds[i]]} );
			std::cout<<zpatches[phyIds[i]].getNFaces()<<" faces in patch "<< phyId2Name[ phyIds[i] ]<<" inserted!"<<std::endl;
		}
	}

	setElemList(elist);
	createNodeDof();
}

}
