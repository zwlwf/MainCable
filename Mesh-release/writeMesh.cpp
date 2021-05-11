#include "feMesh.h"
#include <stdlib.h>

namespace ZMesh{

/* write out to a directory like foam
 * it main write X and T and boundary parts to pathName,
 * the node in them are all from 1, for matlab to use
 * and compressed to continuous number, with max number = nNodes
 */
//note:: no / at the end of pathName
void feMesh::write2Foam(std::string pathName)
{
    std::cout << "feMesh: write2Foam(): writing mesh to directory( no / at the end of pathName) " << pathName << std::endl; 
    std::string myCommand("test -d " + pathName + " || mkdir ");
    system( (myCommand+=pathName).c_str() );  // if pathName is not exist create it. // you have to run it under cygwin?
    
	std::unordered_map<int,int> nodeNumCompressed;
		
	int i=1; // set start index to 1
	for(std::set<ZNode*, ZEntityLessThen>::const_iterator itn = _nlist.begin();
			itn != _nlist.end();
			++itn )
	{
		ZNode * np = *itn;
		nodeNumCompressed[ np->getNum() ] = i++ ;
	}

    std::ofstream outf((pathName+"/nodes").c_str());
    outf << std::setprecision(MESHDATA_DOUBLE_PRECISION);
    //std::cout << "nodes\n";
	for(std::set<ZNode*, ZEntityLessThen>::const_iterator itn = _nlist.begin();
			itn != _nlist.end();
			++itn )
	{
		ZNode * np = *itn;
		outf  << np->x()<< " "
			  << np->y()<< " "
			  << np->z()<< "\n";
	}

    outf.close();
    outf.open((pathName+"/elems").c_str(),std::ofstream::out ) ;
    //std::cout << "elements\n";
	for(std::set<ZElem*, ZEntityLessThen>::const_iterator ite = _elist.begin();
			ite != _elist.end();
			++ite )
	{
		ZElem * ep = *ite;
        for(int i = 0; i < ep->getNNode(); i++)
            outf  << nodeNumCompressed[ (*ep)[i]->getNum() ]<< " ";
        outf << '\n';
    }
    outf.close();

	for( str2ZPatch2Map::iterator it = patches2.begin();
			it != patches2.end();
			it++ )
	{
		outf.open((pathName+"/"+it->first).c_str(),std::ofstream::out ) ;
		ZPatch2 & zp2 = it->second;
		for( ZElem* ep : zp2 )
		{
			for(int i=0; i<ep->getNNode(); i++)
				outf<<ep->getNode(i)->getNum()<<" ";
			outf<<std::endl;
		}
		outf.close();
	}

}

    void feMesh::write2vtk(const std::string fname)
	{
		std::ofstream os( fname.c_str() );
		os
		<<"# vtk DataFile Version 2.0 "<<std::endl
		<<"Volume example"<<std::endl
		<<"ASCII"<<std::endl
		<<"DATASET UNSTRUCTURED_GRID"<<std::endl
		<<"POINTS "<<getNNode()<< " double"<<std::endl;
		for( std::set< ZNode*, ZEntityLessThen>::iterator itn = _nlist.begin();
			 itn != _nlist.end();
			 itn++ )
		{
			ZNode * np = *itn;
			os<< np->x() <<" "<< np->y() <<" "<< np->z()<< std::endl;
		}

		os<<"CELLS "<<getNElem()<<" ";
		int nOut=0;
		long pos = os.tellp();
		os<<"                                              "<<std::endl;
		std::ostringstream os2;
		os2<<"CELL_TYPES "<<getNElem()<<std::endl;
		for (std::set< ZElem*, ZEntityLessThen>::iterator ite = _elist.begin();
				ite != _elist.end();
				ite++)
		{
			ZElem * ep = *ite;
			int etype = ep->getType();
			std::vector<int> dofs = ep->getAllNodeDofs();
			switch ( etype )
			{
				case ( FEM_ZPoint ):
					os<<1<<" "<<dofs[0]<<std::endl;
					os2<<1<<std::endl;
					nOut += 2;
					break;
				case (FEM_ZLink2D):
				case (FEM_ZLink3D):
				case (FEM_ZBeam2D):
				case (FEM_ZBeam3D):
				case (FEM_ZLine2):
				case (FEM_ZLine3Q):
					os<<2<<" "<<dofs[0]<<" "<<dofs[1]<<std::endl;
					os2<<3<<std::endl;
					nOut += 3;
					break;
				case ( FEM_ZTri2D ):
				case ( FEM_ZTri2D1V):
				case ( FEM_ZTri3D ) :
					os<<3<<" "<<dofs[0]<<" "<<dofs[1]<<" "<<dofs[2]<<std::endl;
					os2<<5<<std::endl;
					nOut += 4;
					break;
				case (FEM_ZTet):
					os<<4<<" "<<dofs[0]<<" "<<dofs[1]<<" "<<dofs[2]<<" "<<dofs[3]<<std::endl;
					os2<<10<<std::endl;
					nOut += 5;
					break;
				case (FEM_ZQuad4Q):
					os<<4<<" "<<dofs[0]<<" "<<dofs[1]<<" "<<dofs[3]<<" "<<dofs[2]<<std::endl;
					os2<<8<<std::endl;
					nOut += 5;
					break;
				case (FEM_ZQuad8Q):
					os<<8<<" "<<dofs[0]<<" "<<dofs[1]<<" "<<dofs[2]<<" "<<dofs[3]<<
				" "<<dofs[4]<<" "<<dofs[5]<<" "<<dofs[6]<<" "<<dofs[7]<<std::endl;
					os2<<23<<std::endl;
					nOut += 9;
					break;
				default:
					std::cout<<"Not recognized element type! at [writMesh.cpp::write2vtk()"<<std::endl;

			}
		}
		os<<std::endl
		<< os2.str();
		os.seekp(pos);
		os<<nOut;
		os.close();
	}

    void feMesh::write2Tikz(const std::string fname)
	{
		std::cout<<" output to tikz start !" <<std::endl;
		std::ofstream os( fname.c_str() );
		os<<"\\begin{tikzpicture}"<<std::endl
		  <<"\\newdimen\\radius "<<std::endl
		  <<"\\radius=0.04cm "<<std::endl;
		os<<"% output nodes "<<std::endl;

		for( std::set< ZNode*, ZEntityLessThen>::iterator itn = _nlist.begin();
			 itn != _nlist.end();
			 itn++ )
		{
			ZNode * np = *itn;
			os<< "\\draw [fill] ("<<np->x() <<","<< np->y() <<") circle (\\radius); " << std::endl;
		}

		os<<"% output elements "<<std::endl;
		for (std::set< ZElem*, ZEntityLessThen>::iterator ite = _elist.begin();
				ite != _elist.end();
				ite++)
		{
			ZElem * ep = *ite;
			int etype = ep->getType();
			std::vector<int> dofs = ep->getAllNodeDofs();
			switch ( etype )
			{
				case (1):
				case (2):
				case (3):
				case (5):
				case (9):
				case (10):
					os
					<<"\\draw ("<< ep->getNode(0)->x()<<","<<ep->getNode(0)->y()<<") -- ("<< ep->getNode(1)->x()<<","<<ep->getNode(1)->y() <<"); "<<std::endl;
					break;
				case (4):
				case (6):
					os
					<<"\\draw ("<< ep->getNode(0)->x()<<","<<ep->getNode(0)->y()<<") -- ("
					<< ep->getNode(1)->x()<<","<<ep->getNode(1)->y() <<")-- ( "
					<< ep->getNode(2)->x()<<","<<ep->getNode(2)->y() <<") -- cycle; "<<std::endl;
					break;
				case (7):
					os
					<<"\\draw ("<< ep->getNode(0)->x()<<","<<ep->getNode(0)->y()<<") -- ("
					<< ep->getNode(1)->x()<<","<<ep->getNode(1)->y() <<")-- ( "
					<< ep->getNode(2)->x()<<","<<ep->getNode(2)->y() <<")-- ( "
					<< ep->getNode(3)->x()<<","<<ep->getNode(3)->y() <<") -- cycle; "<<std::endl;
					break;
				case (8):
					os
					<<"\\draw ("<< ep->getNode(0)->x()<<","<<ep->getNode(0)->y()<<") -- ("
					<< ep->getNode(1)->x()<<","<<ep->getNode(1)->y() <<")-- ( "
					<< ep->getNode(3)->x()<<","<<ep->getNode(3)->y() <<")-- ( "
					<< ep->getNode(2)->x()<<","<<ep->getNode(2)->y() <<") -- cycle; "<<std::endl;
					break;
				case (11):
					os
					<<"\\draw ("<< ep->getNode(0)->x()<<","<<ep->getNode(0)->y()<<") -- ("
					<< ep->getNode(1)->x()<<","<<ep->getNode(1)->y() <<")-- ( "
					<< ep->getNode(2)->x()<<","<<ep->getNode(2)->y() <<")-- ( "
					<< ep->getNode(3)->x()<<","<<ep->getNode(3)->y() <<")-- ( "
					<< ep->getNode(4)->x()<<","<<ep->getNode(4)->y() <<")-- ( "
					<< ep->getNode(5)->x()<<","<<ep->getNode(5)->y() <<")-- ( "
					<< ep->getNode(6)->x()<<","<<ep->getNode(6)->y() <<")-- ( "
					<< ep->getNode(7)->x()<<","<<ep->getNode(7)->y() <<") -- cycle; "<<std::endl;
					break;
				default:
					std::cout<<"Not recognized element type!"<<std::endl;

			}
		}
		
		// loop over boudary patch element type
		for ( str2ZPatch2Map::iterator itp = patches2.begin();
				itp != patches2.end();
				itp ++ )
		{
			os<<"% output boundary "<<itp->first<<std::endl;
			ZPatch2 & zp = itp->second;
			for ( std::list<ZElem* >::iterator it = zp.begin();
					it != zp.end();
					it ++ )
			{
				ZElem *ep = *it;
				int etype = ep->getType();
				switch ( etype )
				{
					case (1):
					case (2):
					case (3):
					case (5):
					case (9):
					case (10):
						os
						<<"\\draw[thick] ("<< ep->getNode(0)->x()<<","<<ep->getNode(0)->y()<<") -- ("<< ep->getNode(1)->x()<<","<<ep->getNode(1)->y() <<"); "<<std::endl;
						break;
					case (4):
					case (6):
						os
						<<"\\draw [thick]("<< ep->getNode(0)->x()<<","<<ep->getNode(0)->y()<<") -- ("
						<< ep->getNode(1)->x()<<","<<ep->getNode(1)->y() <<")-- ( "
						<< ep->getNode(2)->x()<<","<<ep->getNode(2)->y() <<") -- cycle; "<<std::endl;
						break;
					case (8):
						os
						<<"\\draw[thick] ("<< ep->getNode(0)->x()<<","<<ep->getNode(0)->y()<<") -- ("
						<< ep->getNode(1)->x()<<","<<ep->getNode(1)->y() <<")-- ( "
						<< ep->getNode(3)->x()<<","<<ep->getNode(3)->y() <<")-- ( "
						<< ep->getNode(2)->x()<<","<<ep->getNode(2)->y() <<") -- cycle; "<<std::endl;
						break;
					default:
						std::cout<<"Not recognized element type!"<<std::endl;

				}
			}
		}
		os<<"\\end{tikzpicture}"<<std::endl;
		os.close();
		std::cout<<" output to tikz done !" <<std::endl;
	}
}
