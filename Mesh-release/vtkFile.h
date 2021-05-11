#include "feMesh.h"

namespace ZMesh
{

class vtkFile : public std::ofstream
{
	public : 
		vtkFile( std::string fname) : std::ofstream(fname.c_str())
		{
			isOpened = false;
		    isPrintPointDataHeader = false;
		    isPrintCellDataHeader = false;
			mfname = fname;
			if ( *this )
			{
				isOpened=true;
			}
			else
				std::cout<<"file " <<fname<<" can not be opened!"<<std::endl;
		}

		~vtkFile()
		{
			if (*this)
			{
				this->close();
			}
		}

		void printHeader()
		{
			(*this)
			<<"# vtk DataFile Version 2.0 "<<std::endl
			<<"Volume example"<<std::endl
			<<"ASCII"<<std::endl
			<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
		}

		void printNodes( std::set<ZNode*,ZEntityLessThen> & nodes )
		{
			// note : set nodes's dof start from 0
			if (!isOpened) 
			{
				std::cout<<"open file at first!"<<std::endl;
				return ;
			}
			int nNodes = nodes.size();
			(*this)
			<<"POINTS "<<nNodes<< " double"<<std::endl; 
			for( ZNode* np : nodes )
			{
				(*this)<< np->x() <<" "<< np->y() <<" "<< np->z()<< std::endl;
			} 
		}

		void printElems( std::set<ZElem*, ZEntityLessThen> & elems )
		{
			if (!isOpened) 
			{
				std::cout<<"open file at first!"<<std::endl;
				return ;
			}
			int nElems = elems.size();
			std::ostringstream os;
			std::ostringstream os2;
			int nOut=0;
			os2<<"CELL_TYPES "<<nElems<<std::endl;
			for ( ZElem* ep : elems)
			{
				int etype = ep->getType();
				os<<ep->getNNode()<<" ";
				for(int ii=0; ii<ep->getNNode(); ii++)
					os<<ep->getNode(ii)->dof()<<" ";
				os<<std::endl;
				switch ( etype )
				{
					case (FEM_ZPoint):
						os2<<1<<std::endl;
						nOut += 2;
						break;
					case (1):
					case (2):
					case (3):
					case (5):
					case (9):
					case (10):
						os2<<3<<std::endl;
						nOut += 3;
						break;
					case (4):
					case (6):
						os2<<5<<std::endl;
						nOut += 4;
						break;
					case (7):
						os2<<10<<std::endl;
						nOut += 5;
						break;
					case (8):
						os2<<8<<std::endl;
						nOut += 5;
						break;
					case (11):
						os2<<23<<std::endl;
						nOut += 9;
						break;
					default:
						std::cout<<"Not recognized element type!"<<std::endl;

				}
			}
			(*this)<<"CELLS "<<nElems<<" "<<nOut<<std::endl
				<<os.str()<<std::endl
				<<os2.str()<<std::endl;
		}

		void printPointScalarData( std::string dataName, const amat::Matrix<double> & V )
		{
			if ( !isPrintPointDataHeader )
				(*this)<<"POINT_DATA "<<V.rows()<<std::endl;
			isPrintPointDataHeader = true;
			(*this)<<"SCALARS "<<dataName<<" double "<<V.cols()<<std::endl
				<<"LOOKUP_TABLE default"<<std::endl;
			for(int i=0; i<V.rows(); i++)
			{
				for(int j=0; j<V.cols(); j++)
					(*this)<<V(i,j)<<" ";
				(*this)<<std::endl;
			}
		}

		void printPointVectorData( std::string dataName, const amat::Matrix<double> & V )
		{
			assert( V.cols() == 3 );
			if ( !isPrintPointDataHeader )
				(*this)<<"POINT_DATA "<<V.rows()<<std::endl;
			isPrintPointDataHeader = true;
			(*this)<<"VECTORS "<<dataName<<" double "<<std::endl;
			for(int i=0; i<V.rows(); i++)
			{
				for(int j=0; j<V.cols(); j++)
					(*this)<<V(i,j)<<" ";
				(*this)<<std::endl;
			}
		}

		void printCellScalarData( std::string dataName, const amat::Matrix<double> & V )
		{
			if ( !isPrintCellDataHeader )
				(*this)<<"CELL_DATA "<<V.rows()<<std::endl;
			isPrintCellDataHeader = true;
			(*this)<<"SCALARS "<<dataName<<" double "<<V.cols()<<std::endl
				<<"LOOKUP_TABLE default"<<std::endl;
			for(int i=0; i<V.rows(); i++)
			{
				for(int j=0; j<V.cols(); j++)
					(*this)<<V(i,j)<<" ";
				(*this)<<std::endl;
			}
		}

		void printCellVectorData( std::string dataName, const amat::Matrix<double> & V )
		{
			assert( V.cols() == 3 );
			if ( !isPrintCellDataHeader )
				(*this)<<"CELL_DATA "<<V.rows()<<std::endl;
			isPrintCellDataHeader = true;
			(*this)<<"VECTORS "<<dataName<<" double "<<std::endl;
			for(int i=0; i<V.rows(); i++)
			{
				for(int j=0; j<V.cols(); j++)
					(*this)<<V(i,j)<<" ";
				(*this)<<std::endl;
			}
		}

	private:
		std::string mfname;
		bool isOpened;
		bool isPrintPointDataHeader;
		bool isPrintCellDataHeader;
};

} // end of namespace ZMesh

