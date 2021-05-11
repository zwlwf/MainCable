/** \file asparsematrix.hpp
 * \brief Sparse storage schemes for matrices (header-only library)
 * \author Aditya Kashi
 * \date 21 July 2015
 */

#ifndef __ASPARSEMATRIX_H
#define __ASPARSEMATRIX_H

// for memcpy() etc
#ifndef _GLIBCXX_CSTRING
#include <cstring>
#endif

#ifdef AMGCL_LIBRARY
#include <amgcl/make_solver.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/solver/gmres.hpp>
#include <amgcl/solver/cg.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#endif

// for malloc()
#ifndef _GLIBCXX_CSTDLIB
#include <cstdlib>
#endif

#ifndef _GLIBCXX_IOSTREAM
#include <iostream>
#endif

#ifndef _GLIBCXX_FSTREAM
#include <fstream>
#endif

#ifndef _GLIBCXX_CMATH
#include <cmath>
#endif

#ifndef _GLIBCXX_VECTOR
#include <vector>
#endif
#include <set>

#ifndef __AMATRIX_H
#include "amatrix.hpp"
#endif

#ifdef _OPENMP
#ifndef OMP_H
#include <omp.h>
#endif
#endif

namespace amat {

typedef Matrix<amc_real> Mat;

/// The type of dense storage to use for vectors in matrix-vector multiplication, etc.
template <typename T>
using DenseMatrix = Matrix<T>;

/// Structure for matrix in compressed row storage
/** \note Set [allocated](@ref allocated) to "true" right after reserving memory for the arrays.
 * Use matrix type "SLU_NR" while using this with SuperLU.
 */
template <class T> struct SMatrixCRS
{
	int nnz;			///< Number of non-zero entries in the matrix
	T* val;				///< Non-zero values of the matrix
	int* col_ind;		///< Column indices of non-zero values of the matrix
	int* row_ptr;		///< Element indices of val and col_ind in which each row begins
	bool allocated;		///< Stores whether the 3 arrays have been allocated

	SMatrixCRS()
	{
		allocated = false;
	}
	
	~SMatrixCRS()
	{
		if(allocated)
		{
			delete [] val;
			delete [] col_ind;
			delete [] row_ptr;
			allocated = false;
		}
	}

	int * getColP() { return col_ind; }
	int * getRowP() { return row_ptr; }
	T * getValP() { return val; }
	int getNnz() { return nnz; }
};


/// Abstract class for sparse storage of matrices
template <class T> class SparseMatrix
{
protected:
	int nnz;
	int nrows;
	int ncols;
public:
	SparseMatrix() {}
	SparseMatrix(int num_rows, int num_cols) : nrows(num_rows), ncols(num_cols)
	{
		nnz = 0;
	}
	int rows() const { return nrows; }
	int cols() const { return ncols; }

	virtual void set(const int x, const int y, const double value) = 0;
	virtual T get(const int x, const int y) const = 0;
	virtual void mprint() const = 0;

	/// Multiplication of this sparse matrix with a row-major vector
	virtual void multiply(const Mat& x, Mat* const a, const char paralel) const = 0;
	virtual void multiply_parts(const Mat* x, const Mat* y, Mat* const ans, const int p) const = 0;
	virtual double getelem_multiply_parts(const int rownum, const Mat* x, const Mat* y, const int p, const double num) const = 0;
};

#define INITIAL_ROW_SIZE 20

/** \brief Implements sparse matrix storage in row-storage format, but with separate arrays for each row of the matrix.
 */
class AMGSolver;

template<class T> class MatrixCRS : public SparseMatrix<T>
{
	protected:
	std::vector<T>* val;
	std::vector<int>* col_ind;			///< col_ind[k][j] contains the column number of val[k][j]
	std::vector<int> rsize;				///< Number of non-zero elements in each row

	using SparseMatrix<T>::nnz;
	using SparseMatrix<T>::nrows;
	using SparseMatrix<T>::ncols;

#ifdef AMGCL_LIBRARY
	typedef amgcl::backend::builtin<double> Backend;
    typedef amgcl::make_solver<
           // Use AMG as preconditioner:
           amgcl::amg<
               Backend,
               amgcl::coarsening::smoothed_aggregation,
               amgcl::relaxation::spai0
               >,
           // And BiCGStab as iterative solver:
           amgcl::solver::bicgstab<Backend>
           //amgcl::solver::cg<Backend>
           //amgcl::solver::gmres<Backend>
           > AMGSolver;
#endif

	AMGSolver *solver;
	bool isAMG; // whether the solver of AMG created

public:
	MatrixCRS() 
	{ 
		val = new std::vector<T>[2];  
		col_ind = new std::vector<int>[2];
		rsize.push_back(0); // by zwl
		rsize.push_back(0);
		isAMG = false;
	}

	MatrixCRS(int num_rows, int num_cols) : SparseMatrix<T>(num_rows,num_cols)
	{
		rsize.resize(nrows,0);
		val = new std::vector<T>[nrows];
		col_ind = new std::vector<int>[nrows];
		for(int i = 0; i < nrows; i++)
		{
			val[i].reserve(INITIAL_ROW_SIZE);
			col_ind[i].reserve(INITIAL_ROW_SIZE);
			rsize[i] = 0; //! add by zwl
		}
		isAMG = false;
		nnz = 0;
	}

	MatrixCRS(const MatrixCRS& other)
	{
		rsize = other.rsize;
		nnz = other.nnz;
		nrows = other.nrows; ncols = other.ncols;

		val = new std::vector<T>[nrows];
		col_ind = new std::vector<int>[nrows];

		for(int i = 0; i < nrows; i++)
			for(int j = 0; j < rsize[i]; j++)
			{
				val[i].push_back(other.val[i][j]);
				col_ind[i].push_back(other.col_ind[i][j]);
			}
		isAMG = false;
	}

	MatrixCRS& operator=(const MatrixCRS& other)
	{
		rsize = other.rsize;
		nnz = other.nnz;
		nrows = other.nrows; ncols = other.ncols;

		delete [] val;
		delete [] col_ind;
		val = new std::vector<T>[nrows];
		col_ind = new std::vector<int>[nrows];

		for(int i = 0; i < nrows; i++)
			for(int j = 0; j < rsize[i]; j++)
			{
				val[i].push_back(other.val[i][j]);
				col_ind[i].push_back(other.col_ind[i][j]);
			}
		return *this;
	}

	void setup(int num_rows, int num_cols)
	{
		delete [] val;
		delete [] col_ind;

		nrows = num_rows; ncols = num_cols;
		rsize.resize(nrows,0); //! by zwl. when nrows <= old one, the elements is not changed to 0
		val = new std::vector<T>[nrows];
		col_ind = new std::vector<int>[nrows];
		for(int i = 0; i < nrows; i++)
		{
			val[i].reserve(INITIAL_ROW_SIZE);
			col_ind[i].reserve(INITIAL_ROW_SIZE);
			rsize[i] = 0; //! add by zwl
		}
		nnz = 0;
	}

	~MatrixCRS()
	{
		delete [] val;
		delete [] col_ind;
	}

	void set(const int x, const int y, const T value)
	{
		if(dabs(double(value)) > ZERO_TOL)
		{
			// search row x for the value
			int pos = -1;
			for(int j = 0; j < rsize[x]; j++)
				if(col_ind[x][j] == y)
					pos = j;
			
			if(pos == -1)
			{
				//#pragma omp critical (omp_setval)
				{
					val[x].push_back(value);
					col_ind[x].push_back(y);
					rsize[x]++;
				}
			}
			else
			{
				val[x][pos] = value;
			}
			nnz++;
		}
	}
	
	int getNnz() { return nnz; }
	// whether this sparse matrix is initialized, by zwl
	bool isAllocated()
	{
		return nrows == rsize.size();
	}

	T get(const int x, const int y) const
	{
		T retval = T(0);
		for(int j = 0; j < rsize[x]; j++)
			if(col_ind[x][j] == y)
				retval = val[x][j];
		return retval;
	}

	/// Combined getter/setter method.
	/// Usage of this method might lead to zeros being stored.
	T& operator()(const int x, const int y)
	{
		T* retval;
		// search row x for the value
		int pos = -1;
		for(int j = 0; j < rsize[x]; j++)
			if(col_ind[x][j] == y)
				pos = j;
		
		if(pos == -1)
		{
			val[x].push_back(T(0));
			col_ind[x].push_back(y);
			rsize[x]++;
			retval = &val[x][rsize[x]-1];
			nnz++;
		}
		else
		{
			retval = &val[x][pos];
		}
		return *retval;
	}

	/** Shrinks allocated memory of val and col_ind to fit the current number of non-zero elements.*/
	void trim()
	{
		for(int i = 0; i < nrows; i++) {
			val[i].shrink_to_fit();
			col_ind[i].shrink_to_fit();
		}
	}


	/*void print_data()
	{
		for(int i = 0; i < nnz; i++)
			std::cout << val[i] << " ";
		std::cout << std::endl;
		for(int i = 0; i < nnz; i++)
			std::cout << row_ind[i] << " ";
		std::cout << std::endl;
		for(int i = 0; i < nnz; i++)
			std::cout << col_ind[i] << " ";
		std::cout << std::endl;
	}*/

	void mprint() const
	{
		std::cout << "\n";
		for(int i = 0; i < nrows; i++)
		{
			for(int j = 0; j < ncols; j++)
				if ( dabs(get(i,j))> ZERO_TOL)
					std::cout <<" ("<<i<<", "<<j<<") "<< get(i,j);
			std::cout << std::endl;
		}
		std::cout<<"ZERO_TOL = "<<ZERO_TOL<<std::endl;
	}

	void fprint(std::ofstream& outfile) const
	{
		//outfile << "\n";
		for(int i = 0; i < nrows; i++)
		{
			for(int j = 0; j < ncols; j++)
				if ( dabs(get(i,j))> ZERO_TOL)
					outfile << " " << get(i,j);
			outfile << std::endl;
		}
	}

	/**	Sorts each row by column index using insertion sort.
	 * We expect each row to have no more than tens of elements, so insertion sort is a good candidate.
	 */
	void sort_rows()
	{
		// insertion sort
		int tempc = 0; T tempv = 0; int j = 0;

		std::vector<int>* rsize = &(this->rsize);
		#ifdef _OPENMP
		std::vector<T>* val = this->val;
		std::vector<int>* col_ind = this->col_ind;
		int nrows = this->nrows;
		int ncols = this->ncols;
		#endif

		for(int l = 0; l < nrows; l++)			// for each row
		{
			for(int i = 1; i < (*rsize)[l]; i++)
			{
				j = i;
				while(j>0 && col_ind[l][j] < col_ind[l][j-1])
				{
					// swap col_ind[l][j] and col_ind[l][j-1], same for vp
					tempc = col_ind[l][j];
					col_ind[l][j] = col_ind[l][j-1];
					col_ind[l][j-1] = tempc;
					tempv = val[l][j];
					val[l][j] = val[l][j-1];
					val[l][j-1] = tempv;
					j--;
				}
			}
		}
	}

	/// Returns product of sparse matrix with vector x and stores it in y_out, for the matrix-vector operation in Spectra,  add by zwl 20190404
	void multiply(const double* x, double *y_out, const char paralel = 'n') const
	{
		for(int i=0; i<nrows; i++)
			y_out[i] = 0;
		
		double temp;
		#ifdef _OPENMP
		const std::vector<T>* val = this->val; const std::vector<int>* col_ind = this->col_ind;
		#endif

 		const std::vector<int>* rsize = &(this->rsize);

		{
			int i;
			//#pragma omp parallel for default(none) private(i,temp) shared(val,col_ind,rsize,y_out,x) if(paralel=='y')
			for(i = 0; i < nrows; i++)
			{
				for(int j = 0; j < (*rsize)[i]; j++)
				{
					temp = val[i][j]*x[col_ind[i][j]];
					y_out[i] += temp;
				}
			}
		}
	}

	/// Returns product of sparse matrix with vector x and stores it in a.
	void multiply(const Mat& x, Mat* const a, const char paralel = 'n') const
	{
		a->zeros();
		
		//std::cout << "MatrixCOO: multiply(): Rows and columns of a: " << a->rows() << " " << a->cols() << std::endl;
		double temp;
		#ifdef _OPENMP
		const std::vector<T>* val = this->val; const std::vector<int>* col_ind = this->col_ind;
		#endif

 		const std::vector<int>* rsize = &(this->rsize);

		for(int k = 0; k < x.cols(); k++)		// for each column of x
		{
			int i;
			//#pragma omp parallel for default(none) private(i,temp) shared(val,col_ind,rsize,k,a,x) if(paralel=='y')
			for(i = 0; i < nrows; i++)
			{
				for(int j = 0; j < (*rsize)[i]; j++)
				{
					temp = val[i][j]*x.get(col_ind[i][j],k);
					(*a)(i,k) += temp;
				}
			}
		}
	}

	Mat operator*(const Mat& x) const
	{
		Mat a( x.rows(), x.cols());
		a.zeros();
		
		//std::cout << "MatrixCOO: multiply(): Rows and columns of a: " << a->rows() << " " << a->cols() << std::endl;
		double temp;
		#ifdef _OPENMP
		const std::vector<T>* val = this->val; const std::vector<int>* col_ind = this->col_ind;
		#endif

 		const std::vector<int>* rsize = &(this->rsize);

		for(int k = 0; k < x.cols(); k++)		// for each column of x
		{
			int i;
			//#pragma omp parallel for default(none) private(i,temp) shared(val,col_ind,rsize,k,a,x) if(paralel=='y')
			for(i = 0; i < nrows; i++)
			{
				for(int j = 0; j < (*rsize)[i]; j++)
				{
					temp = val[i][j]*x.get(col_ind[i][j],k);
					a(i,k) += temp;
				}
			}
		}
		return a;
	}

	/// Returns product of two sparse matrix, by zwl
	MatrixCRS multiply2(const MatrixCRS<T> & R ,const char paralel = 'n') const
	{
		MatrixCRS<T> a(nrows, R.ncols);
		
		//std::cout << "MatrixCOO: multiply(): Rows and columns of a: " << a->rows() << " " << a->cols() << std::endl;
		double temp;
		#ifdef _OPENMP
		const std::vector<T>* val = this->val; const std::vector<int>* col_ind = this->col_ind;
		#endif

		for(int i = 0; i < nrows; i++)
		{
			for(int k = 0; k < R.cols(); k++)	
			{
				for(int j = 0; j < rsize[i]; j++)
				{
					temp = R.get(col_ind[i][j],k);
					if(dabs(double(temp)) > ZERO_TOL)
						a(i,k) = a.get(i,k) + val[i][j]*temp;
				}
			}
		}
		return a;
	}

	/// Like the multiply() method, except the argument matrix is considered x for the first p-1 rows, 0 in the pth row and y for the remaining rows.
	void multiply_parts(const Mat* x, const Mat* y, Mat* const ans, const int p) const
	{
		ans->zeros();
		#ifdef _OPENMP
		const std::vector<T>* val = MatrixCRS::val; const std::vector<int>* col_ind = MatrixCRS::col_ind;
		#endif
		const std::vector<int>* rsize = &this->rsize;

		for(int k = 0; k < x->cols(); k++)		// for each column of x
		{
			int i;
			//#pragma omp parallel for default(none) private(i) shared(val,row_ind,col_ind,nnz,k,ans,x,y,p)
			for(i = 0; i < nrows; i++)
			{
				for(int j = 0; j < (*rsize)[i]; j++)
				{
					if(col_ind[i][j] < p)
						(*ans)(i,k) += val[i][j] * x->get(col_ind[i][j],k);
					else if(col_ind[i][j] > p)
						(*ans)(i,k) += val[i][j] * y->get(col_ind[i][j],k);
				}
			}
		}
	}

	/** Returns dot product of rownum'th row of this matrix with a certain std::vector. 
	 * This vector is composed of x for for the first p-1 elements and y for p+1'th element onwards, with num as pth element.
	 * \note NOTE: Make sure dimensions of x and y are same.
	 */
	T getelem_multiply_parts(const int rownum, const Mat* x, const Mat* y, const int p, const T num) const
	{
		T ans = 0;
		//if(sorted == false) std::cout << "! MatrixCOO: multiply(): Matrix not sorted yet!\n";

		for(int k = 0; k < x->cols(); k++)		// for each column of x
		{
			#ifdef _OPENMP
			const std::vector<T>* val = MatrixCRS::val; const std::vector<int>* col_ind = MatrixCRS::col_ind;
			#endif
			const std::vector<int>* rsize = &this->rsize;
			int j;
			//#pragma omp parallel for default(none) private(j) shared(val,col_ind,rsize,k,rownum,p,x,y,num) reduction(+:ans)
			for(j = 0; j < (*rsize)[rownum]; j++)
			{
					if(col_ind[rownum][j] < p)
						ans += val[rownum][j] * x->get(col_ind[rownum][j],k);
					else if(col_ind[rownum][j] > p)
						ans += val[rownum][j] * y->get(col_ind[rownum][j],k);
					else
						ans += val[rownum][j] * num;
			}
		}
		return ans;
	}

	// use rownum'th row to dot x vector, return a T
	amat::Matrix<T> row_multiply(const int rownum, const Mat* x) const
	{
		amat::Matrix<T> ans(1,x->cols()); // edit by zwl
		ans.zeros();
		for(int k = 0; k < x->cols(); k++)		// for each column of x
		{
			const std::vector<int>* rsize = &this->rsize;
			int j;
			//#pragma omp parallel for default(none) private(j) shared(val,col_ind,rsize,k,rownum,p,x,y,num) reduction(+:ans)
			for(j = 0; j < (*rsize)[rownum]; j++)
			{
						ans(0,k) += val[rownum][j] * x->get(col_ind[rownum][j],k);
			}
		}
		return ans;
	}

	/// D is returned as a column-vector containing diagonal elements of this sparse matrix.
	void get_diagonal(Mat* D) const
	{
		D->zeros();
		for(int i = 0; i < nrows; i++)
		{
			for(int j = 0; j < rsize[i]; j++)	
				if(i == col_ind[i][j])
				{
					//(*D)(row_ind[i]) = val[i];
					D->set(i,0, val[i][j]);
				}
		}
	}
	
	void get_lower_triangle(MatrixCRS& L)
	{
		delete [] L.val;
		delete [] L.col_ind;
		L.val = new std::vector<T>[nrows];
		L.col_ind = new std::vector<int>[nrows];
		
		L.nrows = nrows;
		L.ncols = ncols;
		L.rsize.resize(nrows,0);
		
		
		for(int i = 0; i < nrows; i++)
		{
		
			for(int j = 0; j < rsize[i]; j++)
				if(i > col_ind[i][j])
				{
					L.val[i].push_back(val[i][j]);
					L.col_ind[i].push_back(col_ind[i][j]);
					L.rsize[i]++;
				}
		}
		
		// shrink to fit
		/*if(L.nnz < nnz)
		{
			for(int i = nnz-1; i >= L.nnz; i--)
			{
				delete L.val+i;
			}
		}*/
	}
	
	void get_upper_triangle(MatrixCRS& L)
	{
		delete [] L.val;
		delete [] L.col_ind;
		L.val = new std::vector<T>[nrows];
		L.col_ind = new std::vector<int>[nrows];
		
		L.nrows = nrows;
		L.ncols = ncols;
		L.rsize.resize(nrows,0);
		
		
		for(int i = 0; i < nrows; i++)
		{
		
			for(int j = 0; j < rsize[i]; j++)
				if(i < col_ind[i][j])
				{
					L.val[i].push_back(val[i][j]);
					L.col_ind[i].push_back(col_ind[i][j]);
					L.rsize[i]++;
				}
		}
	}

	MatrixCRS<T> transpose()
	{
		std::cout << "MatrixCRS: transpose(): Transposing matrix" << std::endl;
		MatrixCRS<T> mat;
		delete [] mat.val;
		delete [] mat.col_ind;
		mat.val = new std::vector<T>[nrows];
		mat.col_ind = new std::vector<int>[nrows];
		mat.rsize.resize(nrows);
		mat.nrows = nrows; mat.ncols = ncols;

		int i;

		for(i = 0; i < nrows; i++)
		{
			for(int j = 0; j < rsize[i]; j++)
			{
				mat.val[col_ind[i][j]].push_back(val[i][j]);
				mat.col_ind[col_ind[i][j]].push_back(i);
				mat.rsize[col_ind[i][j]]++;
			}
		}
		// not sorted

		return mat;
	}

	void operator+=( const MatrixCRS<T> & B) // by zwl, 2017/05/15
	{
		if ( B.nrows != nrows || B.ncols != ncols)
		{
			std::cout<<"these two matrixes are not the same size!\n";
			return ;
		}
		for (int i=0; i<nrows; i++)
		{
			for (int j=0; j<B.rsize[i]; j++)
			this->operator()(i, B.col_ind[i][j]) += B.val[i][j];
		}
	}

	void operator *=( T lambda) // by zwl, 2017/06/05
	{
		for (int i=0; i<nrows; i++)
		{
			for (int j=0; j<rsize[i]; j++)
			this->operator()(i, col_ind[i][j]) *= lambda;
		}
	}

	void operator-=( MatrixCRS<T> & B) // by zwl, 2017/05/15
	{
		if ( B.nrows != nrows || B.ncols != ncols)
		{
			std::cout<<"these two matrixes are not the same size!\n";
			return ;
		}
		for (int i=0; i<nrows; i++)
		{
			for (int j=0; j<B.rsize[i]; j++)
			this->operator()(i, B.col_ind[i][j]) -= B.val[i][j];
		}
	}

	void AddEqual( T alpha, MatrixCRS<T> & B) // by zwl, 2018/11/26
	{
		if ( B.nrows != nrows || B.ncols != ncols)
		{
			std::cout<<"these two matrixes are not the same size!\n";
			return ;
		}
		for (int i=0; i<nrows; i++)
		{
			for (int j=0; j<B.rsize[i]; j++)
			this->operator()(i, B.col_ind[i][j]) += alpha*B.val[i][j];
		}
	}

	void MinusEqual( T alpha, MatrixCRS<T> & B) // by zwl, 2018/11/26
	{
		if ( B.nrows != nrows || B.ncols != ncols)
		{
			std::cout<<"these two matrixes are not the same size!\n";
			return ;
		}
		for (int i=0; i<nrows; i++)
		{
			for (int j=0; j<B.rsize[i]; j++)
			this->operator()(i, B.col_ind[i][j]) -= alpha*B.val[i][j];
		}
	}

	// setRow , by zwl
	void setRow(int i, std::vector<T> val_i, std::vector<int> col_i)
	{
		nnz = nnz - rsize[i] + col_i.size();
		rsize[i] = col_i.size();
		val[i] = val_i;
		col_ind[i] = col_i;
		val[i].shrink_to_fit();
		col_ind[i].shrink_to_fit();
	}

	/// Combines A, B, C, and D (4 n-by-n matrices) into one 2n-by-2n matrix. 
	/** CAUTION: All 4 input matrices must have same size! */
	void combine_sparse_matrices(const MatrixCRS& A11, const MatrixCRS& A12, const MatrixCRS& A21, const MatrixCRS& A22)
	{
		std::cout << "MatrixCOO: combine_sparse_matrices(): Combining matrices" << std::endl;
		//std::cout << "MatrixCOO: combine_sparse_matrices(): A11 dimensions: " << A11.nrows << " " << A11.ncols << std::endl;

		nrows = A11.nrows*2; ncols = A11.ncols*2;
		delete [] val; delete [] col_ind;
		val = new std::vector<T>[nrows];
		col_ind = new std::vector<int>[nrows];

		// reserve some amount of space
		int maxrowsize = 6;
		for(int ir = 0; ir < nrows; ir++)
			if(A11.rsize[ir] > maxrowsize) maxrowsize = A11.rsize[ir];

		for(int ir = 0; ir < nrows; ir++) {
			val[ir].reserve(INITIAL_ROW_SIZE*2);
			col_ind[ir].reserve(INITIAL_ROW_SIZE*2);
		}

		rsize.resize(nrows,0);

		// start filling up the matrix

		int i;
		for(i = 0; i < A11.nrows; i++)
		{
			for(int j = 0; j < A11.rsize[i]; j++)
			{
				val[i].push_back(A11.val[i][j]);
				col_ind[i].push_back(A11.col_ind[i][j]);
				rsize[i]++;
			}
		}

		for(i = 0; i < A12.nrows; i++)
		{
			for(int j = 0; j < A12.rsize[i]; j++)
			{
				val[i].push_back(A12.val[i][j]);
				col_ind[i].push_back(A12.col_ind[i][j]+nrows/2);
				rsize[i]++;
			}
		}

		for(i = 0; i < A21.nrows; i++)
		{
			for(int j = 0; j < A21.rsize[i]; j++)
			{
				val[i+nrows/2].push_back(A21.val[i][j]);
				col_ind[i+nrows/2].push_back(A21.col_ind[i][j]);
				rsize[i+nrows/2]++;
			}
		}

		for(i = 0; i < A22.nrows; i++)
		{
			for(int j = 0; j < A22.rsize[i]; j++)
			{
				val[i+nrows/2].push_back(A22.val[i][j]);
				col_ind[i+nrows/2].push_back(A22.col_ind[i][j]+nrows/2);
				rsize[i+nrows/2]++;
			}
		}
	}

	///	Converts to dense matrix.
	Mat toDense() const
	{
		Matrix<T> dense(nrows, ncols);
		dense.zeros();
		for(int i = 0; i < nrows; i++)
			for(int j = 0; j < rsize[i]; j++)
				dense(i,col_ind[i][j]) = val[i][j];
		return dense;
	}

	/// Returns the sparse matrix in compressed row format
	void get_CRS_matrix(SMatrixCRS<T>& rmat) const
	{
		rmat.nnz = 0;
		for(int i = 0; i < nrows; i++)
			rmat.nnz += rsize[i];

		if(rmat.allocated)
		{
			delete [] rmat.val;
			delete [] rmat.col_ind;
			delete [] rmat.row_ptr;
		}

		rmat.val = new T[rmat.nnz];
		rmat.col_ind = new int[rmat.nnz];
		rmat.row_ptr = new int[nrows+1];
		rmat.allocated = true;
		rmat.row_ptr[0] = 0;

		int i, j, k = 0;
		for(i = 0; i < nrows; i++)
		{
			for(j = 0; j < rsize[i]; j++)
			{
				rmat.val[k] = val[i][j];
				rmat.col_ind[k] = col_ind[i][j];
				k++;
			}
			rmat.row_ptr[i+1] = k;
		}
	}

	// get the matrix that do not include dofs in excludeDof
	// add by zwl on 2019/04/07
	void get_CRS_matrix_excludeDof(SMatrixCRS<T>& rmat, const std::set<int> & excludeDof) const
	{
		// call sort_rows() at first if not done ever.
		rmat.nnz = nnz;
		if ( rmat.allocated )
		{
			delete [] rmat.val;
			delete [] rmat.col_ind;
			delete [] rmat.row_ptr;
		}
		rmat.val = new T[rmat.nnz];
		rmat.col_ind = new int[rmat.nnz];
		rmat.row_ptr = new int[nrows+1];
		rmat.row_ptr[0] = 0;
		int k=0;
		int rowTrue=0;
		int *shiftDof = new int[nrows]; 
		if ( excludeDof.find(0) != excludeDof.end() )
			shiftDof[0] = -1;
		else
			shiftDof[0] = 0;
		for(int i=1; i<nrows; i++)
		{
			if ( excludeDof.find(i) != excludeDof.end() )
				shiftDof[i] = shiftDof[i-1];
			else
				shiftDof[i] = shiftDof[i-1]+1;
		}
		for(int i=0; i<nrows; i++)
		{
			if ( excludeDof.find(i) != excludeDof.end() )
			{}
			else
			{
				for(int j=0; j<rsize[i]; j++)
				{
					if( excludeDof.find(col_ind[i][j]) != excludeDof.end() )
					{
					}
					else
					{
						rmat.val[k] = val[i][j];
						rmat.col_ind[k] = shiftDof[col_ind[i][j]];
						k++;
					}
				}
				rowTrue++;
				rmat.row_ptr[rowTrue] = k; 
			}
		}
		rmat.nnz = k;
		delete [] shiftDof;
	}

	void addMatrix( Matrix<T> subA, int* indx, int* indy = NULL )
	{
		if ( indy == NULL)
			indy = indx;
		for(int i =0; i<subA.rows(); i++)
			for(int j=0; j<subA.cols(); j++)
				this->operator()(indx[i], indy[j]) = this->operator()(indx[i], indy[j]) + subA(i,j);
	}

	void couple( std::vector<std::pair<int,int> >& coupleDofs )
	{
		for( std::pair<int,int> cp : coupleDofs )
		{
			for(int j=0; j<rsize[cp.second]; j++)
			{
				this->set( cp.first, col_ind[cp.second][j],
						val[cp.second][j] + this->get(cp.first, col_ind[cp.second][j]) );
			}
			val[cp.second].clear();
			col_ind[cp.second].clear();
			rsize[cp.second] = 1;

			bool foundTwoInd=false;
			for(int i=0; i<this->rows(); i++)
			{
				for(int j=0; j<col_ind[i].size(); j++)
					if( col_ind[i][j] == cp.first || col_ind[i][j] == cp.second )
					{
						foundTwoInd = false;
						if ( col_ind[i][j] == cp.first )
						{
							for(int j1=j+1; j1<col_ind[i].size(); j1++)
								if(col_ind[i][j1] == cp.second)
								{
									val[i][j] += val[i][j1];
									val[i][j1] = 0.0;
									foundTwoInd=true;
									break;
								}
						}
						else
						{
							for(int j1=j+1; j1<col_ind[i].size(); j1++)
								if ( col_ind[i][j1] == cp.first )
								{
									foundTwoInd = true;
									val[i][j1] += val[i][j];
									val[i][j] = 0.0;
									break;
								}
							if ( !foundTwoInd )
								col_ind[i][j]=cp.first;
						}
						break;
					}
			} // end of sum col slaver to master

		} 
		for ( std::pair<int,int> cp : coupleDofs )
		{
			rsize[cp.second] = 2;
			col_ind[cp.second].push_back(cp.first);
			col_ind[cp.second].push_back(cp.second);
			val[cp.second].push_back(-1.0);
			val[cp.second].push_back(1.0);
			col_ind[cp.second].shrink_to_fit();
			val[cp.second].shrink_to_fit();
		}
	}
#ifdef AMGCL_LIBRARY
	void createAMGSolver()
	{
		if ( isAMG )
		{	
			std::cout<<"AMG solver has been created!"<<std::endl;
			return;
		}

		if ( ncols != nrows )
		{	
			std::cout<<"AMG solver only for square matrix!"<<std::endl;
			return;
		}

		SMatrixCRS<T> CRS_tmp;
		sort_rows();
	    get_CRS_matrix(CRS_tmp);

		solver = new AMGSolver( boost::make_tuple(
		nrows,
		boost::make_iterator_range(CRS_tmp.getRowP(), CRS_tmp.getRowP() + nrows+1),
		boost::make_iterator_range(CRS_tmp.getColP(), CRS_tmp.getColP() + CRS_tmp.getNnz()),
		boost::make_iterator_range(CRS_tmp.getValP(), CRS_tmp.getValP() + CRS_tmp.getNnz())
		) );
		//Solver solve( boost::tie(n, ptr, col, val) );

		isAMG = true;
	}	

	void updateAMGSolver()
	{
		if ( isAMG )
			delete solver;
		SMatrixCRS<T> CRS_tmp;
	    get_CRS_matrix(CRS_tmp);

		solver = new AMGSolver( boost::make_tuple(
		nrows,
		boost::make_iterator_range(CRS_tmp.getRowP(), CRS_tmp.getRowP() + ncols+1),
		boost::make_iterator_range(CRS_tmp.getColP(), CRS_tmp.getColP() + CRS_tmp.getNnz()),
		boost::make_iterator_range(CRS_tmp.getValP(), CRS_tmp.getValP() + CRS_tmp.getNnz())
		) );
		isAMG = true;
	}

	void solve(const amat::Matrix<double> & rhs, amat::Matrix<double> & x, 
		     amat::Matrix<int> & iters, amat::Matrix<double> & error )
	{
		if ( rhs.rows() != nrows )
		{
			std::cout<<"rhs is wrong size!"<<std::endl;
			return;
		}
		
		if ( !isAMG )
		{
			createAMGSolver();
		}

		int itr; double err;
	    for (int i=0; i<rhs.cols(); i++)
		{
			std::vector<double> rhs_l(rhs.rows(),0.0);
			std::vector<double> x_l(rhs.rows(),0.0);
			for( int j=0; j<rhs.rows(); j++)
			{
				rhs_l[j] = rhs(j,i);
				x_l[j] = x(j,i);
			}
			boost::tie(itr, err) = (*solver)(rhs_l, x_l);
			iters(i) = itr;
			error(i) = err;
			for( int j=0; j<rhs.rows(); j++)
				x(j,i) = x_l[j];
		}
	}
#endif

#	ifdef COMPILE_STUPID_STUFF
	/// Creates an Eigen3 sparse matrix in row major format
	void get_Eigen3_rowmajor_matrix( Eigen::MappedSparseMatrix<T, Eigen::RowMajor>& A ) const
	{
		SMatrixCRS<T> rmat;
		rmat.nnz = 0;
		for(int i = 0; i < nrows; i++)
			rmat.nnz += rsize[i];

		rmat.val = new T[rmat.nnz];
		rmat.col_ind = new int[rmat.nnz];
		rmat.row_ptr = new int[nrows+1];
		rmat.allocated = true;
		rmat.row_ptr[0] = 0;

		int i, j, k = 0;
		for(i = 0; i < nrows; i++)
		{
			for(j = 0; j < rsize[i]; j++)
			{
				rmat.val[k] = val[i][j];
				rmat.col_ind[k] = col_ind[i][j];
				k++;
			}
			rmat.row_ptr[i+1] = k;
		}

		/*int i, j, k = 0;
		for(i = 0; i < nrows; i++)
		{
			for(j = 0; j < rsize[i]; j++)
			{
				A.valuePtr()[k] = val[i][j];
				A.innerIndexPtr()[k] = col_ind[i][j];
				k++;
			}
			A.outerIndexPtr()[i+1] = k;
		}*/
	}
#	endif

	/* Computes the LU factorization with partial pivoting.
	 * L is the unit lower triangular matrix, U is the upper triangular matrix.
	 */
	/*void LUfactor(MatrixCRS<double>& L, MatrixCRS<double>& U, MatrixCRS<int>& P)
	{
		// copy this matrix (A) into U
		delete [] U.val;
		delete [] U.col_ind;
		U.val = new vector<double>[nrows];
		U.col_ind = new vector<int>[ncols];
		U.rsize = rsize;
		for(int i = 0; i < nrows; i++)
		{
			U.val[i].resize(rsize[i]);
			U.col_ind[i].resize(rsize[i]);
			for(int j = 0; j < rsize[i]; j++)
			{
				U.val[i][j] = val[i][j];
				U.col_ind[i][j] = col_ind[i][j];
			}
		}

		// Make L and P identity matrices
		delete [] L.val;
		delete [] L.col_ind;
		delete [] P.val;
		delete [] P.col_ind;
		L.val = new vector<double>[nrows];
		L.col_ind = new vector<int>[nrows];
		P.val = new vector<int>[nrows];
		P.col_ind = new vector<int>[nrows];
		for(int i = 0; i < nrows; i++)
		{
			L.rsize.push_back(1);
			P.rsize.push_back(1);
			L.val[i].resize(1);
			L.col_ind[i].resize(1);
			P.val[i].resize(1);
			P.col_ind[i].resize(1);

			L.val[i].reserve(10);
			L.col_ind[i].reserve(10);

			for(int j = 0; j < rsize[i]; j++)
			{
				L.val[i][j] = 1.0;
				L.col_ind[i][j] = i;
				P.val[i][j] = 1;
				P.col_ind[i][j] = i;
			}
		}

		// We can now start Gaussian elimination with partial pivoting
		double max = 0;
		int maxi = 0;
		for(int k = 0; k < nrows-1; k++)
		{
			// first find maximum element in the kth column.
			max = U.val[k][0];
			maxi = 0;
			for(int i = k+1; i < nrows; i++)
			{
				// we need a binary search of the ith row to look for col_ind k.
			}
		}
	}*/
};

class MatrixCRS_traditional : public SparseMatrix<double>
{
	double* val;
	int* col_ind;
	int* row_ptr;
	using SparseMatrix::nnz;
	using SparseMatrix::nrows;
	using SparseMatrix::ncols;
	int valsize;		// total size allocated to array val, and therefore also to array col_ind

public:
	MatrixCRS_traditional(int num_rows, int num_cols) : SparseMatrix<double>(num_rows, num_cols)
	{
		val = new double[num_rows*2];
		col_ind = new int[num_rows*2];
		valsize = num_rows*2;
		for(int i = 0; i < valsize; i++)
		{
			val[i] = 0.0;
			col_ind[i] = -1;
		}
		row_ptr = new int[nrows+1];
		for(int i = 0; i < nrows; i++) row_ptr[i] = -1;
		row_ptr[nrows+1] = 0;
		nnz = 0;
	}
	~MatrixCRS_traditional()
	{
		delete [] val;
		delete [] col_ind;
		delete [] row_ptr;
	}

	void set(int x, int y, double value) // 这个到底是什么鬼！！！
	{
		if(col_ind[valsize-1] != -1)		// re-allocate if no space is available
		{
			double* oldval = val;
			int* oldcol = col_ind;
			val = new double[2*valsize];
			col_ind = new int[2*valsize];
			for(int i = 0; i < valsize; i++)		// copy data over
			{
				val[i] = oldval[i];
				col_ind[i] = oldcol[i];
			}
			for(int i = valsize; i < 2*valsize; i++)		// initialize extra space with dummy data
			{
				val[i] = 0;
				col_ind[i] = -1;
			}
			delete [] oldval;
			delete [] oldcol;
		}

		// if value is not zero

		int ind = row_ptr[x];
		int nind = row_ptr[x+1];
		if(ind == nind)				// only zeros in this row
		{
			for(int i = nind; i < nnz; i++)
			{
				val[i+1] = val[i];				// create space for this element
				col_ind[i+1] = col_ind[i];
			}
		}
	}

	double get(const int x, const int y)
	{
		int ind = row_ptr[x];
		//if(col_ind[ind] > y) return 0.0;
		int nind = row_ptr[x+1];
		for(int i = ind; i <= nind-1; i++)
		{
			if(col_ind[i]==y) return val[i];
		}
		return 0.0;
	}
};

#ifndef SPARSE_MATRIX_TO_USE

/** We use MatrixCRS as the default sparse matrix implementation.
 *
 * \note NOTE that the default sparse matrix type stores floating-point entries of type amc_real
 */
typedef MatrixCRS<amc::amc_real> SpMatrix;

// by zwl
class ZSparseMatProd
{
	protected:
	amat::SpMatrix *M;
public:
	ZSparseMatProd( amat::SpMatrix * mat )
	{
		M = mat ;
	}
    int rows() const { return M->rows(); }
    int cols() const { return M->cols(); }
    // y_out = M * x_in
    void perform_op(const double *x_in, double *y_out)
    {
		M->multiply( x_in, y_out );
    }
};

class SMatrixCRSProd
{
	amat::SMatrixCRS<double> *m;
	int nrows;
	int ncols;
	public:
	SMatrixCRSProd( int r, int c, amat::SMatrixCRS<double> *mm )
	{
		m = mm;
		nrows = r;
		ncols = c;
	}

	int rows() const { return nrows; }

	int cols() const { return ncols; }

	void perform_op ( const double* x_in, double *y_out)
	{
		for(int i=0; i<nrows; i++)
		{
			y_out[i] = 0.0;
			for(int j=m->row_ptr[i]; j<m->row_ptr[i+1]; j++ )
			{
				y_out[i] += m->val[j]*x_in[m->col_ind[j]];
			}
		}
	}
};

#define SPARSE_MATRIX_TO_USE 1
#endif

}	// end namespace
#endif
