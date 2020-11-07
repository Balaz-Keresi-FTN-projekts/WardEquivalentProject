#include "matrix.h"
#include <fstream>
#include <complex>
#include <iostream>
#include <iomanip>
#include <mkl.h>

using namespace std;


// ***** MATRIX METHODS *********
Matrix::Matrix()
{
	ColNum = 1;
	RowNum = 1;
	body = new complex<double>(0,0);
}

Matrix::Matrix(const char* filename, int hr, int hc, int lr, int lc)
{
	ColNum = 0;
	RowNum = 0;
	if (lr <= hr || lc <= hc)
	{	
		cerr<<"Matrix constructor error: The starting and ending are not given correctly"<<endl;
		system("pause");
		exit(EXIT_SUCCESS);
	}

	ifstream is(filename);
	int i=0;
	int j=0;
	while (!is.eof())
	{

		char ch;
		is>>ch;
		if (ch == '[')
		{
			j = 0;
			ColNum = 0;
			if(i >= hr && i<=lr)
			{
				RowNum++;
			}
			i++;
		}
		if (isdigit(ch) || ch == '-' || ch == '.' )
		{
			double d;
			is.putback(ch);
			is>>d;
			if( j >= hc && j <= lc)
				ColNum++;
			j++;
		}
	}

	is.close();
	body = new complex<double>[ColNum*RowNum];
	
	i=-1;
	j=0;
	int bodyind = 0;
	is.open(filename);
	while (!is.eof())
	{

		char ch;
		is>>ch;
		if (ch == '[')
		{
			j=0;
			i++;
		}
		if (isdigit(ch) || ch == '-' || ch == '.' )
		{
			double d;
			is.putback(ch);
			is>>d;
			if( j >= hc && j <= lc)
				if( i >= hr && i <= lr )
				{
					body[bodyind++] = complex<double>(0,d);
				}
			j++;
		}
	}
	is.close();
}

Matrix::Matrix(int n, int m=0)
{
	if (n<=0) n=1;
	if (m<=0) m=n;
	ColNum = m;
	RowNum = n;
	body = new complex<double>[m*n];
	for (int i=0; i<m*n; i++)
		body[i] = 0;
}

Matrix::Matrix(const Matrix& M)
{
	ColNum = M.getColNum();
	RowNum = M.getRowNum();
	int N = M.getColNum()*M.getRowNum();
	body = new complex<double>[N];
	
	for (int i=0; i<N; i++)
	{
		body[i] = M.body[i];
	}
		
}

Matrix::~Matrix()
{
	delete[] body;
}


void Matrix::toArrey(complex<double>* arr)const
{
	for(int i=0; i<RowNum; i++)
		for(int j=0; j<ColNum; j++)
		{
			*(arr + i*ColNum + j) = *(body + i*ColNum + j);
		}
}

void Matrix::setFromArrey(const complex<double>* arr)
{
	for(int i=0; i<RowNum*ColNum; i++)
		*(body+i) = *(arr+i);
}


bool Matrix::isNullMatrix()const
{
	for(int i=0; i<RowNum*ColNum; i++)
		if(body[i] != complex<double>(0,0) )
			return false;
	return true;
}

void Matrix::reduceTo(int hr, int hc, int lr, int lc)
{
	if( hr > lr || hc > lc )
	{
		cerr<<"Error: Wrong parameters at the matrix reduction"<<endl;
		system("pause");
		exit(EXIT_SUCCESS);
	}
	if (hr < 0 || hc < 0 || lr >= RowNum || lc >= ColNum)
	{
		cerr<<"Error: Wrong parameters at the matrix reduction"<<endl;
		system("pause");
		exit(EXIT_SUCCESS);
	}
	int RowNum2 = lr-hr+1;
	int ColNum2 = lc-hc+1;
	complex<double>* body2 = new complex<double>[RowNum2*ColNum2];
	for (int i=hr; i<=lr; i++)
		for(int j=hc; j<=lc; j++)
		{
			*(body2 + (i-hr)*ColNum2 + j-hc) = *(body + i*ColNum + j);
		}
	delete[] body;
	body = body2;
	RowNum = RowNum2;
	ColNum = ColNum2;
}

void Matrix::swap(int r1, int r2)
{
	complex<double> obj;
	for (int i=0; i<RowNum; i++)
	{
		obj = *(body + ColNum*i + r1);
		*(body + ColNum*i + r1) = *(body + ColNum*i + r2);
		*(body + ColNum*i + r2) = obj;
	}
	for (int i=0; i<ColNum; i++)
	{
		obj = *(body + ColNum*r1 + i);
		*(body + ColNum*r1 + i) = *(body + ColNum*r2 + i);
		*(body + ColNum*r2 + i) = obj;
	}
}


Matrix& Matrix::operator=(const Matrix& M)
{
	ColNum = M.getColNum();
	RowNum = M.getRowNum();
	int N = M.getColNum()*M.getRowNum();
	body = new complex<double>[N];
	
	for (int i=0; i<N; i++)
	{
		body[i] = M.body[i];
	}
	return *this;
}

Matrix operator+(const Matrix &M1, const Matrix &M2)
{
	if ((M1.getRowNum() != M2.getRowNum()) || (M1.getColNum() != M2.getColNum())) 
	{
		cerr<<"There is a problem with the matrix addition. The sizes are diferrent"<<endl;
		cerr<<"\tThe program cannot continue... "<<endl;
		system("pause");
		exit(EXIT_SUCCESS);
	}
	Matrix M(M1);
	for(int i=0; i<M1.getRowNum(); i++)
		for(int j=0; j<M1.getColNum(); j++)
		{
			M.body[i*M1.getColNum()+j] = M1.body[i*M1.getColNum()+j] + M2.body[i*M1.getColNum()+j];	
			cout<<i*M1.getColNum()+j<<endl;
		}

		
	
	return M;
}

Matrix operator-(const Matrix &M1, const Matrix &M2)
{
	if ((M1.getRowNum() != M2.getRowNum()) || (M1.getColNum() != M2.getColNum())) 
	{
		cerr<<"There is a problem at matrix subtraction. The sizes are diferrent"<<endl;
		cerr<<"\tThe program cannot continue... "<<endl;
		system("pause");
		exit(EXIT_SUCCESS);
	}
	Matrix M(M1);
	for(int i=0; i<M1.getRowNum(); i++)
		for(int j=0; j<M1.getColNum(); j++)
		{
			M.body[i*M1.getColNum()+j] = M1.body[i*M1.getColNum()+j] - M2.body[i*M1.getColNum()+j];	
			
		}

		
	
	return M;
}

Matrix operator*(const Matrix &M1, const Matrix &M2)
{
	if (M1.getColNum() != M2.getRowNum()) 
	{
		cerr<<"There is a problem at matrix multiplication. Sizes mismatch"<<endl;
		cerr<<"\tThe program cannot continue... "<<endl;
		system("pause");
		exit(EXIT_SUCCESS);
	}
	const unsigned int n = M2.getColNum();
	const unsigned int m = M1.getRowNum();
	const unsigned int k = M1.getColNum();
	
	Matrix M(m,n);
	

	complex<double> alfa(1,0);
	complex<double> beta(0,0);

	cblas_zgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, &alfa, 
		M1.body, k, M2.body, n, &beta, M.body, n );
	return M;
}

ostream& operator<<(ostream & out, const Matrix& M)
{
	for(int i=0; i<M.RowNum; i++)
	{
		out<<" [ ";
		for(int j=0; j<M.ColNum; j++)
			out<<setw(6)<<(M.body + i*M.ColNum + j)->imag()<<"\t";
		out<<" ]"<<endl;
	}
	return out;
}

void Matrix::printMatrixRow(ostream& out, int i) const
{
	for(int j=0; j<ColNum; j++)
		out<<setw(6)<<(body + i*ColNum + j)->imag()<<",\t";
}

Matrix Matrix::Inverse()const
{
	if (getColNum() != getRowNum()) 
	{
		cerr<<"Error: Cannot invert a non-square matrix"<<endl;
		cerr<<"\tThe program cannot continue... "<<endl;
		system("pause");
		exit(EXIT_SUCCESS);
	}
	int n=getRowNum();
	Matrix inv(*this);
	int* piv = new int[getRowNum()];

	LAPACKE_zgetrf(LAPACK_ROW_MAJOR, n, n, (MKL_Complex16*)inv.body, n, piv);
	LAPACKE_zgetri( LAPACK_ROW_MAJOR, n, (MKL_Complex16*)inv.body, n, piv);

	delete[] piv;
	return inv;
}


// ********* VECTOR METHODS **********

Vector_col::Vector_col()
{
	length = 1;
	body = new complex<double>(0,0);
	slack = -1;
}

Vector_col::Vector_col (const char* filename, int hp, int lp)
{
	length = 0;
	slack = -1;
	if (lp <= hp)
	{	
		cerr<<"Vector constructor error: The starting and ending are not given correctly"<<endl;
		system("pause");
		exit(EXIT_SUCCESS);
	}

	ifstream is(filename);
	int i=0;
	while (!is.eof())
	{
		char ch;
		is>>ch;
		if (ch == ';')
		{
			
			if(i >= hp && i<=lp)
				length++;
			i++;
		}
	}
	is.close();
	body = new complex<double>[length];
	
	i=0;
	int bodyind = 0;
	is.open(filename);
	bool real = true;
	while (!is.eof())
	{
		char ch;
		is>>ch;
		if (ch == ';')
		{
			i++;
		}
		if (isdigit(ch) || ch == '-' || ch == '.' )
		{
			double d;
			is.putback(ch);
			is>>d;
			if( i >= hp && i <= lp)
			{
				if(real)
					body[bodyind] = complex<double>(d,0);
				else
				{
					body[bodyind] = complex<double>(body[bodyind].real(),d);
					bodyind++;
				}
				real = !real;
			}
		}
		if( i >= hp && i <= lp)
			if( ch == '?')
			{
				body[bodyind] = complex<double>(0,0);
				slack = bodyind++;
			}
		}
		is.close();
}

Vector_col::Vector_col(int size)
{
	length = size;
	body = new complex<double>[size];
	for (int i=0; i<length; i++)
		body[i] = 0;
	slack = -1;
}

Vector_col::Vector_col(const Vector_col& V)
{
	length = V.getLength();
	body = new complex<double>[length];
	for (int i=0; i<length; i++)
	{
		body[i] = V.body[i];
	}
	slack = V.slack;
}

Vector_col::~Vector_col()
{
	delete[] body;
}


void Vector_col::toArrey(complex<double>* arr)const
{
	for(int i=0; i<length; i++)
	{
		*(arr + i) = *(body + i);
	}
}

void Vector_col::setFromArrey(const complex<double>* arr)
{
	for(int i=0; i<length; i++)
		*(body+i) = *(arr+i);
}


void Vector_col::printVectorElem(ostream& out, int i)const
{
	if (i == slack)
		out<<"?"<<"\t\t\t;"<<endl;
	else
		out<<setprecision(4)<<body[i].real()<<"\t\t"<<body[i].imag()<<"\t;"<<endl;
}

void Vector_col::reduceTo(int hp, int lp)
{
	if( hp > lp )
	{
		cerr<<"Error: Wrong parameters at the vector reduction"<<endl;
		system("pause");
		exit(EXIT_SUCCESS);
	}
	if (hp < 0 || lp >= length )
	{
		cerr<<"Error: Wrong parameters at the vector reduction"<<endl;
		system("pause");
		exit(EXIT_SUCCESS);
	}
	int length2 = lp-hp+1;
	complex<double>* body2 = new complex<double>[length2];
	for (int i=hp; i<=lp; i++)
	{
		*(body2 + (i-hp) ) = *(body + i);
	}
	delete[] body;
	body = body2;
	length = length2;
	slack = slack - hp;
	if (slack < -1)
		slack = -1;
}

void Vector_col::swap(int r1, int r2)
{
	complex<double> obj;;
	obj = *(body + r1);
	*(body + r1) = *(body + r2);
	*(body + r2) = obj;
}

Vector_col& Vector_col::operator=(const Vector_col& V)
{
	length = V.getLength();
	body = new complex<double>[length];
	for (int i=0; i<length; i++)
	{
		body[i] = V.body[i];
	}
	slack = V.slack;
	return *this;
}

Vector_col operator+ (const Vector_col& V1, const Vector_col& V2)
{
	if ((V1.getLength() != V2.getLength()) ) 
	{
		cerr<<"There is a problem at vector addition. The sizes are diferrent"<<endl;
		cerr<<"\tThe program cannot continue... "<<endl;
		system("pause");
		exit(EXIT_SUCCESS);
	}
	Vector_col V(V1);
	for(int i=0; i<V1.getLength(); i++)
	{	
			V.body[i] = V1.body[i] + V2.body[i];	
	}
	V.slack = V1.slack + V2.slack + 1;
	return V;
}

Vector_col operator- (const Vector_col& V1, const Vector_col& V2)
{
	if ((V1.getLength() != V2.getLength()) ) 
	{
		cerr<<"There is a problem at vector subtraction. The sizes are diferrent"<<endl;
		cerr<<"\tThe program cannot continue... "<<endl;
		system("pause");
		exit(EXIT_SUCCESS);
	}
	Vector_col V(V1);
	for(int i=0; i<V1.getLength(); i++)
	{	
			V.body[i] = V1.body[i] - V2.body[i];	
	}
	V.slack = V1.slack + V2.slack + 1;
	return V;
}

Vector_col operator*(const Matrix& M, const Vector_col& V)
{
	if (M.getColNum() != V.getLength()) 
	{
		cerr<<"There is with the matrix/vector multiplication. Sizes mismatch"<<endl;
		cerr<<"\tThe program cannot continue... "<<endl;
		system("pause");
		exit(EXIT_SUCCESS);
	}
	const unsigned int n = M.getColNum();
	const unsigned int m = M.getRowNum();
	
	Vector_col Vout(m);
	

	complex<double> alfa(1,0);
	complex<double> beta(0,0);

	cblas_zgemv (CblasRowMajor, CblasNoTrans, m, n, &alfa, 
		M.body, n, V.body, 1, &beta, Vout.body, 1);

	Vout.slack = V.slack;
	return Vout;
}

complex<double> Vector_col::operator[] (int i) const
{
	if (i >= length)
	{
		cerr<<"Error: Vector indexing is beyond the vector size"<<endl;
		system("pause");
		exit(EXIT_SUCCESS);
	}
	if (i < 0)
	{
		cerr<<"Error: Vector negative indexing"<<endl;
		system("pause");
		exit(EXIT_SUCCESS);
	}
	return body[i]; 
}

complex<double>& Vector_col::operator[] (int i)
{
	if (i >= length)
	{
		cerr<<"Error: Vector indexing is beyond the vector size"<<endl;
		system("pause");
		exit(EXIT_SUCCESS);
	}
	if (i < 0)
	{
		cerr<<"Error: Vector negative indexing"<<endl;
		system("pause");
		exit(EXIT_SUCCESS);
	}
	return *(body+i);
}


ostream& operator<< (ostream& out, const Vector_col& V)
{
	for (int i=0; i<V.length; i++)
	{
		if (i == V.slack)
			out<<"?"<<"\t\t\t;"<<endl;
		else
			out<<setprecision(4)<<V.body[i].real()<<"\t\t"<<V.body[i].imag()<<"\t;"<<endl;
	}
	return out;
}