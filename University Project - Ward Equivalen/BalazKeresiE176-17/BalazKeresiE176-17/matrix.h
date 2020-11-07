#ifndef BALAZ_MATRIX_H
#define BALAZ_MATRIX_H



#include <complex>

using namespace std;

class Matrix;
class Vector_col;


class Matrix
{
		complex<double>* body;
		int RowNum;
		int ColNum;
	public:
		Matrix();				  // Podrazumevani konstruktor
		Matrix(const char* filename, int hr, int hc, int lr, int lc); // Konstruktor, koji ocitava matricu iz fajla, na osnovu dva coska
																      // (1) Gornji desni cosak
																	  // (2) Donji levi cosak
		Matrix(int n, int m);	  // Kreira praznu matricu dimenzije nxm
		Matrix(const Matrix& M);  // Konstruktor kopije
		~Matrix();

		int getRowNum()const {return RowNum;}
		int getColNum()const {return ColNum;}
		complex<double> getElement(int i, int j) const {return *(body + i*ColNum + j);}
		void setElement(const complex<double> val, int i, int j) {*(body + i*ColNum + j) = val;} 
		void toArrey(complex<double>*) const;
		void setFromArrey (const complex<double>* );
		
		bool isNullMatrix()const;
		void reduceTo(int hr, int hc, int lr, int lc);
		void swap(int, int);

		Matrix& operator= (const Matrix& M2);
		friend Matrix operator+(const Matrix& M1, const Matrix& M2);
		friend Matrix operator-(const Matrix& M1, const Matrix& M2);
		friend Matrix operator*(const Matrix& M1, const Matrix& M2);
		
		friend Vector_col operator*(const Matrix& M, const Vector_col& V);
		Matrix Inverse()const;

		friend ostream& operator<<(ostream&, const Matrix& M);
		void printMatrixRow(ostream&, int i)const; // Ispise i-tu vrstu matrice
		
};
 
class Vector_col {
		complex<double>* body;
		int length;
		int slack;
	public:
		Vector_col();	// Podrazumevani konstruktor
		Vector_col(const char* filename, int hp, int lp); //konstruktor, koji kreira podvektor na vektora sacuvano u fajlu
		Vector_col(int n);								  //konstruktor, koji kreira prazan vektor dimenzije n
		Vector_col(const Vector_col&);					  //konstruktor kopije
		~Vector_col();

		int getLength()const {return length;}
		int getSlack()const {return slack;}
		complex<double> getElement(int i) const {return *(body + i); }
		void setElement(const complex<double> val, int i) { *(body + i) = val; }

		void toArrey(complex<double>*)const;
		void setFromArrey (const complex<double>* );
		
		
		void reduceTo(int hp, int lp);
		void swap(int, int);

		Vector_col& operator= (const Vector_col&);
		friend Vector_col operator+(const Vector_col& V1, const Vector_col& V2);
		friend Vector_col operator-(const Vector_col& V1, const Vector_col& V2);
		friend Vector_col operator*(const Matrix& M, const Vector_col& V);
		complex<double> operator[] (int i) const;
		complex<double>& operator[] (int i);

		friend ostream& operator<<(ostream&, const Vector_col&V);
		void printVectorElem(ostream& out, int i)const;   // Ispisi jedan element vektora
};


#endif