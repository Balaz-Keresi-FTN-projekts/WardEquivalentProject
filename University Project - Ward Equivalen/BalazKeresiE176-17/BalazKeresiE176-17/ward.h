#ifndef BALAZ_WARD_H
#define BALAZ_WARD_H

#include <ctime>
#include "matrix.h"

//Klase za vardove razlikuju u nacinu skladenje matrice i po nacinu resenje zadatka
//    Zbog toga ovde su realizovani abstraktnom klasom

class WardEQ
{
	protected:
		bool isSolved;
		string YbusFile;
		string JinjFile;
		int setParams(const char* Yfilename);
		int numONodes;
		int gf, gl;
		string type;
	public:
		WardEQ();
		WardEQ(WardEQ& );
		string getYbusFile() {return YbusFile;}
		string getJinjFile() {return JinjFile;}
		int getNumONodes() {return numONodes; }
		string getType() {return type; }
		virtual double solve()=0;
		virtual void writeOut()=0;
		virtual void exportToMatlab()=0;
};

// Ekvivalent dobijen direktnom inverzijom matrice Yss

class WardInverseEQ : public WardEQ
{
	
		Matrix Yss;
		Matrix Ysg;
		Matrix Ysi;
		Matrix Ygs;
		Matrix Ygg;
		Matrix Ygi;
		Matrix Yis;
		Matrix Yig;
		Matrix Yii;
		Matrix Yggp;
		Vector_col Js;
		Vector_col Jg;
		Vector_col Ji;
		Vector_col Jgp;
		
	public:
		WardInverseEQ();
		WardInverseEQ(WardInverseEQ &);
		WardInverseEQ(const char* Yfilename, const char* Jfilename, int, int) ;
	    double solve();
		void writeOut();
		void exportToMatlab() {}
};


// Elvivalent dobijen gausovom eliminacijom
class WardGaussEQ : public WardEQ
{
		Matrix Ysss; // Velika matrica spoljasnjeg sistema
		Matrix Ysi;
		Matrix Ygi;
		Matrix Yis;
		Matrix Yig;
		Matrix Yii;
		Matrix Yggp;
		Vector_col Jss; // Veliki vektor spoljasnjeg sistema
		Vector_col Ji;
		Vector_col Jgp;
		int argmax(int k);
		void swapSystem (int, int);
	public:
		WardGaussEQ ();
		WardGaussEQ(WardGaussEQ &);
		WardGaussEQ(const char* Yfilename, const char* Jfilename, int, int);
		double solve();
		void writeOut();
		void exportToMatlab() {}
};


// Ekvivalent dobijen LU faktorizacijom
class WardLUEQ : public WardEQ
{
		Matrix Yss;
		Matrix Ysg;
		Matrix Ysi;
		Matrix Ygs;
		Matrix Ygg;
		Matrix Ygi;
		Matrix Yis;
		Matrix Yig;
		Matrix Yii;
		Matrix Yggp;
		Vector_col Js;
		Vector_col Jg;
		Vector_col Ji;
		Vector_col Jgp;
	public:
		WardLUEQ();
		WardLUEQ(WardLUEQ &);
		WardLUEQ(const char* Yfilename, const char* Jfilename, int, int);
	    double solve();
		void writeOut();
		void exportToMatlab();
};


#endif