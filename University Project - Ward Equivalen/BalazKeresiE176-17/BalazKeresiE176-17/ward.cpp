#include "ward.h"
#include <iostream>
#include <fstream>
#include <mkl.h>

using namespace std;

/********** BASE CLASS **********/
int WardEQ::setParams(const char* Yfilename)
{
	int Ybusdim;
	Matrix Test;
	int ColNum = 0;
	int RowNum = 0;
	// Funkcija iskoristi cinjenicu, da prilikom parametrizacije matrica, konstruktor ne moze da ocitava vise redove od toga koliki se nalazi u fajlu
	//ni u slucaju, kada vrednost argumenta veca od broja redova matrice.
	// Matrica koju ocitavamo teorijski moze da bude ne-kvadratna, zato posebno proverimo uslove za redove i colone. (Ybus je obavezno kvadratna)
	while (Test.getColNum() > ColNum || Test.getRowNum() > RowNum)
	{
		ColNum++ ;
		RowNum++ ;
		Test = Matrix(Yfilename, 0, 0, ColNum+1, RowNum+1);
	}
	
	if (Test.getColNum() != Test.getRowNum())
	{
		cerr<<"Error: The given matrix is not square"<<endl;
		cerr<<"\tProgram finishes..."<<endl;
		system("pause");
		exit(EXIT_SUCCESS);
	}
	Ybusdim = Test.getColNum() - 1; //Zero-based indexing
	return Ybusdim;
} 

WardEQ::WardEQ()
{
	 YbusFile = "";
	 JinjFile = "";
	 isSolved = false;
	 numONodes = 0;
	 gf = 0;
	 gl = 0;
	 type = "";
}

WardEQ::WardEQ(WardEQ& weq)
{
	YbusFile = weq.YbusFile;
	JinjFile = weq.JinjFile;
	isSolved = weq.isSolved;
	numONodes = weq.numONodes;
	gf = weq.gf;
	gl = weq.gl;
	type = "";
}

/********** EQUIVALENT BY INVERSION CLASS *************/

WardInverseEQ::WardInverseEQ() : Yss(), Ysg(), Ysi(), Ygs(), Ygg(), Ygi(), Yis(), Yig(), Yii(), Yggp(),
							    Js(), Jg(), Ji(), Jgp()							   
{
	type = "Inverse";
}

WardInverseEQ::WardInverseEQ(WardInverseEQ &wieq): Yss(wieq.Yss), Ysg(wieq.Ysg), Ysi(wieq.Ysi), 
												   Ygs(wieq.Ygs), Ygg(wieq.Ygg), Ygi(wieq.Ygi),
												  Yis(wieq.Yis), Yig(wieq.Yig), Yii(wieq.Yii), Yggp(wieq.Yggp),
												  Js(wieq.Js), Jg(wieq.Jg), Ji(wieq.Ji), Jgp(wieq.Jgp)
{
	type = "Inverse";
}


WardInverseEQ::WardInverseEQ(const char* Yfilename, const char* Jfilename, int gf, int gl) 
{
	numONodes = setParams(Yfilename);
	isSolved = false;
	type = "Inverse";

	Yss = Matrix(Yfilename, 0, 0, gf-1, gf-1);
	Ysg = Matrix(Yfilename, 0, gf, gf-1, gl-1);
	Ysi = Matrix(Yfilename, 0, gl, gf-1, numONodes);
	Ygs = Matrix(Yfilename, gf, 0, gl-1, gf-1);
	Ygg = Matrix(Yfilename, gf, gf, gl-1, gl-1);
	Ygi = Matrix(Yfilename, gf, gl, gl-1, numONodes);
	Yis = Matrix(Yfilename, gl, 0, numONodes, gf-1);
	Yig = Matrix(Yfilename, gl, gf, numONodes, gl-1);
	Yii = Matrix(Yfilename, gl, gl, numONodes, numONodes);
	Js = Vector_col(Jfilename, 0, gf-1);
	Jg = Vector_col(Jfilename, gf, gl-1);
	Ji = Vector_col(Jfilename, gl, numONodes);

	this->gf = gf;
	this->gl = gl;

	Yggp = Matrix();
	Jgp = Vector_col();

	YbusFile = Yfilename;
	JinjFile = Jfilename;

	if ( !Ysi.isNullMatrix() || !Yis.isNullMatrix() ) //Matrice sadrze samo 0-e kao elementi
	{
		cerr<<"Error: The given nodes are not make a border"<<endl;
		cerr<<"\tProgram finishes..."<<endl;
		system("pause");
		exit(EXIT_SUCCESS);
	}
}

double WardInverseEQ::solve()
{
	double begin = clock();
	Matrix Yssinv = Yss.Inverse();
	Yggp = Ygg - Ygs*Yssinv*Ysg;
	Jgp = Jg - Ygs*Yssinv*Js;
	double end = clock();

	isSolved = true;
	//cout<<"Inverse: "<<"Ward matrix succesfully computed by matrix inversion"<<endl;
	return end-begin;
}

void WardInverseEQ::writeOut()
{
	// Output
	if (isSolved)
	{
		string Jout(JinjFile);
		string Yout(YbusFile);

		int pos = Jout.find_last_of('.');
		Jout.insert(pos, "_Ward_wInverse");
		pos = Yout.find_last_of('.');
		Yout.insert(pos, "_Ward_wInverse");

		// Ispis struje ekvivalenta
		ofstream os;
		os.open(Jout);
		os<<"";
		os.close();
		os.open(Jout, std::ofstream::app);
		os<<Jgp;
		os<<Ji;
		os.close();

		// Ispis Ybus ekvivalenta
		os.open(Yout);
		os<<"";
		os.close();
		os.open(Yout, std::ofstream::app);
		for(int i=0; i<Ygg.getRowNum(); i++)
		{
			os<<"[ ";
			Yggp.printMatrixRow(os,i); Ygi.printMatrixRow(os, i);
			os<<"\t]"<<endl;
		}
		for(int i=0; i<Yig.getRowNum(); i++)
		{
			os<<"[ ";
			Yig.printMatrixRow(os, i); Yii.printMatrixRow(os, i);
			os<<"\t]"<<endl;
		}
		os.close();

		cout<<"Inverse: "<<"The Ward equivalent matrix was succesfully writen into the output file"<<endl;
	}
	else
		cout<<"Inverse: "<<"Can't write an output from this file because it is not made the Ward equivalent yet"<<endl; 
}


/********** EQUIVALENT BY GAUSSIAN ELIMINATION **********/

WardGaussEQ::WardGaussEQ() : Ysss(), Ysi(), Ygi(), Yis(), Yig(), Yii(), Yggp(),
								 Jss(), Ji(), Jgp()							   
{
	 type = "Gauss";
}

WardGaussEQ::WardGaussEQ(WardGaussEQ &wgeq): Ysss(wgeq.Ysss),  Ysi(wgeq.Ysi), Ygi(wgeq.Ygi),
												  Yis(wgeq.Yis), Yig(wgeq.Yig), Yii(wgeq.Yii), Yggp(wgeq.Yggp),
												  Jss(wgeq.Jss), Ji(wgeq.Ji), Jgp(wgeq.Jgp)
{
	type = "Gauss";
}

WardGaussEQ::WardGaussEQ(const char* Yfilename, const char* Jfilename, int gf, int gl) 
{
	numONodes = setParams(Yfilename); 
	isSolved = false;
	type = "Gauss";

	Ysss = Matrix(Yfilename, 0, 0, gl-1, gl-1); // Velika matrica
	Ysi = Matrix(Yfilename, 0, gl, gf-1, numONodes);
	Ygi = Matrix(Yfilename, gf, gl, gl-1, numONodes);
	Yis = Matrix(Yfilename, gl, 0, numONodes, gf-1);
	Yig = Matrix(Yfilename, gl, gf, numONodes, gl-1);
	Yii = Matrix(Yfilename, gl, gl, numONodes, numONodes);
	Jss = Vector_col(Jfilename, 0, gl-1);
	Ji = Vector_col(Jfilename, gl, numONodes);

	this->gf = gf;
	this->gl = gl;

	Yggp = Ysss;
	Jgp = Jss;

	
	YbusFile = Yfilename;
	JinjFile = Jfilename;

	if ( !Ysi.isNullMatrix() || !Yis.isNullMatrix() ) //Matrice sadrze samo 0-e kao elementi
	{
		cerr<<"Error: The given nodes are not make a border"<<endl;
		cerr<<"\tProgram finishes..."<<endl;
		system("pause");
		exit(EXIT_SUCCESS);
	}
}

int WardGaussEQ::argmax(int k)
{
	int i_max = 0;
	double i_maxval = Ysss.getElement(0,k).imag();
	for (int i=0; i<(Ysss.getRowNum()-(gl-gf)); i++)
		if ( Ysss.getElement(i, k).imag() > i_maxval)
		{
			i_max = i;
			i_maxval = Ysss.getElement(i,k).imag();
		}
	return i_max;
}

void WardGaussEQ::swapSystem(int r1, int r2)
{
	Yggp.swap(r1, r2);
	Jgp.swap( r1, r2);
}

double WardGaussEQ::solve()
{
	Yggp = Ysss;
	Jgp = Jss;
	int Yggpcol = Yggp.getColNum();
	int Yggprow = Yggp.getRowNum();
	int Jgplength = Jgp.getLength();
	
	double begin1 = clock();
	complex<double> nullkomp = complex<double>(0,0);
	complex<double> faktor;
	//int i_max;
	for(int k=0; k<gf; k++)
	{
		//i_max = argmax(k);
		//swapSystem(k, i_max);
		for (int i=k+1; i<Yggpcol; i++)
		{
			if( Yggp.getElement(i,k) == nullkomp)
				continue;
			faktor = Yggp.getElement(i,k) / Yggp.getElement(k,k);
			for(int j=k+1; j<Yggprow; j++)
			{
				if (Yggp.getElement(k,j) == nullkomp) 
					continue;
				Yggp.setElement( Yggp.getElement(i,j) - faktor*Yggp.getElement(k,j), i, j);
			}
			Jgp.setElement( Jgp.getElement(i) - faktor * Jgp.getElement(k), i);
			Yggp.setElement(nullkomp , i, k);
		}
	}
	double end1 = clock();

	Yggp.reduceTo(gf, gf, gl-1, gl-1);
	Jgp.reduceTo(gf, gl-1);
	//cout<<"Gauss: Ward matrix is succesfully computed by gaussian elmineation"<<endl;
	isSolved = true;
	return end1 - begin1;
}

void WardGaussEQ::writeOut()
{
	// Output
	if (isSolved)
	{
		string Jout(JinjFile);
		string Yout(YbusFile);

		int pos = Jout.find_last_of('.');
		Jout.insert(pos, "_Ward_wGauss");
		pos = Yout.find_last_of('.');
		Yout.insert(pos, "_Ward_wGauss");

		// Ispis struje ekvivalenta
		ofstream os;
		os.open(Jout);
		os<<"";
		os.close();
		os.open(Jout, std::ofstream::app);
		os<<Jgp;
		os<<Ji;
		os.close();

		// Ispis Ybus ekvivalenta
		os.open(Yout);
		os<<"";
		os.close();
		os.open(Yout, std::ofstream::app);
		for(int i=0; i<Yggp.getRowNum(); i++)
		{
			os<<"[ ";
			Yggp.printMatrixRow(os,i); Ygi.printMatrixRow(os, i);
			os<<"\t]"<<endl;
		}
		for(int i=0; i<Yig.getRowNum(); i++)
		{
			os<<"[ ";
			Yig.printMatrixRow(os, i); Yii.printMatrixRow(os, i);
			os<<"\t]"<<endl;
		}
		os.close();

		
		cout<<"Gauss: "<<"The Ward equivalent matrix was succesfully writen into the output file"<<endl;
	}
	else
		cout<<"Gauss: "<<"Can't write an output from this file because it is not made the Ward equivalent yet"<<endl; 
}


/********** EQUIVALENT BY LU FACTORIZATION ************/

WardLUEQ::WardLUEQ() : Yss(), Ysg(), Ysi(), Ygs(), Ygg(), Ygi(), Yis(), Yig(), Yii(), Yggp(),
							    Js(), Jg(), Ji(), Jgp()							   
{
	type =  "LU factorization";
}

WardLUEQ::WardLUEQ(WardLUEQ &wlueq): Yss(wlueq.Yss), Ysg(wlueq.Ysg), Ysi(wlueq.Ysi), 
									     Ygs(wlueq.Ygs), Ygg(wlueq.Ygg), Ygi(wlueq.Ygi),
									     Yis(wlueq.Yis), Yig(wlueq.Yig), Yii(wlueq.Yii), Yggp(wlueq.Yggp),
									     Js(wlueq.Js), Jg(wlueq.Jg), Ji(wlueq.Ji), Jgp(wlueq.Jgp)
{
	type =  "LU factorization";
}

WardLUEQ::WardLUEQ(const char* Yfilename, const char* Jfilename, int gf, int gl) 
{
	numONodes = setParams(Yfilename);
	isSolved = false;
	type =  "LU factorization";

	Yss = Matrix(Yfilename, 0, 0, gf-1, gf-1);
	Ysg = Matrix(Yfilename, 0, gf, gf-1, gl-1);
	Ysi = Matrix(Yfilename, 0, gl, gf-1, numONodes);
	Ygs = Matrix(Yfilename, gf, 0, gl-1, gf-1);
	Ygg = Matrix(Yfilename, gf, gf, gl-1, gl-1);
	Ygi = Matrix(Yfilename, gf, gl, gl-1, numONodes);
	Yis = Matrix(Yfilename, gl, 0, numONodes, gf-1);
	Yig = Matrix(Yfilename, gl, gf, numONodes, gl-1);
	Yii = Matrix(Yfilename, gl, gl, numONodes, numONodes);
	Js = Vector_col(Jfilename, 0, gf-1);
	Jg = Vector_col(Jfilename, gf, gl-1);
	Ji = Vector_col(Jfilename, gl, numONodes);

	this->gf = gf;
	this->gl = gl;

	Yggp = Matrix();
	Jgp = Vector_col();

	YbusFile = Yfilename;
	JinjFile = Jfilename;

	if ( !Ysi.isNullMatrix() || !Yis.isNullMatrix() ) //Matrice sadrze samo 0-e kao elementi
	{
		cerr<<"Error: The given nodes are not make a border"<<endl;
		cerr<<"\tProgram finishes..."<<endl;
		system("pause");
		exit(EXIT_SUCCESS);
	}
}





double WardLUEQ::solve()
{
	complex<double>* Yss_arrey = new complex<double>[Yss.getRowNum()*Yss.getColNum() ];
	Yss.toArrey(Yss_arrey);
	complex<double>* Ysg_arrey = new complex<double>[Ysg.getColNum()*Ysg.getRowNum()];
	Ysg.toArrey(Ysg_arrey);
	int pivsize = (Ysg.getColNum() > Ysg.getRowNum()) ? (Ysg.getColNum()) : (Ysg.getRowNum());
	int* ipiv = new lapack_int[pivsize];

	double begin1 = clock();
	LAPACKE_zgetrf ( LAPACK_ROW_MAJOR, Yss.getRowNum(), Yss.getColNum(),  (MKL_Complex16*)Yss_arrey, Yss.getColNum(), ipiv);

	LAPACKE_zgetrs ( LAPACK_ROW_MAJOR, 'N', Yss.getRowNum(), Ysg.getColNum(), (MKL_Complex16*)Yss_arrey, Yss.getRowNum(), ipiv, (MKL_Complex16*)Ysg_arrey,
		Ysg.getColNum());
	double end1 = clock();

	Matrix W(Ysg.getRowNum(), Ysg.getColNum() );
	W.setFromArrey(Ysg_arrey);
	delete[] Ysg_arrey;

	double begin2 = clock();
	Yggp = Ygg - Ygs*W;
	double end2 = clock();

	// Yss je vec faktorizovana zato moze nju iskorisiti za kalkulaciju struje
	complex<double>* Js_arrey = new complex<double>[Js.getLength() ];
	Js.toArrey(Js_arrey);

	double begin3 = clock();
	LAPACKE_zgetrs ( LAPACK_ROW_MAJOR, 'N', Yss.getRowNum(), 1, (MKL_Complex16*)Yss_arrey, Yss.getRowNum(), ipiv, (MKL_Complex16*)Js_arrey, 1);
	double end3 = clock();

	Vector_col Z(Js.getLength() );
	Z.setFromArrey(Js_arrey);
	delete[] Js_arrey;

	double begin4 = clock();
	Jgp = Jg - Ygs*Z;
	double end4 = clock();

	delete[] Yss_arrey;
	delete[] ipiv;

	//cout<<"LU factorization: Ward matrix is succesfully computed by LU factorization"<<endl;
	isSolved = true;
	return (end1-begin1) + (end2-begin2) + (end3-begin3) + (end4-begin4);
}



void WardLUEQ::writeOut()
{
	// Output
	if (isSolved)
	{
		string Jout(JinjFile);
		string Yout(YbusFile);

		int pos = Jout.find_last_of('.');
		Jout.insert(pos, "_Ward_wLUfactorization");
		pos = Yout.find_last_of('.');
		Yout.insert(pos, "_Ward_wLUfactorization");

		// Ispis struje ekvivalenta
		ofstream os;
		os.open(Jout);
		os<<"";
		os.close();
		os.open(Jout, std::ofstream::app);
		os<<Jgp;
		os<<Ji;
		os.close();

		// Ispis Ybus ekvivalenta
		os.open(Yout);
		os<<"";
		os.close();
		os.open(Yout, std::ofstream::app);
		for(int i=0; i<Yggp.getRowNum(); i++)
		{
			os<<"[ ";
			Yggp.printMatrixRow(os,i); Ygi.printMatrixRow(os, i);
			os<<"\t]"<<endl;
		}
		for(int i=0; i<Yig.getRowNum(); i++)
		{
			os<<"[ ";
			Yig.printMatrixRow(os, i); Yii.printMatrixRow(os, i);
			os<<"\t]"<<endl;
		}
		os.close();

		cout<<"LU factorization: "<<"The Ward equivalent matrix was succesfully writen into the output file"<<endl;
	}
	else
		cout<<"LU factorization: "<<"Can't write an output from this file because it is not made the Ward equivalent yet"<<endl; 
}

void WardLUEQ::exportToMatlab()
{
	if (isSolved)
	{
		string Jout(JinjFile);
		string Yout(YbusFile);
		string outputname = YbusFile;
		int pos = outputname.find_last_of('.');
		outputname.insert(pos, "_MATLAB");
		pos = outputname.find_last_of('.');
		outputname.erase(pos, pos+3);
		outputname.append(".m");
		
		// Print out original

		ofstream os;
		os.open(outputname, ofstream::out);
		os<<"Ybusdom = [..."<<endl;

		for (int i=0; i<numONodes+1; i++)
		{
			if (i==gf+Jg.getSlack())
				continue;
			for (int j=0; j<numONodes+1; j++)
			{
				if (j == gf+Jg.getSlack())
					continue;
				if (i < gf && j <gf)
					os<<Yss.getElement(i,j).imag()<<"   ";
				else if (i < gf && j >= gf && j < gl)
					os<<Ysg.getElement(i,j-gf).imag()<<"   ";
				else if (i < gf && j >= gl)
					os<<Ysi.getElement(i,j-gl).imag()<<"   ";

				else if (i >= gf && i < gl && j < gf)
					os<<Ygs.getElement(i-gf,j).imag()<<"   ";
				else if (i >= gf && i < gl && j >= gf && j < gl)
					os<<Ygg.getElement(i-gf,j-gf).imag()<<"   ";
				else if (i >= gf && i < gl && j >= gl )
					os<<Ygi.getElement(i-gf,j-gl).imag()<<"   ";

				else if (i >= gl && j < gf )
					os<<Yis.getElement(i-gl,j).imag()<<"   ";
				else if (i >= gl && j >= gf && j < gl)
					os<<Yig.getElement(i-gl,j-gf).imag()<<"   ";
				else if (i >= gl && j >= gl)
					os<<Yii.getElement(i-gl,j-gl).imag()<<"   ";				
			}
			if (i == numONodes )
				os<<"]";
			os<<";"<<endl;
		}
		
		os<<endl;

		os<<"Ywbusdom = [..."<<endl;
		int g = gl - gf;
		for (int i=0; i<numONodes-gf+1; i++)
		{
			if (i==Jgp.getSlack())
				continue;
			for (int j=0; j<numONodes-gf+1; j++)
			{
				if (j == Jg.getSlack())
					continue;

				if (i < g && j < g )
					os<<Yggp.getElement(i,j).imag()<<"   ";
				else if (i < g && j >= g )
					os<<Ygi.getElement(i,j-g).imag()<<"   ";

				else if (i >= g  && j < g)
					os<<Yig.getElement(i-g,j).imag()<<"   ";
				else if (i >= g && j >= g )
					os<<Yii.getElement(i-g,j-g).imag()<<"   ";				
			}
			if (i== numONodes -gf )
				os<<"]";
			os<<";"<<endl;
		}

		os<<endl;
		os<<endl;

		os<<"Ybusdom = complex(0, Ybusdom);"<<endl;
		os<<"Ywbusdom = complex(0, Ywbusdom);"<<endl;

		os<<endl;
		os<<"Jinjdom = [..."<<endl;

		for (int j=0; j<numONodes+1; j++)
		{
			if (j == gf+Jg.getSlack())
				continue;
			if (j <gf)
				os<<"\tcomplex"<<Js[j] - complex<double>(110,0)*Ysg.getElement(j,Jg.getSlack())<<"\t" ;
			else if (j >= gf && j < gl)
				os<<"\tcomplex"<<Jg[j-gf] - complex<double>(110,0)*Ygg.getElement(j-gf,Jg.getSlack())<<"\t" ;
			else if (j >= gl )
				os<<"\tcomplex"<<Ji[j-gl] - complex<double>(110,0)*Yig.getElement(j-gl,Jg.getSlack())<<"\t" ;
			if (j==numONodes)
				os<<"];"<<endl;
			else
				os<<";"<<endl;
		}

		os<<endl;
		os<<"Jwinjdom = [..."<<endl;
		g = gl - gf;
		for (int j=0; j<numONodes-gf+1; j++)
		{
			if (j == Jg.getSlack())
				continue;
			if (j < g)
				os<<"\tcomplex"<<Jgp[j] - complex<double>(110,0)*Yggp.getElement(j,Jgp.getSlack())<<"\t" ;
			else if (j >= g)
				os<<"\tcomplex"<<Ji[j-g] - complex<double>(110,0)*Yig.getElement(j-g,Jgp.getSlack())<<"\t" ;
			if (j==numONodes-gf)
				os<<"];"<<endl;
			else
				os<<";"<<endl;
		}

	os<<endl;
	os<<"Yslack = [..."<<endl;
	for (int i=0; i<numONodes+1; i++)
		{
			if (i == Jg.getSlack()+gf)
				continue;
			if (i < gf)
				os<<Ygs.getElement(Jg.getSlack(),i).imag()<<"  " ;
			else if (i >= gf && i < gl)
				os<<Ygg.getElement(Jg.getSlack(),i-gf).imag()<<"  " ;
			else if (i >= gl )
				os<<Ygi.getElement(Jg.getSlack(),i-gl).imag()<<"  " ;
		}
	os<<"];"<<endl;
		
	os<<endl;
	os<<"Ywslack = [..."<<endl;
	g = gl - gf;
	for (int i=0; i<numONodes-gf+1; i++)
		{
			if (i == Jg.getSlack())
				continue;
			if (i < g)
				os<<Yggp.getElement(Jg.getSlack(),i).imag()<<"  " ;
			else if (i >= g)
				os<<Ygi.getElement(Jg.getSlack(),i-g).imag()<<"  " ;
		}
	os<<"];"<<endl;
	
	os<<"Yslack = complex(0, Yslack);"<<endl;
	os<<"Ywslack = complex(0, Ywslack);"<<endl;


	os<<endl;
	os<<"Vorig = linsolve(Ybusdom, Jinjdom);"<<endl;
	os<<"Vw = linsolve(Ywbusdom, Jwinjdom);"<<endl;

	os<<endl;
	os<<"slack = "<<Jg.getSlack() + 1<<";"<<endl;
	os<<"islack_w = Ywslack*Vw + Ywbusdom(slack,slack) * 110;"<<endl;
	os<<"slack = "<<Jg.getSlack() + gf + 1<<";"<<endl;
	os<<"islack_orig = Yslack*Vorig + Ybusdom(slack,slack) * 110;"<<endl;
	os<<"length = "<<numONodes <<";"<<endl;

	os<<endl;
	os<<"Verr = (abs(Vorig(slack:length))-abs(Vw))./abs(Vorig(slack:length));"<<endl;
	os<<"Fierr = (angle(Vorig(slack:length))-angle(Vw));"<<endl;
	os<<"islack_err = (abs(islack_orig) - abs(islack_w))./abs(islack_orig);"<<endl;
	os<<"islack_ang_err = (angle(islack_orig) - angle(islack_w));"<<endl;
	}
	else
		cout<<"LU factorization: "<<"Can't write an output from this file because it is not made the Ward equivalent yet"<<endl; 

}
