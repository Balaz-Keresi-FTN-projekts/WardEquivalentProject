//  ************* DOMACI ZADATAK 2 *******************
//
// Balaz Keresi E176/2017
//
// Program na osnovu fajlova inicijalizacije (Ybus.txt, Iinj.txt) ocitava parametre mreze: Matricu admitansi i struje injektiranja cvorova
//  Na osnovu ocitavane podatke napravi Ward-ov ekvivalent mreze (detaljnije objasnjen u dokumentu) i izpise model tog ekvivalenta u jednom fajlu.
//  Za proracun koristio sam 2 klase: Matrix i Vector_col, kojima sam realizovao postupke za matricne i vektorske proracune.


#include <iostream>
#include <fstream>
#include "ward.h"
#include <ctime>

using namespace std;


int main()
{
	cout<<"Initialization..."<<endl;
	// Pozicija granicinih cvorova - zero based indexing
	int gf = 5; // Ybigbus: 21, Ybus: 5
	int gl = 8; // Ybigbus: 24, Ybus: 8
	const char YbusFile[] = "Ybus.txt";
	const char JinjFile[] = "Jinj.txt";
	int gf2 = 21; // Ybigbus: 21, Ybus: 5
	int gl2 = 24; // Ybigbus: 24, Ybus: 8
	const char YbusFile2[] = "Ybigbus.txt";
	const char JinjFile2[] = "Jbiginj.txt";
	
	WardEQ** weq = new WardEQ*[6];
	weq[0] = new WardInverseEQ(YbusFile, JinjFile, gf, gl);
	weq[1] = new WardGaussEQ(YbusFile, JinjFile, gf, gl);
	weq[2] = new WardLUEQ(YbusFile, JinjFile, gf, gl);
	weq[3] = new WardInverseEQ(YbusFile2, JinjFile2, gf2, gl2);
	weq[4] = new WardGaussEQ(YbusFile2, JinjFile2, gf2, gl2);
	weq[5] = new WardLUEQ(YbusFile2, JinjFile2, gf2, gl2);



	cout<<"Speed test begins"<<endl;
	cout<<"n = "<<weq[0]->getNumONodes() + 1<<endl;

	for (int i=0; i<6; i++)
	{		
		double s=0;
		if (i==3)
			cout<<"n = "<<weq[3]->getNumONodes()+1<<endl;
		for (int j=0; j<1000; j++)
			{
				s += weq[i]->solve();
			}
		cout<<weq[i]->getType()<<": time needed for 1000 calculations: "<<s<<" ms"<<endl;
		weq[i]->writeOut();
	}
	weq[2]->exportToMatlab();
	weq[5]->exportToMatlab();
	
	
	delete[] weq;

	system("pause");
	return 0;
}
