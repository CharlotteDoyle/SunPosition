#include "sunpos.h"
#include <iostream>
#include <fstream>
using namespace std;

SunCoord Sun;

double UT(double BT)
{
	return BT - 8;
}

void main()
{
	auto BT = 9.0;

	fstream fout;

	fout.open("d:/data/data.csv");
	Sun = { UT(BT), 22, 10, 2015, 0, 0.6965124321, 2.031412957, 1, 20 };

	while (BT <= 15.0)
	{
		Sun.SetCoord(UT(BT), 22, 10, 2015, 0, 0.6965124321, 2.031412957, 1, 20);
		Sun.Calculate();

		cout << BT << '\t' << Sun.Zenith << endl;
		fout << BT << ',' << Sun.Zenith << endl;
		BT += 0.1;
	}

	fout.close();
	system("pause");
}