#include <iostream>
#include <conio.h> // for getch_
#include <locale> // for overladed character functions and getline
#include <fstream>
#include <cmath>
#include "Copenfile.h" // header file of the Copenfile class

using namespace std;

fstream myfile;
Copenfile array_data[4608];//+1 to the array size
int num_items = 0;

double FStrain(double bottom, double top)
{
	double sum_mrr = 0, sum_mtt = 0, sum_mpp = 0, sum_mrt = 0, sum_mrp = 0, sum_mtp = 0, sum_m0 = 0;
	for (int i = 0;i < 4607;i++)  // looping all earthquakes the find the ones that fit the range
	{
		if (bottom < array_data[i].depth  && array_data[i].depth < top) //the ragne of depths the user is interested in
		{
			sum_mrr += array_data[i].mrr *pow(10, array_data[i].iexp); //This is the Kostrov summation
			sum_mtt += array_data[i].mtt *pow(10, array_data[i].iexp);
			sum_mpp += array_data[i].mpp *pow(10, array_data[i].iexp);
			sum_mrt += array_data[i].mrt *pow(10, array_data[i].iexp);
			sum_mrp += array_data[i].mrp *pow(10, array_data[i].iexp);
			sum_mtp += array_data[i].mtp *pow(10, array_data[i].iexp);
		}
	}

	double **A, **I;

	I = new double*[3];
	for (int i = 0;i < 3;i++) I[i] = new double[3];

	for (int i = 0; i < 3;i++) //initialise Identity  matrix
	{
		for (int j = 0; j < 3;j++)
		{
			if (i == j) I[i][j] = 1.;
			else I[i][j] = 0;
		}
	}

	A = new double*[3]; //initialise matrix
	for (int i = 0;i < 3;i++) A[i] = new double[3];

	A[0][0] = sum_mtt ;	A[0][1] = sum_mtp;	A[0][2] = sum_mrt;

	A[1][0] = sum_mtp;	A[1][1] = sum_mpp;	A[1][2] = sum_mrp;

	A[2][0] = sum_mrt;	A[2][1] = sum_mrp;	A[2][2] = sum_mrr;

	

	for (int k = 0;k < 2;k++) // upper triangle - looping over pivot row
	{
		for (int i = (k + 1);i < 3; i++) //looping over the bellow pivot raw
		{
			double s = (A[i][k] / A[k][k]); //scalling factor
			for (int j = 0;j < 3;j++)
			{
				A[i][j] = A[i][j] - s* A[k][j];
				I[i][j] = I[i][j] - s* I[k][j];
			}
		}
	}

	for (int k = 2;k > -1;--k) //Find the Lower triangular matrix
	{
	
		for (int i = (k - 1); i > -1;--i) //Looping over pivor row
		{
			double s = (A[i][k] / A[k][k]); //Looping below pivor row
			for (int j = 0;j<3;j++)
			{
				A[i][j] = A[i][j] - s* A[k][j];
				I[i][j] = I[i][j] - s* I[k][j];
			}
		}
	}
	
	/*The matrix can be printed to confirm that it is diagonalised
	cout << "matrix A  after lower" << endl;
	for (int i = 0; i < 3;i++) 
	{
		for (int j = 0; j < 3;j++)
			cout << A[i][j] << "\t";
		cout << endl;
	}
	*/
	double miu = 3.3e10;
	double Length = 5.6e6;
	double Depth = (35e3);
	double Width = 1.7e6;
	double Volume = Length *Depth *Width;
	double time = 49 * 365 * 24 * 60 * 60;

	double str_devisor = 1 / (2 * miu*time*Volume);

	double pressure_axis = str_devisor* A[0][0];
	double tension_axis = str_devisor* A[2][2];
	double sec_inv= sqrt(pow(pressure_axis, 2) + pow(tension_axis, 2));
	return sec_inv *(365 * 24 * 60 * 60);
}

/*This function calcutates the plate velocity in SI units
and converts it to mm yr-1. The input comes form the summation of the
Seismic moments */
double Fvelocity(double bottom, double top)
{
	double sum_m0 = 0;
	for (int i = 0;i < 4607;i++) //loop over the list to find earthquakes in the range
	{
		if (bottom < array_data[i].depth  && array_data[i].depth < top)
		{
			 sum_m0 += array_data[i].m0; //sum the M0s
		}
	}

	double miu = 3.3e10;
	double Length = 5.6e6;
	double Width = (35e3)*22.5; // Depth* sin dip
	double time = 49 * 365 * 24 * 60 * 60;
	double vel_in_ms = (sum_m0 *1e-7) / (miu*Length*Width*time);// in m/s
	return (vel_in_ms * 1000 * (365 * 24 * 60 * 60));
	
}

void load(void) // Function that load the columns into the memmory for further calculations

{
	while (!myfile.eof())
	{
		myfile >> array_data[num_items].lon;
		myfile >> array_data[num_items].lat;
		myfile >> array_data[num_items].depth;
		myfile >> array_data[num_items].mrr;
		myfile >> array_data[num_items].mtt;
		myfile >> array_data[num_items].mpp;
		myfile >> array_data[num_items].mrt;
		myfile >> array_data[num_items].mrp;
		myfile >> array_data[num_items].mtp;
		myfile >> array_data[num_items].iexp;
		myfile >> array_data[num_items].name_x;
		myfile >> array_data[num_items].name_y;
		myfile >> array_data[num_items].name;
		myfile >> array_data[num_items].m0;
		myfile >> array_data[num_items].mw;

		if (!myfile.eof())
			num_items++;
	}
}

void save(int bottom, int top) // Function that saves the results 

{
	fstream outfile;
	outfile.open("H:\\tempsavedeq.txt", fstream::out);
	if (outfile.fail())
	{
		cout << "Failed to open file for writting " << endl;
		
	}
	
	outfile << "For earthquakes that happened from " << bottom << "km to " << top << "km deep:" << endl;
	outfile <<"The Second Invariant of strain is: " << FStrain(bottom, top) << " / yr" << endl;
	outfile << "According to Kreemer 2003, the second invariant of strain rate in the NW boundary of the pacific is 250 / yr. " << endl;;
	outfile << "This is a difference of " << ((FStrain(bottom, top) - 250)/ 250)*100<<"%" << endl<<endl<<endl;
	
	outfile << "The velocity is : " << Fvelocity(bottom, top) <<"mm / yr"<< endl;
	outfile << "According to Stern, 2002, the velocity that was obtained geodetically is 90 mm/ yr" << endl;
	outfile << "This is a difference of " << ((Fvelocity(bottom, top) - 90) / 90) * 100 << "%" << endl;
	outfile << "A large difference in the seismically and geodedically acquired plate velocities is because the deformation occures seismically." << endl;
	outfile.close();
			
}


void main(void)
{
	myfile.open("H:\\Visual Studio 2015\\Projects\\Lecture1\\Lecture1\\japaneqm0mw.txt", fstream::in);
	if (myfile.fail())
	{
		cout << "Failed to open file" << endl;
		exit(0);
	}
	
	load(); //load the file
	char ans;
	do // loo to allow multiple calcuations
	{
		double *dep_s, *dep_e;
		dep_s = new double;
		dep_e = new double;

		cout << "Plase enter the lower depth of the range the earthquakes you are interested in?" << endl;
		cin >> *dep_s;

		cout << "Plase enter the upper depth of the range earthquakes you are interested in?" << endl;
		cin >> *dep_e;

		cout << "Velocity is : " <<  Fvelocity( *dep_s, *dep_e) << " mm/yr" << endl;
		cout << "the second invariant of strain is: " << FStrain(*dep_s, *dep_e) << endl;
				
		save(*dep_s, *dep_e);

		cout << "Do you want to perform an other calculation (y//n) ?" << endl;
		ans = _getch();
		cin >> ans;
	} while (toupper(ans) == 'Y');
	 if (toupper(ans) == 'N') exit(0);
	
	myfile.close();

}

