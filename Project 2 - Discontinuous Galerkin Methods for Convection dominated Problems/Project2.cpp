// 20:21 07
// ==============================================     PROJECT 2     =========================================================== 
// ================           Discontinuous Galerkin Methods for Convection-dominated problems           ======================
// ----------------------------------------------------------------------------------------------------------------------------

#include<iostream>
#include<chrono>
#include<fstream>
#include<cmath>
#include<vector>
#include<iomanip>
#include<sstream>
#include<cassert>
#include<cstdlib>
#include<ctime>
#include <algorithm>
#ifndef MVECTOR_H // the 'include guard'
#define MVECTOR_H // see C++ Primer Sec. 2.9.2
// Class that represents a mathematical vector
class MVector
{
public:
	// constructors
	MVector() {}
	explicit MVector(int n) : v(n) {}
	MVector(int n, double x) : v(n, x) {}
	MVector(std::initializer_list<double> l) : v(l) {}

	// access element (lvalue) (see example sheet 5, q5.6)
	double& operator[](int index)
	{
		return v[index];
	}

	// access element (rvalue) (see example sheet 5, q5.7)
	double operator[](int index) const {
		return v[index];
	}

	int size() const { return v.size(); } // number of elements
	// Define addition for vectors of size two
	
	MVector operator+(MVector w)
	{
		MVector returnVal = { v[0] + w[0],v[1] + w[1] };
		return returnVal;
	}
	// Defining scalar multiplication for vector (2 by 2, only)
	friend MVector operator*(double f, const MVector v)
	{
		MVector returnVal = {f*v[0],f*v[1]};
		return returnVal;
	}

private:
	std::vector<double> v;
};
#endif

#ifndef MMATRIX_H // the 'include guard'
#define MMATRIX_H
// Class that represents a mathematical matrix
class MMatrix
{
public:
	// constructors
	MMatrix() : nRows(0), nCols(0) {}
	MMatrix(int n, int m, double x = 0) : nRows(n), nCols(m), A(n* m, x) {}

	// set all matrix entries equal to a double
	MMatrix& operator=(double x)
	{
		for (unsigned i = 0; i < nRows * nCols; i++) A[i] = x;
		return *this;
	}

	// access element, indexed by (row, column) [rvalue]
	double operator()(int i, int j) const
	{
		return A[j + i * nCols];
	}

	// access element, indexed by (row, column) [lvalue]
	double& operator()(int i, int j)
	{
		return A[j + i * nCols];
	}
	// Calculate an inverse for a 2 by 2 matrix 
	MMatrix invert(MMatrix A)
	{
		MMatrix returnMatrix(2,2,0.0);
		double det = A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1);
		if (det == 0)
		{
			return A;
		}
		returnMatrix(0,0) = (1.0 / det) * A(0, 0);
		returnMatrix(0, 1) = (-1.0 / det) * A(0, 1);
		returnMatrix(1, 0) = (-1.0 / det) * A(1, 0);
		returnMatrix(1, 1) = (1.0 / det) * A(1, 1);
		return returnMatrix;
	}
	// Matrix vector multiplication (For 2 by 2 matrices and 2 by 1 vectors)
	friend MVector operator*(MMatrix A, MVector V)
	{
		MVector returnVec(2);
		returnVec[0] = A(0, 0) * V[0] + A(0, 1) * V[1];
		returnVec[1] = A(1, 0) * V[0] + A(1, 1) * V[1];
		return returnVec;
	}

	// size of matrix
	int Rows() const { return nRows; }
	int Cols() const { return nCols; }

private:
unsigned int nRows, nCols;
std::vector<double> A;
};
#endif
//===================================================================================================
void display_vector(std::vector<double>& v)
{
	for (int i = 0; i < v.size(); i++)
	{
		std::cout << v[i] << " ";
	}
	std::cout << std::endl;
}
void display_MVector(MVector v)
{
	for (int j = 0; j < v.size(); j++)
	{
		std::cout << v[j] << " ";
	}
	std::cout << std::endl;
}
void display_MMatrix(MMatrix A)
{}
//===================================================================================================

// ==========================================================  Inintialisation  ==========================================================
/**/
class AdvectionElement
{
public:
	// Pointer to the left neighbour
	AdvectionElement* Left_neighbour_pt;
	// Pointer to the right neighbour
	AdvectionElement* Right_neighbour_pt;
	// Storage for the coordinates
	std::vector<double> X;
	// Storage for the unknowns
	std::vector<double> U;
	// Constructor: initialise the vectors to hold two entries.
	AdvectionElement()
	{
		// Resize the vectors to hold two entries each
		X.resize(2);
		U.resize(2);
	}
	// Return the value of the coordinate at local coordinate s using
	// equation (1.2)
	double interpolated_x(double s)
	{
		double psi_0 = 0.5 * (1 - s);
		double psi_1 = 0.5 * (1 + s);
		return X[0] * psi_0 + X[1] * psi_1;
	}
	// Return the value of the unknown at local coordinate s using
	// equation (1.4)
	double interpolated_u(double s)
	{
		return U[0] * 0.5 * (1 - s) + U[1] * 0.5 * (1 + s);
	}
	// ================  Timestepping loop   ===================
	virtual double flux(double u)  // This is a virtual fn so it can be overloaded later
	{
		if (u == interpolated_u(-1.0)) { u = -1.0; }
		else if (u == interpolated_u(1.0)) { u = 1.0; }
		else {
			if (u>U[0] && u<U[1])
			{
				double s = (2 * u - (U[0] + U[1])) / (U[1] - U[0]); // This doesn't work if using a discontinuous IC
				//std::cout << "s= " << s << std::endl;
				u=s;
			}
			else {
				return u; // "inbetween elements"
			}
		}
		double returnVal = interpolated_u(u); // f(u)=u
		return returnVal;
	}
	// A function that returns the integral of the flux function over the element using two-point Gauss rule
	double integrate_flux()
	{
		// Using two-point Gauss rule
		return flux(interpolated_u(-1.0 / sqrt(3))) + flux(interpolated_u(1.0 / sqrt(3)));
	}
	// A function h(a,b) that returns the numerical flux. Local Lax-Friedrichs flux
	virtual double h(double a, double b)
	{
		// Lax-Friedrichs flux
		double max = 1.0;
		// Return the quantity 
		return 0.5 * (flux(a) + flux(b)) - 0.5 * max * (b - a);
	}
	MVector U_updated;
	// A function that calculates the updated values of the unknowns U using first-order scheme
	void timestep(double dt)
	{
		MVector U_timestep(2),Flux(2),F_e(2),U_;
		MMatrix M_e(2,2,0.0);
		double coef = (X[1] - X[0]) / 6;
		M_e(0, 0) = coef * 2.0;
		M_e(0, 1) = M_e(1, 0) = coef * 1.0;
		M_e(1, 1) = coef * 2.0;
		F_e[0] = -0.5 * integrate_flux();
		F_e[1] = 0.5 * integrate_flux();
		U_ = { U[0],U[1] };
		if ((*Left_neighbour_pt).U[1] == U[0])
		{
			// Continuous
			Flux[0] = flux(U_[0]);
		}
		else
		{
			// Discontinuous
			Flux[0] = h((*Left_neighbour_pt).U[1], U[0]);
		}
		if ((*Right_neighbour_pt).U[0] == U[1])
		{
			// Continuous
			Flux[1] = -flux(U[1]);
		}
		else
		{
			// Discontinuous
			Flux[1] = -h(U[1],(*Right_neighbour_pt).U[0]);
		}
		// Vector equation which calculates updates
		U_timestep = U_ + dt *( M_e.invert(M_e) * (F_e+Flux));
		U_updated=U_timestep;
	}
	
}; //End of the class definition
// =======================================================   Derived class   ===================================================================
// Creating a new class called BurgersElement, which is required to solve Burger's equation using discontinuous Galerkin
// methods. The approriate flux and numerical flux function need to be overloaded for this to work!
class BurgersElement : public AdvectionElement
{
	// This BurgersElement inherits from the AdvectionElements class
public:
	// Add existing members function from base class
		// Pointer to the left neighbour
	BurgersElement* Left_neighbour_pt;
	// Pointer to the right neighbour
	BurgersElement* Right_neighbour_pt;
	// Storage for the coordinates
	std::vector<double> X;
	// Storage for the unknowns
	std::vector<double> U;
	// Constructor: initialise the vectors to hold two entries.
	BurgersElement()
	{
		// Resize the vectors to hold two entries each
		X.resize(2);
		U.resize(2);
	}
	// Return the value of the coordinate at local coordinate s using
	// equation (1.2)
	double interpolated_x(double s)
	{
		double psi_0 = 0.5 * (1 - s);
		double psi_1 = 0.5 * (1 + s);
		return X[0] * psi_0 + X[1] * psi_1;
	}
	// Return the value of the unknown at local coordinate s using
	// equation (1.4)
	double interpolated_u(double s)
	{
		return U[0] * 0.5 * (1 - s) + U[1] * 0.5 * (1 + s);
	}
	// ================  Timestepping loop   ===================
	virtual double flux(double u) 
	{
		if (u == interpolated_u(-1.0)) { u = -1.0; }
		else if (u == interpolated_u(1.0)) { u=1.0; }
		else { return 0.5 * pow(u, 2); 

		}
		double returnVal = 0.5*pow(interpolated_u(u), 2);
		return returnVal;
	}
	// A function that returns the integral of the fluz function over the element using two-point Gauss rule
	double integrate_flux()
	{
		// Using two-point Gauss rule
		return flux(interpolated_u(-1.0 / sqrt(3.0))) + flux(interpolated_u(1.0 / sqrt(3.0)));
	}
	// A function h(a,b) that returns the numerical flux. Local Lax-Friedrichs flux 
	virtual double h(double a, double b)
	{
		std::vector<double> v = { std::abs(a), std::abs(b) };
		double max = 0.0;
		for (int i = 0; i < v.size(); i++)
		{
			if (v[i] > max)
			{
				max = v[i];
			}
		}
		// Return the quantity std::cout << 0.5 * (flux(a) + flux(b)) - 0.5 * max * (b - a) << std::endl;
		return 0.5 * (flux(a) + flux(b)) - 0.5 * max * (b - a);
	}
	MVector U_updated;
	// A function that calculates the updated values of the unknowns U using first-order scheme
	void timestep(double dt)
	{
		MVector U_timestep(2), Flux(2), F_e(2), U_;
		MMatrix M_e(2, 2, 0.0);
		double coef = (X[1] - X[0]) / 6;
		M_e(0, 0) = coef * 2.0;
		M_e(0, 1) = M_e(1, 0) = coef * 1.0;
		M_e(1, 1) = coef * 2.0;
		F_e[0] = -0.5 * integrate_flux();
		F_e[1] = 0.5 * integrate_flux();
		U_ = { U[0], U[1] };
		if ((*Left_neighbour_pt).U[1] == U[0])
		{
			// Continuous
			Flux[0] = flux(U_[0]);
		}
		else
		{
			// Discontinuous
			Flux[0] = h((*Left_neighbour_pt).U[1], U[0]);
		}
		if ((*Right_neighbour_pt).U[0] == U[1])
		{
			// Continuous
			Flux[1] = -flux(U_[1]);
		}
		else
		{
			// Discontinuous
			Flux[1] = -h(U[1],(*Right_neighbour_pt).U[0]);
		}
		// Vector equation which calculates updates
		U_timestep = U_ + dt * (M_e.invert(M_e) * (F_e + Flux));
		//display_MVector(dt * (M_e.invert(M_e) * (F_e + Flux)));
		U_updated = U_timestep;
	}
}; // End of class defninition


// ===================================================================================================================================================================================================
int main()
{
	// =======================================       Question 3  -  Discontinuous profile       =============================================
	// Discontinuous initial profile
	int N = 100;
	double x_start = 0.0, x_end = 2 * std::acos(-1.0);
	std::vector<BurgersElement> elements(N);
	for (int j = 0; j < N; j++)
	{
		// Loop over vector to fill out member data X, U
		// Initialise X values
		elements[j].X[0] = x_start + j * (x_end - x_start) / N;
		elements[j].X[1] = x_start + (j + 1) * (x_end - x_start) / N;
		// Using square-wave profile
		if (elements[j].X[0] >= 0.0 && elements[j].X[0] <= 1.0)
		{
			elements[j].U[0] = 1.0;
			elements[j].U[1] = 1.0;
		}
		else
		{
			//std::cout << "Discontinuity occurs at " << j << std::endl;
			elements[j].U[0] = 0.0;
			elements[j].U[1] = 0.0;
		}

	}
	//elements[N - 1].U[1] = elements[0].U[0];  // Connecting the first and last elements U^0_0 = U^{N-1}_1
	for (int i = 0; i < N; i++)
	{
		// Setting neighbour pointers for each element
		if (i == 0)
		{
			// set neighbour pointer connecting first and last element
			elements[i].Left_neighbour_pt = &elements[N - 1];
			elements[i].Right_neighbour_pt = &elements[i + 1];
		}
		else if (i == N - 1)
		{
			// set neighbour pointer connecting first and last element
			elements[i].Right_neighbour_pt = &elements[0];
			elements[i].Left_neighbour_pt = &elements[i - 1];
		}
		else
		{
			elements[i].Left_neighbour_pt = &elements[i - 1];// Left pointer
			elements[i].Right_neighbour_pt = &elements[i + 1];// Right pointer
		}
	}

	// Some testing
	//std::cout << elements[15].U[1] << std::endl; // =1
	//std::cout << elements[16].U[0] << std::endl; // =0

	//std::cout << elements[16].Left_neighbour_pt->U[1] << std::endl;
	//std::cout << elements[16].U[0] << std::endl;
	//std::cout << elements[0].flux((elements[16].U[0] + elements[15].Left_neighbour_pt->U[1])/2) << std::endl;
	//std::cout << elements[0].h(elements[15].Left_neighbour_pt->U[1], elements[16].U[0]) << std::endl;
	//elements[16].timestep(0.001);
	//std::cout << elements[16].U_updated[0] << std::endl;


	double dt = 0.0005, total_time=0.0;
	for (int k = 0; k < 4000; k++)
	{
		if (k == 0) //t=0
		{
			std::cout << "=========   " << total_time << "   ============" << std::endl;
			std::ofstream Q3square0;
			Q3square0.open("Q3square0.txt");
			if (!Q3square0)
			{
				std::cout << "Couldn't open file!!" << std::endl;
				return 1;
			}
			Q3square0 << elements[0].U[0] << std::endl;
			for (int i = 0; i < N; i++)// Write to file
			{
				//Q3square0 << elements[i].U[0] << std::endl;
				Q3square0 << elements[i].U[1] << std::endl;
			}
			Q3square0.close();
		}
		for (int i = 0; i < N; i++)
		{
			elements[i].timestep(dt);
		}
		for (int j = 0; j < N; j++)
		{
			elements[j].U[0] = elements[j].U_updated[0];
			elements[j].U[1] = elements[j].U_updated[1];
		}
		total_time = total_time + dt;

		/*
		double u_int = 0.0;
		for (int l = 0; l < N - 1; l++)
		{
			u_int = 0.5 * (elements[l].U[1] + elements[l + 1].U[0]);
			elements[l].U[1] = u_int;
			elements[l + 1].U[0] = u_int;
		}*/
		if (k == 499) //t=0.25
		{
			std::cout << "=========   " << total_time << "   ============" << std::endl;
			std::ofstream Q3square025;
			Q3square025.open("Q3square025.txt");
			if (!Q3square025)
			{
				std::cout << "Couldn't open file!!" << std::endl;
				return 1;
			}
			Q3square025 << elements[0].U[0] << std::endl;
			for (int i = 0; i < N; i++)// Write to file
			{
				//Q3square025 << elements[i].U[0] << std::endl;
				Q3square025 << elements[i].U[1] << std::endl;
			}
			Q3square025.close();
		}
		if (k == 1499) //t=0.75
		{
			std::cout << "=========   " << total_time << "   ============" << std::endl;
			std::ofstream Q3square075;
			Q3square075.open("Q3square075.txt");
			if (!Q3square075)
			{
				std::cout << "Couldn't open file!!" << std::endl;
				return 1;
			}
			Q3square075 << elements[0].U[0] << std::endl;
			for (int i = 0; i < N; i++)// Write to file
			{
				//Q3square075 << elements[i].U[0] << std::endl;
				Q3square075 << elements[i].U[1] << std::endl;
			}
			Q3square075.close();
		}
		if (k == 1999) //t=1
		{
			std::cout << "=========   " << total_time << "   ============" << std::endl;
			std::ofstream Q3square1;
			Q3square1.open("Q3square1.txt");
			if (!Q3square1)
			{
				std::cout << "Couldn't open file!!" << std::endl;
				return 1;
			}
			Q3square1 << elements[0].U[0] << std::endl;
			for (int i = 0; i < N; i++)// Write to file
			{
				//Q3square1 << elements[i].U[0] << std::endl;
				Q3square1 << elements[i].U[1] << std::endl;
			}
			Q3square1.close();
		}

		if (k == 2499) //t=1.25
		{
			std::cout << "=========   " << total_time << "   ============" << std::endl;
			std::ofstream Q3square125;
			Q3square125.open("Q3square125.txt");
			if (!Q3square125)
			{
				std::cout << "Couldn't open file!!" << std::endl;
				return 1;
			}
			Q3square125 << elements[0].U[0] << std::endl;
			for (int i = 0; i < N; i++)// Write to file
			{
				//Q3square125 << elements[i].U[0] << std::endl;
				Q3square125 << elements[i].U[1] << std::endl;
			}
			Q3square125.close();
		}

		if (k == 3499) //t=1.75
		{
			std::cout << "=========   " << total_time << "   ============" << std::endl;
			std::ofstream Q3square175;
			Q3square175.open("Q3square175.txt");
			if (!Q3square175)
			{
				std::cout << "Couldn't open file!!" << std::endl;
				return 1;
			}
			Q3square175 << elements[0].U[0] << std::endl;
			for (int i = 0; i < N; i++)// Write to file
			{
				//Q3square1 << elements[i].U[0] << std::endl;
				Q3square175 << elements[i].U[1] << std::endl;
			}
			Q3square175.close();
		}

		if (k == 3999) //t=2
		{
			std::cout << "=========   " << total_time << "   ============" << std::endl;
			std::ofstream Q3square2;
			Q3square2.open("Q3square2.txt");
			if (!Q3square2)
			{
				std::cout << "Couldn't open file!!" << std::endl;
				return 1;
			}
			Q3square2 << elements[0].U[0] << std::endl;
			for (int i = 0; i < N; i++)// Write to file
			{
				//Q3square2 << elements[i].U[0] << std::endl;
				Q3square2 << elements[i].U[1] << std::endl;
			}
			Q3square2.close();
		}

	}
	





	//std::cout << "working" << std::endl;
	return 0;
}



//  ===========================================    Testing Initialisation Q1-3   ================================================
// Testing of class AdvectionElement
	/*
	int N=10;
	std::vector<AdvectionElement> elements(N);
	double x_start = 0.0; const double x_end = 2.0 * 3.14159265358979323846;

	for (int j = 0; j < N; j++)
	{
		// Loop over vector to fill out member data X, U
		elements[j].X[0] = x_start + j * (x_end - x_start) / N;
		elements[j].X[1] = x_start + (j + 1) * (x_end - x_start) / N;
		elements[j].U[0] = 1.5 + std::sin(elements[j].X[0]);
		elements[j].U[1] = 1.5 + std::sin(elements[j].X[1]);
	}
	//std::cout << elements[0].X[0] << std::endl;
	//std::cout << elements[0].X[1] << std::endl;
	for (int i = 0; i < N; i++)
	{
		// Setting neighbour pointers for each element
		if (i == 0)
		{
			// set neighbour pointer connecting first and last element
			elements[i].Left_neighbour_pt = &elements[N - 1];
			elements[i].Right_neighbour_pt = &elements[i + 1];
		}
		else if (i == N - 1)
		{
			// set neighbour pointer connecting first and last element
			elements[i].Right_neighbour_pt = &elements[0];
			elements[i].Left_neighbour_pt = &elements[i - 1];
		}
		else
		{
			elements[i].Left_neighbour_pt = &elements[i - 1];// Left pointer
			elements[i].Right_neighbour_pt = &elements[i + 1];// Right pointer
		}
	}
	//std::cout << (*elements[N-2].Right_neighbour_pt).X[1] << std::endl;
	//std::cout << elements[2].X[1] << std::endl;
	//std::cout << elements[2].U[1] << std::endl;
	*/

	// Program that gives approximation of function for N=10,100,200
	/*
	std::ofstream approximation_u10;
	approximation_u10.open("approximation_u10.txt");
	if (!approximation_u10)
	{
		std::cout << "Couldn't open file!!" << std::endl;
		return 1;
	}
	// Testing of class AdvectionElement
	int N = 10;
	std::vector<AdvectionElement> elements(N);
	double x_start = 0.0; const double x_end = 2.0 * 3.14159265358979323846;
	for (int j = 0; j < N; j++)
	{
		// Loop over vector to fill out member data X, U
		elements[j].X[0] = x_start + j * (x_end - x_start) / N;
		elements[j].X[1] = x_start + (j + 1) * (x_end - x_start) / N;
		elements[j].U[0] = 1.5 + std::sin(elements[j].X[0]);
		elements[j].U[1] = 1.5 + std::sin(elements[j].X[1]);

		approximation_u10 << elements[j].U[0] << "  " << elements[j].X[0];
		approximation_u10 << std::endl;
		approximation_u10 << elements[j].interpolated_u(0.0) << "  " << elements[j].interpolated_x(0.0) ; // Unknown at centre of each element
		approximation_u10 << std::endl;
		if (j == N - 1)
		{
			approximation_u10 << elements[j].U[1] << "  " << elements[j].X[1];
			approximation_u10 << std::endl;
		}
	}
	approximation_u10.close();
	*/

	/*
	std::ofstream approximation_u100;
	approximation_u100.open("approximation_u100.txt");
	if (!approximation_u100)
	{
		std::cout << "Couldn't open file!!" << std::endl;
		return 1;
	}
	// Testing of class AdvectionElement
	int N = 100;
	std::vector<AdvectionElement> elements(N);
	double x_start = 0.0; const double x_end = 2.0 * 3.14159265358979323846;
	for (int j = 0; j < N; j++)
	{
		// Loop over vector to fill out member data X, U
		elements[j].X[0] = x_start + j * (x_end - x_start) / N;
		elements[j].X[1] = x_start + (j + 1) * (x_end - x_start) / N;
		elements[j].U[0] = 1.5 + std::sin(elements[j].X[0]);
		elements[j].U[1] = 1.5 + std::sin(elements[j].X[1]);

		approximation_u100 << elements[j].U[0] << "  " << elements[j].X[0];
		approximation_u100 << std::endl;
		approximation_u100 << elements[j].interpolated_u(0.0) << "  " << elements[j].interpolated_x(0.0); // Unknown at centre of each element
		approximation_u100 << std::endl;
		if (j == N - 1)
		{
			approximation_u100 << elements[j].U[1] << "  " << elements[j].X[1];
			approximation_u100 << std::endl;
		}
	}
	approximation_u100.close();
	*/

	/*
	std::ofstream approximation_u200;
	approximation_u200.open("approximation_u200.txt");
	if (!approximation_u200)
	{
		std::cout << "Couldn't open file!!" << std::endl;
		return 1;
	}
	// Testing of class AdvectionElement
	int N = 200;
	std::vector<AdvectionElement> elements(N);
	double x_start = 0.0; const double x_end = 2.0 * 3.14159265358979323846;
	for (int j = 0; j < N; j++)
	{
		// Loop over vector to fill out member data X, U
		elements[j].X[0] = x_start + j * (x_end - x_start) / N;
		elements[j].X[1] = x_start + (j + 1) * (x_end - x_start) / N;
		elements[j].U[0] = 1.5 + std::sin(elements[j].X[0]);
		elements[j].U[1] = 1.5 + std::sin(elements[j].X[1]);

		approximation_u200 << elements[j].U[0] << "  " << elements[j].X[0];
		approximation_u200 << std::endl;
		approximation_u200 << elements[j].interpolated_u(0.0) << "  " << elements[j].interpolated_x(0.0); // Unknown at centre of each element
		approximation_u200 << std::endl;
		if (j == N - 1)
		{
			approximation_u200 << elements[j].U[1] << "  " << elements[j].X[1];
			approximation_u200 << std::endl;
		}
	}
	approximation_u200.close();
	*/


	//=========================================================================================================================================================================
		//                                                                   Question 1
		// Defining the domain on which we solve the equation on
	/*double x_start = 0.0, x_end = 2 * std::acos(-1.0);
	// ==========================   Solving the advection equation   ==================================
	// Initial sine-wave profile u = 1.5 + std::sin(x)
	// Again begin by initialising the variables
	int N = 250;
	std::vector<AdvectionElement> elements(N);
	for (int j = 0; j < N; j++)
	{
		// Loop over vector to fill out member data X, U
		elements[j].X[0] = x_start + j * (x_end - x_start) / N;
		elements[j].X[1] = x_start + (j + 1) * (x_end - x_start) / N;
		elements[j].U[0] = 1.5 + std::sin(elements[j].X[0]);
		elements[j].U[1] = 1.5 + std::sin(elements[j].X[1]);
	}
	elements[N - 1].U[1] = elements[0].U[0];  // Connecting the first and last elements U^0_0 = U^{N-1}_1
	for (int i = 0; i < N; i++)
	{
		// Setting neighbour pointers for each element
		if (i == 0)
		{
			// set neighbour pointer connecting first and last element
			elements[i].Left_neighbour_pt = &elements[N - 1];
			elements[i].Right_neighbour_pt = &elements[i + 1];
		}
		else if (i == N - 1)
		{
			// set neighbour pointer connecting first and last element
			elements[i].Right_neighbour_pt = &elements[0];
			elements[i].Left_neighbour_pt = &elements[i - 1];
		}
		else
		{
			elements[i].Left_neighbour_pt = &elements[i - 1];// Left pointer
			elements[i].Right_neighbour_pt = &elements[i + 1];// Right pointer
		}
	}
	//========================================================================================================================================================================================================================
	// Initialise the solution for t=0

	std::cout << elements[0].U[0] << std::endl;
	elements[0].timestep(0.025);
	std::cout << "flux(1.5) = " << elements[0].flux(1.5) << std::endl;
	std::cout << elements[0].U_updated[0] << std::endl;
	*/


	// =====================================   Data for t=0   =============================================
	/*
	std::ofstream time0;
	time0.open("time0.txt");
	if (!time0)
	{
		std::cout << "Couldn't open file!!" << std::endl;
		return 1;
	}
	time0 << elements[0].U[0] << std::endl;
	for (int i = 0; i < N; i++)// Write to file
	{
		time0 << elements[i].U[1] << std::endl;
	}
	time0.close();
	// TIMESTEPPING	and setting resolution
	double dt = 0.025, total = 0.0;
	for (int k = 0; k < 40; k++)
	{
		// Loop over all elements
		for (int i = 0; i < N; i++)
		{
			elements[i].timestep(dt);
		}
		// Mechanism for updating all the values once timestep(dt) has been called for every element.
		for (int j = 0; j < N; j++)
		{
			elements[j].U[0] = elements[j].U_updated[0];
			elements[j].U[1] = elements[j].U_updated[1];
		}
		// Enforce constraint that adjacent nodes must be equal "smoothing"
		double u_int = 0.5 * (elements[0].U[0] + elements[N - 1].U[1]);
		elements[0].U[0] = u_int;
		elements[N - 1].U[1] = u_int;  // Connecting the first and last elements U^0_0 = U^{N-1}_1
		for (int l = 0; l < N - 1; l++)
		{
			u_int = 0.5 * (elements[l].U[1] + elements[l + 1].U[0]);
			elements[l].U[1] = u_int;
			elements[l + 1].U[0] = u_int;
		}
		total = total + dt;
		// =====================================   Data for t=0.25   =============================================
		if (k == 9) // time = 0.25
		{
			// Write the solution at this time step to a txt file
			std::ofstream time025;
			time025.open("time025.txt");
			if (!time025)
			{
				std::cout << "Couldn't open file!!" << std::endl;
				return 1;
			}
			time025 << elements[0].U[0] << std::endl;
			for (int i = 0; i < N; i++)// Write to file
			{
				//time025 << elements[i].U[0] << std::endl;
				time025 << elements[i].U[1] << std::endl;
			}
			time025.close();
			std::cout << total << "====================" << std::endl;
			for (int m = 0; m < N; m++)
			{
				display_vector(elements[m].U);
			}
		}
		// =====================================   Data for t=0.5   =============================================
		if (k == 19) // time = 0.5
		{
			// Write the solution at this time step to a txt file
			std::ofstream time05;
			time05.open("time05.txt");
			if (!time05)
			{
				std::cout << "Couldn't open file!!" << std::endl;
				return 1;
			}
			time05 << elements[0].U[0] << std::endl;
			for (int i = 0; i < N; i++)// Write to file
			{
				time05 << elements[i].U[1] << std::endl;
			}
			time05.close();
			std::cout << total << "====================" << std::endl;
			for (int m = 0; m < N; m++)
			{
				display_vector(elements[m].U);
			}
		}
		// =====================================   Data for t=1.0   =============================================
		if (k == 39) // time = 1.0
		{
			// Write the solution at this timestep to a txt file
			std::ofstream time1;
			time1.open("time1.txt");
			if (!time1)
			{
				std::cout << "Couldn't open file!!" << std::endl;
				return 1;
			}
			time1 << elements[0].U[0] << std::endl;
			for (int i = 0; i < N; i++)// Write to file
			{
				time1 << elements[i].U[1] << std::endl;
			}
			time1.close();
			std::cout << total << "====================" << std::endl;
			for (int m = 0; m < N; m++)
			{
				display_vector(elements[m].U);
			}

		}
	}
	*/


/*
	//  ============================================================  Question 2  =============================================================================================
	// Discontinuous initial profile
	int N = 200;
	double x_start = 0.0, x_end = 2 * std::acos(-1.0);
	std::vector<AdvectionElement> elements(N);
	for (int j = 0; j < N; j++)
	{
		// Loop over vector to fill out member data X, U
		// Initialise X values
		elements[j].X[0] = x_start + j * (x_end - x_start) / N;
		elements[j].X[1] = x_start + (j + 1) * (x_end - x_start) / N;
		// Using square-wave profile
		if (elements[j].X[0] >= 0.0 && elements[j].X[0] <= 1.0)
		{
			elements[j].U[0] = 1.0;
			elements[j].U[1] = 1.0;
		}
		else
		{
			elements[j].U[0] = 0.0;
			elements[j].U[1] = 0.0;
		}

	}
	//elements[N - 1].U[1] = elements[0].U[0];  // Connecting the first and last elements U^0_0 = U^{N-1}_1
	for (int i = 0; i < N; i++)
	{
		// Setting neighbour pointers for each element
		if (i == 0)
		{
			// set neighbour pointer connecting first and last element
			elements[i].Left_neighbour_pt = &elements[N - 1];
			elements[i].Right_neighbour_pt = &elements[i + 1];
		}
		else if (i == N - 1)
		{
			// set neighbour pointer connecting first and last element
			elements[i].Right_neighbour_pt = &elements[0];
			elements[i].Left_neighbour_pt = &elements[i - 1];
		}
		else
		{
			elements[i].Left_neighbour_pt = &elements[i - 1];// Left pointer
			elements[i].Right_neighbour_pt = &elements[i + 1];// Right pointer
		}
	}

	//std::cout << "Numerical flux = " << elements[0].h(elements[0].U[1], elements[0].Right_neighbour_pt->U[0]) << std::endl;

	//std::cout << "Flux(U[1]) = " << elements[0].flux(elements[0].U[1]) << std::endl;
	/*
	for (int i = 0; i < N; i++)
	{
		display_vector(elements[i].U);
	}

	// Perform timestepping
	//
	//double dt = 0.025;
	//std::cout << elements[0].U[0] << std::endl;
	//elements[0].timestep(dt);
	//std::cout << elements[0].U_updated[0] << std::endl;
	//std::cout << elements[15].U[0] << std::endl;

	double dt = 0.001, total_time = 0.0;
	for (int k = 0; k < 1000; k++)
	{
		for (int i = 0; i < N; i++)
		{
			elements[i].timestep(dt);
		}
		for (int j = 0; j < N; j++)
		{
			elements[j].U[0] = elements[j].U_updated[0];
			elements[j].U[1] = elements[j].U_updated[1];
		}
		total_time = total_time + dt;
	}*/


	// =====================================  TIMESTEPPING  ==================================================
/*
double dt = 0.001, total = 0.0;
for (int k = 0; k < 1000; k++)
{
	if (k == 0)
	{
		std::cout << "============  " << total << "   ============" << std::endl;
		for (int i = 0; i < N; i++)
		{
			display_vector(elements[i].U);
		}
		std::ofstream Q2data;
		Q2data.open("Q2data.txt");
		if (!Q2data)
		{
			std::cout << "Couldn't open file!!" << std::endl;
			return 1;
		}
		//Q2data << elements[0].U[0] << std::endl;
		for (int l = 0; l < N; l++)
		{
			Q2data << elements[l].U[0] << std::endl;
			Q2data << elements[l].U[1] << std::endl;
		}
		Q2data.close();
	}
	for (int i = 0; i < N; i++)
	{
		elements[i].timestep(dt);
	}
	for (int j = 0; j < N; j++)
	{
		elements[j].U[0] = elements[j].U_updated[0];
		elements[j].U[1] = elements[j].U_updated[1];
	}
	//elements[0].U[0] = elements[N - 1].U[1];  // Connecting the first and last elements U^0_0 = U^{N-1}_1
	/*
	double u_int = 0.5 * (elements[0].U[0] + elements[N - 1].U[1]);
	elements[0].U[0] = u_int;
	elements[N - 1].U[1] = u_int;  // Connecting the first and last elements U^0_0 = U^{N-1}_1
	for (int l = 0; l < N - 1; l++)
	{
		u_int = 0.5 * (elements[l].U[1] + elements[l + 1].U[0]);
		elements[l].U[1] = u_int;
		elements[l + 1].U[0] = u_int;
	}*/
/*
	total = total + dt;
	if (k == 249)
	{
		std::cout << "============  " << total << "   ============" << std::endl;
		for (int i = 0; i < N; i++)
		{
			display_vector(elements[i].U);
		}
		std::ofstream Q2data0;
		Q2data0.open("Q2data0.txt");
		if (!Q2data0)
		{
			std::cout << "Couldn't open file!!" << std::endl;
			return 1;
		}
		//Q2data0 << elements[0].U[0] << std::endl;
		for (int l = 0; l < N; l++)
		{
			Q2data0 << elements[l].U[0] << std::endl;
			Q2data0 << elements[l].U[1] << std::endl;
		}

		Q2data0.close();
	}

	if (k == 499)
	{
		std::cout << "============  " << total << "   ============" << std::endl;
		for (int i = 0; i < N; i++)
		{
			display_vector(elements[i].U);
		}

		std::ofstream Q2data1;
		Q2data1.open("Q2data1.txt");
		if (!Q2data1)
		{
			std::cout << "Couldn't open file!!" << std::endl;
			return 1;
		}
		//Q2data1 << elements[0].U[0] << std::endl;
		for (int l = 0; l < N; l++)
		{
			Q2data1 << elements[l].U[0] << std::endl;
			Q2data1 << elements[l].U[1] << std::endl;
		}
		Q2data1.close();
	}

	if (k == 999)
	{
		std::cout << "============  " << total << "   ============" << std::endl;
		for (int i = 0; i < N; i++)
		{
			display_vector(elements[i].U);
		}
		std::ofstream Q2data2;
		Q2data2.open("Q2data2.txt");
		if (!Q2data2)
		{
			std::cout << "Couldn't open file!!" << std::endl;
			return 1;
		}
		//Q2data2 << elements[0].U[0] << std::endl;
		for (int l = 0; l < N; l++)
		{
			Q2data2 << elements[l].U[0] << std::endl;
			Q2data2 << elements[l].U[1] << std::endl;
		}
		Q2data2.close();
	}
}
*/





// ==================================================================================================================================================


	//  ============================================================  Question 3  =============================================================================================
	// Initialise the sine-wave initial velocity profile in the same we did in Q1
/*
	double x_start = 0.0, x_end = 2 * std::acos(-1.0);
	int N = 200;
	std::vector<BurgersElement> elements(N);  // Note the BurgersElement
	for (int j = 0; j < N; j++)
	{
		// Loop over vector to fill out member data X, U
		elements[j].X[0] = x_start + j * (x_end - x_start) / N;
		elements[j].X[1] = x_start + (j + 1) * (x_end - x_start) / N;
		elements[j].U[0] = 1.5 + std::sin(elements[j].X[0]);
		elements[j].U[1] = 1.5 + std::sin(elements[j].X[1]);
	}
	//elements[N - 1].U[1] = elements[0].U[0];  // Connecting the first and last elements U^0_0 = U^{N-1}_1
	for (int i = 0; i < N; i++)
	{
		// Setting neighbour pointers for each element
		if (i == 0)
		{
			// set neighbour pointer connecting first and last element
			elements[i].Left_neighbour_pt = &elements[N - 1];
			elements[i].Right_neighbour_pt = &elements[i + 1];
		}
		else if (i == N - 1)
		{
			// set neighbour pointer connecting first and last element
			elements[i].Right_neighbour_pt = &elements[0];
			elements[i].Left_neighbour_pt = &elements[i - 1];
		}
		else
		{
			elements[i].Left_neighbour_pt = &elements[i - 1];// Left pointer
			elements[i].Right_neighbour_pt = &elements[i + 1];// Right pointer
		}
	}

	// Testing the derived class
	/*
	std::cout << elements[0].U[0] << std::endl;
	std::cout << elements[0].U[1] << std::endl;
	std::cout << elements[0].flux(-1.0/sqrt(3.0)) << "+" << elements[0].flux(1.0/sqrt(3.0)) << " = "  << elements[0].integrate_flux() << std::endl;
	std::cout << 0.5*elements[0].interpolated_u(-1.0/sqrt(3.0))* elements[0].interpolated_u(-1.0 / sqrt(3.0)) << std::endl;
	*/
	//std::cout << elements[0].integrate_flux() << std::endl;
/*
double dt = 0.0005, total = 0.0;
for (int i = 0; i < N; i++)
{
	display_vector(elements[i].U);
}

for (int k = 0; k < 4000; k++)
{
	// Initial profile
	if (k == 0)
	{
		std::ofstream Q3sine_0;
		Q3sine_0.open("Q3sine_0.txt");
		if (!Q3sine_0)
		{
			std::cout << "Couldn't open file!" << std::endl; return 1;
		}
		Q3sine_0 << elements[0].U[0] << std::endl;
		for (int j = 0; j < N; j++)
		{
			//Q3sine_0 << elements[j].U[0] << std::endl;
			Q3sine_0 << elements[j].U[1] << std::endl;
		}
		Q3sine_0.close();
	}
	for (int i = 0; i < N; i++)
	{
		elements[i].timestep(dt);
	}
	for (int j = 0; j < N; j++)
	{
		elements[j].U[0] = elements[j].U_updated[0];
		elements[j].U[1] = elements[j].U_updated[1];
	}
	total = total + dt;
	if (k == 999) // t = 0.5
	{
		std::cout << "============  " << total << "  ============" << std::endl;
		for (int i = 0; i < N; i++)
		{
			display_vector(elements[i].U);
		}
		std::ofstream Q3sine_05;
		Q3sine_05.open("Q3sine_05.txt");
		if (!Q3sine_05)
		{
			std::cout << "Couldn't open file!" << std::endl; return 1;
		}
		Q3sine_05 << elements[0].U[0] << std::endl;
		for (int j = 0; j < N; j++)
		{
			//Q3sine_05 << elements[j].U[0] << std::endl;
			Q3sine_05 << elements[j].U[1] << std::endl;
		}
		Q3sine_05.close();
	}
	if (k == 1799) //t =0.9
	{
		std::cout << "============  " << total << "  ============" << std::endl;
		for (int i = 0; i < N; i++)
		{
			display_vector(elements[i].U);
		}
		std::ofstream Q3sine_09;
		Q3sine_09.open("Q3sine_09.txt");
		if (!Q3sine_09)
		{
			std::cout << "Couldn't open file!" << std::endl; return 1;
		}
		Q3sine_09 << elements[0].U[0] << std::endl;
		for (int j = 0; j < N; j++)
		{
			//Q3sine_09 << elements[j].U[0] << std::endl;
			Q3sine_09 << elements[j].U[1] << std::endl;
		}
		Q3sine_09.close();
	}

	if (k == 1999) // t = 1
	{
		std::cout << "============  " << total << "  ============" << std::endl;
		for (int i = 0; i < N; i++)
		{
			display_vector(elements[i].U);
		}
		std::ofstream Q3sine_1;
		Q3sine_1.open("Q3sine_1.txt");
		if (!Q3sine_1)
		{
			std::cout << "Couldn't open file!" << std::endl; return 1;
		}
		Q3sine_1 << elements[0].U[0] << std::endl;
		for (int j = 0; j < N; j++)
		{
			//Q3sine_1 << elements[j].U[0] << std::endl;
			Q3sine_1 << elements[j].U[1] << std::endl;
		}
		Q3sine_1.close();
	}

	if (k == 2199) // t=1.1
	{
		std::cout << "============  " << total << "  ============" << std::endl;
		for (int i = 0; i < N; i++)
		{
			display_vector(elements[i].U);
		}
		std::ofstream Q3sine_11;
		Q3sine_11.open("Q3sine_11.txt");
		if (!Q3sine_11)
		{
			std::cout << "Couldn't open file!" << std::endl; return 1;
		}
		Q3sine_11 << elements[0].U[0] << std::endl;
		for (int j = 0; j < N; j++)
		{
			//Q3sine_11 << elements[j].U[0] << std::endl;
			Q3sine_11 << elements[j].U[1] << std::endl;
		}
		Q3sine_11.close();
	}

	if (k == 2499) // t=1.25
	{
		std::cout << "============  " << total << "  ============" << std::endl;
		for (int i = 0; i < N; i++)
		{
			display_vector(elements[i].U);
		}
		std::ofstream Q3sine_125;
		Q3sine_125.open("Q3sine_125.txt");
		if (!Q3sine_125)
		{
			std::cout << "Couldn't open file!" << std::endl; return 1;
		}
		Q3sine_125 << elements[0].U[0] << std::endl;
		for (int j = 0; j < N; j++)
		{
			//Q3sine_125 << elements[j].U[0] << std::endl;
			Q3sine_125 << elements[j].U[1] << std::endl;
		}
		Q3sine_125.close();
	}

	if (k == 2999) // t=1.5
	{
		std::cout << "============  " << total << "  ============" << std::endl;
		for (int i = 0; i < N; i++)
		{
			display_vector(elements[i].U);
		}
		std::ofstream Q3sine_15;
		Q3sine_15.open("Q3sine_15.txt");
		if (!Q3sine_15)
		{
			std::cout << "Couldn't open file!" << std::endl; return 1;
		}
		Q3sine_15 << elements[0].U[0] << std::endl;
		for (int j = 0; j < N; j++)
		{
			//Q3sine_15 << elements[j].U[0] << std::endl;
			Q3sine_15 << elements[j].U[1] << std::endl;
		}
		Q3sine_15.close();
	}

	if (k == 3499) // t=1.75
	{
		std::cout << "============  " << total << "  ============" << std::endl;
		for (int i = 0; i < N; i++)
		{
			display_vector(elements[i].U);
		}
		std::ofstream Q3sine_175;
		Q3sine_175.open("Q3sine_175.txt");
		if (!Q3sine_175)
		{
			std::cout << "Couldn't open file!" << std::endl; return 1;
		}
		Q3sine_175 << elements[0].U[0] << std::endl;
		for (int j = 0; j < N; j++)
		{
			//Q3sine_175 << elements[j].U[0] << std::endl;
			Q3sine_175 << elements[j].U[1] << std::endl;
		}
		Q3sine_175.close();
	}

	if (k == 3999)
	{
		std::cout << "============  " << total << "  ============" << std::endl;
		for (int i = 0; i < N; i++)
		{
			display_vector(elements[i].U);
		}
		std::ofstream Q3sine_2;
		Q3sine_2.open("Q3sine_2.txt");
		if (!Q3sine_2)
		{
			std::cout << "Couldn't open file!" << std::endl; return 1;
		}
		Q3sine_2 << elements[0].U[0] << std::endl;
		for (int j = 0; j < N; j++)
		{
			//Q3sine_2 << elements[j].U[0] << std::endl;
			Q3sine_2 << elements[j].U[1] << std::endl;
		}
		Q3sine_2.close();
	}
}

}*/

/*
	// ========================================================================================================================================================
	// =======================================       Question 3  -  Discontinuous profile       =============================================
	// Discontinuous initial profile
int N = 100;
double x_start = 0.0, x_end = 2 * std::acos(-1.0);
std::vector<BurgersElement> elements(N);
for (int j = 0; j < N; j++)
{
	// Loop over vector to fill out member data X, U
	// Initialise X values
	elements[j].X[0] = x_start + j * (x_end - x_start) / N;
	elements[j].X[1] = x_start + (j + 1) * (x_end - x_start) / N;
	// Using square-wave profile
	if (elements[j].X[0] >= 0.0 && elements[j].X[0] <= 1.0)
	{
		elements[j].U[0] = 1.0;
		elements[j].U[1] = 1.0;
	}
	else
	{
		//std::cout << "Discontinuity occurs at " << j << std::endl;
		elements[j].U[0] = 0.0;
		elements[j].U[1] = 0.0;
	}

}
//elements[N - 1].U[1] = elements[0].U[0];  // Connecting the first and last elements U^0_0 = U^{N-1}_1
for (int i = 0; i < N; i++)
{
	// Setting neighbour pointers for each element
	if (i == 0)
	{
		// set neighbour pointer connecting first and last element
		elements[i].Left_neighbour_pt = &elements[N - 1];
		elements[i].Right_neighbour_pt = &elements[i + 1];
	}
	else if (i == N - 1)
	{
		// set neighbour pointer connecting first and last element
		elements[i].Right_neighbour_pt = &elements[0];
		elements[i].Left_neighbour_pt = &elements[i - 1];
	}
	else
	{
		elements[i].Left_neighbour_pt = &elements[i - 1];// Left pointer
		elements[i].Right_neighbour_pt = &elements[i + 1];// Right pointer
	}
}
std::cout << elements[15].U[1] << std::endl; // =1
std::cout << elements[16].U[0] << std::endl; // =0

std::cout << elements[16].Left_neighbour_pt->U[1] << std::endl;
std::cout << elements[16].U[0] << std::endl;
//std::cout << elements[0].flux((elements[16].U[0] + elements[15].Left_neighbour_pt->U[1])/2) << std::endl;
//std::cout << elements[0].h(elements[15].Left_neighbour_pt->U[1], elements[16].U[0]) << std::endl;
//elements[16].timestep(0.001);
//std::cout << elements[16].U_updated[0] << std::endl;


double dt = 0.001;
for (int k = 0; k < 100; k++)
{
	for (int i = 0; i < N; i++)
	{
		elements[i].timestep(dt);
	}
	for (int j = 0; j < N; j++)
	{
		elements[j].U[0] = elements[j].U_updated[0];
		elements[j].U[1] = elements[j].U_updated[1];
	}

	/*
	double u_int = 0.0;
	for (int l = 0; l < N - 1; l++)
	{
		u_int = 0.5 * (elements[l].U[1] + elements[l + 1].U[0]);
		elements[l].U[1] = u_int;
		elements[l + 1].U[0] = u_int;
	}*/
/*
}
std::cout << "  =====================   " << std::endl;
for (int k = 0; k < N; k++)
{
	display_vector(elements[k].U);
}

*/










// CRAP
/*
// Compute updated values for this particular timestep
for (int i = 0; i < U.size(); i++)
{
	double diff = (elements[i].X[1] - elements[i].X[0]) / 18.0;
	U_tmp.push_back(elements[i].U[0] + dt * diff * (-1.5 * integrate_flux() + 2 * h(elements[i - 1].U[1], elements[i].U[0]) + h(elements[i].U[1], elements[i + 1].U[0])));
	U_tmp.push_back(elements[i].U[1] + dt * diff * (1.5 * integrate_flux() - 2 * h(elements[i].U[1], elements[i + 1].U[0]) - h(elements[i - 1].U[1], elements[i].U[0])));
}
*/

/* 08/01/21 18:53
	std::cout << elements[0].U[0] << std::endl;
	for (int i = 0; i < N; i++)
	{
		elements[i].timestep(1.0);
	}
	for (int j = 0; j < N; j++)
	{
		elements[j].U[0] = elements[j].U_updated[0];
		elements[j].U[1] = elements[j].U_updated[1];
	}
	elements[N - 1].U[1] = elements[0].U[0];  // Connecting the first and last elements U^0_0 = U^{N-1}_1
	std::cout << elements[0].U[0] << std::endl;
	double dt = 0.1;
	for (int i = 0; i < N; i++)
	{
		elements[i].timestep(dt);
	}
	for (int j = 0; j < N; j++)
	{
		elements[j].U[0] = elements[j].U_updated[0];
		elements[j].U[1] = elements[j].U_updated[1];
	}
	std::cout << elements[0].U[0] << std::endl;
*/

/*
	//double diff = (elements[0].X[1] - elements[0].X[0]) / 18.0;
	//std::cout << (std::acos(-1)/3.0)*1.0* (-1.5 * elements[0].integrate_flux() + 2 * elements[0].flux(elements[0].U[0]) + elements[0].flux(elements[0].U[1])) << std::endl;





	std::cout << "U_0= " << elements[0].U[0] << "   U_1= " << elements[0].U[1] << std::endl;
	for (int k = 0; k < 3; k++)
	{
		double dt = k * 1.0 / 10.0;
		for (int i = 0; i < N; i++)
		{
			//std::cout << dt << std::endl;

			elements[i].timestep(dt);
		}
		for (int j = 0; j < N; j++)
		{
			elements[j].U[0] = elements[j].U_updated[0];
			elements[j].U[1] = elements[j].U_updated[1];
		}
		elements[N - 1].U[1] = elements[0].U[0];  // Connecting the first and last elements U^0_0 = U^{N-1}_1

		//std::cout << "dt= " << dt << std::endl;
		//std::cout << "U_0= " << elements[0].U[0] << "   U_1= " << elements[0].U[1] << std::endl;
		std::cout << "dt= " << dt << std::endl;
		std::cout << "U_0^1 = " << elements[0].U[1] << std::endl;
		std::cout << "U_1^0 = " << elements[1].U[0] << std::endl;


		std::cout << elements[1].interpolated_u(0.0) << std::endl;
	}
*/


// 09/01/21 9:59
// Define the updated values for the unknowns
/*
double U_0_updated, U_1_updated;
double diff = (X[1] - X[0]) / 18.0;
if (U[0] == (*Left_neighbour_pt).U[1])
{
	// The flux is continuous
	U_0_updated = U[0] + dt * diff * (-1.5 * integrate_flux() + 2 * flux(U[0]) + flux(U[1]));
	U_1_updated = U[1] + dt * diff * (1.5 * integrate_flux() - 2 * flux(U[1]) - flux(U[0]));
	//std::cout << "timestep = " <<  dt * diff * (-1.5 * integrate_flux() + 2 * flux(U[0]) + flux(U[1])) << std::endl;
}
else
{
	// Flux is discontinuous
	U_0_updated = U[0] + dt * diff * (-1.5 * integrate_flux() + 2 * h((*Left_neighbour_pt).U[1], U[0]) + h(U[1], (*Right_neighbour_pt).U[0]));
	U_1_updated = U[1] + dt * diff * (1.5 * integrate_flux() - 2 * h(U[1], (*Right_neighbour_pt).U[0]) - h((*Left_neighbour_pt).U[1], U[0]));
}
*/




/* THINK THIS WORKS!
int i = 1;
MVector U_timestep(2), Flux(2), F_e(2), U,update;
MMatrix M_e(2, 2, 0.0),M_e_inv;
double coef = (elements[i].X[1] - elements[i].X[0]) / 6.0,dt=0.5;
M_e(0, 0) = coef * 2.0;
M_e(0, 1) = M_e(1, 0) = coef * 1.0;
M_e(1, 1) = coef * 2.0;
M_e_inv = M_e.invert(M_e);
F_e[0] = -0.5 * elements[i].integrate_flux();
F_e[1] = 0.5 * elements[i].integrate_flux();

Flux[0] = elements[i].flux(elements[i].U[0]);
Flux[1] = -elements[i].flux(elements[i].U[1]);
U = { elements[i].U[0],elements[i].U[1] };
// Vector equation which calculates updates

display_MVector(F_e+Flux);
display_MVector(Flux);
//std::cout << M_e_inv(0,1) << std::endl;

U_timestep = U + dt * (M_e.invert(M_e) * (F_e + Flux));
MVector w = dt * (M_e.invert(M_e) * (F_e + Flux));
display_MVector(w);
display_MVector(U_timestep);
*/

