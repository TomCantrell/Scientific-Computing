// 3/12/2020
// C++ file for Project 1
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

#ifndef MVECTOR_H // the 'include guard'
#define MVECTOR_H // see C++ Primer Sec. 2.9.2
#include <vector>
#include <cassert>

// Class that represents a mathematical vector
class MVector
{
public:
	// constructors
	MVector() {}
	explicit MVector(int n) : v(n) {}
	MVector(int n, double x) : v(n, x) {}
	MVector(std::initializer_list<double> l) : v(l) {}
	// access element (lvalue) // not decalred as const allowing for modification(see example sheet 5, q5.6)
	double& operator[](int index) 
	{
		// Assert Error
		if (index <= v.size()-1 && index >= 0) { return v[index]; }
		else { assert(index>=v.size(),"not accessible"); }// std::cout << "not accessible" << std::endl;  exit(1);
	}
	// access element (rvalue) (see example sheet 5, q5.7)
	double operator[](int index) const {
		//Assert Errors
		if (index <= v.size() - 1 && index >= 0) { return v[index]; }
		else { assert("not accessible"); }// std::cout << "not accessible" << std::endl;  exit(1);
	}
	int size() const { return v.size(); } // number of elements

	void resize(int x){ v.resize(x); }
	void push_back(double x) { v.push_back(x); }

	// Member function to swap elements in MVector
	void swap(int i, int j)
	{
		double tmp = v[i];
		v[i] = v[j];
		v[j] = tmp;
	}
	// Creating inital data
	void initialise_random(double xmin, double xmax)
	{
		size_t s = v.size();
		for (size_t i = 0; i < s; i++)
		{
			v[i] = xmin + (xmax - xmin) * rand() / static_cast<double>(RAND_MAX);
		}
	}
private:
	std::vector<double> v;
};
#endif


// Display vector elements
void display_vector(MVector v)
{
	for (int j = 0; j < v.size(); j++)
	{
		std::cout << v[j] << " ";
	}
	std::cout << std::endl;
}
void Display_vector(std::vector<double> v)
{
	for (int j = 0; j < v.size(); j++)
	{
		std::cout << v[j] << " ";
	}
	std::cout << std::endl;
}

// ==========================================    Sorting Algorithms    =================================================

// ==============================================  Bubblesort  ==============================================================
// function which takes a MVector and sorts it in ascending order using bubblesort algorithm
void bubble(MVector& v)
{
	int stop = 1;
	while (stop < v.size()) // loop over the vector until stopping condition is reached
	{
		stop = 1;
		for (int i = 0; i < (v.size() - 1); i++) 
		{
			if (v[i] > v[i + 1]) // Comparing adjacent elements
			{
				v.swap(i, i + 1); 
				stop--;
			}
			else { stop++; }
		}
	}
}

// ===============================================  Quicksort  ================================================================
// function to perform quicksort on section of MVector between indices
void quick_recursive(MVector& v, int start, int end)
{
	int x = start + (end - start) * static_cast<double>(rand()) / static_cast<double>(RAND_MAX); // random pivot
	MVector s1, s2, s3;
	double pivotvalue = v[x];
	for (int j = 0; j <= end; j++) // create vectors s1, s2 and s3
	{
		if (v[j] < v[x]) { s1.push_back(v[j]); }
		else if (v[j] == v[x]) { s2.push_back(v[j]); }
		else { s3.push_back(v[j]);}
	}
	// Call fn recursively
	if (s3.size() > 1)
	{
		quick_recursive(s3, 0, s3.size()-1);
	}
	if (s1.size() > 1)
	{
		quick_recursive(s1, 0, s1.size()-1);
	}
	for (int k = start; k < start + s1.size(); k++) // add the sorted sub-vectors back to main vector
	{
		v[k] = s1[k - start];
	}
	for (int k = start + s1.size(); k < start + s1.size() + s2.size(); k++)
	{
		v[k] = s2[k - s1.size() - start];
	}
	for (int k = start + s1.size() + s2.size(); k < start + s1.size() + s2.size() + s3.size(); k++)
	{
		v[k] = s3[k - s1.size() - s2.size() - start];
	}
}
// wrapper to recursive function 
void quick(MVector& v) { quick_recursive(v, 0, v.size()-1); }

// ==========================================================  Heapsort  ==========================================================
/* heap_from_root makes a heap with i-th vertex as the root, assuming first n values stored in the vector make up the heap data */
void heap_from_root(MVector& v, int i, int n)
{
	// just for 3 nodes, parent and two children make a heap
	int l_child = 2 * i + 1, r_child = 2 * i + 2;
	if (l_child > n) { return; }
	int count1 = 0, count2 = 0;		
	if (2 * i + 2 == n)
		{
		if (v[(2 * i) + 1] > v[i])			
		{
			v.swap(i, 2 * i + 1);
		}
	}
	else if (2 * i + 1 == n)
	{
		if (v[2 * i + 1] > v[i]){v.swap(i, 2 * i + 1);}
	}
	else {
		// Check which leaf is the largest                     
		double max_leaf = std::max(v[2 * i + 2], v[2 * i + 1]);
		if (max_leaf > v[i]) // Only execute if child value > parent value
		{
			if (max_leaf == v[2 * i + 2])
			{
				v.swap(i, 2 * i + 2);
				count2++;
			}
			else {
				v.swap(i, 2 * i + 1);
				count1++;
			}
		}
	}
	// recursively, if our swap counter is positive then call the fn again
	if (count1 > 0)
	{
		if (2 * (2 * i + 1) + 1 < n) // if the swapped node has children
		{
			heap_from_root(v, 2 * i + 1, n);
		}
	}
	else if (count2 > 0)
	{
		if (2 * (2 * i + 2) + 1 < n) // if the swapped node has children
		{
			heap_from_root(v, 2 * i + 2, n);
		}
	}
	else {return; }
}

/* Implement the heapsort algorithm using heap_from_root to build the 
heap and then to perform sorting */
void heap(MVector& v)
{
	int j = v.size();
	// first non-leaf element
	while(2 * j > v.size() - 2) // This loop locates where the first parent vertices in the vector
	{
		j--;
	}
	// loop over the vector using heap_from_root to create a heap 
	for (int i = j; i >= 0; i--)
	{
		heap_from_root(v, i, v.size());
	}
	v.swap(0, v.size() - 1);
	// implementing algorithm
	for (int k = 1; k < v.size(); k++)
	{	
		heap_from_root(v, 0, v.size() - (k));
		v.swap(0, v.size() - (k+1));
		if (v.size() - (k + 1) == 1){return;}
	}
}
// ======================================================================================================================

int main()
{
	std::srand(std::time(NULL)); // // Seeds generator, avoiding same seq of numbers

	// Task 1
	// Simple program that creates an MVector of size 10
	/*MVector v(10, 1.0);
	// Changing MVector v
	for (int i = 1; i < v.size(); i++)
	{
		v[i] += v[i - 1];
	}

	// Output data to file
	std::ofstream Output_1;
	Output_1.open("Output_1.txt");
	if (!Output_1)
	{
		std::cout << "Couldn't open file!!" << std::endl;
		return 1;
	}
	for (int j = 0; j < v.size(); j++)
	{
		Output_1.width(5);  Output_1 << v[j];
	}
	Output_1 << std::endl;
	Output_1.close(); //Close file
	*/

	// Testing assertions in [] operator member functions
	/*
	MVector v(10);
	std::cout << v[11] << std::endl;
	*/




	// A program that computes the average sort time for a randomly initialised vector of length n
	// ======================================================================================================================================
	//Program that computes av. run time of the bubble sort algorithm, outputting this data to a file
	/*
	auto strt = std::chrono::steady_clock::now();
	
	std::ofstream Bubble_sort_time;
	Bubble_sort_time.open("Bubble_sort_time.txt");
	if (!Bubble_sort_time)
	{
		std::cout << "Couldn't open file!" << std::endl;
		return 1;
	}
	for (int i=1; i< 14; i++)
	{
		MVector average_bubble;
		// average sort time of 10 vectors
		for (int j = 0; j < 10; j++)
		{
			MVector v(pow(2, i));
			v.initialise_random(-100000.0, 100000.0);
			auto start = std::chrono::steady_clock::now();
			bubble(v);
			auto end = std::chrono::steady_clock::now();
			auto diff = end - start;
			average_bubble.push_back(std::chrono::duration<double, std::milli>(diff).count());
		}
		double sum=0.0;
		for (int k = 0; k < average_bubble.size(); k++)
		{
			sum = sum + average_bubble[k];
		}
		double average_sort_time_bubble = sum / average_bubble.size();
		//std::cout << "time taken = " << std::chrono::duration<double, std::milli>(diff).count() << "ms" << std::endl;// Using milli = sec * 10^-3
		// Write this to txt file
		Bubble_sort_time.width(10); Bubble_sort_time << pow(2,i);
		Bubble_sort_time.width(20); Bubble_sort_time << average_sort_time_bubble;
		Bubble_sort_time << std::endl;
	}
	Bubble_sort_time.close();
	
	std::cout << "quicksort" << std::endl;
	// Program that computes average run time of quicksort algorithm
	std::ofstream Quick_sort_time;
	Quick_sort_time.open("Quick_sort_time.txt");
	if (!Quick_sort_time)
	{
		std::cout << "Couldn't open file!" << std::endl;
		return 1;
	}
	
	for (int i = 1; i < 14; i++)
	{
		MVector average_quick;
		for (int j = 0; j < 100; j++)
		{
			//std::cout << "n = " << n + i * 10 << std::endl;
			MVector v(pow(2, i));
			v.initialise_random(-100000.0, 100000.0);
			auto start = std::chrono::steady_clock::now();
			quick(v);
			auto end = std::chrono::steady_clock::now();
			auto diff = end - start;
			average_quick.push_back(std::chrono::duration<double, std::milli>(diff).count());
		}
		double sum = 0.0;
		for (int k = 0; k < average_quick.size(); k++)
		{	
			//std::cout << average_quick[k] << std::endl;
			sum = sum + average_quick[k];
		}
		
		double average_sort_time_quick = sum / average_quick.size();

		//std::cout << "time taken = " << std::chrono::duration<double, std::milli>(diff).count() << "ms" << std::endl;// Using milli = sec * 10^-3
		// Write this to txt file
		Quick_sort_time.width(10); Quick_sort_time << pow(2,i);
		Quick_sort_time.width(20); Quick_sort_time << average_sort_time_quick;
		Quick_sort_time << std::endl;
	}
	Quick_sort_time.close();
	
	std::cout << "heapsort" << std::endl;
	// Program that computes average run time of heapsort algorithm
	std::ofstream Heap_sort_time;
	Heap_sort_time.open("Heap_sort_time.txt");
	if (!Heap_sort_time)
	{
		std::cout << "Couldn't open file!" << std::endl;
		return 1;
	}
	for (int i = 1; i < 14; i++)
	{
		//std::cout << "n = " << n + i * 10 << std::endl;
		MVector average_heap;
		for (int j = 0; j < 100; j++)
		{
			MVector v(pow(2, i));
			v.initialise_random(-100000.0, 100000.0);
			auto start = std::chrono::steady_clock::now();
			heap(v);
			auto end = std::chrono::steady_clock::now();
			auto diff = end - start;
			average_heap.push_back(std::chrono::duration<double, std::milli>(diff).count());
		}
		double sum = 0.0;
		std::cout << "average_heap.size()= "<< average_heap.size() << std::endl;
		for (int k = 0; k < average_heap.size(); k++)
		{
			sum = sum + average_heap[k];
		}
		double average_sort_time_heap = sum / average_heap.size();

		//std::cout << "time taken = " << std::chrono::duration<double, std::milli>(diff).count() << "ms" << std::endl;// Using milli = sec * 10^-3
		// Write this to txt file
		Heap_sort_time.width(10); Heap_sort_time << pow(2,i);
		Heap_sort_time.width(20); Heap_sort_time << average_sort_time_heap;
		Heap_sort_time << std::endl;
	}
	Heap_sort_time.close();

	auto ennd = std::chrono::steady_clock::now();
	auto final_diff = ennd - strt;

	std::cout << "done in " << std::chrono::duration<double, std::milli>(final_diff).count() << "ms" << std::endl;
	*/

	/*
	MVector w = { 12.0,4.0,30.8,50.9,0.01,0.02,10000.0,0.02 };
	display_vector(w);
	quick(w);
	display_vector(w);
	*/
	// =====================================================================================================================================
	//display_vector(v);




	// ==========   Testing sorting algorithms  =========== 

	// Testing Bubblesort
	/*MVector v(3);
	v[0] = 5.5; v[1] = 2.0; v[2] = 1.0;
	display_vector(v);
	//v.swap(1, 2);
	bubble(v);
	display_vector(v);
	int pivot = (v.size() + 1) / 2;
	std::cout << pivot << std::endl;
	*/

	// Testing Quicksort
	// Testing initialise_random(,)
	/*MVector w(20);
	w.initialise_random(0.0, 100.0);
	//MVector v = { 4.0 ,0.01 ,12.0 ,18.0 ,15.0 ,14.0 ,21.0 ,19.0 ,30.0 ,17.0, 31.1 ,40.0,40.0,0.002,1000.0,25.0 };
	display_vector(w);
	quick(w);
	display_vector(w);
	

	// Testing Heapsort
	//MVector v = { 2.0,8.0,5.0,3.0,9.0,1.0,4.0,4.0 };
	//MVector v = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 0.002, 25000 };
	//MVector v = { 3, 10, 14, 8, 9, 12, 13, 7, 1, 0, 4, 11, 5, 2, 6, 15, 0.5, 60, -1.0, -1.1 };
	//MVector v(16834);
	//v.initialise_random(-10000.0, 10000.0);
	//std::cout << "Start:" << std::endl;
	//display_vector(v);
	auto start1 = std::chrono::steady_clock::now();
	//heap_from_root(v,0,v.size());
	bubble(v);
	auto end1 = std::chrono::steady_clock::now();
	auto time1_diff = end1 - start1;
	//std::cout << "End:" << std::endl;
	//display_vector(v);

	std::cout << std::endl;  // std::chrono::duration<double, std::milli>(time_diff).count() << "ms"
	std::cout << "time taken: " << std::chrono::duration<double, std::milli>(time1_diff).count() << "ms" << std::endl;
	*/

	std::ofstream test_random;
	test_random.open("test_random.txt");
	if (!test_random)
	{
		std::cout << "Couldn't open file!" << std::endl;
		return 1;
	}
	MVector v(10000);
	v.initialise_random(-1000000.0, 1000000.0);
	for (int j = 0; j < v.size(); j++)
	{
		test_random << v[j];
		test_random << std::endl;
	}	
	test_random.close();













	/*
	int n=v.size(), j = n; // set to j, don't want n to change value
	while (2 * j > n - 2) // This loop locates where the first parent vertices in the vecto
	{
		j--;
	}
	int i = 0;
	std::cout << "i=" << i << std::endl;
	for (int k = j; k >= i; k--)
	{
		std::cout << "k=" << k << std::endl;
		int count1 = 0, count2 = 0;
		if (2 * k + 2 == n)
		{
			if (v[(2 * k) + 1] > v[k])
			{
				//display_vector(v);
				v.swap(k, 2 * k + 1);
				//counter++;
				std::cout << "a" << std::endl;
				display_vector(v);
			}
		}
		else {
			// Check which leaf is the largest                    -->This works<-- 
			double max_leaf = std::max(v[2 * k + 2], v[2 * k + 1]);
			//std::cout << "max_leaf = " << max_leaf << std::endl;
			if (max_leaf > v[k]) // Only execute if child value > parent value
			{
				if (max_leaf == v[2 * k + 2])
				{
					//display_vector(v);
					v.swap(k, 2 * k + 2);
					count2++;
					//std::cout << "k=" << k << std::endl;
					display_vector(v);
				}
				else {
					//display_vector(v);
					v.swap(k, 2 * k + 1);
					count1++;
					//std::cout << "k=" << k << std::endl;
					display_vector(v);
				}
			}
		}
		//std::cout << "k=" << k << std::endl;
		
		if (k == i)
		{
			std::cout << "count1= " << count1 << " count2= " << count2 << std::endl;
			
			if (count1 == 1)
			{
				// create temporary vector
				MVector tempvec1;
				tempvec1.push_back(v[2 * k + 1]);
				for (int m = i; pow(2, m) + 1 <= v.size(); m++)
				{
					tempvec1.push_back(v[pow(2, m) + 1]);
					tempvec1.push_back(v[pow(2, m) + 2]);
				}
				//heap_from_root(tempvec1, 0, tempvec1.size());
				// Add the temporary vector back
				

			}
			
			if (count2 == 1)
			{
				// create temporary vector
				MVector tempvec2;
				tempvec2.push_back(v[pow(2, k + 1)]);
				for (int p = i; pow(2, p) + 1 <= v.size()/2; p++) //i always 0
				{	
					int l = p;
					std::cout <<"p="<< p << std::endl;
					while (l < pow(2, p))
					{
						std::cout << l << std::endl;
						//tempvec2.push_back(v[pow(2, p) + l]);
						l++;
					}
					
					//tempvec2.push_back(v[pow(2, p) + 1]);
					//tempvec2.push_back(v[pow(2, p) + 2]);
				}
				display_vector(tempvec2);
				//heap_from_root(tempvec2, 0, tempvec2.size());
				// Add the temporary vector back
			}
		}
	
	}
	
	*/
	// Testing all sorting algorithms
	/*
	MVector w(4);
	w.initialise_random(-1000000.0, 1000000.0);
	display_vector(w);
	std::cout << "Using bubblesort algorithm: " << std::endl;
	bubble(w);
	display_vector(w);
	std::cout << std::endl;
	MVector x(10);
	x.initialise_random(-100000.0, 100000.0);
	display_vector(x);
	std::cout << " Using quicksort algorithm: " << std::endl;
	quick(x);
	display_vector(x);
	std::cout << std::endl;
	MVector y(10);
	y.initialise_random(-100000.0, 100000.0);
	display_vector(y);
	std::cout << " Using heapsort algorithm: " << std::endl;
	heap(y);
	display_vector(y);
	*/
	/*
	// Testing initalise random 
	// Testing randomness of data
	std::ofstream test_randomness;
	test_randomness.open("test_randomness.txt");
	if (!test_randomness)
	{
		std::cout << "Couldn't open file!" << std::endl;
		return 1;
	}
	// Write this to txt file
	for (int k = 0; k < 10; k++)
	{
		MVector rand_vec(2048);
		rand_vec.initialise_random(-1000000.0, 1000000.0);
		for (int j = 0; j < rand_vec.size(); j++)
		{
			test_randomness.width(20); test_randomness << rand_vec[j];
		}
		test_randomness << std::endl;
	}
	test_randomness.close();

	
	MVector a(10,1), b(30),c(100); 
	a.initialise_random(0.0, 121.0);
	b.initialise_random(10.0, 20.0);
	c.initialise_random(1.0, 100.0);
	//std::cout << "Vector a; " << std::endl; display_vector(a);
	//std::cout << "Vector b; " << std::endl; display_vector(b);
	//std::cout << "Vector c; " << std::endl; display_vector(c);
	// Testing the swap function 
	*/
	/*
	display_vector(v);
	v.swap(3, 4);
	display_vector(v);
	*/

	return 0;
}






/*  =========  4.2.1  ==========
// Tasks 1 - 4
#ifndef MVECTOR_H // the 'include guard'
#define MVECTOR_H // see C++ Primer Sec. 2.9.2
#include <vector>
#include <cassert>
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
	double &operator[](int index)
	{
		// Assert Errors

		return v[index];
	}
	// access element (rvalue) (see example sheet 5, q5.7)
	double operator[](int index) const {
		//Assert Errors
		//assert(index>=v.size());
		return v[index];
	}
	int size() const { return v.size(); } // number of elements
	void swap(int i, int j)
	{
		double tmp = v[i];
		v[i] = v[j];
		v[j] = tmp;
	}
private:
	std::vector<double> v;
};
#endif




int main()
{
	// Task 1
	// Simple program that creates an MVector of size 10
	MVector v(10, 1.0);
	// Changing MVector v
	for (int i = 1; i < v.size(); i++)
	{
		v[i] += v[i - 1];
	}
	// Output data to file
	std::ofstream Output_1;
	Output_1.open("Output_1.txt");
	if (!Output_1)
	{
		std::cout << "Couldn't open file!!" << std::endl;
		return 1;
	}
	for (int j = 0; j < v.size(); j++)
	{
		Output_1.width(5);  Output_1 << v[j];
	}
	Output_1 << std::endl;
	Output_1.close(); //Close file


	// Testing the swap function
	std::cout << "v[3]= " << v[3] << ", v[4]= " << v[4] << std::endl;
	v.swap(3, 4);
	std::cout << "v[3]= " << v[3] << ", v[4]= " << v[4] << std::endl;

	//std::cout << a[3] << " " << a[0] << std::endl;
	//std::cout << a[a.size()] << std::endl;
	return 0;
}



*/


