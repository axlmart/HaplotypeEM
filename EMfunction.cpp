#include <iostream> 
#include <fstream> 
#include <map> 
#include <vector> 
#include <string>
#include "Matrix615.h"
#include <cmath>
#include <set>
using namespace std;


vector<int> EMalgoDC (vector<int> M1, vector<int> M2);
//set <int> EM_Merge (set <int> List1, set <int> List2);


int main(int argc, char** argv){
	Matrix615<int> my_matrix;
	my_matrix.readFromFile(argv[1]);


// Divide the matrix in submatrices of 2 columns including all rows.
// Call outside EM function on each subproblem. Hi Axel!

// Note for now (at least) the entries are intergers (because changing the Matrix.h 
// to allow for strings seemed complicated). Thus each value represents a different genotype :
// integer value --> genotype 

// 0 --> 0/0
// 1 --> 0/1
// 10 --> 1/0
// 11 --> 1/1
//
//
// Note we will need to modify this loop to account for the possibility of odd # of markers.

vector<int> M1_1;
vector<int> M2_1;
vector<int> M1_2;
vector<int> M2_2;
//vector<vector<int, int>, int> Haplolist;
for ( int i = 0 ; i < my_matrix.colNums(); i=i+4) {
	M1_1 = my_matrix.GetCol(i);
	M2_1 = my_matrix.GetCol(i+1);
	EMalgoDC(M1_1,M2_1);

        M1_2 = my_matrix.GetCol(i+2);
        M2_2 = my_matrix.GetCol(i+3);

	EMalgoDC(M1_2,M2_2);
}

return 0;
}

vector<int> EMalgoDC (vector<int> M1, vector<int> M2) {

// First step : Find unambiguous haplotypes

	// Set initial counts to zero for all haplotypes :
	int n00 = 0, n01 = 0, n10 = 0 , n11 = 0 ;

	for (int j =0 ; j < M1.size() ; j++ ) {
		//all +2 cases
		if (M1[j] == 0 && M2[j] == 0) {
			n00 = n00 +2;}
		if (M1[j] == 11 && M2[j] == 11) {
                        n11 = n11 +2;}
                if (M1[j] == 0 && M2[j] == 11) {
                        n01 = n01 +2;}
		if (M1[j] == 11 && M2[j] == 0) {
                        n10 = n10 +2;}		
		// all +1 cases
		if (M1[j] == 0 && M2[j] == 1) {
                        n01 = n01 +1;
			n00 = n00 +1;}
		if (M1[j] == 1 && M2[j] == 0) {
                        n10 = n10 +1;
			n00 = n00 +1;}
		if (M1[j] == 11 && M2[j] == 1) {
                        n10 = n10 +1;
                        n11 = n11 +1;}
                if (M1[j] == 1 && M2[j] == 11) {
                        n01 = n01 +1;
                        n11 = n11 +1;}
		if (M1[j] == 10 && M2[j] == 0) {
                        n10 = n10 +1;
                        n00 = n00 +1;}
		if (M1[j] == 0 && M2[j] == 10) {
                        n01 = n01 +1;
                        n00 = n00 +1;}
		if (M1[j] == 10 && M2[j] == 11) {
                        n01 = n01 +1;
                        n11 = n11 +1;}
		if (M1[j] == 11 && M2[j] == 10) {
                        n10 = n10 +1;
                        n11 = n11 +1;}
	}

	// We can merge some of the if statements 
	//
	//
	// Now that all unambiguous counts have been generated we need the counts of all ambiguous genotypes :
	// all remaining genotypes all give the possibility of haveing all 4 haplotypes thus we can do a simple substraction to get the count

int Nam = 2*M1.size() - n00 - n01 -n10 - n11;
cout << "Number of ambiguous haplotypes : " << Nam << endl;
// Ready to start EM ! We need a thrshol bound on precision for this not to run for ever and to set initial values for the probabilities

double p00 = 0.25 , p01 = 0.25 , p10 = 0.25 , p11 = 0.25;
double p00old = 0, p01old = 0 , p10old = 0 , p11old = 0;
int it = 0;
while( abs(p00-p00old)+abs(p01-p01old)+abs(p10-p10old)+abs(p11-p11old) > pow(10,-6)  ) {

it++;

double denom = 2*p01*p10+2*p00*p11;
double En00 = Nam*p00*p11/denom + n00;
double En01 = Nam*p01*p10/denom + n01;
double En10 = Nam*p01*p10/denom + n10;
double En11 = Nam*p00*p11/denom + n11;

p00old = p00;
p01old = p01;
p10old = p10;
p11old = p11;

p00 = En00/(2*M1.size());
p01 = En01/(2*M1.size());
p10 = En10/(2*M1.size());
p11 = En11/(2*M1.size());

cout << "Iteration : " << it << ", p00 : " << p00 << ", p01 : " << p01 << ", p10 : " << p10 << ", p11 : " << p11 << endl;
}

// Get list of possible haplotypes from final haplotype probs :

vector<int> HaplotypeList;
if ( p00 > pow(10,-10)) {HaplotypeList.push_back(0);}
if ( p01 > pow(10,-10)) {HaplotypeList.push_back(1);}
if ( p10 > pow(10,-10)) {HaplotypeList.push_back(10);}
if ( p11 > pow(10,-10)) {HaplotypeList.push_back(11);}


return HaplotypeList;	

}

// Need merging function here : 
//
// Takes Input : Genotypes from 4 columns, list of possible haplotype on these columns, set equal frequency to each (within function.
// outputs possible haplotypes at these positions (perhaps keep track of frequencies for later use).
 
/*
vector<int> EM_Merge(set<int> List1 , set<int> List2) {

// First step : Read in values from both list and create all possible combinations of haplotypes based on these.
// Note there is a maximum of 2^4 haplotypes (ie if our list contain all contain the previous 4 possible haplotype we will have gained nothing.
//
// We have recreate all possible haplotypes from the given ones : find combinatorical arrangements functions in c++ ? /do it ourselves ?
//
// After that we will have to generate flexible structure based on the above results in order to solve the greater problem.
//
//
// Need to think about how to deal with this when confronted to large data bases (more than 4 markers). --> or make larger and larger stacks of two until we solve the overall pb ?


set<int>::iterator iter;
cout << " First haplotype list : " ;
for (iter=List1.begin(); iter !=List1.end() ; ++iter) {
 cout << *iter << "    ";
}

cout << endl << endl;

cout << " Second haplotype list : " ;
for (iter=List2.begin(); iter !=List2.end() ; ++iter) {
 cout << *iter << "    ";
}

cout << endl << endl;

vector<int> PosHaplo;

if(std::find(M1.begin(), M1.end(), 0) != M1.end() && std::find(M2.begin(), M2.end(), 0) != M2.end()) {
PosHaplo.push_back(0);
}

if(std::find(M1.begin(), M1.end(), 0) != M1.end() && std::find(M2.begin(), M2.end(), 1) != M2.end()) {
PosHaplo.push_back(1);
}

if(std::find(M1.begin(), M1.end(), 0) != M1.end() && std::find(M2.begin(), M2.end(), 10) != M2.end()) {
PosHaplo.push_back(10);
}

if(std::find(M1.begin(), M1.end(), 0) != M1.end() && std::find(M2.begin(), M2.end(), 11) != M2.end()) {
PosHaplo.push_back(11);
}

if(std::find(M1.begin(), M1.end(), 1) != M1.end() && std::find(M2.begin(), M2.end(), 0) != M2.end()) {
PosHaplo.push_back(0100);
}

if(std::find(M1.begin(), M1.end(), 1) != M1.end() && std::find(M2.begin(), M2.end(), 1) != M2.end()) {
PosHaplo.push_back(0101);
}

if(std::find(M1.begin(), M1.end(), 1) != M1.end() && std::find(M2.begin(), M2.end(), 10) != M2.end()) {
PosHaplo.push_back(0110);
}

if(std::find(M1.begin(), M1.end(), 1) != M1.end() && std::find(M2.begin(), M2.end(), 11) != M2.end()) {
PosHaplo.push_back(0111);
}

if(std::find(M1.begin(), M1.end(), 10) != M1.end() && std::find(M2.begin(), M2.end(), 0) != M2.end()) {
PosHaplo.push_back(1000);
}

if(std::find(M1.begin(), M1.end(), 10) != M1.end() && std::find(M2.begin(), M2.end(), 1) != M2.end()) {
PosHaplo.push_back(1001);
}

if(std::find(M1.begin(), M1.end(), 10) != M1.end() && std::find(M2.begin(), M2.end(), 10) != M2.end()) {
PosHaplo.push_back(1010);
}

if(std::find(M1.begin(), M1.end(), 10) != M1.end() && std::find(M2.begin(), M2.end(), 11) != M2.end()) {
PosHaplo.push_back(1011);
}

if(std::find(M1.begin(), M1.end(), 11) != M1.end() && std::find(M2.begin(), M2.end(), 0) != M2.end()) {
PosHaplo.push_back(1100);
}

if(std::find(M1.begin(), M1.end(), 11) != M1.end() && std::find(M2.begin(), M2.end(), 1) != M2.end()) {
PosHaplo.push_back(1101);
}

if(std::find(M1.begin(), M1.end(), 11) != M1.end() && std::find(M2.begin(), M2.end(), 10) != M2.end()) {
PosHaplo.push_back(1110);
}

if(std::find(M1.begin(), M1.end(), 1111) != M1.end() && std::find(M2.begin(), M2.end(), 1111) != M2.end()) {
PosHaplo.push_back(1111);
}


// Now vector PosHaplo contains all possbile haplotypes at these 4 markers. Need to find all unambiguous genotypes again ... Things are going to get messy here. Need ideas

}*/


