#include <iostream> 
#include <fstream> 
#include <map> 
#include <vector> 
#include <string>
#include "Matrix615.h"
#include <cmath>
#include <set>
#include <map>
#include <utility>
#include <functional>
#include <algorithm>
using namespace std;

void init_EM( map<int, std::set <int> > pos_miniHaps_byPerson , map<int, int> obs_miniHaps , Matrix615<int>&
	my_matrix , map <int , vector<int> >& postEM1_miniHaps_byLocus , int locus_num, map<int, vector<set<int> > >& pos_miniHaps_allLoci_byPerson, map<int , vector<int> >& indicatorMap_byLoci){
	
	//pos_miniHaps_byPerson is the map of possible haplotypes
	//obs_miniHaps counts the number of observed haplotypes
	//my_matrix is the data file
	//postEM1_miniHaps_byLocus gets updated with possible haplotypes at each 2 or 3 bp "locus"
	//locus_num is the set of 2 or 3 bps we are looking at. Locus num ranges from o to floor(#markers)
	
	int locus_size = 2; //start out assuming we are looking at 2 markers
	if ((my_matrix.colNums()-locus_num*2)==3){ locus_size=3;} //if we are looking at the last locus
								//and it is of size 3, account for it
	int index_adj=(locus_size==3 ? 1 : 0)*4;
	
	//should allocate and delete this memory
	map<int, double> hap_prob;
	for (int i=0+index_adj ; i<pow(2,locus_size)+index_adj; ++i){
		hap_prob[i]= 1 / double(pow(2,locus_size));
	}
	
	// Getting diplotypes 
	map<int, set<pair<int,int> > > Markers3_dip ;
	
	for (int n=0; n<my_matrix.rowNums(); ++n){//for N rows in locus matrix 
		if (locus_size==2){// if num markers in the locus =2
			for(int m=0; m<locus_size; ++m){ //for each marker in the locus of size 2
				if(m==0){
					if(my_matrix.data[n][locus_num*2 + m]==0){//if we have two 0s at the first marker
						pos_miniHaps_byPerson[n].erase(2); pos_miniHaps_byPerson[n].erase(3); }//remove 2,3 from pos_miniHaps_byPerson[n]
					else if(my_matrix.data[n][locus_num*2 + m]==11){//if we have two ones at the first marker
						pos_miniHaps_byPerson[n].erase(0); pos_miniHaps_byPerson[n].erase(1); } //remove 0,1
				} else{
					if(my_matrix.data[n][locus_num*2 + m]==0){//if we have two  0s at the second marker
						pos_miniHaps_byPerson[n].erase(1); pos_miniHaps_byPerson[n].erase(3); } //remove 1 and 3
					else if(my_matrix.data[n][locus_num*2 + m]==11){ //if we have two 1s at the second marker
						pos_miniHaps_byPerson[n].erase(0); pos_miniHaps_byPerson[n].erase(2); } //remove 0 and 2
				}
			}
			
			std::set<int>::iterator itDip1;
			std::set<int>::iterator itDip2;
			itDip1=pos_miniHaps_byPerson[n].begin();
			itDip2= pos_miniHaps_byPerson[n].end();
			--itDip2;
			if(pos_miniHaps_byPerson[n].size() == 1) {Markers3_dip[n].insert(std::make_pair(*itDip1,*itDip1) ) ; }
			if(pos_miniHaps_byPerson[n].size() == 2) {Markers3_dip[n].insert(std::make_pair(*itDip1,*itDip2) ); }
			if(pos_miniHaps_byPerson[n].size() == 4) {Markers3_dip[n].insert(std::make_pair(0,3) ); 
													  Markers3_dip[n].insert(std::make_pair(1,2) );} 
		}
		
		else if (locus_size==3){// if num markers in the locus = 3
			for (int m=0; m<locus_size; ++m){ //for each marker in the locus of size 3
				if (m==0){
					if (my_matrix.data[n][locus_num*2 + m]==0){
						pos_miniHaps_byPerson[n].erase(8); pos_miniHaps_byPerson[n].erase(9);
						pos_miniHaps_byPerson[n].erase(10); pos_miniHaps_byPerson[n].erase(11);}
					else if(my_matrix.data[n][locus_num*2 + m]==11){
						pos_miniHaps_byPerson[n].erase(4); pos_miniHaps_byPerson[n].erase(5);
						pos_miniHaps_byPerson[n].erase(6); pos_miniHaps_byPerson[n].erase(7);}
				} else if (m==1){
					if (my_matrix.data[n][locus_num*2 + m]==0){
						pos_miniHaps_byPerson[n].erase(6); pos_miniHaps_byPerson[n].erase(7);
						pos_miniHaps_byPerson[n].erase(10); pos_miniHaps_byPerson[n].erase(11);}
					else if(my_matrix.data[n][locus_num*2 + m]==11){
						pos_miniHaps_byPerson[n].erase(4); pos_miniHaps_byPerson[n].erase(5);
						pos_miniHaps_byPerson[n].erase(8); pos_miniHaps_byPerson[n].erase(9);}
				} else{
					if (my_matrix.data[n][locus_num*2 + m]==0){
						pos_miniHaps_byPerson[n].erase(5); pos_miniHaps_byPerson[n].erase(7);
						pos_miniHaps_byPerson[n].erase(9); pos_miniHaps_byPerson[n].erase(11);}
					else if(my_matrix.data[n][locus_num*2 + m]==11){
						pos_miniHaps_byPerson[n].erase(4); pos_miniHaps_byPerson[n].erase(6);
						pos_miniHaps_byPerson[n].erase(8); pos_miniHaps_byPerson[n].erase(10);}
				}
			}
			
			std::set<int>::iterator itDip1;
			std::set<int>::iterator itDip2;
			itDip1=pos_miniHaps_byPerson[n].begin();
			itDip2= pos_miniHaps_byPerson[n].end();
			--itDip2;
			if(pos_miniHaps_byPerson[n].size() == 1) {Markers3_dip[n].insert(std::make_pair(*itDip1,*itDip1) ) ; }
			if(pos_miniHaps_byPerson[n].size() == 2) {Markers3_dip[n].insert(std::make_pair(*itDip1,*itDip2) ); }
			if(pos_miniHaps_byPerson[n].size() > 2) {
				while((int)*itDip1 <= (int)*itDip2){
					Markers3_dip[n].insert(std::make_pair(*itDip1,*itDip2) );
					--itDip2;
					++itDip1;
				}
			}
		}
		else{}
		
		std::set<int>::iterator it;
		if (pos_miniHaps_byPerson[n].size()==2){//if size pos_miniHaps_byPerson[n]==2 iterate through set using *iter and  ++obs_miniHaps[*it]
			for(it=pos_miniHaps_byPerson[n].begin(); it!=pos_miniHaps_byPerson[n].end(); ++it){
				++obs_miniHaps[*it];
			}
		}
		else if (pos_miniHaps_byPerson[n].size()==1){ //if size pos_miniHaps_byPerson[n]==1 obs_miniHaps[*it]=obs_miniHaps[*it]+2
			it=pos_miniHaps_byPerson[n].begin();
			obs_miniHaps[*it]=obs_miniHaps[*it]+2;
		}
		else{}
		
		pos_miniHaps_allLoci_byPerson[n].push_back(pos_miniHaps_byPerson[n]);
		
	}
			

	//Run EM udating observed map and hap_prob
	// First get the diplotypes 
	
	//cout<< Markers3_dip.size();
	
	std::set<pair<int,int> >::iterator it6;
	for (int i=0; i<Markers3_dip.size();i++){
		for(it6=Markers3_dip[i].begin();it6!=Markers3_dip[i].end();++it6){
		}
	}
	
		//initialize prob :
	for (int i=0+index_adj; i<pow(2,locus_size)+index_adj; ++i){ 
		hap_prob[i] = (double)(1/pow(2,locus_size));
	}

	// create temp storage for prob starting with values = 0 
	map<int, double> hap_prob_temp;
	for (int i=0+index_adj; i<pow(2,locus_size)+index_adj; ++i){ 
		hap_prob_temp[i] = (double)(0);
	}
	
 	int pos_dip;
	std::set<pair<int,int> >::iterator it3M;
	int iterations = 0;
	double dif = 0;
	double Phats;
	double denom;
	int breaker = 0 ;
	
	// EM Algo
	// First loop until convergence
		// second goes through each haplotype
			// third goes each individual
				// fourth goes through all possible diplotypes of that person
	
	
	while (breaker == 0) { 
                
		dif = 0;// will be used to break	
		// set our new values from prev iterations to "old" values if we have already run the first iteration only
		if (iterations != 0 ) {
			for (int i=0+index_adj; i<pow(2,locus_size)+index_adj; ++i){ 
				hap_prob[i] = hap_prob_temp[i];
				}
		}
		iterations++;
	
		for (int k=0+index_adj; k<pow(2,locus_size)+index_adj; ++k) { // or up to length of possible haplotypes 
		
			double obs_counts =  obs_miniHaps[k];
			for (int i = 0 ; i < my_matrix.rowNums() ; ++i) {
			
				pos_dip = Markers3_dip[i].size() ; // get length or number of possible diplotypes for indiv i 
				Phats = 0; //initialize numerator
				denom = 0 ;//initialize denominator
				
				if ( pos_dip == 1) {} // do nothing as we have already found these as observable counts 
				
				else {
					for(it3M=Markers3_dip[i].begin();it3M!=Markers3_dip[i].end();++it3M){
						denom = denom + hap_prob[(*it3M).first] * hap_prob[(*it3M).second] ; // multiply together all paired estimates of prob and sum over
						
						// We also want to keep track of the numerator that is the specific pair of interest
						if ( (*it3M).first == k || (*it3M).second == k) {
							Phats = hap_prob[(*it3M).first] * hap_prob[(*it3M).second];
						}
					}
					obs_counts = obs_counts + Phats/denom; // will be updated to new count every iteration, always simply overwritten
				}
			}
			
			hap_prob_temp[k] = obs_counts / (2*my_matrix.rowNums()); //once obscounts found simply divide by # of haplos
		}
		
		// to compensate now having statement in while loop
		for (int i=0+index_adj; i<pow(2,locus_size)+index_adj; ++i){ 
				dif = dif + abs(hap_prob[i] - hap_prob_temp[i]);}
		if(dif < pow(10,-5)) {breaker = 1;} // if small enough get out of the while loop
	}
	
	cout << "EM algorithm on locus " << locus_num+1 << " yields : " << endl ; 
	// now just have to push back temp values into original map
	for (int i=0+index_adj; i<pow(2,locus_size)+index_adj; ++i){ 
		hap_prob[i] = hap_prob_temp[i];
		cout << "Probability haplotype " << i << " = " << hap_prob[i] << endl; ;
		}
		cout << endl << endl;
	
	//push back hap probs >0 into postEM1_miniHaps_byLocus[locus_num]
	
	for (int i=0+index_adj; i<pow(2,locus_size)+index_adj; ++i){
		if (hap_prob[i] > 0.000001){ //change to fit whatever we decide
			postEM1_miniHaps_byLocus[locus_num].push_back(i);
			indicatorMap_byLoci[locus_num].push_back(1);
		}
		else{indicatorMap_byLoci[locus_num].push_back(0);}
	}
	
	
	
	
}


//add element of one vector to end of another, preserving both vectors
 vector<int> append_element(vector<int>& vec1 , int& hap){
 	vector<int> vec2;
	vec2.reserve(vec1.size() +1 );
	vec2.insert( vec2.end(), vec1.begin(), vec1.end());
	vec2.push_back(hap);
	return vec2;
} 
	


//trying to make this work for the first step, see if I can make it recursive.
//for each vector in the vector of vectors, add each element of the locus_miniHaps, and return a new vector of vectors
//this works needs to make it recursive
vector<vector<int> > combine_haps(vector<vector<int> >& fullHaps_preEM2, vector<int>& locus_miniHaps){
	vector<vector<int> > new_fullHaps_preEM2;
	for (int i=0; i<fullHaps_preEM2.size(); ++i){
		for (int j=0; j<locus_miniHaps.size(); ++j){
			new_fullHaps_preEM2.push_back(append_element(fullHaps_preEM2[i],locus_miniHaps[j]));
		}
	}
	
	return new_fullHaps_preEM2;
}

vector<vector<int> > combine_haps_recursive(map<int, vector<int> >& postEM1_miniHaps_byLocus){
	vector<vector<int> > new_fullHaps_preEM2;
	vector<int> dummy;
	for (int h=0; h< postEM1_miniHaps_byLocus[0].size(); ++h){
		new_fullHaps_preEM2.push_back(append_element(dummy,postEM1_miniHaps_byLocus[0][h]));
	}
	for (int i=1; i<postEM1_miniHaps_byLocus.size(); ++i){
		new_fullHaps_preEM2=combine_haps(new_fullHaps_preEM2, postEM1_miniHaps_byLocus[i]);
	}
	
	return new_fullHaps_preEM2;
}

std::set <int> remove_impossible_haps(std::set<int> fullHaps_enum , vector<set<int> >& pos_miniHaps_allLoci_personn , map <int , vector<int> >& postEM1_miniHaps_byLocus, 
	vector<int>& num_fullHaps_vec, map<int , vector<int> >& indicatorMap_byLoci){
	
	
	for( int i=0; i<postEM1_miniHaps_byLocus.size(); ++i){
		int index_adj=(indicatorMap_byLoci[i].size()==8 ? 1 : 0)*4; //is this a 3bp locus?
		for(int j=0; j<postEM1_miniHaps_byLocus[i].size(); ++j){
			int curr_hap=postEM1_miniHaps_byLocus[i][j];
			int is_poss_for_n=(pos_miniHaps_allLoci_personn[i].find(curr_hap) != pos_miniHaps_allLoci_personn[i].end() ? 1 : 0);
			if (is_poss_for_n==0){ //if the person cannot have that hapotype at the locus, remove all all_locus haps with that hap at the locus
				int skip1=0;
				for (int m=0+index_adj; m<curr_hap; m++){ //potential problem site!!!!
					skip1 = skip1 + num_fullHaps_vec[i+1]*indicatorMap_byLoci[i][m-index_adj];
				}
				int skip2=0;
				for (int m=curr_hap+1; m<indicatorMap_byLoci[i].size()+index_adj; ++m){
					skip2= skip2 + num_fullHaps_vec[i+1]*indicatorMap_byLoci[i][m-index_adj];
				}
				int l=0;
				while(l < num_fullHaps_vec[0]){
					l=l+skip1;
					int curr_l=l;
					while(l < curr_l+num_fullHaps_vec[i+1]){
						fullHaps_enum.erase(l);
						++l;
					}
					l=l+skip2;
					//cout<<l << " is less than" << num_fullHaps_vec[0] <<endl;
				}
			}
		}
	}
	
	std::set<int>::iterator it7;
	for(it7=fullHaps_enum.begin();it7!=fullHaps_enum.end();++it7){
		//cout<< (*it7) << " ";
		}
	//cout<<endl;
	return fullHaps_enum;
}

map<int, vector<int> > expand_coded_haps(vector<vector<int> >& full_hap_by_loci){
	map<int, vector<int> > expanded_haps;
	for (int i=0; i<full_hap_by_loci.size(); ++i){
		for (int j=0; j<full_hap_by_loci[i].size(); ++j){
			if(full_hap_by_loci[i][j]==0){expanded_haps[i].push_back(0); expanded_haps[i].push_back(0);}
			else if(full_hap_by_loci[i][j]==1){expanded_haps[i].push_back(0); expanded_haps[i].push_back(1);}
			else if(full_hap_by_loci[i][j]==2){expanded_haps[i].push_back(1); expanded_haps[i].push_back(0);}
			else if(full_hap_by_loci[i][j]==3){expanded_haps[i].push_back(1); expanded_haps[i].push_back(1);}
			else if(full_hap_by_loci[i][j]==4){expanded_haps[i].push_back(0); expanded_haps[i].push_back(0); expanded_haps[i].push_back(0);}
			else if(full_hap_by_loci[i][j]==5){expanded_haps[i].push_back(0); expanded_haps[i].push_back(0); expanded_haps[i].push_back(1);}
			else if(full_hap_by_loci[i][j]==6){expanded_haps[i].push_back(0); expanded_haps[i].push_back(1); expanded_haps[i].push_back(0);}
			else if(full_hap_by_loci[i][j]==7){expanded_haps[i].push_back(0); expanded_haps[i].push_back(1); expanded_haps[i].push_back(1);}
			else if(full_hap_by_loci[i][j]==8){expanded_haps[i].push_back(1); expanded_haps[i].push_back(0); expanded_haps[i].push_back(0);}
			else if(full_hap_by_loci[i][j]==9){expanded_haps[i].push_back(1); expanded_haps[i].push_back(0); expanded_haps[i].push_back(1);}
			else if(full_hap_by_loci[i][j]==10){expanded_haps[i].push_back(1); expanded_haps[i].push_back(1); expanded_haps[i].push_back(0);}
			else {expanded_haps[i].push_back(1); expanded_haps[i].push_back(1); expanded_haps[i].push_back(1);}
			
		}
	}
	return expanded_haps;
}


			
map<int, set<pair<int,int> > > diplotyper2(map<int, set<int> >& pos_fullHaps_enum_byPerson, Matrix615<int>& my_matrix, map<int, vector<int> >& Expanded_haps){
	map<int, set<pair<int,int> > > diplotype_map;
	for(int n=0; n<pos_fullHaps_enum_byPerson.size(); ++n){
		
		set<int>::iterator itup;
		set<int>::iterator itdown;
		itup=pos_fullHaps_enum_byPerson[n].begin();
		itdown=pos_fullHaps_enum_byPerson[n].end();
		--itdown;
		if(pos_fullHaps_enum_byPerson[n].size()<=2){
		
			diplotype_map[n].insert(std::make_pair(*itup,*itdown));

		} else{
			while((int)*itup < int(*itdown)){
				bool Match = true;
				for (int m = 0; m < my_matrix.colNums() ; m++ ) {
					
					std::pair <int,int> a_diplotype = std::make_pair(*itup,*itdown);
					int G = my_matrix.data[n][m];
					int h1 = Expanded_haps[a_diplotype.first][m];
					int h2 = Expanded_haps[a_diplotype.second][m];
						if ( ( h1 == 0 && h2 == 0 && G != 0) || (h1==1 && h2==1 && G!=11) || (h1==0 && h2==1 && G!=1) || (h1==1 && h2==0 && G!=1) ) {
							Match = false;	
							break;	
						}			
				}	
				
				if (Match == true ) {diplotype_map[n].insert(std::make_pair(*itup,*itdown));}
								
				++itup; --itdown;
			}
		}
	}
	cout << endl << endl ;
	std::set<pair<int,int> >::iterator it12;
	for (int i=0; i<diplotype_map.size();i++){
		cout<< "Person " << i+1 << " diplotypes: ";
		for(it12=diplotype_map[i].begin();it12!=diplotype_map[i].end();++it12){
			cout<<"(" <<(*it12).first << "," << (*it12).second <<") ";
		}
		cout<<endl;
	}
	
	cout << endl<< endl;
	
	return diplotype_map;
}
						
		


void final_EM(Matrix615<int>& my_matrix , map <int , vector<int> >& postEM1_miniHaps_byLocus, map<int, double>& hap_freqs, map<int, vector<set<int> > >& pos_miniHaps_allLoci_byPerson,
			 map<int , vector<int> >& indicatorMap_byLoci ){
	
	int num_locus=postEM1_miniHaps_byLocus.size();
	int num_locus_plus[num_locus]; //the number of haplotypes that can be consturcted using all subsets to the right of and including the current subset
						//for example, if we found 3,2,4,and 2 possible haploytpes subsets 0-3, when looking subset 1, we would have 16=4*2*2
						//this is used for enumeration purposes later. The first elelment won't end up being used, but we keep it for interpritability
	
	num_locus_plus[num_locus]=1;
	int j=num_locus-1;
	while(j>=0){
		num_locus_plus[j]=num_locus_plus[j+1]*postEM1_miniHaps_byLocus[j].size();
		j=j-1;
	}
	
	
	int total_haps=num_locus_plus[0];
	
	vector<int> num_fullHaps_vec (num_locus_plus, num_locus_plus+num_locus+1);

	
	//create hap frequency map
	
	for(int i=0; i<num_locus_plus[0]; i++) { //num_locus_plus[0] is total number of possible haplotypes
		hap_freqs[i]=1/double(num_locus_plus[0]);
	}
	
	//need a map that keeps track of the possible haplotypes for each person. each person will have up to num_locus_plus[0] possible haplotypes
	//go through each person, systematically add or remove possible haplotypes
	//at each d
	int hap_enumeration[total_haps];
	for(int i=0; i<total_haps; ++i){
		hap_enumeration[i]=i;
	}
	
	std::set<int> fullHaps_enum (hap_enumeration, hap_enumeration+total_haps);
	
	//create a function that removes haplotypes from hap_enumeration if a person cant have that haplotype
	//call this function for each person and put results into an possible haplotypes table and update obs_counts just like in init_EM
	//can use postEM1_miniHaps_byLocus to create indicator for if a hap is possible at a loci, which is needed to figure out which enumerations to remove
	
	map<int, set<int> > pos_fullHaps_enum_byPerson;
	
	for(int i=0; i<my_matrix.rowNums(); ++i){
		pos_fullHaps_enum_byPerson[i]=remove_impossible_haps(fullHaps_enum, pos_miniHaps_allLoci_byPerson[i], postEM1_miniHaps_byLocus, num_fullHaps_vec, indicatorMap_byLoci);} 
	
		
	std::set<int>::iterator it4;

	std::set<int>::iterator it9;
		
	map<int,int> obs_fullHaps;
	for(int p=0; p<pos_fullHaps_enum_byPerson.size(); ++p){
		for(it9=pos_fullHaps_enum_byPerson[p].begin();it9!=pos_fullHaps_enum_byPerson[p].end();++it9){	
			if(pos_fullHaps_enum_byPerson[p].size()==2){
				++obs_fullHaps[*it9];
			}
			else if(pos_fullHaps_enum_byPerson[p].size()==1){
				obs_fullHaps[*it9]=obs_fullHaps[*it9]+2;
			}
			else{}
		}
	}
	
	
	//initizialize hap probs
	
	map<int, double> fullHap_prob;
	for (int i=0; i<total_haps; ++i){
		fullHap_prob[i]= 1 / double(total_haps);
	}
	
	vector<vector<int> > full_hap_by_loci=combine_haps_recursive(postEM1_miniHaps_byLocus);
	
	//define function that eturns 2b or 3bp loci codes back into 1s and 0s
	
	map<int,vector<int> > Expanded_haps=expand_coded_haps(full_hap_by_loci);
	
	map<int, set<pair<int,int> > >diplotypes = diplotyper2(pos_fullHaps_enum_byPerson, my_matrix, Expanded_haps);
	
	//runEM
	
	// create temp storage for prob starting with values = 0 
	map<int, double> hap_prob_temp;
	for (int i=0; i<total_haps; ++i){ 
		hap_prob_temp[i] = (double)(0);
		//cout << "initial Hap Prob temp = " << hap_prob_temp[i] << endl;
		}
	
	//cout << endl<< endl;
	
	int pos_dip;
	std::set<pair<int,int> >::iterator it3M;
	std::set<pair<int,int> >::iterator it4M;
	int iterations = 0;
	double dif = 0;
	double Phats;
	double denom;
	int breaker = 0 ;
	// First loop until convergence
		// second goes through each haplotype
			// third goes each individual
				// fourth goes through all possible diplotypes of that person
				

	std::vector<double> MaxPhats(my_matrix.rowNums(), 0);
	map<int, set<pair<int,int> > > Phase;
	
	while (breaker == 0) { 
	
		dif = 0;// will be used to break
				// set our new values from prev iterations to "old" values if we have already run the first iteration only
		if (iterations != 0 ) {
			for (int i=0; i<total_haps; ++i){ 
				fullHap_prob[i] = hap_prob_temp[i];
			}
		}
		iterations++;
	
		for (int k=0; k<total_haps; ++k) { // or up to length of possible haplotypes 
			
			double obs_counts =  obs_fullHaps[k];

			for (int i = 0 ; i < my_matrix.rowNums() ; ++i) {
			
				pos_dip = diplotypes[i].size() ; // get length or number of possible diplotypes for indiv i 
				Phats = 0; //initialize numerator
				denom = 0 ;//initialize denominator
				
				if ( pos_dip == 1) {
				std::set<pair<int,int> >::iterator it4M;
				it4M=diplotypes[i].begin();
				Phase[i].insert(std::make_pair((*it4M).first,(*it4M).second) );
				MaxPhats[i] = 1;
				} // do nothing as we have already found these as observable counts 
							
				
				
				else {
					for(it3M=diplotypes[i].begin();it3M!=diplotypes[i].end();++it3M){
						denom = denom + 2*fullHap_prob[(*it3M).first] * fullHap_prob[(*it3M).second] ; // multiply together all paired estimates of prob and sum over

						// We also want to keep track of the numerator that is the specific pair of interest
						if ( (*it3M).first == k || (*it3M).second == k) {
							Phats = 2*fullHap_prob[(*it3M).first] * fullHap_prob[(*it3M).second];
							
							if ( MaxPhats[i] < Phats ) {
								Phase[i].clear();
								Phase[i].insert(std::make_pair((*it3M).first,(*it3M).second) );
							}
							
						}
						
					}
					
					if ( MaxPhats[i] < Phats ) {
						MaxPhats[i] = Phats/denom;
					}
					
					
					obs_counts = obs_counts + Phats/denom; // will be updated to new count every iteration, always simply overwritten
				}
			}
			
			hap_prob_temp[k] = obs_counts / (2*my_matrix.rowNums()); //once obscounts found simply divide by # of haplos
		}
		
		// to compensate now having statement in while loop
		for (int i=0; i<total_haps; ++i){ 
				dif = dif + abs(fullHap_prob[i] - hap_prob_temp[i]);}
		if(dif < pow(10,-5) ) {breaker = 1;} // if small enough get out of the while loop
	}

	
	// now just have to push back temp values into original map
	for (int i=0; i<total_haps; ++i){ 
		fullHap_prob[i] = hap_prob_temp[i];
		cout<< "Hap " << i << " prob=" << fullHap_prob[i]<<endl;}
		

	// Phaser : 
	cout << endl << endl << "Phaser " << endl << endl ;
	
	std::set<pair<int,int> >::iterator itP1;
	std::set<pair<int,int> >::iterator itP2;
	for (int people = 0 ; people < MaxPhats.size();people++){
		for(itP1=Phase[people].begin();itP1!=Phase[people].end();++itP1){
			cout << "Individual " << people + 1 << " : (" << (*itP1).first <<","<<(*itP1).second << ")	" ;
		}
		itP2=Phase[people].begin();
		cout << "Diplotype Probability = " << MaxPhats[people] << endl;
		
		// Get the full diplotypes for each person 
	
		cout << "First Haplotype :\t" ;
		for (int i = 0; i < my_matrix.colNums(); i++) {
			cout << Expanded_haps[(*itP2).first][i] << " ";
		}
			
		cout << endl;
	
		cout << "Second Haplotype :\t" ;
		for (int i = 0; i < my_matrix.colNums(); i++) {
			cout << Expanded_haps[(*itP2).second][i] << " ";
		}
	
		cout << endl << endl ;	
	}	
	
		cout << "Iterations needed  : " << iterations << endl << endl;
	
}
					
				
//////////////////////
	
	
	


int main(int argc, char** argv){
	Matrix615<int> my_matrix;
	my_matrix.readFromFile(argv[1]);


// Divide the matrix in submatrices of 2 columns including all rows.
// Call outside EM function on each subproblem. Hi Axel!

// Note for now (at least) the entries are intergers (because changing the Matrix.h 
// to allow for strings seemed complicated). Thus each value represents a different genotype :
// integer value --> genotype 

// 0 --> 00
// 1 --> 01
// 11 --> 11
//
//
// Note we will need to modify this loop to account for the possibility of odd # of markers.

	int N = my_matrix.rowNums();
	int M = my_matrix.colNums();

/*create a map of possible haplotypes where each key is a person. 
0 corresponds to <0,0>
1 corresponds to <0,1>
2 = <1,0>
3 = <1,1>
4 = <0,0,0>
5 = <0,0,1>
6 = <0,1,0>
7 = <0,1,1>
8 = <1,0,0>
9 = <1,0,1>
10 = <1,1,0>
11 = <1,1,1>
*/
 
	int hap_nums_2bp[]={0,1,2,3}; //arrays to be turned into sets
	int hap_nums_3bp[]={4,5,6,7,8,9,10,11};
	std::set<int> hap_nums_2bp_set (hap_nums_2bp, hap_nums_2bp+4);
	std::set<int> hap_nums_3bp_set (hap_nums_3bp, hap_nums_3bp+8);

	map<int, std::set <int> > master_hapmap_2bp;

	for(int n=0; n<N; ++n){
		master_hapmap_2bp[n] = hap_nums_2bp_set;
	}

	map<int , std::set <int> > master_hapmap_3bp;

	for( int n=0; n<N; ++n){
		master_hapmap_3bp[n] = hap_nums_3bp_set;
	}
	
	map<int, int> observed_map_2bp; //these keep count of the number of haplotypes we observed
	map<int, int> observed_map_3bp;
	
	map<int, vector<int> > postEM1_miniHaps_byLocus;
	map<int, vector<std::set<int> > > pos_miniHaps_allLoci_byPerson;
	map<int, vector<int> > indicatorMap_byLoci;
	
	
	//add for loop to go through every 2/3 bp hap
	cout << endl;
	cout << "First, divide problems into subsets : " << endl << endl;
	for (int i=0; i< floor((double)M/2); i++){
		if ((M-i*2)==3){
			init_EM(master_hapmap_3bp, observed_map_3bp, my_matrix,
				postEM1_miniHaps_byLocus, i, pos_miniHaps_allLoci_byPerson, indicatorMap_byLoci);
		}
		else{
			init_EM(master_hapmap_2bp, observed_map_2bp, my_matrix,
				postEM1_miniHaps_byLocus, i, pos_miniHaps_allLoci_byPerson, indicatorMap_byLoci);
		}
	}
	
	
	map<int, double> hap_freqs;
	
	cout<<"Now merge all the subproblems (Conquer step) : "<<endl;
	
	final_EM(my_matrix , postEM1_miniHaps_byLocus, hap_freqs, pos_miniHaps_allLoci_byPerson, indicatorMap_byLoci );

	return 0;
	
		

}