int pseudo( input matrix file ) {

  // for possible haplotype counts, initialize Px(H+1) matrix with 0s (final column is "total" possible haplotype counts)
  // vector of observed haplotype counts, to be updated 

  for (marker 2m and 2m+1) { // 2 markers at a time
    for (person p) {
      for (possible haplotype h) {
        if (haplotype h not possible given person p's genotype) 
          {add 0 to person p, haplotype h in count matrix;} // have a matrix for possible haplotypes and total possible haplotype count for each person 
        else { add 1 to person p, haplotype h; // adds 1 to possible haplotype count
          add 1 to person p, total possible haplotypes; // adds 1 to total possible haplotype count
          if (total=3) {go to next person (break)}
      }
      if (total=1) { update corresponding observed haplotype count by 2 ;}
      else if (total=2) {update each corresponding observed haplotype by 1 ;}
      else {do nothing;}
    }
    
    Using observed counts, do EM to get reduced set of haplotypes;
    Add this set to list of possible haplotypes by marker pairing;
     // So we only have to keep the possible haplotypes list? //
    
  }

}

/* Ideas for improving efficiency
- Currently, time is: m/2*n*4^(m/2) --> m*n*2^(m-1)
- Keep list of genotypes recorded so if someone matches, we don't have to run the code again
- Say we have 8 markers, what is more efficient?
  12  34  56  78 -> EM on 2 markers           12  34  56   78 -> EM on 2 markers 
    \/      \/                                  \/    /    /
   1234    5678 -> EM on 4 markers             1234  /    /   -> EM 4 markers
        \/                                         \/    /
     12345678 -> EM 8 markers                    123456 /  -> EM 6 markers
                                                      \/
                                                   12345678  -> EM 8 markers
*/
