int pseudo( input matrix file ) {

  for (marker 2m and 2m+1 // two markers at a time // ) {
    for (person p) {
      for (possible haplotype h) {
        if (haplotype h not possible given person p's genotype) 
          {add 0 to person p, haplotype h in count matrix}
        else { add 1 to person p, haplotype h;
          add 1 to person p, total possible haplotypes;
          if (total=3) {go to next person (break)}
      }
      if (total=1) { update corresponding observed haplotype count by 2 }
      else if (total=2) {update each corresponding observed haplotype by 1}
      else {do nothing}
    }
    
    Using observed counts, do EM to get reduced set of ambiguous haplotypes;
    Add this set to list of possible haplotypes by marker pairing;
    
  }

}
