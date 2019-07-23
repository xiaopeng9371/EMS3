#ifndef __EMS2_H_
#define __EMS2_H_

#include "ems.hpp"
#include "motif_tree.hpp"
#include "motif_tree_fast.hpp"


template<class T>
struct ems { static const char * name; };
template<> const char *ems<MotifTreeFast>::name = "ems2";
template<> const char *ems<MotifTreeSlow>::name = "ems2m";

template <class MotifTree>
class Ems2 : public MotifFinder< Ems2<MotifTree> > {
  using MotifFinder< Ems2<MotifTree> >::motifs;
  using MotifFinder< Ems2<MotifTree> >::l;
  using MotifFinder< Ems2<MotifTree> >::d;
  using MotifFinder< Ems2<MotifTree> >::reads;
  using MotifFinder< Ems2<MotifTree> >::domain;
  using MotifFinder< Ems2<MotifTree> >::domain_size;
  MotifTree main_tree;
  MotifTree *curr_tree;
  size_t count = 0;
  int leftmost;
  int rightmost;
  std::string kmer;
  std::vector<std::string> nbrs;
  int stars;
  double gen_time=0.0;

  void gen_nbrhood3(int start, int alpha) {
  int len = kmer.size();
    if (alpha > 0) {
    int j;
    for (j=start; j<len; j++) {
      if (kmer[j] == '*') continue;
      if (kmer[j] == '-') continue;
      if ((j>0) && (kmer[j-1] == '-')) continue;
      kmer.insert(j,1,'*');
        gen_nbrhood3(j+1, alpha-1); 
      kmer.erase(j,1);
    }
    if (rightmost && (kmer[j-1] != '-')) {
      kmer.insert(j,1,'*');
      gen_nbrhood3(j+1, alpha-1); 
      kmer.erase(j,1);
      }
    } else {
      std::string t(kmer);
      int i=t.size()-1;
      while (i>=0) {
        if (t[i] == '-') t.erase(i,1);
        if (t[i] == '*') t[i] = domain_size;
        i--;
      }
      curr_tree->insert(t);
    }
  }

  void gen_nbrhood2(int start, int sigma, int alpha) {
  int len = kmer.size();
    if (sigma > 0) {
    for (int j=start; j<len; j++) {
      if (kmer[j] == '-') { while (kmer[++j] == '-'); continue; }
      if (!rightmost && (stars+1 == j) && (kmer[j+1] == '-')) continue;
      char t = kmer[j];
      kmer[j] = '*';
      int old_stars = stars;
      if (stars+1 == j) stars++;
      gen_nbrhood2(j+1, sigma-1, alpha); 
      kmer[j] = t;
      stars = old_stars;
      }
    } else {
    int new_start = leftmost ? 0 : stars+2;
    gen_nbrhood3(new_start, alpha);
    }
  }

void gen_nbrhood(int start, int delta, int sigma, int alpha) {
  int len = kmer.size();
    if (delta > 0) {
    for (int j=start; j<len; j++) {
      char t = kmer[j];
      //kmer.erase(j,1);
      kmer[j] = '-';
      gen_nbrhood(j+1, delta-1,sigma,alpha);
      //kmer.insert(j,1,t);
      kmer[j] = t;
      }
    } else {
    gen_nbrhood2(0, sigma, alpha);
    }
  }

  void gen_all(std::string& seq, int l, int d) {
    int m = seq.size();
    for (int q=-d; q<=+d; q++) {
      int k = l+q;
      for (int delta = std::max(0,q); delta <= (d+q)/2; delta++) {
        int alpha = delta - q;
        int sigma = d - alpha - delta;
        for (int i=0; i<m-k+1; i++) {
        kmer = seq.substr(i, k);
        int start = 1; stars = -1;
        leftmost = rightmost = 0;
        if (i == 0) { leftmost = 1; }
        if (i+k == m) { start = 0; rightmost = 1; }
        gen_nbrhood(start, delta, sigma, alpha);
        }
      }
    }
  }

  public:

  Ems2<MotifTree>(const std::string &input, int l, int d, Params &params):
    MotifFinder<Ems2<MotifTree> >(ems<MotifTree>::name, input, l, d, params), main_tree(l, motifs, "main"), curr_tree(&main_tree) {
      main_tree.setDomain(domain);
    }

  void search() {
    //std::cout << "Processing sequence " << 0 << "..." << std::endl;
    std::cout.flush();
    gen_all(reads[0], l, d);
    for (size_t i=1; i<reads.size(); i++) {
      //std::cout << "Processing sequence " << i << "..." << std::endl;
      std::cout.flush();
      Motifs tmp;
      MotifTree tmp_tree(l, tmp, "tmp");
      tmp_tree.setDomain(domain);
      curr_tree = &tmp_tree;
      gen_all(reads[i], l, d);
      main_tree.intersect(curr_tree);
    }
    main_tree.traverse();
  }

};

#endif // __EMS2_H_

