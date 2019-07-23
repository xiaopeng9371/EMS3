#ifndef __EMS1_H_
#define __EMS1_H_

#include "ems.hpp"

class Ems1 : public MotifFinder<Ems1> {

  uint8_t *A1, *A2;
  uint8_t current_seq;
  uint64_t num_lmer;

  void pr(string &s) {
    for (size_t i=0; i<s.size(); i++) {
      cout << (int)s[i];
    }
    cout << endl;
  }

  string to_str(string &s) {
    string t;
    for (size_t i=0; i<s.size(); ++i) {
      t.push_back(domain[s[i]]);
    }
    return t;
  }

  std::string to_str(uint64_t id, int l) {
    std::string s;
    do {
      int c = id & 0x3;
      if (c==0) { s.insert(0,1,'A'); }
      else if (c==1) { s.insert(0,1,'C'); }
      else if (c==2) { s.insert(0,1,'G'); }
      else if (c==3) { s.insert(0,1,'T'); }
      id >>= 2; l--;
    } while (id>0);
    while (l>0) {l--; s.insert(0,1,'A');}
    return s;
  }
  
  uint64_t to_int(std::string s) {
    uint64_t idx = 0;
    int len = s.size();
    for(int i=0;i<len;i++) {
          idx = idx<<2;
          idx += s[i];
    }
    return idx;
  }

  void gen_nbrhood(std::string& x, int errors) {
    int len = x.size();
    if (errors > 0) {
      for (int j=0; j<len; j++) {
        char t = x[j];
        x.erase(j,1);
        gen_nbrhood(x, errors-1);
        x.insert(j,1,t);
      }
      for (int j=0; j<len; j++) {
        for (size_t k=0; k<domain.size(); k++) {
          if (x[j] == (int)k) continue;
          char t = x[j];
          x[j] = k;
          gen_nbrhood(x, errors-1); 
          x[j] = t;
        }
      }
      for (int j=0; j<len+1; j++) {
        for (size_t k=0; k<domain.size(); k++) {
          x.insert(j,1,k);
          gen_nbrhood(x, errors-1); 
          x.erase(j,1);
        }
      }
    } else if (len == l) {
      //pr(x);
      uint64_t idx = to_int(x);
      //cout << to_str(x) << "," << idx << "," << to_str(idx, l) << endl;
      if (A1[idx] != current_seq) {
        A1[idx] = current_seq;
        A2[idx]++;
      }
    }
  }

  void gen_all(std::string& x) {
    int m = x.size();
    for (int q=-d; q<=+d; q++) {
      int k = l+q;
      for (int i=0; i<m-k+1; i++) {
        std::string kmer = x.substr(i, k);
        //pr(kmer);
        //cout << domain << endl;
        gen_nbrhood(kmer, d);
      }
    }
  }

  public:

  Ems1(const std::string &input, int l, int d, Params &params):
    MotifFinder("ems1", input, l, d, params), A1(0), A2(0), current_seq(0), num_lmer(0) {
      num_lmer = (uint64_t)1 << (2*l); 
      std::cout << "Memory required = " << 2*num_lmer << " bytes." << "\n";
      try {
        A1 = new uint8_t[num_lmer];
        A2 = new uint8_t[num_lmer];
      } catch (std::bad_alloc& exc) {
        std::cerr << "Could not allocate enough memory!!!" << std::endl;
        exit(0);
      }
      memset(A1, 0, num_lmer*sizeof(uint8_t));
      memset(A2, 0, num_lmer*sizeof(uint8_t));
    }

  void search() {
    size_t n = reads.size();
    for (size_t i=0; i<n; i++) {
      std::cout << "Processing sequence " << i << "...";
      std::cout.flush();
      current_seq = i+1;
      gen_all(reads[i]);
      std::cout << " Done." << std::endl;
    }

    for (uint64_t i=0; i<num_lmer; i++) {
      if (A2[i] == n) {
        motifs.push_back(to_str(i,l));
      }
    }
  }

};

#endif // __EMS1A_H_

