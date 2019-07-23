#ifndef __EMS2P_H
#define __EMS2P_H_

#include "ems.hpp"
#include<iostream>
#include<string>
#include<vector>
#include<map>
#include<thread>
#include<utility>
using namespace std;

#include <string.h>
#include <stdio.h>
#include <cassert>

#include "motif.hpp"
#include "motif_set.hpp"

string pr(const string& x, uchar c) {
   char rev_code[] = "ACGT";
   string s;
   for (size_t i=0; i<x.size(); i++) {
     if (x[i] == c) s.push_back('*');
     else s.push_back(rev_code[(int)x[i]]);
   }
   return s;
}

bool compare(Motif &m, Auxif& x, string &s, uchar c) {
for (size_t i=0; i<s.size(); i++) {
  if (s[i] == c) {
    if (m.get_2bits(2*i) || !x.get_bit(i)) return false;
  } else {
    if ((m.get_2bits(2*i) != s[i]) || x.get_bit(i)) return false;
  }
}
return true;
    
}

class NbdGenerator {
  private:
    size_t domain_size;
    vector<Motif> &curr_array;
    vector<Auxif> &curr_aux_array;
    std::string x;
    std::string mark;
    int l, d, k;
    bool leftmost, rightmost;
    size_t expanded_count;

    void gen_nbrhood3(int start, int alpha, int count) {
      int len = x.size();
      if (alpha > 0) {
        for (int j=start; j<len+1; j++) {
          if ((mark[j] == 'Y') || (mark[j] == 'Z')) continue;
          x.insert(j,1,domain_size);
          mark.insert(j,1,'N');
          gen_nbrhood3(j+1, alpha-1, 4*count);
          x.erase(j,1);
          mark.erase(j,1);
        }
      } else {
        curr_array.push_back(Motif(x));
        curr_aux_array.push_back(Auxif(x, domain_size));
        //cout << pr(x, domain_size) << endl;
        //cout << curr_array[curr_array.size()-1].get_kmer(l) << endl;
        //assert(curr_array.size() == curr_aux_array.size());
        //assert(compare(curr_array[curr_array.size()-1], curr_aux_array[curr_aux_array.size()-1], x, domain_size));
        expanded_count += count;
      }
    }

    void gen_nbrhood2(int start, bool stars_before, int sigma, int alpha, int count) {
      int len = x.size();
      if (sigma > 0) {
        for (int j=start; j<len; j++, stars_before=false) {
          if ((mark[j] == 'Y')) continue;
          char rr = mark[j+1];
          if (!rightmost && stars_before && (rr == 'Y')) continue;
          char t = x[j];
          char r = mark[j];
          x[j] = domain_size;
          mark[j] = 'Y';
          if (!leftmost && stars_before && (rr == 'N')) mark[j+1] = 'Z';
          gen_nbrhood2(stars_before, j+1, sigma-1, alpha, 4*count);
          x[j] = t;
          mark[j] = r;
          mark[j+1] = rr;
        }
      } else {
        gen_nbrhood3(0, alpha, count);
      }
    }

    void gen_nbrhood(int start, int end, int delta, int sigma, int alpha, int count) {
      if (delta > 0) {
        for (int j=start; j<end; j++) {
          char t = x[j];
          char rr = mark[j];
          x.erase(j,1);
          mark.erase(j,1);
          char r = mark[j];
          mark[j] = 'Y';
          gen_nbrhood(j, end-1, delta-1, sigma, alpha, count);
          x.insert(j,1,t);
          mark[j] = r;
          mark.insert(j,1, rr);
        }
      } else {
        gen_nbrhood2(0, true, sigma, alpha, count);
      }
    }

  public:

    NbdGenerator(size_t _ds, vector<Motif>& _ca, vector<Auxif>& _caa, std::string& _x, int _l, int _d, bool _lm, bool _rm) : 
      domain_size(_ds), curr_array(_ca), curr_aux_array(_caa), x(_x), l(_l), d(_d), k(x.size()), leftmost(_lm), rightmost(_rm) {
      }

    size_t generate()
    {
      expanded_count = 0;
      for (int delta = std::max(0,k-l); delta <= (d+k-l)/2; delta++) {
        int q = k-l;
        int alpha = delta - q;
        int sigma = d - alpha - delta;
        mark = std::string(k+1, 'N');
        int start;
        if (!leftmost) { mark[0] = 'Z'; }
        if (!rightmost) { mark[k] = 'Z'; start = 1;} else { start = 0;}
        gen_nbrhood(start, k, delta, sigma, alpha, 1);
      }
      return expanded_count;
    }
};



class Worker
{
private:
    size_t domain_size;
    const string &seq;
    const vector<Motif> &main_array;
    int l, d, m;
    vector<pair<int, int> > works;

    void radix_sort(vector<Motif> & vCmbItem, vector<Auxif> & vCmbItemAux, vector<Motif> & vCmbTmp, vector<Auxif> & vCmbTmpAux, uint64 compact_count)
    {
      uint64 bucket_count[4]; // use domain_size
      for(int i = 0; i < l; i++)
      {
        for (size_t bucket =0; bucket < domain_size; bucket++)
          bucket_count[bucket] = 0;
        for(uint64 j = 0; j < compact_count;  j++ )
        {
          uchar y = vCmbItemAux[j].get_bit(i);
          if (y) {
            for (size_t bucket =0; bucket < domain_size; bucket++)
              bucket_count[bucket]++;
          } else {
            uchar x = vCmbItem[j].get_2bits(2*i);
            bucket_count[x]++;
          }
        }
        for (size_t bucket=1; bucket < domain_size; bucket++)
          bucket_count[bucket] += bucket_count[bucket-1];
        uint64 j = compact_count;
        compact_count = bucket_count[domain_size-1];
        while (j--)
        {
          uchar y = vCmbItemAux[j].get_bit(i);
          if (y) {
            for (size_t bucket=0; bucket < domain_size; bucket++) {
              vCmbTmp[--bucket_count[bucket]].set(vCmbItem[j]);
              vCmbTmp[bucket_count[bucket]].set_2bits(bucket, 2*i);
              vCmbTmpAux[bucket_count[bucket]].set(vCmbItemAux[j]);
            }
          } else {
	      uchar x = vCmbItem[j].get_2bits(2*i);
              vCmbTmp[--bucket_count[x]].set(vCmbItem[j]);
              vCmbTmpAux[bucket_count[x]].set(vCmbItemAux[j]);
          }
        }
        swap(vCmbItem, vCmbTmp);
        swap(vCmbItemAux, vCmbTmpAux);
      }
    }

    void intersect(const vector<Motif> & str1,  vector<Motif> & str2, vector<Motif>& strTmp)
    {
      auto first1 = str1.begin(), last1 = str1.end();
      auto first2 = str2.begin(), last2 = str2.end();
      auto result = first2;
      int count = 0;
      while ((first1 != last1) && (first2 != last2)) {
        if (*first1 < *first2) {
          first1++;
        } else if (*first1 > *first2) {
          first2++;
        } else {
          *result = *first2;
          first1++;
          first2++;
          count = 1;
          break;
        }
      }
      while ((first1 != last1) && (first2 != last2)) {
        if (*first1 < *first2) {
          first1++;
        } else if (*first1 > *first2) {
          first2++;
        } else {
          if (*result != *first2)
            *(++result) = std::move(*first2);
          first1++;
          first2++;
        }
      }

      str2.erase(result+count, last2);
    }

public:
  Worker(size_t ds, const string &_seq, const vector<Motif> &_ma, int _l, int _d): 
    domain_size(ds), seq(_seq), main_array(_ma), l(_l), d(_d), m(seq.size())
  {
    for (int k=l-d; k<=l+d; k++) {
      for (int i=0; i<m-k+1; i++) {
        works.push_back(make_pair(i,k));
      }
    }
    uint64 wsize = works.size();
    for (uint64 i=0; i<wsize; i++) {
        uint64 j = i+rand()%(wsize-i);
        std::swap(works[i], works[j]);
    }
  }
  Worker(const Worker &w): domain_size(w.domain_size), seq(w.seq), main_array(w.main_array) {
  }
  ~Worker() { }
  uint64 get_load() { return works.size(); }
  void process(vector<Motif>& curr_array, vector<Auxif>& curr_aux_array, vector<Motif> &tmp_array,vector<Auxif> &tmp_aux_array, size_t start, size_t end)
  {
    curr_array.clear();
    curr_aux_array.clear();
    size_t expanded_count = 0;
    for (size_t j=start; j<end; j++) {
      int i = works[j].first;
      int k = works[j].second;
      string x = seq.substr(i, k);
      NbdGenerator generator(domain_size, curr_array, curr_aux_array, x, l, d, (i<=0), (i+k>=m));
      expanded_count += generator.generate();
    }
    uint64 compact_count = curr_array.size();
    curr_array.resize(expanded_count);
    curr_aux_array.resize(expanded_count);
    tmp_array.resize(expanded_count);
    tmp_aux_array.resize(expanded_count);
    radix_sort(curr_array, curr_aux_array, tmp_array, tmp_aux_array, compact_count);
    if (main_array.size()) {
      intersect(main_array, curr_array, tmp_array);
    } else {
      curr_array.erase( unique( curr_array.begin(), curr_array.end() ), curr_array.end() );
    }
  }
};

/*
void merge_motifs(const vector<Motif> &a, const vector<Motif> &b, vector<Motif> &c) {
    MotifSet motif_set;
    motif_set.init_add(a.data(), 0, a.size());
    motif_set.init_add(b.data(), 0, b.size());
    c.clear();
    Motif motif, next_motif;
    if (motif_set.get_min(motif)) {
        while (motif_set.get_min(next_motif)) {
            if (!(motif == next_motif)) {
                c.push_back(motif);
                motif = next_motif;
            }
        }
        c.push_back(motif);
    }
}
*/
void merge_motifs(const vector<Motif> &a, const vector<Motif> &b, vector<Motif> &c) {
    auto i = a.begin(), m = a.end();
    auto j = b.begin(), n = b.end();
    c.clear();
    // ASSUMPTION: a, b each has at least 1 elements
    if (*i < *j) {
      c.push_back(*i++);
    } else {
      if (*i == *j) ++i;
      c.push_back(*j++);
    }
    while ((i != m) && (j != n)) {
      if (*i < *j) {
         if (*i != c.back()) c.push_back(*i);
         ++i;
      } else {
	 if (*i == *j) ++i;
         if (*j != c.back()) c.push_back(*j);
         ++j;
      }
    }
    c.insert(end(c), i, m);
    c.insert(end(c), j, n);
}

class Ems2p : public MotifFinder<Ems2p> {
  vector<Motif> main_array;
  vector<Motif> *tmp_array;
  vector<Motif> *curr_array;
  vector<Auxif> *tmp_aux_array;
  vector<Auxif> *curr_aux_array;

  void gen_all(std::string& seq, int l, int d) {
	  vector<thread> work_threads;
    Worker worker(domain_size, seq, main_array, l, d);
    uint64 total_load = worker.get_load();
    uint64 ind_load = (total_load + params.num_threads - 1)/params.num_threads;
    uint64 start = 0, pos;
    for (int i=0; i<params.num_threads; ++i) {
      pos = min(start+ind_load, total_load);
      work_threads.push_back(thread(&Worker::process, &worker, std::ref(curr_array[i]), std::ref(curr_aux_array[i]), std::ref(tmp_array[i]), std::ref(tmp_aux_array[i]), start, pos));
      start = pos;
    }
    for (auto& p : work_threads)
      p.join();
#if 0
    main_array.clear();
    MotifSet motif_set;
    Motif motif, next_motif;
    for (int i=0; i<params.num_threads; ++i) {
        size_t end_pos = curr_array[i].size();
        if (end_pos) motif_set.init_add(curr_array[i].data(), 0, end_pos);
    }
    if (motif_set.get_min(motif)) {
        while (motif_set.get_min(next_motif)) {
            if (!(motif == next_motif)) {
                main_array.push_back(motif);
                motif = next_motif;
            }
        }
        main_array.push_back(motif);
    }
#else
    vector<int> todo;
    size_t total = 0;
    for (int i=0; i<params.num_threads; ++i) {
        size_t size = curr_array[i].size();
        if (size) {
          todo.push_back(i);
          total += size;
        }
    }
    int todo_size = todo.size();
    while (todo_size>2) {
      vector<thread> merge_threads;
      for (int i=todo_size-1; i > 0; i -= 2) {
          merge_threads.push_back(thread(merge_motifs, std::ref(curr_array[todo[i-1]]), std::ref(curr_array[todo[i]]), std::ref(tmp_array[todo[i-1]])));
      }
      for (auto &p : merge_threads)
          p.join();
      for (int i=todo_size-1; i > 0; i -= 2) {
          swap(curr_array[todo[i-1]], tmp_array[todo[i-1]]);
          todo.erase(begin(todo)+i);
      }
      todo_size = todo.size();
    } 
    if (todo_size == 1) {
        main_array.clear();
        main_array.insert(end(main_array), begin(curr_array[todo[0]]), end(curr_array[todo[0]]));
    } else {
        merge_motifs(curr_array[todo[0]], curr_array[todo[1]], main_array);
    } 
#endif
  }

  public:

  Ems2p(const std::string &input, int l, int d, Params &params):
    MotifFinder("ems2p", input, l, d, params) {
      tmp_array = new vector<Motif>[params.num_threads];
      tmp_aux_array = new vector<Auxif>[params.num_threads];
      curr_array = new vector<Motif>[params.num_threads];
      curr_aux_array = new vector<Auxif>[params.num_threads];
    }
  ~Ems2p() {
    delete [] tmp_array;
    delete [] tmp_aux_array;
    delete [] curr_array;
    delete [] curr_aux_array;
  }

  void search() {
    for (size_t i=0; i<reads.size(); i++) {
      std::cout << "Processing sequence " << i << "..." << std::endl;
      std::cout.flush();
      gen_all(reads[i], l, d);
//      cout << "candidates remaining " << main_array.size() << endl;
//      for (auto& motif : main_array) {
//          cout << pr(motif.get_kmer(l), domain_size) << endl;
//      }
      if (!main_array.size()) {
        motifs.clear();
        return;
      }
    }
    motifs.resize(main_array.size());
    for (size_t i=0; i<main_array.size(); ++i) {
      motifs[i] = main_array[i].get_kmer(l);
      for (int j=0; j<l; ++j) {
        motifs[i][j] = domain[motifs[i][j]];
      }
    }
    cout << "Num threads = " << params.num_threads << endl;
  }
};

#endif // __EMS2_H_

