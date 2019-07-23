#ifndef __EMS_HPP__
#define __EMS_HPP__

#include <sys/time.h>
#include <sys/resource.h>
#include "utils.h"
#include "omp.h"
#include <iostream>
#include <fstream>

struct Params {
int l;
int d;
int num_threads;
};

template <typename DerivedMotifFinder>
class MotifFinder {

protected:
  std::string name;
  std::string domain;
  size_t domain_size;
  std::string input;
  Reads reads;
  int l, d;
  Params &params;
  Motifs motifs;

public:

  MotifFinder(const std::string& _name, const Reads &_reads, int _l, int _d, Params &_params) : name(_name), reads(_reads), l(_l), d(_d), params(_params) { }
  MotifFinder(const std::string& _name, const std::string &_input, int _l, int _d, Params &_params) : name(_name), input(_input), l(_l), d(_d), params(_params) {
    read_file(input.c_str(), reads);
    domain = getAlphabet(reads);
    std::cout << "Domain is " << domain << std::endl;
    domain_size = domain.size();
    std::cout << "Domain size is " << domain_size << std::endl;
    encodeStrings(reads, domain);
  }
  ~MotifFinder() { }

  void search() { std::cout << "Please implement search function in the derived class." << std::endl; }

  Motifs& searchGetMotifs() {
    static_cast<DerivedMotifFinder*>(this)->search();
    return motifs;
  }

  void searchWriteMotifs(Params &params) {
    std::string output = get_out_file(input, l, d, name);
    //std::cout << "l      = " << l << ", d = " << d << std::endl;
    //std::cout << "input  = " << input << std::endl;
    //std::cout << "output = " << output << std::endl;
   // clock_t begin=clock();
 //  omp_set_num_threads(4);
    //ofstream myFile;
    //myFile.open("emsTimeMemory",ios::app);
    double begin = omp_get_wtime();
    static_cast<DerivedMotifFinder*>(this)->search();
   // clock_t end=clock();
    double end = omp_get_wtime();
    //double elapsed = diffclock(end,begin);
    double elapsed = end-begin;
    struct rusage usage;
    getrusage( RUSAGE_SELF, &usage );
    std::ofstream out(output);
    std::cout << "\r" << name << ": (" << l << "," << d << ") Edited Motifs found (in "<< elapsed << " sec, using " << (size_t)usage.ru_maxrss << " KB): " << motifs.size() << "       " << std::endl;
    //myFile << name << ": (" << l << "," << d << ") Edited Motifs found using " << params.num_threads << " threads:(in "<< elapsed << " sec, using " << (size_t)usage.ru_maxrss << " KB): " << motifs.size() << "       " << std::endl;
    //myFile.close();
    for (size_t i=0; i<motifs.size(); ++i) {
      out << motifs[i] << std::endl;
    }
    out.close();
  }
};

#endif // __EMS_HPP__



