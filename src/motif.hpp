#ifndef __MOTIFP_H_
#define __MOTIFP_H_

#include "utils.h"

struct Motif {
  unsigned long long data;
  typedef unsigned long long data_t;
  Motif() {
    data = 0;
  }
  Motif(const string &x) {
    set(x);
  }
  inline void set(const Motif &x) {
    data = x.data;
  }
  inline void set(const std::string &x) {
    data = 0;
    int len = x.size();
    for (int i=0, j=len-1; i<len; ++i, --j) {
      set_2bits(x[i] & 3, 2*j);
    }
  }
  inline std::string get_kmer(int k) {
    std::string x(k, ' ');
    uint64_t t = data;
    for (int i=k-1; i>=0; --i) {
      x[i] = t & 3;
      t >>= 2;
    }
    return x;
  }
  inline void set_2bits(const uint64_t x, const uint32_t p) {
    data |= x << p;
  }
  inline unsigned char get_2bits(const uint32_t p) {
    return (data >> p) & 3;
  }
//  inline void insert_3bits(const uint64_t x, const uint32_t p) {
//    data = ((data & ~mask[p]) << 3) | (x << p) | (data & mask[p]);
//  }
//  inline void delete_3bits(const uint32_t p) {
//    data = ((data & ~mask[p+3]) >> 3) | (data & mask[p]);
//  }
//  inline void substitute_3bits(const uint64_t x, const uint32_t p) {
//    data = (data & (~mask[p+3] | mask[p])) | (x << p);
//  }
  inline void SHL_insert_2bits(const uint64_t x) {
    data = (data << 2) | x;
  }
  inline void SHR_insert_2bits(const uint64_t x, const uint32_t p) {
    data >>= 2;
    data |= x << p;
  }
  inline void clear(void) {
    data = 0;
  }
  inline bool operator<(const Motif &x) const { return data < x.data; }
  inline bool operator>(const Motif &x) const { return data > x.data; }
  inline bool operator==(const Motif &x) const { return data == x.data; }
  inline bool operator!=(const Motif &x) const { return data != x.data; }
};

struct Auxif {
  uint32 data;
  Auxif() {
    data = 0;
  }
  Auxif(const string &x, uchar c) {
    set(x, c);
  }
  inline void set(const Auxif &x) {
    data = x.data;
  }
  inline void set(const std::string &x, uchar c) {
    data = 0;
    int len = x.size();
    for (int i=0, j=len-1; i<len; ++i, --j) {
      if (x[i] == c) {
        set_bit(1, j);
      }
    }
  }
  inline void set_bit(const uint32_t x, const uint32_t p) {
    data |= x << p;
  }
  inline unsigned char get_bit(const uint32_t p) {
    return (data >> p) & 1;
  }
//  inline void insert_3bits(const uint64_t x, const uint32_t p) {
//    data = ((data & ~mask[p]) << 3) | (x << p) | (data & mask[p]);
//  }
//  inline void delete_3bits(const uint32_t p) {
//    data = ((data & ~mask[p+3]) >> 3) | (data & mask[p]);
//  }
//  inline void substitute_3bits(const uint64_t x, const uint32_t p) {
//    data = (data & (~mask[p+3] | mask[p])) | (x << p);
//  }
  inline void SHL_insert_bit(const uint64_t x) {
    data = (data << 1) | x;
  }
  inline void SHR_insert_bit(const uint64_t x, const uint32_t p) {
    data >>= 1;
    data |= x << p;
  }
  inline void clear(void) {
    data = 0;
  }
};


#endif //__MOTIFP_H_

