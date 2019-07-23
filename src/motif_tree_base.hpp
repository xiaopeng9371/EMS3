#ifndef __MOTIF_TREE_BASE__
#define __MOTIF_TREE_BASE__

#include <iostream>
#include <string>
#include <bitset>
#include <vector>
#include "utils.h"

#define POOL_SIZE 1024

//char codes[] = {
// -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
// -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
// -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,100, -1, -1, -1, -1, -1,
// -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
// -1,  0, -1,  1, -1, -1, -1,  2, -1, -1, -1, -1, -1, -1, -1, -1,
// -1, -1, -1, -1,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
// -1,  0, -1,  1, -1, -1, -1,  2, -1, -1, -1, -1, -1, -1, -1, -1,
// -1, -1, -1, -1,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
// -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
// -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
// -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
// -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
// -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
// -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
// -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
// -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
//};

template <typename DerivedMotifTree, typename TreeNode>
class MotifTreeBase {
protected:
  string domain;
  size_t domain_size;
  TreeNode* root;
  uint32_t max_depth;
  std::string x;
  uint32_t mask;
  std::string name;
  Motifs &motifs;

public:
  void traverseRecursive(TreeNode* node, uint32_t depth) {}
  void insertRecursive(TreeNode* node, const std::string& motif, uint32_t depth) {}
  void intersectRecursive(TreeNode* to, const TreeNode* from, uint32_t depth) {}
  void setDomain(const string& _domain) { domain = _domain; domain_size= domain.size(); mask = (1 << domain_size) - 1; }
  
  MotifTreeBase(uint32_t _max_depth, Motifs& _motifs, const char* _name) : domain_size(domain.size()), root(TreeNode::allocateNode()), max_depth(_max_depth), x(max_depth, ' '), mask((1<<domain_size) - 1), name(_name), motifs(_motifs) {}
  ~MotifTreeBase() {
  }
  void traverse() {
    motifs.clear();
    static_cast<DerivedMotifTree*>(this)->traverseRecursive(root, 0);
  }
  void traverseOut() {
    std::cout << "######## Traverse " << name << " #############" << std::endl;
    static_cast<DerivedMotifTree*>(this)->traverse();
    std::cout << "count=" << motifs.size() << std::endl;
    for (uint32_t i=0; i<motifs.size(); i++) {
      std::cout << motifs[i] << std::endl;
    }
  }
  void insert(const std::string& motif) {
    static_cast<DerivedMotifTree*>(this)->insertRecursive(root, motif, 0);
  }
  void intersect(const DerivedMotifTree* other) {
    static_cast<DerivedMotifTree*>(this)->intersectRecursive(root, other->root, 0);
  }
};

#endif // __MOTIF_TREE_BASE__

