#include <iostream>
#include <string>
#include <bitset>
#include <vector>
#include "utils.h"

const std::string alphabet="*ACGT";
const size_t alphabet_size = alphabet.size();


inline size_t getTrieCode(char base) {
   switch (base) {
     case '*':
       return 0;
     case 'A':
       return 1;
     case 'C':
       return 2;
     case 'G':
       return 3;
     case 'T':
       return 4;
   }
   return alphabet_size+1;
}

union TrieNode {
  struct {
  TrieNode* children[5];  // TODO make dynamic
  };
  TrieNode* next;
};

void print_node(TrieNode *t) {
  for (size_t j=0; j<alphabet_size; j++) {
    std::cout << " child " << j << " ptr=" << t->children[j] << std::endl;
  }
}

#define POOL_SIZE 1024
TrieNode* free_trie_head = 0;

TrieNode* allocateTrieNode() {
  if (!free_trie_head) {
    free_trie_head = new TrieNode[POOL_SIZE];
    TrieNode* head = free_trie_head;
    for (size_t i = 1; i < POOL_SIZE; i++) {
      head->next = head+1;
      head++;
    }
    head->next = 0;
  }
  TrieNode* head = free_trie_head;
  free_trie_head = head->next;
  memset((char*)head, 0, sizeof(TrieNode));
  return head;
}

void deallocateTrieNode(TrieNode* node) {
  node->next = free_trie_head;
  free_trie_head = node;
}

void valid(const TrieNode *t) {
  assert(t);
  assert(t != t->children[0]);
  assert(t != t->children[1]);
  assert(t != t->children[2]);
  assert(t != t->children[3]);
  assert(t != t->children[4]);
}

void valid_recursive(const TrieNode *t) {
  valid(t);
  for (size_t j=0; j<alphabet_size; j++) {
      if (t->children[j]) {
        valid_recursive(t->children[j]);
      }
  }
}



class MotifTrie {
  TrieNode* root;
  size_t max_depth;
  std::string x;
  std::string name;
  Motifs &motifs;
  size_t count;
  void traverseRecursive(TrieNode* node, size_t depth) {
    //std::string tmp(6*depth, '#');
    //std::cout << tmp << depth << std::endl;
    if (depth == max_depth) {
      //count++;
      //std::cout << x << std::endl;
      motifs.push_back(x);
      return;
    }
    for (size_t i=1; i<alphabet_size; ++i) {
      x[depth] = alphabet[i];
      if (node->children[0]) {
        traverseRecursive(node->children[0], depth+1);
      }
      if (node->children[i]) {
        traverseRecursive(node->children[i], depth+1);
      }
    }
  }
  void print_recursive(const TrieNode *t) {
    static size_t depth = 0;
    depth++;
    std::string tmp(5*depth, '+');
    //  std::cout << tmp << " " << (t!=0) 
    //                               << ", " << (t->children[4] ? "x" : "0")
    //                               << ", " << (t->children[3] ? "x" : "0")
    //                               << ", " << (t->children[2] ? "x" : "0")
    //                               << ", " << (t->children[1] ? "x" : "0")
    //                               << ", " << (t->children[0] ? "x" : "0")
    //                               << std::endl;
    if (depth==max_depth+1) {
      std::cout << x << std::endl;
      return;
    }
    for (size_t j=0; j<alphabet_size; j++) {
      x[depth-1] = alphabet[j];
      if (t->children[j]) {
        print_recursive(t->children[j]);
      }
    }
    depth--;
  }
  void insertRecursive(TrieNode* node, const std::string& motif, size_t depth) {
    //std::cout << "motif=" << motif << ", depth=" << depth << std::endl;
    if (depth >= max_depth) {
      return;
    }
    size_t child = getTrieCode(motif[depth]);
    if (!node->children[child]) {
      node->children[child] = allocateTrieNode();
    }
    insertRecursive(node->children[child], motif, depth+1);
  }
  void createChildIfNotHave(TrieNode* node, size_t j) {
    if (!node->children[j]) {
      node->children[j] = allocateTrieNode();
    }
  }
  void cleanUpUnused(TrieNode* node, size_t j) {
    TrieNode* t = node->children[j];
    for (size_t k=0; k<alphabet_size; k++) {
      if (t->children[k]) return;
    }
    deallocateTrieNode(t);
    node->children[j] = 0;
  }
  void deleteNode(TrieNode* t) {
    for (size_t k=0; k<alphabet_size; k++) {
      if (t->children[k]) deleteNode(t->children[k]);
    }
    deallocateTrieNode(t);
  }

  void insertCommonRecursive(TrieNode* to, const TrieNode* from, const std::string& motif, size_t depth) {
    if (depth >= max_depth) {
      return;
    }
    //std::string tmp(depth, '#');
    //std::cout << "motif=" << motif << ", depth=" << depth << ", char=" << motif[depth] << std::endl;
    size_t child = getTrieCode(motif[depth]);

    if (0 == child) {
      for (size_t j=0; j<alphabet_size; j++) {
        if (from->children[j]) {
          //std::cout << "id=" << child << ", going to recurse with j=" << j << std::endl;
          createChildIfNotHave(to, j);
          insertCommonRecursive(to->children[j], from->children[j], motif, depth+1);
          if (depth < max_depth-1) cleanUpUnused(to, j);
        }
      }
    } else {
        if (from->children[0]) {
          //std::cout << "id=" << child  << ", going to recurse with 0" << std::endl;
          createChildIfNotHave(to, child);
          insertCommonRecursive(to->children[child], from->children[0], motif, depth+1);
          if (depth < max_depth-1) cleanUpUnused(to, child);
        }
        if (from->children[child]) {
          //std::cout << "id=" << child  << ", going to recurse with j=" << child  << std::endl;
          createChildIfNotHave(to, child);
          insertCommonRecursive(to->children[child], from->children[child], motif, depth+1);
          if (depth < max_depth-1) cleanUpUnused(to, child);
        }
    }
  }
  void intersectRecursive(TrieNode* to, const TrieNode* x, const TrieNode* y, size_t depth) {
    if (x->children[0]) {
      for (size_t j=0; j<alphabet_size; j++) {
        if (y->children[j]) {
          //std::cout << "id=" << child << ", going to recurse with j=" << j << std::endl;
          createChildIfNotHave(to, j);
          intersectRecursive(to->children[j], x->children[0], y->children[j], depth+1);
          if (depth < max_depth-1) cleanUpUnused(to, j);
        }
      }
    }
    for (size_t j=1; j<alphabet_size; j++) {
      if (x->children[j]) {
        if (y->children[0]) {
          createChildIfNotHave(to, j);
          intersectRecursive(to->children[j], x->children[j], y->children[0], depth+1);
          if (depth < max_depth-1) cleanUpUnused(to, j);
        }
        if (y->children[j]) {
          createChildIfNotHave(to, j);
          intersectRecursive(to->children[j], x->children[j], y->children[j], depth+1);
          if (depth < max_depth-1) cleanUpUnused(to, j);
        }
      }
    }
  }
  public:
  MotifTrie(size_t _max_depth, Motifs& _motifs, const char* _name) : root(allocateTrieNode()), max_depth(_max_depth), x(max_depth, ' '), name(_name), motifs(_motifs) {}
  MotifTrie(const MotifTrie* a, const MotifTrie* b, size_t _max_depth, Motifs& _motifs, const char* _name="") : root(allocateTrieNode()), max_depth(_max_depth), x(max_depth, ' '), name(_name), motifs(_motifs) {
    intersectRecursive(root, a->root, b->root, 0);
  }
  ~MotifTrie() {
    //std::cout << "In deallocator " << std::endl;
    deleteNode(root);
  }
  void traverse() {
//    std::cout << "######## Traverse " << name << " #############" << std::endl;
//    count = 0;
    motifs.clear();
    traverseRecursive(root, 0);
    sort( motifs.begin(), motifs.end() );
    motifs.erase( unique( motifs.begin(), motifs.end() ), motifs.end() );
//    std::cout << "count=" << motifs.size() << std::endl;
//    for (size_t i=0; i<motifs.size(); i++) {
//      std::cout << motifs[i] << std::endl;
//    }
  }
  void insert(const std::string& motif) {
    insertRecursive(root, motif, 0);
  }
  void insert(const MotifTrie* other, const std::string& motif) {
    insertCommonRecursive(root, other->root, motif, 0);
  }
  void print() {
    std::cout << "root=" << root << ", name=" << name << std::endl;
    print_recursive(root);
  }
};

#if MAIN_TREE

int main() {
  Motifs m1, m2;
  MotifTrie mt(4, m1, "tree1");
  MotifTrie mt2(4, m2, "tree2");
  mt.insert("**GA");
  mt.insert("*G*A");
  mt.insert("*GA*");
  mt.traverse();
  mt2.insert(&mt, "*GGA");
  mt2.insert(&mt, "A*GT");
  mt2.print();
  mt2.traverse();
}

#endif 

