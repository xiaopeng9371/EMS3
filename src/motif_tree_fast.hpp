#ifndef __MOTIF_TREE_FAST__
#define __MOTIF_TREE_FAST__

#include "motif_tree_base.hpp"


struct TreeNodeFast {
  union {
    struct {
      TreeNodeFast* children[20];  // TODO make dynamic
      size_t sharing_info;
    };
    TreeNodeFast* next;
  };
  static TreeNodeFast* free_head;
  static TreeNodeFast* allocateNode() {
    if (!free_head) {
      free_head = new TreeNodeFast[POOL_SIZE];
      TreeNodeFast* head = free_head;
      for (uint32_t i = 1; i < POOL_SIZE; i++) {
        head->next = head+1;
        head++;
      }
      head->next = 0;
    }
    TreeNodeFast* head = free_head;
    free_head = head->next;
    memset((char*)head, 0, sizeof(TreeNodeFast));
    return head;
  }
  static void deallocateNode(TreeNodeFast* node) {
    node->next = free_head;
    free_head = node;
  }
  void print() {
    std::cout << "Node " << this << ", sharing info=" << sharing_info << std::endl;
//    for (size_t j=0; j<domain_size; j++) {
//      std::cout << " child " << j << " ptr=" << children[j] << std::endl;
//    }
  }
};
TreeNodeFast* TreeNodeFast::free_head = 0;


class MotifTreeFast : public MotifTreeBase<MotifTreeFast, TreeNodeFast> {
  void print() {
    std::cout << "root=" << root << ", name=" << name << std::endl;
    print_recursive(root);
  }
  void print_recursive(const TreeNodeFast *t) {
    static int depth = 0;
    std::string tmp(5*depth, '+');
    depth++;
    std::cout << tmp << " " << (t!=0) << ", " << std::bitset<4>(t->sharing_info).to_string() 
      << ", " << (t->children[3] ? "x" : "0")
      << ", " << (t->children[2] ? "x" : "0")
      << ", " << (t->children[1] ? "x" : "0")
      << ", " << (t->children[0] ? "x" : "0")
      << std::endl;
    size_t current_info = 0;
    size_t current_mask = 1;
    for (size_t j=0; j<domain_size; j++) {
      if (!(current_info & current_mask)) {
        if (t->children[j]) {
          print_recursive(t->children[j]);
          current_info |= t->children[j]->sharing_info;
        }
      }
      current_mask <<= 1;
    }
    depth--;
  }
  void deleteNode(TreeNodeFast* node) {
    assert(node);
    size_t current_info = mask;
    size_t current_mask = 1;
    for (size_t j=0; j<domain_size; j++) {
      if ((current_info & current_mask) && (node->children[j]) ) {
        current_info &= ~(node->children[j]->sharing_info);
        deleteNode(node->children[j]);
      }
      current_mask <<= 1;
    }
    TreeNodeFast::deallocateNode(node);
  }
  bool emptyNode(const TreeNodeFast* node) {
    for (size_t j=0; j<domain_size; j++) {
      if (node->children[j]) return false;
    }
    return true;
  }
  void copy(TreeNodeFast* to, const TreeNodeFast* from);
  void copyUnion(TreeNodeFast* to, const TreeNodeFast* from, const std::string& motif, size_t depth);
  void copyIntersection(TreeNodeFast* to, const TreeNodeFast* from, const TreeNodeFast* other, size_t depth);
  void insertCommonRecursive(TreeNodeFast* to, const TreeNodeFast* from, const std::string& motif, size_t depth);
  bool hasIntersect(TreeNodeFast* node, const std::string& motif, size_t depth);

  public:
  MotifTreeFast(uint32_t _max_depth, Motifs& _motifs, const char* _name) : MotifTreeBase(_max_depth, _motifs, _name) {}
  ~MotifTreeFast() {
    deleteNode(root);
  }
  void traverseRecursive(TreeNodeFast* node, size_t depth);
  void insertRecursive(TreeNodeFast* node, const std::string& motif, size_t depth);
  void insertRecursiveNew(TreeNodeFast* node, const std::string& motif, size_t depth);
  void intersectRecursive(TreeNodeFast* to, const TreeNodeFast* from, size_t depth);
};


void MotifTreeFast::traverseRecursive(TreeNodeFast* node, size_t depth) {
  //std::string tmp(6*depth, '#');
  //std::cout << tmp << depth << "," << std::bitset<4>(node->sharing_info).to_string() << std::endl;
  if (depth == max_depth) {
    motifs.push_back(x);
    return;
  }
  for (size_t i=0; i<domain_size; ++i) {
    x[depth] = domain[i];
    if (node->children[i]) {
      traverseRecursive(node->children[i], depth+1);
    }
  }
}

void MotifTreeFast::copy(TreeNodeFast* to, const TreeNodeFast* from) {
  assert(from);
  size_t current_info = 0;
  size_t current_mask = 1;
  for (size_t j=0; j<domain_size; j++) {
    if (!(current_info & current_mask)) {
      if (from->children[j]) {
        TreeNodeFast *new_node  = TreeNodeFast::allocateNode();
        new_node->sharing_info = from->children[j]->sharing_info;
        copy(new_node, from->children[j]);
        for (size_t k=0, new_mask=1; k<domain_size; k++, new_mask <<= 1) {
          if (new_node->sharing_info & new_mask) {
            to->children[k] = new_node;
          }
        }
        current_info |= from->children[j]->sharing_info;
      }
    }
    current_mask <<= 1;
  }
}

void MotifTreeFast::copyUnion(TreeNodeFast* to, const TreeNodeFast* from, const std::string& motif, size_t depth) {
  //std::cout << "Copy Union: motif=" << motif << ", depth=" << depth << std::endl;
  if (depth == max_depth) return;
  size_t child = motif[depth];
  if (child >= domain_size) {
    size_t remaining = mask;
    size_t current_mask = 1;
    for (size_t j=0; j<domain_size; j++, current_mask <<= 1) {
      if ((remaining & current_mask) && from->children[j]) {
        TreeNodeFast *new_node  = TreeNodeFast::allocateNode();
        new_node->sharing_info = from->children[j]->sharing_info;
        for (size_t k=0, new_mask=1; k<domain_size; k++, new_mask <<= 1) {
          if (new_node->sharing_info & new_mask) to->children[k] = new_node;
        }
        copyUnion(new_node, from->children[j], motif, depth+1);
        remaining &= ~from->children[j]->sharing_info;
      }
    }
    if (remaining) {
      TreeNodeFast *new_node  = TreeNodeFast::allocateNode();
      new_node->sharing_info = remaining;
      for (size_t k=0, new_mask=1; k<domain_size; k++, new_mask <<= 1) {
        if (remaining & new_mask) to->children[k] = new_node;
      }
      insertRecursive(new_node, motif, depth+1);
    }
  } else {
    size_t child_mask = (1<<child);
    size_t current_info = mask & (~child_mask);
    size_t current_mask = 1;
    for (size_t j=0; j<domain_size; j++, current_mask <<= 1) {
      if ((current_info & current_mask) && from->children[j]) {
        TreeNodeFast *new_node  = TreeNodeFast::allocateNode();
        new_node->sharing_info = from->children[j]->sharing_info & current_info;
        for (size_t k=0, new_mask=1; k<domain_size; k++, new_mask <<= 1) {
          if (new_node->sharing_info & new_mask) to->children[k] = new_node;
        }
        copy(new_node, from->children[j]);
        current_info &= ~(new_node->sharing_info);
      }
    }
    to->children[child] = TreeNodeFast::allocateNode();
    to->children[child]->sharing_info = child_mask;
    if (from->children[child]) {
      copyUnion(to->children[child], from->children[child], motif, depth+1);
    } else {
      insertRecursive(to->children[child], motif, depth+1);
    }
  }
}

void MotifTreeFast::copyIntersection(TreeNodeFast* to, const TreeNodeFast* from, const TreeNodeFast* other, size_t depth) {
  size_t current_info = mask;
  size_t current_mask = 1;
  for (size_t j=0; j<domain_size; j++) {
    if ((current_info & current_mask) && from->children[j]) {
      size_t from_info = from->children[j]->sharing_info;
      size_t other_info = (other && other->children[j]) ? other->children[j]->sharing_info : 0;
      size_t common = other_info & from_info;
      if (common) {
        TreeNodeFast* new_node = TreeNodeFast::allocateNode();
        new_node->sharing_info = common;
        copyIntersection(new_node, from->children[j], other->children[j], depth+1);
        if ((depth < max_depth-1) && emptyNode(new_node)) {
          TreeNodeFast::deallocateNode(new_node);
          new_node = 0;
        }
        for (size_t k=0, new_mask=1; k<domain_size; k++, new_mask <<= 1) {
          if (common & new_mask) to->children[k] = new_node;
        }
        current_info &= ~common;
      }
    }
    current_mask <<= 1;
  }
}

void MotifTreeFast::insertRecursive(TreeNodeFast* node, const std::string& motif, size_t depth) {
  //std::cout << "INSERT RECURSIVE: motif=" << motif << ", depth=" << depth << std::endl;
  if (depth >= max_depth) {
    return;
  }
  size_t child = motif[depth];
  if (child >= domain_size) {
    size_t current_info = 0;
    size_t current_mask = 1;
    for (size_t j=0; j<domain_size; j++) {
      if (!(current_info & current_mask)) {
        if (node->children[j]) {
          insertRecursive(node->children[j], motif, depth+1);
          current_info |= node->children[j]->sharing_info;
        }
      }
      current_mask <<= 1;
    }
    // std::cout << "current_info = " << std::bitset<4>(current_info).to_string() << std::endl;
    // std::cout << "current_mask = " << std::bitset<4>(current_mask).to_string() << std::endl;
    // std::cout << "        mask = " << std::bitset<4>(mask).to_string() << std::endl;
    size_t remaining = (~current_info) & mask;
    // std::cout << "remaining = " << std::bitset<4>(remaining).to_string() << std::endl;
    if (!remaining) return;
    TreeNodeFast* new_node = TreeNodeFast::allocateNode();
    new_node->sharing_info = remaining;
    current_mask = 1;
    for (size_t j=0; j<domain_size; j++) {
      if (remaining & current_mask) {
        node->children[j] = new_node;
      } 
      current_mask <<= 1;
    }
    insertRecursive(new_node, motif, depth+1);
  } else {
    size_t child_mask = (1<<child);
    if (!node->children[child]) {
      node->children[child] = TreeNodeFast::allocateNode();
      node->children[child]->sharing_info = child_mask;
      insertRecursive(node->children[child], motif, depth+1);
    } else {
      TreeNodeFast* old_node = node->children[child];
      if (node->children[child]->sharing_info & (~child_mask)) {
        node->children[child] = TreeNodeFast::allocateNode();
        old_node->sharing_info &= ~child_mask;
        node->children[child]->sharing_info = child_mask;
        copyUnion(node->children[child], old_node, motif, depth+1);
      } else {
        insertRecursive(node->children[child], motif, depth+1);
      }
    }
  }
}

void MotifTreeFast::insertCommonRecursive(TreeNodeFast* to, const TreeNodeFast* from, const std::string& motif, size_t depth) {
  //std::cout << "motif=" << motif << ", depth=" << depth << std::endl;
  if (depth >= max_depth) {
    return;
  }
  std::string tmp(depth, '#');
  size_t min_j, max_j;
  size_t child = motif[depth];
  if (child >= domain_size) {
    min_j = 0; max_j = domain_size-1;
  } else {
    min_j = max_j = child;
  }
  size_t current_mask = (1 << min_j);
  size_t current_info = mask;
  for (size_t j=min_j; j<=max_j; j++) {
    //std::cout << tmp << " j=" << j << ", min_j=" << min_j << ", max_j=" << max_j << std::endl;
    if ((current_info & current_mask) && from->children[j]) {
      size_t from_info = (child >= domain_size) ? from->children[j]->sharing_info : current_mask;
      size_t to_info = (to && to->children[j]) ? to->children[j]->sharing_info : 0;
      size_t common = to_info & from_info;
      //std::cout << "ch=" << child << ", frm=" << from_info << ", to=" << to_info << ", common=" << common << std::endl;
      if (common) {
        if (to_info != common) {
          TreeNodeFast* old_node = to->children[j];
          TreeNodeFast* new_node = TreeNodeFast::allocateNode();
          old_node->sharing_info = to_info & ~common;
          new_node->sharing_info = common;
          for (size_t k=0, new_mask=1; k<domain_size; k++, new_mask <<= 1) {
            if (common & new_mask) to->children[k] = new_node;
          }
          copy(new_node, old_node);
        }
        current_info &= ~common;
      } else {
        TreeNodeFast* new_node = TreeNodeFast::allocateNode();
        new_node->sharing_info = from_info;
        for (size_t k=0, new_mask=1; k<domain_size; k++, new_mask <<= 1) {
          if (from_info & new_mask) {
            to->children[k] = new_node;
          }
        }
        current_info &= ~from_info;
      }
      insertCommonRecursive(to->children[j], from->children[j], motif, depth+1);
    }
    current_mask <<= 1;
  }
}

void MotifTreeFast::intersectRecursive(TreeNodeFast* to, const TreeNodeFast* from, size_t depth) {
  size_t current_info = mask;
  size_t current_mask = 1;
  for (size_t j=0; j<domain_size; j++) {
    if ((current_info & current_mask) && to->children[j]) {
      size_t to_info = to->children[j]->sharing_info;
      size_t from_info = (from && from->children[j]) ? from->children[j]->sharing_info : 0;
      size_t common = to_info & from_info;
      if (common) {
        if (to_info != common) {
          TreeNodeFast* old_node = to->children[j];
          old_node->sharing_info = to_info & ~common;
          TreeNodeFast* new_node = TreeNodeFast::allocateNode();
          new_node->sharing_info = common;
          copyIntersection(new_node, old_node, from->children[j], depth+1);
          if ((depth < max_depth-1) && emptyNode(new_node)) {
            TreeNodeFast::deallocateNode(new_node);
            new_node = 0;
          }
          for (size_t k=0, new_mask=1; k<domain_size; k++, new_mask <<= 1) {
            if (common & new_mask) to->children[k] = new_node;
          }
        } else {
          intersectRecursive(to->children[j], from->children[j], depth+1);
        }
        current_info &= ~common;
      } else {
        to->children[j]->sharing_info &= ~current_mask;
        if (!to->children[j]->sharing_info) {
          deleteNode(to->children[j]);
        }
        to->children[j] = 0;
      }
    }
    current_mask <<= 1;
  }
}

bool MotifTreeFast::hasIntersect(TreeNodeFast* node, const std::string& motif, size_t depth) {
  if (node) {
    size_t child = motif[depth];
    if (child >= domain_size) {
      size_t current_info = 0;
      size_t current_mask = 1;
      for (size_t j=0; j<domain_size; j++) {
        if (!(current_info & current_mask)) {
          if (node->children[j]) {
            if (hasIntersect(node->children[j], motif, depth+1)) {
              return true;
            }
            current_info |= node->children[j]->sharing_info;
          }
        }
        current_mask <<= 1;
      }
    } else {
      if (!node->children[child]) return false;
      return hasIntersect(node->children[child], motif, depth+1);
    }
  }
  return true;
}


#endif // __MOTIF_TREE_FAST__

