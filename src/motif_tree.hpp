#ifndef __MOTIF_TREE__
#define __MOTIF_TREE__

#include "motif_tree_base.hpp"

struct TreeNodeSlow {
  union {
    struct {
      TreeNodeSlow* left_child;
      TreeNodeSlow* right_sibling;
      uint32_t sharing_info;
    };
    TreeNodeSlow* next;
  };
  static TreeNodeSlow* free_head;
  static TreeNodeSlow* allocateNode() {
    if (!free_head) {
      free_head = new TreeNodeSlow[POOL_SIZE];
      TreeNodeSlow* head = free_head;
      for (uint32_t i = 1; i < POOL_SIZE; i++) {
        head->next = head+1;
        head++;
      }
      head->next = 0;
    }
    TreeNodeSlow* head = free_head;
    free_head = head->next;
    memset((char*)head, 0, sizeof(TreeNodeSlow));
    return head;
  }
  static void deallocateNode(TreeNodeSlow* node) {
    node->next = free_head;
    free_head = node;
  }
  void print() {
    std::cout << "Node " << this << ", sharing info=" << sharing_info << std::endl;
    std::cout << " left child " << left_child << std::endl;
    std::cout << " right sibling " << right_sibling << std::endl;
  }
};

TreeNodeSlow* TreeNodeSlow::free_head = 0;


class MotifTreeSlow : public MotifTreeBase<MotifTreeSlow, TreeNodeSlow> {
  void print() {
    std::cout << "root=" << root << ", name=" << name << std::endl;
    print_recursive(root);
  }
  void print_recursive(const TreeNodeSlow *t) {
    static int depth = 0;
    std::string tmp(5*depth, '+');
    depth++;
    std::cout << tmp << " " << (t!=0) << ", " << std::bitset<4>(t->sharing_info).to_string() << std::endl;
    for (TreeNodeSlow* child = t->left_child; child; child = child->right_sibling) {
      print_recursive(child);
    }
    depth--;
  }
  void copy(TreeNodeSlow* to, const TreeNodeSlow* from);
  void copyUnion(TreeNodeSlow* to, const TreeNodeSlow* from, const std::string& motif, uint32_t depth);
  void copyIntersection(TreeNodeSlow* to, const TreeNodeSlow* from, const TreeNodeSlow* other, uint32_t depth);
  void deleteNode(TreeNodeSlow* node);

public:
  MotifTreeSlow(uint32_t _max_depth, Motifs& _motifs, const char* _name) : MotifTreeBase(_max_depth, _motifs, _name) {}
  ~MotifTreeSlow() {
    deleteNode(root);
  }
  void traverseRecursive(TreeNodeSlow* node, uint32_t depth);
  void insertRecursive(TreeNodeSlow* node, const std::string& motif, uint32_t depth);
  void insertRecursiveNew(TreeNodeSlow* node, const std::string& motif, uint32_t depth);
  void intersectRecursive(TreeNodeSlow* to, const TreeNodeSlow* from, uint32_t depth);
};

void MotifTreeSlow::traverseRecursive(TreeNodeSlow* node, uint32_t depth) {
  //std::string tmp(6*depth, '#');
  //std::cout << tmp << depth << "," << std::bitset<4>(node->sharing_info).to_string() << std::endl;
  if (depth == max_depth) {
    motifs.push_back(x);
    return;
  }
  TreeNodeSlow* children[4] = {0}; // TODO make domain_size
  for (TreeNodeSlow *child = node->left_child; child; child = child->right_sibling) {
    uint32_t info = child->sharing_info;
    for (uint32_t j=0, mask=1; j<domain_size; j++, mask <<= 1) {
      if (info & mask) children[j] = child;
    }
  }
  for (uint32_t i=0; i<domain_size; ++i) {
    x[depth] = domain[i];
    if (children[i]) {
      traverseRecursive(children[i], depth+1);
    }
  }
}

void MotifTreeSlow::copy(TreeNodeSlow* to, const TreeNodeSlow* from) {
  assert(from);
  TreeNodeSlow *from_child = from->left_child, *to_child = 0;
  while (from_child) {
    to_child = TreeNodeSlow::allocateNode();
    to_child->sharing_info = from_child->sharing_info;
    to_child->right_sibling = to->left_child;
    to->left_child = to_child;
    copy(to_child, from_child);
    from_child = from_child->right_sibling;
  }
}

void MotifTreeSlow::copyUnion(TreeNodeSlow* to, const TreeNodeSlow* from, const std::string& motif, uint32_t depth) {
  //std::cout << "Copy Union: motif=" << motif << ", depth=" << depth << std::endl;
  if (depth == max_depth) return;
  uint32_t ch = motif[depth];
  //std::cout << "ch=" << ch << std::endl;
  uint32_t remaining = (ch >= domain_size) ? mask : (1 << ch);
  TreeNodeSlow *from_child = from->left_child, *to_child = 0;
  while (from_child) {
    to_child = TreeNodeSlow::allocateNode();
    to_child->sharing_info = from_child->sharing_info;
    to_child->right_sibling = to->left_child;
    to->left_child = to_child;
    uint32_t common = from_child->sharing_info & remaining;
    if (common) {
      if (common == from_child->sharing_info) {
        copyUnion(to_child, from_child, motif, depth+1);
      } else {
        to_child->sharing_info &= ~common;
        copy(to_child, from_child);
        to_child = TreeNodeSlow::allocateNode();
        to_child->sharing_info = common;
        to_child->right_sibling = to->left_child;
        to->left_child = to_child;
        copyUnion(to_child, from_child, motif, depth+1);
      }
    } else {
      copy(to_child, from_child);
    }
    remaining &= ~from_child->sharing_info;
    from_child = from_child->right_sibling;
  }
  if (remaining) {
    to_child = TreeNodeSlow::allocateNode();
    to_child->sharing_info = remaining;
    to_child->right_sibling = to->left_child;
    to->left_child = to_child;
    insertRecursive(to_child, motif, depth+1);
  }
}

void MotifTreeSlow::copyIntersection(TreeNodeSlow* to, const TreeNodeSlow* from, const TreeNodeSlow* other, uint32_t depth) {
  //std::cout << "Copy intersection : to=" << to << ", from=" << from << ", depth=" << depth << std::endl;
  TreeNodeSlow *from_child = from->left_child;
  to->left_child = 0;
  while (from_child) {
    for (TreeNodeSlow *other_child = other->left_child; other_child; other_child = other_child->right_sibling) {
      uint32_t common = from_child->sharing_info & other_child->sharing_info;
      if (common) {
        TreeNodeSlow *new_node = TreeNodeSlow::allocateNode();
        copyIntersection(new_node, from_child, other_child, depth+1);
        if ((depth < max_depth-1) && (!new_node->left_child)) {
          TreeNodeSlow::deallocateNode(new_node);
        } else {
          new_node->sharing_info = common;
          new_node->right_sibling = to->left_child;
          to->left_child = new_node;
        }
      }
    }
    from_child = from_child->right_sibling;
  }
}

void MotifTreeSlow::deleteNode(TreeNodeSlow* node) {
  TreeNodeSlow *child = node->left_child;
  while (child) {
    TreeNodeSlow* tmp = child;
    child = child->right_sibling;
    deleteNode(tmp);
  }
  TreeNodeSlow::deallocateNode(node);
}

void MotifTreeSlow::insertRecursive(TreeNodeSlow* node, const std::string& motif, uint32_t depth) {
  //std::cout << "INSERT RECURSIVE: motif=" << motif << ", depth=" << depth << std::endl;
  if (depth >= max_depth) {
    return;
  }
  uint32_t ch = motif[depth];
  //std::cout << "ch=" << ch << std::endl;
  if (ch >= domain_size) {
    uint32_t info = 0;
    TreeNodeSlow *child = node->left_child;
    while (child) {
      info |= child->sharing_info;
      insertRecursive(child, motif, depth+1);
      child = child->right_sibling;
    }
    uint32_t remaining = (~info) & mask;
    // std::cout << "remaining = " << std::bitset<4>(remaining).to_string() << std::endl;
    if (!remaining) return;
    TreeNodeSlow* new_node = TreeNodeSlow::allocateNode();
    new_node->sharing_info = remaining;
    new_node->right_sibling = node->left_child;
    node->left_child = new_node;
    insertRecursive(new_node, motif, depth+1);
  } else {
    uint32_t child_mask = (1<<ch);
    TreeNodeSlow *child = node->left_child;
    while (child) {
      uint32_t info = child->sharing_info;
      if (info & child_mask) break;
      child = child->right_sibling;
    }
    if (!child) {
      TreeNodeSlow* new_node = TreeNodeSlow::allocateNode();
      new_node->sharing_info = child_mask;
      new_node->right_sibling = node->left_child;
      node->left_child = new_node;
      insertRecursive(new_node, motif, depth+1);
    } else {
      if (child->sharing_info & (~child_mask)) {
        TreeNodeSlow* new_node = TreeNodeSlow::allocateNode();
        new_node->sharing_info = child_mask;
        child->sharing_info &= ~child_mask;
        new_node->right_sibling = node->left_child;
        node->left_child = new_node;
        copyUnion(new_node, child, motif, depth+1);
      } else {
        insertRecursive(child, motif, depth+1);
      }
    }
  }
}

void MotifTreeSlow::insertRecursiveNew(TreeNodeSlow* node, const std::string& motif, uint32_t depth) {
  //std::cout << "INSERT RECURSIVE: motif=" << motif << ", depth=" << depth << std::endl;
  if (depth >= max_depth) {
    return;
  }
  uint32_t ch = motif[depth];
  uint32_t remaining = (ch >= domain_size) ? mask : (1 << ch);
  TreeNodeSlow *child = node->left_child;
  while (remaining && child) {
    uint32_t common = child->sharing_info & remaining;
    remaining &= ~child->sharing_info; 
    if (common) {
      if (common != child->sharing_info) {
        child->sharing_info &= ~common;
        TreeNodeSlow* new_node = TreeNodeSlow::allocateNode();
        new_node->sharing_info = common;
        new_node->right_sibling = node->left_child;
        node->left_child = new_node;
        copyUnion(new_node, child, motif, depth+1);
      } else {
        insertRecursive(child, motif, depth+1);
      }
    }
    child = child->right_sibling;
  }
  // std::cout << "remaining = " << std::bitset<4>(remaining).to_string() << std::endl;
  if (!remaining) return;
  TreeNodeSlow* new_node = TreeNodeSlow::allocateNode();
  new_node->sharing_info = remaining;
  new_node->right_sibling = node->left_child;
  node->left_child = new_node;
  insertRecursive(new_node, motif, depth+1);
}

void MotifTreeSlow::intersectRecursive(TreeNodeSlow* to, const TreeNodeSlow* from, uint32_t depth) {
  //std::cout << "Intersect recursive : to=" << to << ", from=" << from << ", depth=" << depth << std::endl;
  TreeNodeSlow *to_child = to->left_child;
  to->left_child = 0;
  while (to_child) {
    uint32_t intersecting = 0;
    for (TreeNodeSlow *from_child = from->left_child; from_child; from_child = from_child->right_sibling) {
      if (from_child->sharing_info & to_child->sharing_info) intersecting++;
    }
    TreeNodeSlow *node = to_child;
    to_child = to_child->right_sibling;
    if (!intersecting) {
      if (depth < max_depth-1) deleteNode(node);
    } else {
      TreeNodeSlow *from_child;
      for (from_child = from->left_child; from_child && (intersecting > 1); from_child = from_child->right_sibling) {
        //std::cout << "  to_child=" << to_child << "(" << to_child->sharing_info << ")" << ", from=" << intersecting[i] << "(" << intersecting[i]->sharing_info << "), depth=" << depth << std::endl;
        TreeNodeSlow *new_node = TreeNodeSlow::allocateNode();
        uint32_t common = from_child->sharing_info & node->sharing_info;
        if (common) {
          copyIntersection(new_node, node, from_child, depth+1);
          if ((depth < max_depth-1) && (!new_node->left_child)) {
            TreeNodeSlow::deallocateNode(new_node);
          } else {
            new_node->sharing_info = common;
            new_node->right_sibling = to->left_child;
            to->left_child = new_node;
          }
          intersecting--;
        }
      }
      for (; from_child; from_child = from_child->right_sibling) {
        uint32_t common = from_child->sharing_info & node->sharing_info;
        if (common) {
          intersectRecursive(node, from_child, depth+1);
          if ((depth < max_depth-1) && (!node->left_child)) {
            TreeNodeSlow::deallocateNode(node);
          } else {
            node->sharing_info = common;
            node->right_sibling = to->left_child;
            to->left_child = node;
          }
          break;
        }
      }
    }
  }
}

#endif // __MOTIF_TREE__

