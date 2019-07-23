cmake_minimum_required(VERSION 3.8)
project(src)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        ems.cpp
        ems.hpp
        ems1.hpp
        ems2.hpp
        ems2p.hpp
        intersect.hpp
        motif.hpp
        motif_set.hpp
        motif_tree.hpp
        motif_tree_base.hpp
        motif_tree_fast.hpp
        motif_trie.hpp
        utils.h)

add_executable(src ${SOURCE_FILES})