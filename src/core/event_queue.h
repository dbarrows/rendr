#pragma once

#include <RcppArmadillo.h>
#include <unordered_map>
#include <functional>
#include "array3.h"

using namespace std;
using namespace arma;

struct node {
    double time;
    uvec3 index;
    node* up = nullptr;
    node* left = nullptr;
    node* right = nullptr;

    uint size();
    node* smaller_twig();
};
uint node::size() {
    return 1 +
        (left != nullptr ? left->size() : 0) +
        (right != nullptr ? right->size() : 0);
}
node* node::smaller_twig() {
    return left == nullptr && right == nullptr ? nullptr :
        left != nullptr && (right == nullptr || left->time < right->time) ? left :
        right;
}


// hashing and equality functions for arma::uvec3.
// Required to use arma::uvec3 as a key in std::unordered_map.
template <>
struct hash<uvec3> {
    std::size_t operator()(const uvec3& k) const {
        size_t res = 17;
        res *= 31 + hash<uint>()(k[0]);
        res *= 31 + hash<uint>()(k[1]);
        res *= 31 + hash<uint>()(k[2]);
        return res;
    }
};
template <>
struct equal_to<uvec3> {
    bool operator()(const uvec3& a, const uvec3& b) const {
        return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
    }
};


class event_queue {
public:
    event_queue(array3<double> times);
    /*event_queue() {
        auto hash = [](const uvec3& v){
            size_t res = 17;
            res *= 31 + std::hash<uint>()(v[0]);
            res *= 31 + std::hash<uint>()(v[1]);
            res *= 31 + std::hash<uint>()(v[2]);
            return res;
        };
        auto equal = [](const uvec3& a, uvec3& b) {
            return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
        };
        node_map = std::unordered_map<
            uvec3,
            node*,
            function<size_t(const uvec3& v)>,
            function<bool(const uvec3& a, uvec3& b)>
        >(hash, equal);
    }*/

    uint size();
    uvec3 next() { return root->index; }
    void push(double time, uvec3 index);


private:
    node* root = nullptr;
    /*std::unordered_map<uvec3,
                       node*, function<size_t(const uvec3& v)>, function<bool(const uvec3& a, uvec3& b)>> node_map;*/
    unordered_map<uvec3, node*> node_map;
    
    void insert_new(node * root, double time, uvec3 index);
    void update(double time, uvec3 index);
    void swap(node* a, node* b);
};
uint event_queue::size() { return root != nullptr ? root->size() : 0; }
event_queue::event_queue(array3<double> times) {

}

void event_queue::insert_new(node * root, double time, uvec3 index) {
    // determine key-value pair to be pushed further
    if (time < root->time) {
        double tmp_time = root->time;
        uvec3 tmp_index = root->index;
        root->time = time;
        root->index = index;
        time = tmp_time;
        index = tmp_index;
        
        node_map.erase(index);
        node_map[root->index] = root;
    }

    // push key-value pair down tree
    if (root->left == nullptr || root->right == nullptr) {
        node n;
        n.time = time;
        n.index = index;
        n.up = root;
        if (root->left == nullptr)
            root->left = &n;
        else
            root->right = &n;
    } else {
        insert_new(root->right->size() < root->right->size() ?
                       root->right :
                       root->left,
                   time,
                   index);
    }
}

void event_queue::update(double time, uvec3 index) {

}

void event_queue::swap(node* a, node* b) {

}

void event_queue::push(double time, uvec3 index) {

}