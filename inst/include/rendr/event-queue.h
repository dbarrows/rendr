#pragma once

#include <RcppArmadillo.h>
#include <map>
#include <functional>
#include <core/array3.h>
#include <core/arma.h>

namespace rendr {

using namespace arma;
using namespace std;
using namespace core;
using uint = unsigned int;

struct node {
    double time;
    uvec3 index;
    node* up = nullptr;
    node* left = nullptr;
    node* right = nullptr;

    node(double time, uvec3 index) : time(time), index(index) {}
    ~node() {
        if (left != nullptr) delete left;
        if (right != nullptr) delete right;
    }

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

class event_queue {
public:
    event_queue() {}
    event_queue(array3<double> times);
    ~event_queue() { /*if (root != nullptr) delete root;*/ }
    uint size();
    pair<double, uvec3> next() { return pair<double, uvec3>(root->time, root->index); }
    void push(double time, uvec3 index);
private:
    node* root = nullptr;
    map<uvec3, node*> node_map;
    
    void insert(node * root, double time, uvec3 index);
    void update(double time, uvec3 index);
    template <typename T> void swap(T* a, T* b);
    void swap_data(node* a, node* b);
};
uint event_queue::size() { return root != nullptr ? root->size() : 0; }
event_queue::event_queue(array3<double> times) {
    for (uint i = 0; i < times.size(); i++)
        push(times[i], times.index3(i));
}

template <typename T>
void event_queue::swap(T* a, T* b) {
    T a_tmp = *a;
    *a = *b;
    *b = a_tmp;
}

void event_queue::insert(node* n, double time, uvec3 index) {
    // determine key-value pair to be pushed further
    if (time < n->time) {
        swap(&(n->time), &time);
        swap(&(n->index), &index);
        
        node_map.erase(index);
        node_map[n->index] = n;
    }

    // push key-value pair down tree
    if (n->left == nullptr || n->right == nullptr) {
        auto new_node = new node(time, index);
        node_map[index] = new_node;
        new_node->up = n;
        if (n->left == nullptr)
            n->left = new_node;
        else
            n->right = new_node;
    } else {
        insert(n->right->size() < n->right->size() ?
                       n->right :
                       n->left,
                   time,
                   index);
    }
}

void event_queue::update(double time, uvec3 index) {
    node* n = node_map[index];
    n->time = time;

    // if new value is smaller than branch, move up
    if (n->up != nullptr && time < n->up->time) {
        swap_data(n, n->up);
        update(time, index);
    // if larger than at least one twig, move down
    } else if ((n->left != nullptr && n->left->time < time) ||
               (n->right != nullptr && n->right->time < time)) {
        node* twig =
            n->left == nullptr ? n->right :
            n->right == nullptr ? n->left :
            n->left->time < n->right->time ? n->left : n->right;
        swap_data(n, twig);
        update(time, index);
    }
}

void event_queue::swap_data(node* a, node* b) {
    swap(&(a->time), &(b->time));
    swap(&(a->index), &(b->index));
    node_map[a->index] = a;
    node_map[b->index] = b;
}

void event_queue::push(double time, uvec3 index) {
    // if this is the first node, assign values to root
    if (root == nullptr) {
        root = new node(time, index);
        node_map[index] = root;
        return;
    }

    // if this is an update, use updating instead of insertion
    if (node_map.find(index) != node_map.end()) {
        update(time, index);
    } else {
        insert(root, time, index);
    }
}

}
