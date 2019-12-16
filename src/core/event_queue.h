#pragma once

#include <RcppArmadillo.h>
#include <unordered_map>
#include <functional>
#include "array3.h"
#include "arma_helpers.h"

using namespace std;
using namespace arma;

struct node {
    double time;
    uvec3 index;
    node* up = nullptr;
    node* left = nullptr;
    node* right = nullptr;

    node(double time, uvec3 index) : time(time), index(index) {}

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
    uint size();
    pair<double, uvec3> next() { return pair<double, uvec3>(root->time, root->index); }
    void push(double time, uvec3 index);
private:
    node* root = nullptr;
    unordered_map<uvec3, node*> node_map;
    
    void insert(node * root, double time, uvec3 index);
    void update(double time, uvec3 index);
    template <typename T> void swap(T* a, T* b);
    void swap(node* a, node* b);
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

void event_queue::insert(node * root, double time, uvec3 index) {
    // if this is the first node, assign values to root
    if (root == nullptr) {
        *root = node(time, index);
        return;
    }

    // determine key-value pair to be pushed further
    if (time < root->time) {
        swap(&(root->time), &time);
        swap(&(root->index), &index);
        
        node_map.erase(index);
        node_map[root->index] = root;
    }

    // push key-value pair down tree
    if (root->left == nullptr || root->right == nullptr) {
        auto n = node(time, index);
        n.up = root;
        if (root->left == nullptr)
            root->left = &n;
        else
            root->right = &n;
    } else {
        insert(root->right->size() < root->right->size() ?
                       root->right :
                       root->left,
                   time,
                   index);
    }
}

void event_queue::update(double time, uvec3 index) {
    node* n = node_map[index];

    // if new value is smaller than branch, move up
    if (n->up != nullptr && time < n->up->time) {
        swap(n, n->up);
        update(time, index);
    // if larger than at least one twig, move down
    } else if ((n->left != nullptr && n->left->time < time) ||
               (n->right != nullptr && n->right->time < time)) {
        node* twig =
            n->left == nullptr ? n->right :
            n->right == nullptr ? n->left :
            n->left->time < n->right->time ? n->left : n->right;
        swap(n, twig);
        update(time, index);
    }
}

void event_queue::swap(node* a, node* b) {
    swap(&(a->time), &(b->time));
    swap(&(a->index), &(b->index));
    node_map[a->index] = a;
    node_map[b->index] = b;
}

void event_queue::push(double time, uvec3 index) {
    // if this is an update, use updating instead of insertion
    if (node_map.find(index) != node_map.end())
        update(time, index);
    else
        insert(root, time, index);
}
