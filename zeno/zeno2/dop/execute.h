#pragma once


#include <zeno2/dop/Node.h>


namespace zeno2::dop {


ztd::zany getval(Input const &input);
ztd::zany resolve(Input const &input, std::set<Node *> &visited);
void touch(Input const &input, std::vector<Node *> &tolink, std::set<Node *> &visited);
void sortexec(Node *root, std::vector<Node *> &tolink, std::set<Node *> &visited);


inline ztd::zany resolve(Input const &input) {
    std::set<Node *> visited;
    return resolve(input, visited);
}


}
