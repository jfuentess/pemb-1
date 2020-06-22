#ifndef NODE_COMP_HPP
#define NODE_COMP_HPP

// Used to represent graphs and trees
class Node {
  
private:
  unsigned int padre; // Position of the first incident edge of a vertex in E
  unsigned int *hijos; // Position of the last incident edge of a vertex in E
 
public:
  Node(){}
  
  ~Node(){}

  unsigned int getFirst() {
    return this->first;
  }

  unsigned int getLast() {
    return this->last;
  }

  void setFirst(unsigned int f) {
    this->first = f;
  }

  void setLast(unsigned int l) {
    this->last = l;
  }
};

#endif
