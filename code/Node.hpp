#ifndef NODE_COMP_HPP
#define NODE_COMP_HPP
#include <bits/stdc++.h>

using namespace std;
class Node {
  
private:
  unsigned int padre; 
  vector<unsigned int> hijos;
  vector <unsigned int> vecinos; 
 
public:
  Node(){
    this->padre = -1;
  }
  
  ~Node(){}

  unsigned int getPadre() {
    return this->padre;
  }

  unsigned int getSizeH() {
    return this->hijos.size();
  }
  unsigned int getSizeV() {
    return this->vecinos.size();
  }
  unsigned int gethijo(unsigned int f) {
    return this->hijos[f];
  }

  vector<unsigned int> getVhijo() {
    return this->hijos;
  }

  unsigned int getvecino(unsigned int f) {
    return this->vecinos[f];
  }

  void setPadre(unsigned int f) {
    this->padre = f;
  }

  void setHijo(unsigned int l) {
    this->hijos.push_back(l);
  }
  void setVecino(unsigned int l) {
    this->vecinos.push_back(l);
  }
};

#endif
