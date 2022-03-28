#ifndef AUXILIAR_COMP_HPP
#define AUXILIAR_COMP_HPP
#include <bits/stdc++.h>
// Used to represent graphs and trees
class auxiliar {
  
private:
  unsigned int levels,ni; // Cantidad de niveles de la jerarquia,ni = nodo inicial
  vector < vector<char> > parentesis;// Position of the last incident edge of a vertex in E
  bool *bits,*cbits;
  int *mapeo,*reverso,*cierres;



public:
  auxiliar(unsigned int niveles,vector < vector<char> > parentesis, bool *bits, int *mapeo ,int *reverso,int *cierres,bool *cbits){
    this->levels = niveles;
    this->parentesis = parentesis;
    this->bits = bits;
        this->cbits = cbits;
    this->mapeo = mapeo;
    this->reverso = reverso;
    this->cierres = cierres;
  }
    auxiliar(){
  }
  
  ~auxiliar(){}

  unsigned int geLevels() {
    return this->levels;
  }

  vector < vector<char> >getParentesis() {
    return this->parentesis;
  }

  bool *getBits() {
    return this->bits;
  }
    bool *getcBits() {
    return this->cbits;
  }
  int *getcierres() {
    return this->cierres;
  }
  int *getMapeo(){
    return this->mapeo;
  }
  int getposMapeo(int pos){
    return this->mapeo[pos];
  }
  int *getReverso(){
    return this->reverso;
  }
  int getposReverso(int pos){
    return this->reverso[pos];
  }
};

#endif
