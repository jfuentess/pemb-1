#ifndef GRAPH_COMP_HPP
#define GRAPH_COMP_HPP

#include <iostream>
#include <stack>
#include <complementary/Tree.hpp>
#include <map>
#include <utility>
#include <vector>
#include "auxiliar.hpp"
#include <sdsl/int_vector.hpp>

using namespace std;
using namespace sdsl;

#define CLOSE_PAR 0
#define OPEN_PAR 1
#define EMPTY 2
#define NONE 3
#define no_visited false
#define no_traversed false
#define visitado true
#define traversed true


class Graph {

private:
  Vertex *V;      // Array of vertices of the tree
  Edge *E;        // Array of edges of the tree. It is the concatenation of the
            // adjacency lists of all nodes 
  unsigned int n; // Number of vertices in the tree
  unsigned int m; // Number of edges in the tree


public:
  Graph () {
    this->n = 0;
    this->m = 0;
  }

  Graph (unsigned int n, unsigned int m) {
    this->n = n;
    this->m = m;
    this->V = new Vertex[n];
    this->E = new Edge[2*m];
  }

  unsigned int vertices() {
    return n;
  }

  unsigned int edges() {
    return m;
  }

  // E[i].src = src and E[i].tgt = tgt
  void setEdge(int i, unsigned int src, unsigned int tgt) {
    this->E[i].setSrc(src);
    this->E[i].setTgt(tgt);
    this->E[i].setCmp(-1);
    this->E[i].setOrientation(1); // Direct
  }

  // V[i].first = first and V[i].last = last
  void setVertex(int i, unsigned int first, unsigned int last) {
    this->V[i].setFirst(first);
    this->V[i].setLast(last);
  }

  // V[i].first = first
  void setVertexFirst(int i, unsigned int first) {
    this->V[i].setFirst(first);
  }

  // V[i].last = last
  void setVertexLast(int i, unsigned int last) {
    this->V[i].setLast(last);
  }

  void setEdgeSrc(int i, unsigned int s) {
    this->E[i].setSrc(s);
  }

  void setEdgeTgt(int i, unsigned int t) {
    this->E[i].setTgt(t);
  }

  void setEdgeCmp(int i, int c) {
    this->E[i].setCmp(c);
  }

  void setEdgeOrientation(int i, bool o) {
    this->E[i].setOrientation(o);
  }

  unsigned int getEdgeSrc(int i) {
    return this->E[i].getSrc();
  }

  unsigned int getEdgeTgt(int i) {
    return this->E[i].getTgt();
  }

  int getEdgeCmp(int i) {
    return this->E[i].getCmp();
  }

  bool getEdgeOrientation(int i) {
    return this->E[i].getOrientation();
  }

  // Return the next edge in counter-clockwise order
  unsigned int getNextEdgeCCW(int i) {

    uint src = this->getEdgeSrc(i);
    uint first = this->getVertexFirst(src);
    uint last = this->getVertexLast(src);
    
    return ((uint)(i+1) > last)? first : i+1;
  }

  // Return the next edge in clockwise order
  unsigned int getNextEdgeCW(int i) {

    uint src = this->getEdgeSrc(i);
    int first = (int)this->getVertexFirst(src);
    int last = (int)this->getVertexLast(src);
    
    return (i-1 < first)? last : i-1;
  }

  Vertex getVertex(int i) {
    return this->V[i];
  }

  Edge getEdge(int i) {
    return this->E[i];
  }

  unsigned int getVertexFirst(int i) {
    return this->V[i].getFirst();
  }

  unsigned int getVertexLast(int i) {
    return this->V[i].getLast();
  }


// niveles es h, la cantidad de granularidades, 
// init es el nodo inicial, nodo 0
// jerarquias es una lista que contiene a que region se mapea cada region perteneciente al nivel mas fino
// total es la cantidad total de regiones existentes
// vector que tiene los limites de cada nivel de granularidad, tomamos como 0 la mas fina, por tanto limites[0] = 0 , limite[1] = n, esto es limite acumulado
// es decir limite[2] = limite[1] + x donde x es la cantidad de regiones a tal nivel
//Retorna estructura auxiiar creada con el fin de solo almacenar informacion calculada

auxiliar dfs_spanning_tree_propio(unsigned int init, unsigned int
        *jerarquias, unsigned int *limites, unsigned int niveles,unsigned int total) {
    unsigned int n = this->vertices();
    unsigned int m = this->edges();    
    bool *visited_vertex = new bool[n](); 
    bool *visited_edge = new bool[m*2]();
    bool *pertenece = new bool[m*2]();
    bool *visited_region = new bool[total](); 
    bool *terminada = new bool[total](); 
    bool *bits = new bool[n*niveles*2](); 
    bool *bits2 = new bool[n*niveles](); 
    int *original = new int[n*niveles]();
    int *reverso = new int[n*niveles]();
    int *cierres = new int[n*niveles]();
    int *posiciones = new int[niveles]();
    int *totalabiertos = new int[niveles]();
    vector <map <pair <int,int>,int > > prueba;
    for(int i = 0; i < niveles; i++)prueba.push_back(map<pair<int,int>,int>());
          vector <map <pair <int,int>,int > > prueba2;
   for(int i = 0; i < niveles; i++)prueba2.push_back(map<pair<int,int>,int>());
    vector <stack <int> > analizados(niveles);
    for(int i = 0; i < niveles; i++)posiciones[i] = 0;
      for(int i = 0; i < niveles; i++)totalabiertos[i] = 0;
    original [0] = init;
    reverso[init] = 0;
    map<pair <int,int>, int > mapa;
        map<pair <int,int>, int > mapa2;
    int cant1 = 0, cant2 = 0,cant3 = 0,cant4=0;
    vector <int> st;
    int a1= 0,a2= 0,a3 = 0,a4 = 0;
    for(int i = 0; i < n ; i++)visited_vertex[i] = no_visited;
    for(int i = 0; i < m*2 ; i++)visited_edge[i] = no_traversed;
    for(int i = 0; i < m*2 ; i++)pertenece[i] = no_traversed;
    for(int i = 0; i < total ; i++)visited_region[i] = no_visited;
    for(int i = 0; i < total ; i++)visited_region[i] = no_visited;
    for(int i = 0; i < n*niveles;i++)cierres[i] = -1;
    vector< vector <char> >secuencia(niveles,vector<char>()); 
    stack <unsigned int> o,e;
    visited_vertex[init] = visitado;
    for(int i = 0; i < niveles; i++){
     int regionaux =  jerarquias[init*niveles + i];
     visited_region[regionaux + limites[i]] = visitado;
    }

    for(int i = this->V[init].getLast(); i >= this->V[init].getFirst();i-- ) {
      o.push(i);
      if(i == this->V[init].getFirst())break;
    }
    for(int i = 0; i < niveles; i++){
    secuencia[i].push_back('[');
    secuencia[i].push_back('(');
    }
    int posbit = 0;
        int posbit2 = 0;
    for(int i = 0; i < niveles; i++){
      bits[posbit + (i* n*2)] = true;
            bits2[posbit2 + (i* n)] = true;
           analizados[i].push(init);
      cierres[posiciones[i] + (i*n)] = posbit2 + (i*n);
    }
    posbit++;
    posbit2++;
    int contprueba = 0;

    while(!o.empty()){
      unsigned int curr = o.top(); o.pop();
      bool auxcheck = visited_vertex[this->E[curr].getTgt()] , auxcheckr = true;
      //Compruebo condicion para que pertenezca al st, que no forme ciclo a niveles superiores.
      for(int i = 1; i < niveles; i++){
        int regionA =  jerarquias[this->E[curr].getSrc() *niveles + i] + limites[i];
        int regionB =  jerarquias[this->E[curr].getTgt() *niveles + i] + limites[i];
        if(regionA == regionB){
          auxcheckr = true;
          break;
        }
        if(visited_region[regionB] == visitado){
          auxcheckr = false;
          break;
        }
      }
      if(!auxcheck && auxcheckr && visited_edge[curr] == no_visited){ //cumple con condicion para visitar vertices y es primera vez que se visita, pertenece al st
        secuencia[0].push_back('(');
        cant1++;
        st.push_back(curr);
        visited_edge[curr] = traversed;
        visited_edge[this->E[curr].getCmp()] = traversed;
        visited_vertex[this->E[curr].getTgt()] = visitado;

        pertenece[curr] = visitado;
        pertenece[this->E[curr].getCmp()] = visitado;
        e.push(curr);
        original[posbit2] = this->E[curr].getTgt();
        reverso[this->E[curr].getTgt()] = posbit2; 
        unsigned int nodotgt = this->E[curr].getTgt();
        for(unsigned int i = this->E[curr].getCmp(),j = 0; j < (this->V[nodotgt].getLast()-this->V[nodotgt].getFirst())+1; i--,j++){
          o.push(i);
          if(i == this->V[nodotgt].getFirst())i = this->V[nodotgt].getLast()+1;
        }
        for(unsigned int i = 0; i < niveles; i++){
          unsigned int p =  jerarquias[this->E[curr].getTgt() *niveles + i] + limites[i];
          unsigned int q =  jerarquias[this->E[curr].getSrc() *niveles + i] + limites[i];

     // cout << (2*i*n) + posbit << endl;
          //Calculo secuencia parentesis para niveles superiores al 0 y calculo sencuencia de bits(secuencia representa primera visita a una region)
          if(visited_region[p] == no_visited && mapa[make_pair(p,q)] == 0 && p!= q ){
 if(i!= 0)secuencia[i].push_back('(');
 if(i == 1)a1++;
          mapa2[make_pair(p,q)]++;
          mapa2[make_pair(q,p)]++;
          mapa[make_pair(p,q)]++;
          mapa[make_pair(q,p)]++;
          if(i != 0)prueba[i][make_pair(q,p)] = jerarquias[this->E[curr].getTgt() *niveles + (i-1) ] + limites[i-1] +1;
          if(i != 0)prueba[i][make_pair(p,q)] = jerarquias[this->E[curr].getSrc() *niveles + (i-1) ] + limites[i-1] +1;
              visited_region[p] = visitado;
              bits[(i*n*2) + posbit] = true;
              bits2[(i*n) + posbit2] = true;
              totalabiertos[i]++;
              analizados[i].push(totalabiertos[i]);
              posiciones[i] = analizados[i].top();
              cierres[posiciones[i] + i*n] = posbit2 + (i*n);
          }
          else {
            posiciones[i] = analizados[i].top();
            cierres[posiciones[i] + i*n] = posbit2 + (i*n);
          }
        }
        posbit++;
        posbit2++;

      }

      else if( visited_edge[curr] == no_visited ){ // segundo caso, primera vez visito arista pero esta no pertenece al st
                          visited_edge[curr] = traversed;
                visited_edge[this->E[curr].getCmp()] = traversed;
        cant2++;
        for(int i = 0; i < niveles;i++){
        unsigned int p =  jerarquias[this->E[curr].getTgt() *niveles + i] + limites[i];
        unsigned int q =  jerarquias[this->E[curr].getSrc() *niveles + i] + limites[i];




bool auxcheckraux = true;
      for(int j = i+1; j < niveles; j++){

        int regionA =  jerarquias[this->E[curr].getSrc() *niveles + j] + limites[j];
        int regionB =  jerarquias[this->E[curr].getTgt() *niveles + j] + limites[j];
        if(regionA == regionB){
          auxcheckraux = true;
          break;
        }
        if(visited_region[regionB] == visitado){
          auxcheckraux = false;
          break;
        }
      }



        if(   !(auxcheckraux && visited_region[p] == no_visited) && p != q && mapa[make_pair(p,q)]==0 ){
                    if(i != 0)prueba2[i][make_pair(q,p)] = jerarquias[this->E[curr].getTgt() *niveles + (i-1) ] + limites[i-1] +1;
                               if(i != 0)prueba2[i][make_pair(p,q)] = jerarquias[this->E[curr].getSrc() *niveles + (i-1) ] + limites[i-1] +1;
          secuencia[i].push_back('[');
           if(i == 1)a2++;
          mapa[make_pair(p,q)]++;
          mapa[make_pair(q,p)]++;
          }
          else if(i == 0)secuencia[i].push_back('[');
        }
      }

//Tercer caso, arista ya fue visitada una vez y esta si pertenecia al st
      else if(pertenece[curr]){
        cant3++;

        for(int i = 0; i < niveles;i++){
        unsigned int p =  jerarquias[this->E[curr].getTgt() *niveles + i] + limites[i];
        unsigned int q =  jerarquias[this->E[curr].getSrc() *niveles + i] + limites[i];
        if(p != q && mapa[make_pair(p,q)]==1 && mapa2[make_pair(p,q)]==1 && (jerarquias[this->E[curr].getSrc() *niveles + (i-1)] + limites[i-1]) == prueba[i][make_pair(p,q)]-1 && i != 0
           && (jerarquias[this->E[curr].getTgt() *niveles + (i-1)] + limites[i-1]) == prueba[i][make_pair(q,p)]-1){
          secuencia[i].push_back(')');
                      bits[(2*i*n) + posbit] = true;
           if(i == 1)a3++;
          if(i != 0){
          analizados[i].pop();
            posiciones[i] = analizados[i].top();
          }
                    mapa[make_pair(p,q)]++;
          mapa[make_pair(q,p)]++;
          }

        }
                    secuencia[0].push_back(')');
            bits[posbit] = true;
       posbit++;
      }

//cuarto caso, arista ya fue visitada 1 vez y esta no pertenece al st.
      else {
        cant4++;
        for(int i = 0; i < niveles;i++){
        unsigned int p =  jerarquias[this->E[curr].getTgt() *niveles + i] + limites[i];
        unsigned int q =  jerarquias[this->E[curr].getSrc() *niveles + i] + limites[i];
        if(p != q && mapa[make_pair(p,q)]==1 && i != 0 && (jerarquias[this->E[curr].getSrc() *niveles + (i-1)] + limites[i-1]) == prueba2[i][make_pair(p,q)]-1 && i != 0
          && (jerarquias[this->E[curr].getTgt() *niveles + (i-1)] + limites[i-1]) == prueba2[i][make_pair(q,p)]-1){
          secuencia[i].push_back(']');
          mapa[make_pair(p,q)]++;
          mapa[make_pair(q,p)]++;
          }
                    else if(i == 0)secuencia[i].push_back(']');
        }
      }


    }


    //agrego para cada nivel parentesis de cierre
        for(int i = 0; i < niveles; i++){
        secuencia[i].push_back(')');
        secuencia[i].push_back(']');
      }
    //      puts("aun vivo");
      //calculo vector auxiliar que sirve para realizar mapeo, esto solo es utilizado para comparar la misma region en la implementacion base.
      vector<int> nivelaux(niveles,0);
    for(int i = 0; i < n; i++){
      int source = original[i];
      for(int j = 1; j < niveles; j++){
        if(bits2[j*n +i] == 1){
          original[j*n +i] = jerarquias[source*niveles + j];
          reverso[jerarquias[source*niveles+j]+n*j]= nivelaux[j];
          nivelaux[j]++;
        }
      }
    }
   // cout << "n = " << n << endl; 
   // puts("imprimire");
    /*
        for(int i = 0; i < niveles; i++)bits[n*i*2 +n*2 -1] = 1;
    for(int i = 0; i < niveles; i++){
      for(int j = 0; j < n*2;j++){
        printf("%d",bits[i*n*2 + j]);
      }
      puts("");
    }
    */
   // puts("termine");
//for(int j = 0; j < secuencia.size();j++)        for(int i = 0; i <secuencia[j].size();i++)printf("%c",secuencia[j][i]);
	for(int i = 0; i < prueba.size();i++){
prueba[i].clear();
prueba2[i].clear();
}	
mapa.clear();
mapa2.clear();
    return auxiliar(niveles,secuencia,bits,original,reverso,cierres);
  }





/*
  void connected_graph() {
    unsigned int n = this->vertices();
    stack <unsigned int> s;
    char *visited = new char[n]();
    int curr = -1;
    unsigned int edge;
    int first = 1;
    int num_vertices = 0;
while(!s.empty() || first) {
      if(first) { // Root
  curr = 0;
  first = 0;
      }
      else {
  edge = s.top(); s.pop();
  curr = this->E[edge].getTgt(); // current
      }
      visited[curr] = 1;
      for(unsigned int i = this->V[curr].getFirst(); i <= this->V[curr].getLast();
     {
  if(!visited[this->E[i].getTgt()])
    s.push(i);      
      }
    }

    for(unsigned int i = 0; i < n; i++)
      if(visited[i] == 0) { // There are unvisited vertices
  num_vertices++;
      }
    
    cout << "unvisited vertices: " << num_vertices << ", visited vertices: " <<
      n-num_vertices << endl;
  }

  int_vector<> ps_tree_encoding() {
    unsigned int n = this->vertices();
    unsigned int m = this->edges();

    int_vector<> S(4*n-5, NONE, 8); // Output string
    char *visited = new char[2*m]();
    /*** Initial setting ***/
    
    // Set the external edges as visited (they are not considered in the
    // traversal)
/*
    uint ext0 = 0; // By definition
    uint ext1 = this->getEdgeTgt(this->getVertexFirst(ext0));
    uint ext2 = this->getEdgeTgt(this->getVertexLast(ext0));
    visited[this->getVertexFirst(ext0)] = 1;
    visited[this->getVertexLast(ext0)] = 1;
    visited[this->getVertexFirst(ext1)] = 1;
    visited[this->getVertexLast(ext1)] = 1;
    visited[this->getVertexFirst(ext2)] = 1;
    visited[this->getVertexLast(ext2)] = 1;

    // For any maximal planar graph, the following entries of S are already defined
    S[0] = OPEN_PAR;       S[4*n-6] = CLOSE_PAR;
    S[1] = OPEN_PAR;       S[4*n-7] = CLOSE_PAR;
    S[2] = OPEN_PAR;       S[4*n-8] = CLOSE_PAR;
    S[3] = EMPTY;
    
    /*** Traversal ***/
  /*  
    uint v = 0; // Starting vertex (by definition)
    uint e = this->getVertexLast(v); // Starting edge (by definition)
    uint end_e = this->getVertexFirst(v); // Stop condition
    uint ee = this->getNextEdgeCW(e);
    uint id = 4;

    while(ee != end_e) {
      if(!visited[ee] && this->getEdgeOrientation(ee)) {
      visited[ee] = 1;
      visited[this->getEdgeCmp(ee)] = 1;
      S[id++] = EMPTY;
      e = ee;
      }
      else if(!visited[ee] && !this->getEdgeOrientation(ee)) {
      visited[ee] = 1;
      visited[this->getEdgeCmp(ee)] = 1;
      S[id++] = OPEN_PAR;
      e = this->getEdgeCmp(ee);
      }
      else if(visited[ee] && this->getEdgeOrientation(ee)) {
      S[id++] = CLOSE_PAR;
      e = this->getEdgeCmp(ee);
      }
      else {
      e = ee;
      }
      ee = this->getNextEdgeCW(e);

    }
/*    
    /* Convert the cw traversal into a ccw traversal */
/*
    int mid = (4*n-5) / 2;

    // Reverse the sequence S
    for(int i=0, j=4*n-6; i < mid; i++, j--) {
      char cp = S[i];
      S[i] = S[j];
      S[j] = cp;
    }

    // Change close parentheses by open parentheses
    for(uint i=0; i < 4*n-5; i++) {
      if(S[i] == OPEN_PAR) S[i] = CLOSE_PAR;
      else if(S[i] == CLOSE_PAR) S[i] = OPEN_PAR;
    }

    return S;
  }
*/  
};

#endif

