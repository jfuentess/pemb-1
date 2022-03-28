#ifndef BASELINE_HPP
#define BASELINE_HPP
#include <bits/stdc++.h>
#include "Node.hpp"
using namespace std;

class baseline{

private:
	Node *nodos;
	unsigned int tam;
	unsigned int level;
public:
	 baseline(){}
	 baseline(unsigned int tam, unsigned int level){
	 	this->tam = tam;
    	this->level = level;
    	this->nodos = new Node[tam];
	 }
	 void setPadre(unsigned int i,unsigned int padre){
	 	this->nodos[i].setPadre(padre);
	 }
	 void setVecino(unsigned int i,unsigned int padre){
	 	this->nodos[i].setVecino(padre);
	 }
	 void setHijo(unsigned int i,unsigned int padre){
	 	this->nodos[i].setHijo(padre);
	 }
	 unsigned int getPadre(unsigned int pos){
	 	return this->nodos[pos].getPadre();
	 }
	 unsigned int getSizeV(unsigned int pos){
	 	return this->nodos[pos].getSizeV();
	 }
	unsigned int getSizeH(unsigned int pos){
	 	return this->nodos[pos].getSizeH();
	 }
	 unsigned int getvecino(unsigned int a, unsigned int b){
	 	return this->nodos[a].getvecino(b);
	 }
	 	 unsigned int gethijo(unsigned int a, unsigned int b){
	 	return this->nodos[a].gethijo(b);
	 }
	 //recibe por entrada las posiciones en mi arreglo de nodos de los nodos a comparar, aca comparo si a esta contenido en b
	 bool contenido(unsigned int a,unsigned int b){
	 	while(a <= b){
	 		if(a == b)return true;
	 		a=getPadre(a);
	 	}
	 	return false;
	 }
	 //recibo como entrada la posicion de a en mi arreglo de nodos, ademas del nivel al cual pertenece y el nivel objetivo
	 // b siempre debe ser menor a a(es decir b es un nivel que representa un detalle mas fino)
	 vector <int> listcontenidos(unsigned int a, unsigned int nivela, unsigned int nivelb){
	 	int obj = nivela-nivelb;
	 	queue <unsigned int> analizar;
	 	vector <int> retorno;
	 	analizar.push(a);
	 	int aux = 1, cont = 0;
	 	while(obj > 0){
	 		cont = 0;
	 		for(int i = 0; i < aux; i++){
	 			int tope = analizar.front();
	 			analizar.pop();
	 			for(int j = 0; j < this->nodos[tope].getSizeH();j++){
	 				cont++;
	 				analizar.push(this->nodos[tope].gethijo(j));
	 			//	printf("estoy pusheando %d\n",this->nodos[tope].gethijo(j));
	 			}
	 			//puts("termino for");

	 		}
	 	aux = cont;
	 	obj--;
	 	}
	 	while(!analizar.empty()){
	 		retorno.push_back(analizar.front());
	 		analizar.pop();
	 	}
	 	return retorno;
	 }

	 bool touches(unsigned int a, unsigned int nivela, unsigned int b,unsigned int nivelb){
	 	if(a > b){
	 		swap(nivelb,nivela);
	 		swap(b,a);
	 	}
	 	bool comprobar = contenido(a,b);
	 	int obj = nivelb - nivela;
	 	for(int i = 0 ; i < getSizeV(a); i++){
	 		int aux = getvecino(a,i);
	 		for(int j = 0; j < obj; j++){
	 			aux = getPadre(aux);
	 		}
	 		if(aux == b && !comprobar)return true;
	 		else if(aux != b && comprobar)return true;
	 	}
	 	return false;
	 }
	 

	 vector<int> touches2(unsigned int a, unsigned int nivela, unsigned int b,unsigned int nivelb){
	 	if(nivela > nivelb){
	 		swap(nivelb,nivela);
	 		swap(b,a);
	 	}
	 	bool comprobar = contenido(a,b);
	 	int obj = nivelb - nivela;
	 	vector <int> retorno;
	 	for(int i = 0 ; i < getSizeV(a); i++){
	 		int aux = getvecino(a,i);
	 		retorno.push_back(aux);
	 	}
	 	return retorno;
	 }




};

#endif