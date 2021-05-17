#include <iostream>
#include <stack>
#include <fstream>
#include <cstdlib>
#include <string.h>
#include "graph.hpp"//
#include "utils.hpp"//
#include "auxiliar.hpp"
#include "pemb-wt.hpp"


using namespace std;
using namespace sdsl;

pemb<> *pe;

int main(int argc, char **argv){
	int total;
	unsigned int n,niveles,aux,mitotal = 0;
	Graph g = read_graph_from_file(argv[1]);
	FILE *fp = fopen(argv[2], "r");
	if (!fp) {
    fprintf(stderr, "Error opening file \"%s\".\n", argv[2]);
    exit(EXIT_FAILURE);
	}
	n = g.vertices();
    fscanf(fp,"%d",&niveles);
	unsigned int *limites = new unsigned int[niveles+1](); 
	unsigned int *totalnivel = new unsigned int[niveles](); 
	unsigned int *jerarquias = new unsigned int[niveles*n](); 
	vector <int> cantidad;
	for(int i = 0; i < niveles; i++){
		fscanf(fp,"%d",&aux);
		cantidad.push_back(aux);
		limites[i] = mitotal;
		mitotal+= aux;
		totalnivel[i] = aux;
	}
	limites[niveles] = mitotal;	
	for(int i = 0; i < n*niveles; i++){
		fscanf(fp,"%d",&aux);
		jerarquias[i] = aux;
	}
	vector<int> regionA,regionB,nivelA,nivelB;
//	    puts("vivo");
    pe = new pemb<>(g,0,jerarquias,limites,niveles,mitotal,n);
   puts("vivo");
	pe->iniresumen(atoi(argv[7]));
	puts("TERMINE INI");
//	pe->frecacumulada();
//	return 0;
/*
//write_structure<HTML_FORMAT>(*pe,argv[5]);
for(int i = 0; i < niveles-1; i++){
	for(int j = 0; j < cantidad[i]; j++){
		for(int k = 0; k < cantidad[i+1]; k++){
  if(pe->Inside(pe->getReverso(j,i),i,pe->getReverso(k,i+1),i+1)){
  	printf("%d\n",k);
  	break;
  }

		}
	}
}
return 0;
*/










	srand (time(NULL));
	FILE *fp3 = fopen(argv[3], "r");
	fscanf(fp3,"%d",&total);
	for(int i = 0; i < total; i++){
		int a1,a2,a3,a4;
		fscanf(fp3,"%d %d %d %d",&a1,&a2,&a3,&a4);
		regionA.push_back(a1);
		nivelA.push_back(a2);
		regionB.push_back(a3);
		nivelB.push_back(a4);
	}


    int divisor = (niveles*(niveles-1))/2;
    divisor *= 200;
    	int antA=0,antB = 1,contadornivel = 0;


	FILE * asalida;
	asalida = fopen (argv[4],"a+");
	FILE *anavecinos;
	anavecinos = fopen (argv[6],"a+");
	unsigned t0, t1,auxt1, auxt2;
	double time,maximo,minimo,tactual,promedio = 0.0;
	bool primero = true;
	
	t0=clock();
	int contador = -1;
	for(int i = 0; i < total;i++){
		if(nivelA[i] == nivelB[i])continue;
		if( (antA != nivelA[i] || antB != nivelB[i]) && nivelA[i] != nivelB[i]){
	primero = true;
	t1=clock();
	time = (double(t1-t0)/1000);	
	promedio +=time/contadornivel;
	fprintf(asalida,"%d %d %lf %lf %lf\n",antA,antB,time/contadornivel,maximo,minimo);		
	antA = nivelA[i];
	antB = nivelB[i];
	t0=clock();
	contadornivel = 0;
		}
//		cout << "consulto por "<< nivelA[i] << " " << nivelB[i] << endl;
	auxt1=clock();
	//cout << pe->getReverso(regionA[i],nivelA[i]) << " " << pe->getReverso(regionB[i],nivelB[i]) << endl;
  if(pe->Inside(pe->getReverso(regionA[i],nivelA[i]),nivelA[i],pe->getReverso(regionB[i],nivelB[i]),nivelB[i]))printf("%d %d %d %d\n",regionB[i], nivelB[i], regionA[i],nivelA[i]);
    auxt2=clock();
    tactual = double(auxt2-auxt1)/1000.0;
    if(primero){
    	primero = false;
    	maximo = tactual;
    	minimo = maximo;
    }
    else{
    	if(tactual > maximo)maximo = tactual;
    	if(tactual < minimo) minimo = tactual;
    }
    contadornivel++;
	}
	t1=clock();
	
	time = (double(t1-t0)/1000.0);
	promedio += time/contadornivel;
	fprintf(asalida,"%d %d %lf %lf %lf\n",antA,antB,time/contadornivel,maximo,minimo);	
	fprintf(asalida, "PROMEDIO = %lf\n",promedio/15.0);
	double tiempofinalh = pe->obtenertiempoh();
	double tiempofinals = pe->obtenertiempos();
        double tiempofinalhm = pe->obtenertiempohm();
        double tiempofinalsm = pe->obtenertiemposm();
//	tiempofinalh /=1000.0;
//	tiempofinals /= 1000.0;
//	        fprintf(asalida, "PROMEDIO H = %lf\n",tiempofinalh);
//        fprintf(asalida, "PROMEDIO S = %lf\n",tiempofinals);
                fprintf(asalida, "PROMEDIO H = %lf\n",tiempofinalhm);
        fprintf(asalida, "PROMEDIO S = %lf\n",tiempofinalsm);
//return 0;

write_structure<HTML_FORMAT>(*pe,argv[5]);
//return 0;
	primero = true;
	t0=clock();
	antA=0,antB = 0,contadornivel = 0;
	for(int i = 0; i < total;i++){

		
		if(antA != nivelA[i] || antB != nivelB[i]){
	primero = true;
	t1=clock();
	time = (double(t1-t0)/1000);	
	fprintf(asalida,"%d %d %lf %lf %lf\n",antA,antB,time/contadornivel,maximo,minimo);		
	antA = nivelA[i];
	antB = nivelB[i];
	t0=clock();
	contadornivel = 0;
		}


		auxt1=clock();
    if(pe->Touches(pe->getReverso(regionA[i],nivelA[i]),nivelA[i],pe->getReverso(regionB[i],nivelB[i]),nivelB[i]))printf("%d\n",regionB[i]);
    auxt2=clock();
        tactual = double(auxt2-auxt1)/1000.0;
    if(primero){
    	primero = false;
    	maximo = tactual;
    	minimo = maximo;
    }
    else{
    	if(tactual > maximo)maximo = tactual;
    	if(tactual < minimo) minimo = tactual;
    }
    contadornivel++;
	}

	t1=clock();
	time = (double(t1-t0)/1000.0);
	fprintf(asalida,"%d %d %lf %lf %lf\n",antA,antB,time/contadornivel,maximo,minimo);	

		regionA.clear();
		nivelA.clear();
		regionB.clear();
		nivelB.clear();




		fscanf(fp3,"%d",&total);
	for(int i = 0; i < total; i++){
		int a1,a2,a3,a4;
		fscanf(fp3,"%d %d %d",&a1,&a2,&a3);
		regionA.push_back(a1);
		nivelA.push_back(a2);
		nivelB.push_back(a3);
	}


	return 0;
	primero = true;
		antA=niveles-1,antB = niveles-2,contadornivel = 0;
			vector <double> acumuladotiempo(9,0.0);
	vector <double> acumuladocantidad(9,0.0);
	time = 0.0;

	vector <double> minimo2(9,100000000.0);
	vector <double> maximo2(9,0.0); 





	for(int i = niveles-1; i > 0; i--){
		for(int j = i-1; j >= 0; j--){
			vector <int> cantvecinos(9,0);
			for(int k = 1; k < totalnivel[i]; k++){
 			pe->Contained(pe->getReverso(k,i),i,j);
			auxt1=clock();
			vector <int> probar = pe->Contained(pe->getReverso(k,i),i,j);
   			 auxt2=clock();
        	 tactual = (double(auxt2-auxt1)/1000.0)/(double)probar.size();
			time += tactual;
            for(int l = 8; l  >= 0; l--){
    			if(pow(10,l) <= probar.size()){
    				cantvecinos[l]++;
    				acumuladocantidad[l]++;
    				acumuladotiempo[l]+= tactual;
    		    	minimo2[l] = min(tactual,minimo2[l]);
    				maximo2[l] = max(tactual,maximo2[l]);
    				break;
    				}
    			}
			if(primero){
    		primero = false;
    		maximo = tactual;
    		minimo = maximo;
    		}
    		else{
    		if(tactual > maximo)maximo = tactual;
    		if(tactual < minimo) minimo = tactual;
    		}


    contadornivel++;
			}
			for(int m = 0; m < 9; m++)fprintf(anavecinos,"%d ",cantvecinos[m]);
			fprintf(anavecinos,"\n");
			fprintf(asalida,"%d %d %lf %lf %lf\n",i,j,time/contadornivel,maximo,minimo);		
			time = 0.0;
			contadornivel = 0;	
			primero =true;		
		}
	}
	//fprintf(asalida,"%d %d %lf %lf %lf\n",antA,antB,time/contadornivel,maximo,minimo);	
		for(int i = 0; i < 9; i++)	fprintf(asalida,"%lf %lf %lf\n",acumuladotiempo[i]/acumuladocantidad[i],maximo2[i],minimo2[i]);	

	return 0;

}
