#include <bits/stdc++.h>

using namespace std;

int main(int argc, char **argv){
	vector <int> limites;
		FILE *fp = fopen(argv[1], "r");
        int niveles,a,totalactual = 0;
        fscanf(fp,"%d",&niveles);
        int limite = 200;
limites.push_back(0); 
        for(int i = 0 ; i < niveles; i++){
        	fscanf(fp,"%d",&a);
        	totalactual+= a;
            limites.push_back(totalactual); 
        }
   //     limites.push_back(0+4761354); 
     //   limites.push_back(0+4761354 + 2233031); 
     //   limites.push_back(0+4761354 + 2233031 + 33804); 
     //   limites.push_back(0+4761354 + 2233031 + 33804 + 11626); 
     //   limites.push_back(0+4761354 + 2233031 + 33804 + 11626 + 595); 
     //   limites.push_back(0+4761354 + 2233031 + 33804 + 11626 + 595 + 9); 
        int total = (niveles * (niveles-1) )/2;

	printf("%d\n",(total+niveles)*limite);
	srand(time(NULL));
	for(int i = 0; i < niveles; i++){
		for(int j = i; j < niveles; j++){
			for(int k = 0; k < limite; k++){
			int a = rand()% (limites[i+1] - limites[i]-1)+1;
			int b = rand()%(limites[j+1] - limites[j]-1)+1;

			printf("%d %d %d %d\n",a,i,b,j);
			}
		}
	}

/*
	printf("%d\n",total*limite);
	for(int i = 0; i < niveles-1; i++){
		for(int j = i+1; j < niveles; j++){
						for(int k = 0; k < limite; k++){
			int a = rand()% (limites[i+1] - limites[i]);
			int b = rand()%(limites[j+1] - limites[j]);
			printf("%d %d %d %d\n",a,i,b,j);
			}
		}
	}
*/
	limite = 2000;
	printf("%d\n",total*limite);

	for(int i = niveles-1; i > 0; i--){
		for(int j = i- 1; j >= 0; j--){
				for(int k = 0; k < limite; k++){
			int a = rand()% (limites[i+1] - limites[i]-1)+1;
			printf("%d %d %d\n",a,i,j);
			}
		}
	}


	return 0;
}