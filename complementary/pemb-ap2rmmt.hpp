#ifndef INCLUDED_SDSL_PEMB
#define INCLUDED_SDSL_PEMB

#include <sdsl/sdsl_concepts.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/rank_support_v.hpp>
#include <sdsl/select_support_mcl.hpp>
#include <sdsl/wt_helper.hpp>
#include "util.hpp"
#include <sdsl/bp_support_sada.hpp>
#include <algorithm> 
#include <stdexcept>
#include <vector>
#include <stack>
#include <utility>
#include "auxiliar2.hpp"
#include "Vertex.hpp"
#include "Edge.hpp"
#include "Tree.hpp"
#include <vector>
#include "graph4.hpp"
#include <sdsl/sd_vector.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/int_vector.hpp>
#include <iostream>
#include <fstream>
#include <cstdio>
//#include <math.h>
 
//! Namespace for the succinct data structure library.
namespace sdsl

{

//! A class that provides support for planar embeddings
/*! This data structure supports the following operations:
 *   - first
 *   - mate
 *   - next
 *   - vertex
 *   - degree
 *   - face
 *   - list_neighbors
 *  An opening parenthesis in the balanced parentheses sequence B and B_star is
 *  represented by a 1 in the bit_vector and a closing parenthesis by a 0.
---- *
 *  \par References
 *      - Leo Ferres, José Fuentes Sepúlveda, Travis Gagie, Meng He and Gonzalo
 *        Navarro:
 *        Fast and Compact Planar Embeddings
  *        WADS 2017
 *      - Leo Ferres, José Fuentes Sepúlveda, Travis Gagie, Meng He and Gonzalo
 *        Navarro:
 *        Fast and Compact Planar Embeddings
 *        Computational Geometry: Theory and Applications
 *
 */
using namespace std;
  template<class t_bitvector   = bit_vector,
  		class n_bitvector = sd_vector <>,
		class intvector =  int_vector<>,
	   class t_succ_tree    = bp_support_sada<>,
	   class t_rank        = typename t_bitvector::rank_1_type,
	   class t_select1     = typename t_bitvector::select_1_type,
	   class t_select0     = typename t_bitvector::select_0_type,
	   class n_rank        = typename n_bitvector::rank_1_type,
	   class n_select1     = typename n_bitvector::select_1_type,
	   class n_select0     = typename n_bitvector::select_0_type>




class pemb
{
    public:
        typedef int_vector<>::size_type              size_type;
        typedef int_vector<>::value_type             value_type;
        typedef random_access_const_iterator<pemb> const_iterator;
        typedef const_iterator                       iterator;
        typedef t_bitvector                          bit_vector_type;
	typedef intvector 			     intv;
        typedef t_rank                               rank_1_type;
        typedef t_select1                            select_1_type;
        typedef t_select0                            select_0_type;
        typedef n_bitvector                          nbit_vector_type;
        typedef n_rank                               nrank_1_type;
        typedef n_select1                            nselect_1_type;
        typedef n_select0                            nselect_0_type;
        typedef t_succ_tree                          succ_tree;

protected:
        size_type         m_vertices  = 0;
        size_type         m_edges  = 0;
        size_type		  m_levels = 1;
        bit_vector_type *m_A;
        nbit_vector_type *B;
        nbit_vector_type *cB;
        bit_vector_type *BP;
        nbit_vector_type *BF;
        nbit_vector_type *B_F;
        bit_vector_type *m_B;
	bit_vector_type *m_BA;
        rank_1_type       *m_A_rank;
        rank_1_type       *m_B_rank;
        nrank_1_type       *B_rank;
        rank_1_type       *BP_rank;
        nrank_1_type       *BF_rank;
        nrank_1_type       *B_F_rank;
        select_1_type     *m_A_select1;
        select_1_type     *m_B_select1;
        nselect_1_type     *B_select1;
        select_1_type     *BP_select1;
        nselect_1_type     *B_F_select1;
        nselect_1_type     *BF_select1;
        nselect_0_type     *BF_select0;
        select_0_type     *m_A_select0;
        select_0_type     *m_B_select0;
        bit_vector_type    *m_B_star;
        succ_tree          *BF_st;
        succ_tree          *m_B_st;
        succ_tree          *m_B_star_st;
        auxiliar			t;
//	int *plano;
	intv plano;
        wt_hutu_int<rrr_vector<63>>      wt,wt2,wtba,cwtba;
        intv  *resumen,*cresumen,*resumen2,*cresumen2;
	int *mapeo, *reverso;
//	intv *rmt_v;
//	intv *rmt_m;
	vector <vector <int> > rmt_v,rmt_m,rmt_r,rmt_r2,rmt_ar,rmt_ar2,rmt_v2,rmt_m2;
	intv *excess;
	double tiempoh =1000.0,tiempohm= 0.0,contadoreh = 0.0;
	double tiempos = 1000.0,tiemposm = 0.0,contadores =0.0, tiempobloque =0.0;
	vector <int> m_height, m_num_internal, m_num_leaves,posinicial,cm_height, cm_num_internal, cm_num_leaves,cposinicial;
	vector <int>m_aheight, m_num_ainternal, m_num_aleaves,posiniciala;
vector <int> m_heightba, m_num_internalba, m_num_leavesba,posinicialba;
vector <int> cm_heightba, cm_num_internalba, cm_num_leavesba,cposinicialba;
vector <int> ultimos;
//        int32_t           m_height  = 0; // height of the rmMt
  //      uint32_t          m_num_internal  = 0; // height of the rmMt
  //      uint32_t          m_num_leaves  = 0; // Leaves of the rmMt
        uint32_t          m_bs = 128;
        uint32_t          posinicial2 = 0,posinicial3 = 0, posinicial4 = 0,posinicial5 = 0;

        void copy(const pemb& p) {
            m_vertices          = p.m_vertices;
            m_edges         = p.m_edges;
            m_levels = p.m_levels;
            m_A       = p.m_A;
            B       = p.B;
            BP       = p.BP;
            BF       = p.BF;
            B_rank  = p.B_rank;
            BP_rank  = p.BP_rank;
            BF_rank  = p.BF_rank;
            B_select1     = p.B_select1;
            BP_select1     = p.BP_select1;
            B_F       = p.B_F;
            B_F_rank  = p.B_F_rank;
            B_F_select1     = p.B_F_select1;
            BF_select1     = p.BF_select1;
            BF_select0     = p.BF_select0;
            m_B_rank  = p.m_B_rank;
            m_B_select1     = p.m_B_select1;
            m_A_rank  = p.m_A_rank;
            m_A_select1     = p.m_A_select1;
            m_A_select0     = p.m_A_select0;
	    m_A_rank.set_vector(&m_A);
	    B_rank.set_vector(&B);
	    BP_rank.set_vector(&BP);
	    BF_rank.set_vector(&B);
	    B_F_rank.set_vector(&B);
	    B_select1.set_vector(&B);
	    BP_select1.set_vector(&BP);
	    BF_select1.set_vector(&BF);
	    BF_select0.set_vector(&BF);
	    B_F_select1.set_vector(&B_F);
	    m_A_select1.set_vector(&m_A);

	    m_A_select0.set_vector(&m_A);
		m_B_select0.set_vector(&m_B);
	    m_B = p.m_B;
	    m_B_star = p.m_B_star;
	    BF_st = p.BF_st;
	    m_B_st = p.m_B_st;
	    m_B_star_st = p.m_B_star_st;
	    BF_st.set_vector(BF);
	    m_B_st.set_vector(m_B);
	    m_B_star_st.set_vector(m_B_star);
        }


    public:

        //! Default constructor
        pemb() {};

     	pemb(Graph g,unsigned int inicial,char **argv) {
	  

    int total;
    unsigned int n,m,niveles,aux,mitotal = 0,tam = 0,a,b,c;
FILE *fp11 = fopen(argv[2], "r");
    if (!fp11) {
    fprintf(stderr, "Error opening file \"%s\".\n", argv[2]);
    exit(EXIT_FAILURE);
    }
    n = g.vertices();
    vector <int> cantidad,acumulado;
    fscanf(fp11,"%d",&niveles);
    unsigned int *limites = new unsigned int[niveles+1](); 
    unsigned int *totalnivel = new unsigned int[niveles](); 
    unsigned int *jerarquias = new unsigned int[niveles*n](); 
//  cantidad.clear();
    for(int i = 0; i < niveles; i++){
        fscanf(fp11,"%d",&aux);
        cantidad.push_back(aux);
        limites[i] = mitotal;
        mitotal+= aux;
        totalnivel[i] = aux;
    }
    limites[niveles] = mitotal; 
    for(int i = 0; i < n*niveles; i++){
        fscanf(fp11,"%d",&aux);
        jerarquias[i] = aux;
    }
    total = mitotal;
    tam = n;







	  //inicializo todo
	  m_vertices = g.vertices();
	  m_edges = g.edges();
	  m_levels = niveles;
     m_A = new bit_vector_type[niveles]();
    resumen = new intv[niveles];
    resumen2 = new intv[niveles]; 
    cresumen2 = new intv[niveles];
    cresumen = new intv[niveles];
rmt_v.resize(niveles);
rmt_m.resize(niveles);
rmt_v2.resize(niveles);
rmt_m2.resize(niveles);
rmt_r.resize(niveles);
rmt_r2.resize(niveles);
rmt_ar.resize(niveles);
rmt_ar2.resize(niveles);
//	rmt_v = new intv[niveles];
//	rmt_m = new intv[niveles];
	excess = new intv[niveles];
     B = new nbit_vector_type[niveles]();
     cB = new nbit_vector_type[niveles]();
BP = new bit_vector_type[niveles]();
     BF = new nbit_vector_type[niveles]();
      B_F = new nbit_vector_type[niveles]();
m_B = new bit_vector_type[niveles]();
m_B_star = new bit_vector_type[niveles]();
m_A_rank = new rank_1_type[niveles]();
m_A_select1 = new select_1_type[niveles]();
B_rank = new nrank_1_type[niveles]();
BP_rank = new rank_1_type[niveles]();
BF_rank = new nrank_1_type[niveles]();
B_select1 = new nselect_1_type[niveles]();
BP_select1 = new select_1_type[niveles]();
B_F_rank = new nrank_1_type[niveles]();
B_F_select1 = new nselect_1_type[niveles]();
BF_select1 = new nselect_1_type[niveles]();
BF_select0 = new nselect_0_type[niveles]();
m_B_rank = new rank_1_type[niveles]();
m_B_select1 = new select_1_type[niveles]();
m_A_rank = new rank_1_type[niveles]();
m_A_select1 = new select_1_type[niveles]();
m_A_select0 = new select_0_type[niveles]();
m_B_select0 = new select_0_type[niveles]();
BF_st = new succ_tree[niveles]();
m_B_st = new succ_tree[niveles]();
//mBA = new  bit_vector_type[1]();
m_B_star_st = new succ_tree[niveles]();
 m_height.resize(niveles,0);
m_num_internal.resize(niveles,0);
 m_num_leaves.resize(niveles,0);
posinicial.resize(niveles,0);

 m_heightba.resize(niveles,0);
m_num_internalba.resize(niveles,0);
 m_num_leavesba.resize(niveles,0);
posinicialba.resize(niveles,0);

 cm_heightba.resize(niveles,0);
cm_num_internalba.resize(niveles,0);
 cm_num_leavesba.resize(niveles,0);
cposinicialba.resize(niveles,0);

ultimos.resize(niveles,0);

 m_aheight.resize(niveles,0);
m_num_ainternal.resize(niveles,0);
 m_num_aleaves.resize(niveles,0);
posiniciala.resize(niveles,0);
 cm_height.resize(niveles,0);
cm_num_internal.resize(niveles,0);
 cm_num_leaves.resize(niveles,0);
cposinicial.resize(niveles,0);
   mapeo = new int[6*m_vertices];
reverso = new int[6*m_vertices];



	  t = g.dfs_spanning_tree_propio(inicial,jerarquias,limites,niveles,total,m_edges-m_vertices+2);
	  vector <int> corch,parente;
	  vector <vector <char> >ParenteAux =t.getParentesis();
	puts("aun vivo");


	  for(int i = 0; i < niveles; i++){
	  	int aux1 = 0, aux2= 0;
	  	for(int j = 0; j < ParenteAux[i].size();j++){
	  		if(ParenteAux[i][j] == '(' || ParenteAux[i][j] == ')')aux1++;
	  		else aux2++;
	  	}
	  	parente.push_back(aux1);
	  	corch.push_back(aux2);
	  }
	  bit_vector_type B_local(2*m_vertices,0);
	  bit_vector_type B_star_local(2*m_edges-2*m_vertices+4,0);

	  unsigned int *marked_edges = new unsigned
	    int[2*m_edges-2*m_vertices+2]();	  
	  unsigned int idx = 0;
	  unsigned int ii = 0;
	  unsigned int pos = 0;
	  	  	  bool *bits_aux = t.getBits();	
bool *cbits_aux = t.getcBits();
	  	  	  int *cierres_aux = t.getcierres();	


	//construyo bitvectors a partir de sencuencia de parentesis y corchetes
	  for(int i = 0; i < niveles; i++){
	  bit_vector_type A_local(t.getParentesis()[i].size(),0);
	  bit_vector_type B_local(parente[i],0);
	  bit_vector_type B_star_local(corch[i],0);
	  	  	int aux1 = 0, aux2= 0,contadorprueba = 0;

	  	for(int j = 0; j < ParenteAux[i].size();j++){
	  		if(ParenteAux[i][j] == '(' || ParenteAux[i][j] == ')')A_local[j] = 1;
	  		if(ParenteAux[i][j] == '(' ){
	  			B_local[aux1] = 1;
	  			aux1++;
	  			contadorprueba++;
	  		}
	  		if(ParenteAux[i][j] == ')' ){
	  			aux1++;	  			
	  		}
	  		if(ParenteAux[i][j] == '[' ){  		
	  			B_star_local[aux2] = 1;
	  			aux2++;	  			
	  		}
	  		if(ParenteAux[i][j] == ']' ){
	  			aux2++;	  			
	  		}
	  	}
	  	m_A[i].swap(A_local);
	    m_B[i].swap(B_local);
	    m_B_star[i].swap(B_star_local);
	  }
//puts("aun no muero");
	  for(int i = 0; i <  niveles; i++){
	  	bit_vector_type B_l(tam*2,0);
	  	bit_vector_type B_lF(tam*2,0);
	  	for(int j = 0; j < tam*2; j++){
	  		if(bits_aux[j + i*tam*2] == true){
	  			B_l[j] = 1;
	  		}
	  	}
	  		  	nbit_vector_type B_l1(B_l);
	  	nbit_vector_type B_lF1(B_lF);
	  	B[i].swap(B_l1);
	  	B_F[i].swap(B_lF1);
	  }

	  	  for(int i = 0; i <  niveles; i++){
	  	bit_vector_type B_l(tam*2,0);
	  	for(int j = 0; j <tam*2 ; j++){
	  		if(bits_aux[j + i*tam*2] == true){
	  			B_l[j] = 1;
	  		}
	  	}
	  	BP[i].swap(B_l);
	  } 
                  for(int i = 0; i <  niveles; i++){
                bit_vector_type B_l(( m_edges-m_vertices+2)*2,0);
                for(int j = 0; j <  ( m_edges-m_vertices+2)*2; j++){
                        if(cbits_aux[j + i*( m_edges-m_vertices+2)*2] == true){
                                B_l[j] = 1;
                        }
                }
B_l[( m_edges-m_vertices+2)*2-1] = 1;
nbit_vector_type B_l1(B_l);
                cB[i].swap(B_l1);
          }

//cB[0][ ( m_edges-m_vertices+2)*2 -1] = 1;
//cB[1][( m_edges-m_vertices+2)*2 -1] = 1;
//cB[2][( m_edges-m_vertices+2)*2 -1] = 1;
//cB[3][( m_edges-m_vertices+2)*2 -1] = 1;
//cB[4][( m_edges-m_vertices+2)*2 -1] = 1;
//cB[5][( m_edges-m_vertices+2)*2 -1] = 1;





/*

ofstream myfile ("ofCB2");
for(int i = 0; i<  6; i++){
myfile << cB[i].size() << endl;
for(int j = 0; j < cB[i].size();j++){
int salida = 0;
if(cB[i][j])salida = 1;
myfile << salida << " ";
}
}
myfile.close();

*/


/*
ifstream myfile10 ("ofmapeo2");
for(int i = 0; i<  6; i++){
 for(int j = 0; j < m_vertices;j++){
int salida = 0;
 myfile10 >> salida;
 mapeo[i*m_vertices + j] = salida;
}                                                                                                                                                                                                                  }
myfile10.close();

ifstream myfile11 ("ofreverso2");
for(int i = 0; i<  6; i++){
 for(int j = 0; j < m_vertices;j++){
int salida = 0;
 myfile11 >> salida;
  reverso[i*m_vertices + j] = salida;
}
}
myfile11.close();


ifstream myfile ("ofB2");
for(int i = 0; i<  6; i++) {
int stam = 0;
myfile >> stam;
cout << stam << endl;
 bit_vector_type A_local(stam,0);
for(int j = 0; j < stam;j++){
int salida = 0;
myfile >> salida;
//if(salida != B[i][j]) cout << salida << " " << B[i][j]  << endl;
//cout << salida <<endl;
A_local[j] = salida;
}                                                                                                                                                                                                                  nbit_vector_type B_l1(A_local);
B[i].swap(B_l1);
}
myfile.close();

ifstream myfile20 ("ofCB2");
for(int i = 0; i<  6; i++) {
int stam = 0;
myfile20 >> stam;
cout << stam << endl;
 bit_vector_type A_local(stam,0);
for(int j = 0; j < stam;j++){
int salida = 0;
myfile20 >> salida;
//if(salida != B[i][j]) cout << salida << " " << B[i][j]  << endl;
//cout << salida <<endl;
A_local[j] = salida;
}                                                                                                                                                                                                                  nbit_vector_type B_l1(A_local);
cB[i].swap(B_l1);
}
myfile20.close();


ifstream myfile2 ("ofBP2");
for(int i = 0; i < 6; i++) {
int stam = 0;
myfile2 >> stam;
 bit_vector_type A_local(stam,0);
for(int j = 0; j < stam;j++) {
int salida = 0;
myfile2 >> salida;
//if(salida != BP[i][j]) cout << salida << " " << BP[i][j]  << endl;
A_local[j] = salida;
}
 BP[i].swap(A_local);
}
myfile2.close();


ifstream myfile3 ("ofm_A2");
for(int i = 0; i < 6; i++) {
int stam = 0;
myfile3 >> stam;
 bit_vector_type A_local(stam,0);
for(int j = 0; j < stam;j++) {
int salida = 0;
myfile3 >> salida;
//if(salida != m_A[i][j]) cout << salida << " " << m_A[i][j]  << endl;
A_local[j] = salida;
}
 m_A[i].swap(A_local);
}
myfile3.close();

ifstream myfile4 ("ofm_B2");
for(int i = 0; i < 6; i++) {
int stam = 0;
myfile4 >> stam;
 bit_vector_type A_local(stam,0);
for(int j = 0; j < stam;j++) {
int salida = 0;
myfile4 >> salida;
//if(salida != m_B[i][j]) cout << salida << " " << m_B[i][j]  << endl;
A_local[j] = salida;                                                                                                                                                                                               }
 m_B[i].swap(A_local);
}
myfile4.close();

ifstream myfile5 ("ofm_B_star2");
for(int i = 0; i < 6; i++) {
int stam = 0;
myfile5 >> stam;
 bit_vector_type A_local(stam,0);
for(int j = 0; j < stam;j++) {
int salida = 0;
myfile5 >> salida;
//if(salida != m_B_star[i][j]) cout << salida << " " << m_B_star[i][j]  << endl;
A_local[j] = salida;
}
 m_B_star[i].swap(A_local);
}
myfile5.close();
*/


/*
for(int i = 0; i < niveles; i++){
	int contpos = 0;
	bit_vector_type BF_l(parente[i],0);
	for(int j = 0; j < tam; j++){
		if(B[i][j]==1){	
			BF_l[contpos] = 1;
			contpos++;
		}
		if(B_F[i][j]==1){
			contpos++;
		}
	}
	nbit_vector_type BF_l1(BF_l);
	BF[i].swap(BF_l1);
}



*/

vector <int>vaux;
for(int i = 0; i < tam*2; i++){
bool entro =false;
for(int j = niveles-1; j > 1; j--){
	if(B[j][i]){
		entro = true;
		vaux.push_back(j);
		break;
	}
}
if(!entro){
	if(BP[1][i])vaux.push_back(1);
	else  vaux.push_back(0);
}
}
//puts("VIVO1");
/*
cout << "FINAL " << cB[0][( m_edges-m_vertices+2)*2 -1 ]  << endl;
for(int i = 0; i < 6; i++){
int sumacor = 0;

for(int j = 0; j < ( m_edges-m_vertices+2)*2; j++){

if(cB[i][j])sumacor++;
}

cout << "NIVEL " << i << " TIENE " << sumacor << " " << m_B_star[i].size() << endl;

}

*/


vector <int>cvaux;
for(int i = 0; i < ( m_edges-m_vertices+2)*2; i++){
bool entro =false;
for(int j = niveles-1; j > 1; j--){
        if(cB[j][i]){
                entro = true;
                cvaux.push_back(j);
                break;
        }
}
if(!entro){
        if(cB[1][i])cvaux.push_back(1);
        else  cvaux.push_back(0);
}

}




/*
	for(int i = 0 ; i < 40; i++)cout << m_B_star[1][i] << " ";
puts("");

        for(int i = 0 ; i < 40; i++)cout << cB[1][i] << " ";
puts("");

        for(int i = 0 ; i < 40; i++)cout << cvaux[i] << " ";
puts("");

cout << m_B_star[0].size()<< "  "<< m_B_star[1].size() << endl;
*/
int_vector <> secuencia(vaux.size());
plano.resize(vaux.size()); //= new vint	[vaux.size()];
for(int i = 0; i < vaux.size();i++){
secuencia[i] = vaux[i];
plano[i] = vaux[i];
}
construct_im(wt,secuencia);

//instancio los rank,select y sst


int_vector <> csecuencia(cvaux.size());
//plano.resize(vaux.size()); //= new vint [vaux.size()];
for(int i = 0; i < cvaux.size();i++){
csecuencia[i] = cvaux[i];
//plano[i] = cvaux[i];
}

/*
int malos = 0;
for(int  i = 0 ; i < csecuencia.size(); i++){
int actual = 1;
for(int j = 0 ; j < 6; j++){
if(actual == 0 && cB[j][i] == 1){
if(malos == 0) cout << "PRIMER MALO = " << i  << " " << j << endl;
malos++;
break;
}
actual = cB[j][i];
}

}
cout << "MALOS = " << malos << endl;
*/
construct_im(wt2,csecuencia);

/*
cout<< "COMIENZO WT2" << endl;
for(int i = 0; i < 20; i++)cout << wt2[i] << " ";
cout << endl;

int totales1 = 0;
for(int i = 0; i < csecuencia.size();i++)if(csecuencia[i] >= 1)totales1++;
cout << "totales 1 = " << totales1 << " " << m_B_star[1].size()  << endl;
totales1 = 0;
for(int i = 0; i < csecuencia.size();i++)if(csecuencia[i] >= 4)totales1++;
cout << "totales 1 = " << totales1 << " " << m_B_star[4].size()  << endl;

totales1 = 0;
for(int i = 0; i < csecuencia.size();i++)if(csecuencia[i] >= 0)totales1++;
cout << "totales 0 = " << totales1 << " " << m_B_star[0].size()  << endl;

*/
//cout << "PESOS TOTALES = "<< secuencia.size() << " " << csecuencia.size() << endl;

/*
for(int i = 0 ; i < niveles; i++){for(int j = 0; j < m_B_star[i].size();j++)cout << m_B_star[i][j] << " ";
cout << endl;

}
*/


int_vector <> secuencia3(wt.size(),0);
vector <int> tempcont(niveles,0);
int pruebaentro = 0;


for(int i = 0; i < wt.size(); i++){
//if(wt[i] >= 4)pruebaentro++;
if(m_B[wt[i]][ tempcont[wt[i]] ] == 1){
 secuencia3[i] = wt[i];
//for(int j = wt[i]; j >=0;j --)tempcont[j]++;
}
for(int j = wt[i]; j >=0;j --)tempcont[j]++;
}



construct_im(wtba,secuencia3);












int_vector <> secuencia4(wt2.size(),0);
vector <int> tempcont2(niveles,0);
int pruebaentro2 = 0;


for(int i = 0; i < wt2.size(); i++){
//if(wt[i] >= 4)pruebaentro++;
if(m_B_star[wt2[i]][ tempcont2[wt2[i]] ] == 1){
 secuencia4[i] = wt2[i];
//for(int j = wt[i]; j >=0;j --)tempcont[j]++;
}
for(int j = wt2[i]; j >=0;j --)tempcont2[j]++;
}



construct_im(cwtba,secuencia4);


//ultimos[i] = m_A_select1[nivel](pos+1)


cout << "total 4 " << tempcont[4] <<" " <<  pruebaentro<< " "  << m_B[4].size()  <<endl;

for(int i = 0; i < niveles; i++){
	  util::init_support(m_A_rank[i], &m_A[i]);
	  util::init_support(m_B_rank[i], &m_B[i]);
	  //  puts("ENTRE");
 	  util::init_support(B_rank[i], &B[i]);
 	  util::init_support(BP_rank[i], &BP[i]);
// 	  util::init_support(BF_rank[i], &BF[i]);
// 	  util::init_support(B_F_rank[i], &B_F[i]);
 	   // puts("ENTRE");
 	  util::init_support(m_A_select1[i], &m_A[i]);
  	  util::init_support(m_B_select1[i], &m_B[i]);
 	  util::init_support(B_select1[i], &B[i]); 
 	  util::init_support(BP_select1[i], &BP[i]); 
 //	  util::init_support(B_F_select1[i], &B_F[i]); 
 //	  util::init_support(BF_select1[i], &BF[i]); 
 //	  util::init_support(BF_select0[i], &BF[i]); 
	  util::init_support(m_A_select0[i], &m_A[i]);
	  util::init_support(m_B_select0[i], &m_B[i]);
	//  succ_tree B_F_local_st(&BF[i]);
	  succ_tree B_local_st(&m_B[i]);
	  succ_tree B_star_local_st(&m_B_star[i]);
	  //BF_st[i].swap(B_F_local_st);
	  //BF_st[i].set_vector(&BF[i]);
	  m_B_st[i].swap(B_local_st);
	  m_B_st[i].set_vector(&m_B[i]);
	  m_B_star_st[i].swap(B_star_local_st);
	  m_B_star_st[i].set_vector(&m_B_star[i]);
}
/*
for(int i = 1; i < niveles; i++){
	for(int j = 0; j < parente[i]; j++){
		cout << BF[i][j] << endl;
	}
	puts("");
}
*/

for(int i = 0; i < niveles; i++){
//      for(int j = wt.size()-1; j >=0 ; j--){
//              if(wt[j]>= i){
//                      ultimos[i] = j;
//                      break;
//              }
//      }
ultimos[i] = m_A_select1[i](m_B[i].size());
}
        }

        //! Copy constructor
        pemb(const pemb& g) {
            copy(g);
        }

        //! Copy constructor
        pemb(pemb&& g) {
            *this = std::move(g);
        }

        //! Assignment operator
        pemb& operator=(const pemb g) {
            if (this != &g) {
                copy(g);
            }
            return *this;
        }

        //! Assignment move operator
        pemb& operator=(pemb&& g) {
            if (this != &g) {
                m_vertices          = g.m_vertices;
                m_edges         = g.m_edges;
                m_levels = g.m_levels;
		
		m_A             = std::move(g.m_A);
		B       = std::move(g.B);
				BP       = std::move(g.BP);
		BF       = std::move(g.BF);
		B_F       = std::move(g.B_F);
		B_rank       = std::move(g.B_rank);
		BF_rank       = std::move(g.BF_rank);
		B_select1       = std::move(g.B_select1);
		B_F_rank       = std::move(g.B_F_rank);
		B_F_select1       = std::move(g.B_F_select1);
		BF_select1       = std::move(g.BF_select1);
		BF_select0       = std::move(g.BF_select0);
		m_A_rank        = std::move(g.m_A_rank);
		m_A_select1     = std::move(g.m_A_select1);
		m_B_rank        = std::move(g.m_B_rank);
		m_B_select1     = std::move(g.m_B_select1);
		m_A_select0     = std::move(g.m_A_select0);
		m_B_select0     = std::move(g.m_B_select0);


		B_rank.set_vector(&B);
		B_select1.set_vector(&B);
		BP_rank.set_vector(&BP);
		BP_select1.set_vector(&BP);
		B_F_rank.set_vector(&BF);
		B_F_select1.set_vector(&BF);
		BF_rank.set_vector(&BF);
		BF_select1.set_vector(&BF);
		BF_select0.set_vector(&BF);
		m_A_rank.set_vector(&m_A);
		m_A_select1.set_vector(&m_A);
		m_B_rank.set_vector(&m_B);
		m_B_select1.set_vector(&m_B);
		m_A_select0.set_vector(&m_A);
		m_B_select0.set_vector(&m_B);

		m_B             = std::move(g.m_B);
		m_B_star        = std::move(g.m_B_star);

		BF_st          = std::move(g.BF_st);
		BF_st.set_vector(&BF);

		m_B_st          = std::move(g.m_B_st);
		m_B_st.set_vector(&m_B);

		m_B_star_st          = std::move(g.m_B_star_st);
		m_B_star_st.set_vector(&m_B_star);
            }
            return *this;
        }

        //! Swap operator
        void swap(pemb& g) {
            if (this != &g) {
                std::swap(m_vertices, g.m_vertices);
                std::swap(m_levels, g.m_levels);
                std::swap(m_edges,  g.m_edges);
		
		m_A.swap(g.m_A);
                util::swap_support(m_A_rank, g.m_A_rank, &m_A, &(g.m_A));
                util::swap_support(m_A_select1, g.m_A_select1, &m_A, &(g.m_A));
                util::swap_support(m_A_select0, g.m_A_select0, &m_A, &(g.m_A));
		B.swap(g.B);
                util::swap_support(B_rank, g.B_rank, &B, &(g.B));
                util::swap_support(B_select1, g.B_select1, &B, &(g.B));

		BP.swap(g.BP);
                util::swap_support(BP_rank, g.BP_rank, &BP, &(g.BP));
                util::swap_support(BP_select1, g.BP_select1, &BP, &(g.BP));  
        BF.swap(g.BF);
        		util::swap_support(BF_rank, g.BF_rank, &BF, &(g.BF));
        		util::swap_support(BF_select1, g.BF_select1, &BF, &(g.BF));
        		util::swap_support(BF_select0, g.BF_select0, &BF, &(g.BF));
        		util::swap_support(BF_st, g.BF_st, &BF, &(g.BF));

		B_F.swap(g.B_F);
                util::swap_support(B_F_rank, g.B_F_rank, &B_F, &(g.B_F));
                util::swap_support(B_F_select1, g.B_F_select1, &B_F, &(g.B_F));

		m_B.swap(g.m_B);
		        util::swap_support(m_B_rank, g.m_B_rank, &m_B, &(g.m_B));
                util::swap_support(m_B_select1, g.m_B_select1, &m_B, &(g.m_B));
                util::swap_support(m_B_st, g.m_B_st, &m_B, &(g.m_B));
                util::swap_support(m_B_select0, g.m_B_select0, &m_B, &(g.m_B));

		m_B_star.swap(g.m_B_star);
                util::swap_support(m_B_star_st, g.m_B_star_st, &m_B_star,
				   &(g.m_B_star));
            }
        }

        //! Returns the number of vertices of the graph
        size_type vertices()const {
            return m_vertices;
        }

        //! Returns the number of edges of the graph
        size_type edges()const {
            return m_edges;
        }
  
        //! Returns a const_iterator to the first vertex.
        const_iterator begin()const {
            return const_iterator(this, 0);
        }

        //! Returns a const_iterator to the element after the last vertex.
        const_iterator end()const {
            return const_iterator(this, vertices());
        }


        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
	  structure_tree_node* child = structure_tree::add_child(v, name,
								 util::class_name(*this));

	  size_type written_bytes = 0;
	  written_bytes += write_member(m_vertices, out, child, "vertices");
	  written_bytes += write_member(m_edges, out, child, "edges");
	  written_bytes += write_member(m_levels, out, child, "levels");
	 double peso = 0.0;
	  for(int i = 0; i < m_levels;i++){
	if(i != 0 && i != m_levels-1)	peso += size_in_mega_bytes(resumen[i]);
	  written_bytes += m_A[i].serialize(out, child, "A");
	  written_bytes += m_A_rank[i].serialize(out, child, "A");
	  written_bytes += m_A_select1[i].serialize(out, child, "A");
	  written_bytes += m_A_select0[i].serialize(out, child, "A");
	  if(i > 1){
	  written_bytes += B[i].serialize(out, child, "B");
	  written_bytes += B_rank[i].serialize(out, child, "B");
	  written_bytes += B_select1[i].serialize(out, child, "B");
	  }
	  else if(i == 1){
	  written_bytes += BP[i].serialize(out, child, "B");
	  written_bytes += BP_rank[i].serialize(out, child, "B");
	  written_bytes += BP_select1[i].serialize(out, child, "B");
          written_bytes += m_B[i-1].serialize(out, child, "mB2");	  
	  written_bytes += m_B_rank[i-1].serialize(out, child, "mB2");
          written_bytes += m_B_select1[i-1].serialize(out, child, "mB2");
  //        written_bytes += m_B_st[i-1].serialize(out, child, "mB2");
          written_bytes += m_B_select0[i-1].serialize(out, child, "mB2");
          written_bytes += m_B_star[i-1].serialize(out, child, "mB2_star");
 written_bytes += m_A[i-1].serialize(out, child, "A0");
//    written_bytes += m_B_star_st[i-1].serialize(out, child, "mB2_star");   
	  }
	  written_bytes += m_B[i].serialize(out, child, "mB");
	  written_bytes += m_B_rank[i].serialize(out, child, "mB");
	  written_bytes += m_B_select1[i].serialize(out, child, "mB");
	  written_bytes += m_B_st[i].serialize(out, child, "mB_succ_tree");
	  written_bytes += m_B_select0[i].serialize(out, child, "mB");
	  
	  written_bytes += m_B_star[i].serialize(out, child, "mB_star");
	  written_bytes += m_B_star_st[i].serialize(out, child, "mB_star_succ_tree");
	  }
//	cout << "PESO TOTAL " << peso << endl;
	cout << "PESO WT = " << size_in_mega_bytes(wt) <<endl;
 cout << "PESO WT2 = " << size_in_mega_bytes(wt2) <<endl;
 cout << "PESO B0 = " << size_in_mega_bytes(BP[0]) <<endl;
 cout << "PESO B0 = " << size_in_mega_bytes(B[0]) <<endl;
	  structure_tree::add_size(child, written_bytes);
	  return written_bytes;
        }

        //! Loads the data structure from the given istream.




        void load(std::istream& in) {
            read_member(m_vertices, in);
            read_member(m_edges, in);
             read_member(m_levels, in);
	    m_A.load(in);
	    m_A_rank.load(in, &m_A);
	    m_A_select1.load(in, &m_A);
	    m_A_select0.load(in, &m_A);
	    B.load(in);
	    	    BP.load(in);
	    BF.load(in);
	    B_rank.load(in, &B);
	    BP_rank.load(in, &BP);
	    BF_rank.load(in, &BF);
	    BF_st.load(in, &BF);
	    BF_select1.load(in, &BF);
	    BF_select0.load(in, &BF);
	    B_select1.load(in, &B);
BP_select1.load(in, &BP);
	   	B_F.load(in);
	    B_F_rank.load(in, &B_F);
	    B_F_select1.load(in, &B_F);


            m_B.load(in);
            m_B_st.load(in, &m_B);
            m_B_select0.load(in, &m_B);
	    
            m_B_star.load(in);
            m_B_star_st.load(in, &m_B_star);
        }

       /* Assuming indices start with 0 */
        //primera ocurrencia del nodo en la secuencia de parentesis y corchetes



size_type mBselect(size_type v, int nivel){

if(nivel == 0)return m_B_select1[0](v);
else return rankwt(selectwtba(nivel,v),nivel);

}


    size_type    first2(size_type v, int nivel) {
	 if(v >= 0) {
	   size_type pos = m_B_select1[nivel](v+1);
	   size_type pos2 = selectmb(nivel,v+1);
	   if(pos != mBselect(v+1,nivel))cout << "select " <<  pos << " " << selectwtba(nivel,v+1) <<" " << selectwt(nivel,v+1)  <<  " "  <<  rankwt(selectwt(nivel,v+2),nivel)   <<
             " " <<  mBselect(v+1,nivel)<< " "  << nivel     <<endl;
	   size_type edge = 1;
	    edge = m_A_select1[nivel](pos+1);
	  //   return edge;
//cout << m_A_rank[nivel](edge+1) << " " << edge << " " << pos << endl;
	    return pos;
	 } else
	   return -1;
       }


         size_type       first(size_type v, int nivel) {
	 if(v >= 0) {
	   size_type pos = m_B_select1[nivel](v+1);
	   size_type edge = 1;
	    edge = m_A_select1[nivel](pos+1);
	     return edge;
	    return pos;
	 } else
	   return -1;
       }

      /* Assuming indices start with 0 */
       //dado una posicion con parentesis que abre me da la posicion donde cierra, es circular

       void iniresumen(int32_t bs){
//cout<< "COMIENZO WT" << endl;
for(int i = 0; i < 20; i++)cout << wt[i] << " ";
cout << endl;
	int valor;
	cout << "tam = " << m_B[0].size()<<endl;
	for(int l = 0; l < m_levels; l++)resumen[l].width(32);
       	m_bs = bs;
       	m_num_leaves[0] = ceil((double)B[0].size()/m_bs);
        m_height[0] = ceil(log(m_num_leaves[0])/log(2));
        m_num_internal[0] =(pow(2,m_height[0]+1)-1);
        for(int l = 0; l < m_levels; l++)resumen[l].resize(m_num_internal[0]);

	for(int l = 0; l < m_levels; l++)cresumen[l].width(32);
       	cm_num_leaves[0] = ceil((double)cB[0].size()/m_bs);
        cm_height[0] = ceil(log(cm_num_leaves[0])/log(2));
        cm_num_internal[0] =(pow(2,cm_height[0]+1)-1);
        for(int l = 0; l < m_levels; l++)cresumen[l].resize(cm_num_internal[0]);



        for(int l = 0; l < m_levels; l++)resumen2[l].width(32);
	m_num_leavesba[0] = ceil((double)B[0].size()/m_bs);
	m_heightba[0] = ceil(log(m_num_leavesba[0])/log(2));
	m_num_internalba[0] =(pow(2,m_heightba[0]+1)-1);
	for(int l = 0; l < m_levels; l++)resumen2[l].resize(m_num_internalba[0]); 



        for(int l = 0; l < m_levels; l++)cresumen2[l].width(32);
        cm_num_leavesba[0] = ceil((double)cB[0].size()/m_bs);
        cm_heightba[0] = ceil(log(cm_num_leavesba[0])/log(2));
cm_num_internalba[0] =(pow(2,cm_heightba[0]+1)-1);
for(int l = 0; l < m_levels; l++)cresumen2[l].resize(cm_num_internalba[0]);


//	cout << "TAMAÑOS RESUMENES = " << m_num_internal[0] << " " << cm_num_internal[0] <<  "  "  <<  B[0].size() << "  " <<  cB[0].size() <<" " << m_B_star[0].size() <<endl;
	for(int l = 0; l < m_levels;l++){
		int actual = 0;
		excess[l].resize(m_B[l].size());
	for(int j = 0; j < m_B[l].size();j++){
			if(m_B[l][j]==1)actual++;
			else actual--;
			excess[l][j] = actual; 
		}

}
//	puts("vivo1");
	for(int l = 0; l < m_levels; l++)util::set_to_value(resumen[l], 0);
for(int l = 0; l < m_levels; l++)util::set_to_value(cresumen[l], 0);
for(int l = 0; l < m_levels; l++)util::set_to_value(resumen2[l], 0);
for(int l = 0; l < m_levels; l++)util::set_to_value(cresumen2[l], 0);
//	cout << m_num_leaves << "  " << m_height << "  " <<  m_num_internal << "  " <<  B[0].size()  <<endl;
	int comprobar = 0;





	for(int l = 0; l < m_levels; l++){
	valor = l;
        for(int i = 0; i < m_num_leaves[0]; i++){
        	int pos = i*m_bs;
//		if(i == m_num_leaves - 1)cout << pos << endl;
       	//	if(i == 0)posinicial = pos;
        	int cont = 0;
        	for(int j = 0; j <  m_bs && pos +j < B[0].size();j++){
        		if(wt[pos+j] >= valor)cont++;
        	}
        	int top_nodes = pow(2, m_height[0])-1;
		if(i == valor)posinicial2 = top_nodes;
//		if(i == m_num_leaves-1)cout << top_nodes+i<<endl;
        	resumen[l][top_nodes +i] = cont;
		comprobar += cont;
        }
//	cout << "comprobar = " << comprobar << endl;
//	puts("vivo2");
        for(int lvl=m_height[0]-1; lvl >= 0 ; lvl--){
	    int num_curr_nodes = pow(2, lvl);
	    int node = 0;
//		cout << "ENTRE A NIVEL " << lvl << endl;
	    	for( ;node < num_curr_nodes  ;node++){
		      int pos = (pow(2,lvl)-1) + node;
		     // if(resumen[pos] != 0)cout<< "SOBREESCRIBI" << endl; 
		      resumen[l][pos] = (resumen[l][2*pos+1] + resumen[l][2*pos+2]);
//			if(lvl == height-2)
	    	}
       	}
}





	for(int l = 0; l < m_levels; l++){
	valor = l;
        for(int i = 0; i < cm_num_leaves[0]; i++){
        	int pos = i*m_bs;
        	int cont = 0;
        	for(int j = 0; j <  m_bs && pos +j < cB[0].size();j++){
        		if(wt2[pos+j] >= valor)cont++;
        	}
        	int top_nodes = pow(2, cm_height[0])-1;
		if(i == valor)posinicial3 = top_nodes;
			if(i == 0 && l == 1)cout << "PRIMER BLOQUE TIENE " << cont << endl;
        	cresumen[l][top_nodes +i] = cont;
		comprobar += cont;
        }
        for(int lvl=cm_height[0]-1; lvl >= 0 ; lvl--){
	    int num_curr_nodes = pow(2, lvl);
	    int node = 0;
	    	for( ;node < num_curr_nodes  ;node++){
		      int pos = (pow(2,lvl)-1) + node;
		      cresumen[l][pos] = (cresumen[l][2*pos+1] + cresumen[l][2*pos+2]);
	    	}
       	}
}






        for(int l = 0; l < m_levels; l++){
        valor = l;
        for(int i = 0; i < m_num_leavesba[0]; i++){
                int pos = i*m_bs;
                int cont = 0;
                for(int j = 0; j <  m_bs && pos +j < B[0].size();j++){
                        if(wtba[pos+j] >= valor)cont++;
                }
                int top_nodes = pow(2, m_heightba[0])-1;
                if(i == valor)posinicial4 = top_nodes;
 //                       if(i == 0 && l == 1)cout << "PRIMER BLOQUE TIENE " << cont << endl;
                resumen2[l][top_nodes +i] = cont;
//		if(l == 4)cout << "total: " << i << " " << cont << endl;
//                comprobar += cont;
        }
                          for(int lvl=m_heightba[0]-1; lvl >= 0 ; lvl--){
            int num_curr_nodes = pow(2, lvl);
            int node = 0;
                for( ;node < num_curr_nodes  ;node++){
                      int pos = (pow(2,lvl)-1) + node;
                      resumen2[l][pos] = (resumen2[l][2*pos+1] + resumen2[l][2*pos+2]);
                }
        }
}





        for(int l = 0; l < m_levels; l++){
        valor = l;
        for(int i = 0; i < cm_num_leavesba[0]; i++){
                int pos = i*m_bs;
                int cont = 0;
                for(int j = 0; j <  m_bs && pos +j < cB[0].size();j++){
                        if(cwtba[pos+j] >= valor)cont++;
                }
                int top_nodes = pow(2, cm_heightba[0])-1;
                if(i == valor)posinicial5 = top_nodes;
 //                       if(i == 0 && l == 1)cout << "PRIMER BLOQUE TIENE " << cont << endl;
                cresumen2[l][top_nodes +i] = cont;
//              if(l == 4)cout << "total: " << i << " " << cont << endl;
//                comprobar += cont;
        }
                          for(int lvl=cm_heightba[0]-1; lvl >= 0 ; lvl--){
            int num_curr_nodes = pow(2, lvl);
            int node = 0;
                for( ;node < num_curr_nodes  ;node++){
                      int pos = (pow(2,lvl)-1) + node;
                      cresumen2[l][pos] = (cresumen2[l][2*pos+1] + cresumen2[l][2*pos+2]);
                }
        }
}





//COMIENZO PARTE DEL RMT, SI FALLA REVISAR TAM DE RMT, REVISAR SI CODIGO ANTERIOR CONSIDERA A LAS HOJAS O NO.





        for(int l = 0; l < m_levels; l++){
        valor = l;
        m_bs = bs;
        cm_num_leaves[l] = ceil((double)m_B_star[l].size()/m_bs);
        cm_height[l] = ceil(log(cm_num_leaves[l])/log(2));
        cm_num_internal[l] =(pow(2,cm_height[l]+1)-1);
	rmt_v2[l].resize(cm_num_internal[l]);
	rmt_m2[l].resize(cm_num_internal[l]);
	for(int ll = 0; ll < cm_num_internal[l];ll++){
	rmt_v2[l][ll] = 0;
	rmt_m2[l][ll] = 0;
}
	int auxcompro= 0;
	int rmb = 0;
        for(int i = 0; i < cm_num_leaves[l]; i++){
                int pos = i*m_bs;
                int cont = 0,vmin = 1,cont2 = 0;
                for(int j = 0; j <  m_bs && pos +j < m_B_star[l].size();j++){
			m_B_star[l][pos+j] == 1? cont++: cont--;
			if(m_B_star[l][pos+j])rmb++;
			if(m_B_star[l][pos+j])cont2++;
			if(vmin > cont) vmin = cont;
                }
                int top_nodes = pow(2, cm_height[l])-1;
                if(i == valor)cposinicial[l] = top_nodes;
		rmt_v2[l][top_nodes + i] = cont;
		rmt_m2[l][top_nodes + i] = vmin;
                comprobar += cont;
        }
        for(int lvl=cm_height[l]-1; lvl >= 0 ; lvl--){
            int num_curr_nodes = pow(2, lvl);
            int node = 0;
                for( ;node < num_curr_nodes  ;node++){
                      int pos = (pow(2,lvl)-1) + node;
			rmt_v2[l][pos] = rmt_v2[l][2*pos +1] + rmt_v2[l][2*pos + 2];
			int auxpos1 = rmt_m2[l][2*pos +1];
			int auxpos2 = rmt_v2[l][2*pos+1]+rmt_m2[l][2*pos+2];
			rmt_m2[l][pos] = min(auxpos1,auxpos2);
                }
        }
}














	cout << m_num_leaves[0] << " contra   " << ceil((double)m_B[0].size()/m_bs) << endl;

        for(int l = 0; l < m_levels; l++){
	cout << "ENTRE CON NIVEL " << l << endl;
        valor = l;
        m_bs = bs;
        m_num_leaves[l] = ceil((double)m_B[l].size()/m_bs);
        m_height[l] = ceil(log(m_num_leaves[l])/log(2));
        m_num_internal[l] =(pow(2,m_height[l]+1)-1);
	rmt_v[l].resize(m_num_internal[l]);
	rmt_m[l].resize(m_num_internal[l]);
	rmt_r[l].resize(m_num_leaves[l]+1);
	rmt_r2[l].resize(m_num_internal[l]);
	for(int ll = 0; ll < m_num_internal[l];ll++){
	rmt_v[l][ll] = 0;
	rmt_m[l][ll] = 0;
	rmt_r2[l][ll] = 0;
//	rmt_r[l][ll] = 0;
}
	int auxcompro= 0;
//	for(int kk = 0; kk < m_B[l].size();kk++)m_B[l][kk] == 1? auxcompro++: auxcompro--;
//	 cout << "COMPROBAR = " << auxcompro << endl;
	int rmb = 0;
//	cout << "NUMEROS = "<< m_num_leaves[l] << "  " << m_B[l].size() << endl;
	cout << "TOP = " << pow(2, m_height[l])-1 <<"  "  <<  m_num_leaves[l] <<endl;
	rmt_r[l][0] = 0;
        for(int i = 0; i < m_num_leaves[l]; i++){
                int pos = i*m_bs;
//              if(i == m_num_leaves - 1)cout << pos << endl;
        //      if(i == 0)posinicial = pos;
                int cont = 0,vmin = 1,cont2 = 0;
                for(int j = 0; j <  m_bs && pos +j < m_B[l].size();j++){
                       // if(excess[l][j] < cont)cont = excess[l][j];
			m_B[l][pos+j] == 1? cont++: cont--;
			if(m_B[l][pos+j])rmb++;
			if(m_B[l][pos+j])cont2++;
 			//m_B[l][pos+j] == 1? auxcompro++: auxcompro--;
			if(vmin > cont) vmin = cont;
                }

		rmt_r[l][i+1] = rmb;
                int top_nodes = pow(2, m_height[l])-1;
		rmt_r2[l][top_nodes + i] = cont2;
                if(i == valor)posinicial[l] = top_nodes;
//              if(i == m_num_leaves-1)cout << top_nodes+i<<endl;
   //             rmt_m[l][top_nodes +i] = cont;
		int resta;
		i == 0? resta = 0: resta = excess[l][ (i-1)*m_bs +m_bs -1 ];
		rmt_v[l][top_nodes + i] = cont;
//		if(top_nodes + i == 65536) cout <<"REVISAR : " << rmt_v[l][top_nodes + i] << "  " <<   cont  <<endl;
//		if(rmt_v[l][top_nodes+ i] != cont) cout <<"EXCESS DISTINTO : " << cont  <<"  " << rmt_v[l][top_nodes+i] <<"  "  << i << "  " <<   m_num_leaves[l] <<endl;
	//	auxcompro += rmt_v[l][top_nodes+i];
		rmt_m[l][top_nodes + i] = vmin;
                comprobar += cont;
        }
//	cout << "COMPROBAR = " << auxcompro << endl;
//      cout << "comprobar = " << comprobar << endl;
//      puts("vivo2");
        for(int lvl=m_height[l]-1; lvl >= 0 ; lvl--){
            int num_curr_nodes = pow(2, lvl);
            int node = 0;
//              cout << "ENTRE A NIVEL " << lvl << endl;
                for( ;node < num_curr_nodes  ;node++){
                      int pos = (pow(2,lvl)-1) + node;
                     // if(resumen[pos] != 0)cout<< "SOBREESCRIBI" << endl;
                     // resumen[l][pos] = (resumen[l][2*pos+1] + resumen[l][2*pos+2]);
//			if(lvl == m_height[l]-1)cout << rmt_v[l][2*pos +1] << "  " << rmt_v[l][2*pos +2] <<  "  "  << 2*pos +1  <<  "  "  << 2*pos +2   <<endl;
			rmt_v[l][pos] = rmt_v[l][2*pos +1] + rmt_v[l][2*pos + 2];
			rmt_r2[l][pos] = rmt_r2[l][2*pos+1] + rmt_r2[l][2*pos + 2];
			int auxpos1 = rmt_m[l][2*pos +1];
			int auxpos2 = rmt_v[l][2*pos+1]+rmt_m[l][2*pos+2];
			rmt_m[l][pos] = min(auxpos1,auxpos2);
//                      if(lvl == height-2)
                }
        }
cout << "PESO = "   << "  "  <<  size_in_mega_bytes(rmt_v[l]) << endl;
cout << "PESO2 = "   << "  "  <<  size_in_mega_bytes(rmt_m[l]) << endl;
cout << "PESO3 = "   << "  "  <<  size_in_mega_bytes(cresumen[l]) << endl;
cout << "PESO4 = "   << "  "  <<  size_in_mega_bytes(resumen[l]) << endl;
cout << "PESO5 = "   << "  "  <<  size_in_mega_bytes(rmt_r[l]) << endl;  
}



cout <<  rmt_v[0][0]  << "  " <<   rmt_v[0][1]  <<  "  "  <<  rmt_v[0][2] << endl;



//COMIENZO RELLENO PARA M_A RANK Y SELECT
	        for(int l = 0; l < m_levels; l++){
	cout << "ENTRE CON NIVEL " << l << endl;
        valor = l;
        m_bs = bs;
        m_num_aleaves[l] = ceil((double)m_A[l].size()/m_bs);
        m_aheight[l] = ceil(log(m_num_aleaves[l])/log(2));
        m_num_ainternal[l] =(pow(2,m_aheight[l]+1)-1);
	rmt_ar[l].resize(m_num_aleaves[l]+1);
	rmt_ar2[l].resize(m_num_ainternal[l]);
	for(int ll = 0; ll < m_num_internal[l];ll++){
	rmt_r2[l][ll] = 0;
}
	int auxcompro= 0;
	int rmb = 0;
	rmt_ar[l][0] = 0;
        for(int i = 0; i < m_num_aleaves[l]; i++){
                int pos = i*m_bs;
                int cont = 0,vmin = 1,cont2 = 0;
                for(int j = 0; j <  m_bs && pos +j < m_A[l].size();j++){
			if(m_A[l][pos+j])rmb++;
			if(m_A[l][pos+j])cont2++;
                }

		rmt_ar[l][i+1] = rmb;
                int top_nodes = pow(2, m_aheight[l])-1;
		rmt_ar2[l][top_nodes + i] = cont2;
                if(i == valor)posiniciala[l] = top_nodes;
        }
        for(int lvl=m_aheight[l]-1; lvl >= 0 ; lvl--){
            int num_curr_nodes = pow(2, lvl);
            int node = 0;
                for( ;node < num_curr_nodes  ;node++){
                      int pos = (pow(2,lvl)-1) + node;
			rmt_ar2[l][pos] = rmt_ar2[l][2*pos+1] + rmt_ar2[l][2*pos + 2];
                }
        }
}


//        m_bs = bs;
//        m_num_leaves = ceil((double)B[0].size()/m_bs);
//        m_height = ceil(log(m_num_leaves)/log(2));
//        m_num_internal =(pow(2,m_height+1)-1);

//	puts("vivo3");

	for(int i = 0; i < m_levels; i++){
//	cout	 rankwt(B[0].size(),i) <<  "  " << resumen[i][0] <<endl;	
	util::bit_compress(resumen[i]);
 util::bit_compress(cresumen[i]);
cout << "PESO3 = "   << "  "  <<  size_in_mega_bytes(cresumen[i]) << endl;
cout << "PESO4 = "   << "  "  <<  size_in_mega_bytes(resumen[i]) << endl;
	}
//	cout << B[0].size() << "  " << resumen[0] << endl;


/*
puts("COMIENZO COMPROBACION REAL");
for(int i = 1; i < 500 ; i++){

//if(m_B_star[1][crankwt(cselectwt(1,i)+1,1)]   != m_B_star[0][cselectwt(1,i)])cout << i-1 <<" "  << cselectwt(1,i)  << " "  << m_B_star[1][i]  <<  "  " << m_B_star[0][cselectwt(1,i)] << " "  << crankwt(cselectwt(1,i)+1,1) <<endl;

if(m_B_star[1][i-1] != m_B_star[0][cselectwt(1,i)]) cout<< "MALO EN COMPROBACION: " << i-1 << " " << cselectwt(1,i) << " " <<  crankwt(cselectwt(1,i),1) << " " << m_B_star[1][i-1] <<" "  << m_B_star[0][cselectwt(1,i)+1] <<" "  <<m_B_star[0][cselectwt(1,i)-1]  <<endl;
//cout << cselectwt(1,i) << endl;
//cout << m_B_star[0][m_A_select0[0](cselectwt(1,i))] << endl;

//if(m_B_star[1][i-1] != m_B_star[0][cselectwt(1,i)]) cout << i-1  << " " <<cselectwt(1,i) << " " <<  m_A_select0[0](cselectwt(1,i)+1) <<  endl;


}


cout<< "COMPRUEBO ARREGLO A" << endl;
int contcerrados = 0;
for(int i = 0; i < 10000; i++){

if(m_A[1][i] == 0){
contcerrados++;
if(m_A[1][i] != m_A[0][m_A_select0[0](cselectwt(1,contcerrados)+1)]) cout << i << " " <<contcerrados << " " << cselectwt(1,contcerrados)  <<endl;

}


}
cout << "TOTAL CERRADOS = " << contcerrados << endl;



cout<< "COMPRUEBO ARREGLO A1" << endl;
int contcerrados1 = 0;
for(int i = 0; i < 1000; i++){

if(m_A[1][i] == 1){
contcerrados1++;
if(m_A[1][i] != m_A[0][m_A_select1[0](selectwt(1,contcerrados1)+1)]) cout << i << " " <<contcerrados1 << " " << selectwt(1,contcerrados1)  <<endl;
// cout<<" COMPROBACION DE CUANTO DEBE SER: " << i << " " <<contcerrados1 << " " << selectwt(1,contcerrados1)  <<" " << m_A_rank[1](i+1)<<endl;  
}


}
cout << "TOTAL ABIERTOS = " << contcerrados1 << endl;


puts("inicio comprobacion");
for(int i = 0, j = 0; i < 500 && j < 30; i++){

if(wt2[i] >= 1){
//if(m_B_star[1][j] != m_B_star[0][i]) cout << "MALO " << i << " " << j << endl;
//cout << m_B_star[1][j] << " " << m_B_star[0][i] << " " << i << " " << j<< endl;
j++;

}


}
*/
 }

	size_type cfopen(size_type i, int nivel){
	int bloque = floor(i/m_bs);
	int top_nodes = pow(2, cm_height[0])-1;
	int tope = bloque*m_bs + m_bs;
	int contador = 0;
	int resp = 0;
	for(int l =  i; l < tope; l++){
	if(m_B_star[0][l] == 1)contador++;
	else contador--;
	
	if(contador == 0){
		resp = l;
		break;
	}

	}
	if(contador == 0) return resp;

	int actual = top_nodes+bloque;
	while(actual != 0){
	if(actual %2 == 0)actual = (actual -2)/2;
	else{
  	if(contador + rmt_m2[0][actual+1] <= 0){
		actual++;
		break;
		}
	contador   += rmt_v2[0][actual+1];
	actual = (actual -1 )/2;
	}
	}

	while(actual < posinicial[0]){
 	int auxactual = actual*2 +1;
 if(contador + rmt_m2[0][auxactual]  <= 0)actual = actual*2 +1;
else {
contador += rmt_v2[0][auxactual];
actual = auxactual+1;
}
	}
int auxpos =  (actual- cposinicial[0])*m_bs;
int entro = 0;

while(contador > 0 && entro < m_bs){
if(m_B_star[0][auxpos] == 1)contador++;
else contador--;
auxpos++;
entro++;
}

return auxpos-1;
}

      size_type mate(size_type i,int nivel) {
	  size_type pos_in_B = m_A_rank[nivel](i); 
	  if(m_A[nivel][i] == 1)pos_in_B++;
	if(m_A[nivel][i] == 1) {
	  size_type match_in_B;
	  if(m_B[nivel][pos_in_B-1] == 1)
	    match_in_B = m_B_st[nivel].find_close(pos_in_B-1);
	  else {
	  	match_in_B = m_B_st[nivel].find_open(pos_in_B-1);
	  }
	  return m_A_select1[nivel](match_in_B+1);
	}
	else
	  {
	    size_type pos_in_B_star = i - pos_in_B+1; 

	    size_type match_in_B_star;
	    if(m_B_star[nivel][pos_in_B_star-1] == 1){
	    match_in_B_star = m_B_star_st[nivel].find_close(pos_in_B_star-1);
	    return m_A_select0[nivel](match_in_B_star+1);
	    }
	    else{
	    match_in_B_star = m_B_star_st[nivel].find_open(pos_in_B_star-1);
	    return m_A_select0[nivel](match_in_B_star+1);	
	    }
	  }  
	return -1;
      }



//dada una posicion y un nivel retorno si el parentesis  es abierto o cerrado
	size_type Eat(size_type posicion, int nivel){
	return  m_B[0][selectwt(nivel,posicion)];
	}
//CREAR NUEVO SELECTWT BASADO EN [] Y CAMBIAR m_B
        size_type Eat2(size_type posicion, int nivel){
        return  m_B_star[0][cselectwt(nivel,posicion)];
        }



	size_type rankmb(size_type pos, int nivel){
   int bloque = floor(pos/m_bs);   	
int inicio = bloque*m_bs;
	int sumar = 0;
	for(int  l = inicio; l  < pos; l++){
	int pos = (selectwt(nivel,l+1));
	if(m_B[0][pos] == 1)sumar++;
		}
	return sumar+rmt_r[nivel][bloque];
	}


	unsigned int selectmb(unsigned int nivel,unsigned int posx){
int pos = 0,vactual = 0,romper = 0;
long long int mov;
//puts("ENTRE AL SELECT");
	for(int i = 0; i < m_height[nivel]; i++){
		mov = pos<<1;
		if(rmt_r2[nivel][mov+1]+vactual >= posx){
			pos = mov+1;
		}
		else{
			vactual+= rmt_r2[nivel][mov+1];
			pos = mov+2;
		}
	}
int actual = vactual;
int posactual =  (pos-posinicial[nivel])*m_bs;
pos = posactual;


    unsigned int r = posactual+ m_bs,l = posactual,m,totalentro = 0;

    	for(; posactual  < r && actual < posx; posactual++){
	int pos2 = (selectwt(nivel,posactual+1));
	if(m_B[0][pos2] == 1)actual++;
		}
//	cout <<"SELECT = " << posactual << endl;
		return posactual-1;
	}


	size_type bsearch(size_type i, int nivel, int obj){
	int bloque = floor(i/m_bs);
	int top_nodes = pow(2, m_height[nivel])-1;
	int tope = bloque*m_bs;
	int contador = 0;
	int resp = 0;
	for(int l =  i; l >= tope; l--){
	int pos = selectwt(nivel,l+1);
	if(m_B[0][pos] == 1)contador--;
	else contador++;
	if(contador == obj){
		resp = l;
		break;
	}

	}
	if(contador == obj) return resp;
	int actual = top_nodes+bloque;
//cout << "PASE BUSQUEDA LINEAL CON " << actual << " " << i <<endl;
	while(actual != 0 && contador - rmt_v[nivel][actual-1] + rmt_m[nivel][actual-1] > obj){
	if(actual %2 != 0)actual = (actual -1)/2;
	else{
/*
  	if(contador + rmt_m[0][actual-1] >= obj){
		actual--;
		break;
		}
*/
	contador   -= rmt_v[nivel][actual-1];
	actual = (actual -2 )/2;
	}

	}


//	cout << "ACTUAL = " << actual << endl;
//AGREGAR CONDICION CON COMPARAR CON RAIZ EN CASO DE QUE VALOR NO SE ENCUENTRA
actual-=1;
	while(actual <  posinicial[nivel]){
 	int auxactual = actual*2 +2;
 if(contador - rmt_v[nivel][auxactual] + rmt_m[nivel][auxactual]  <= obj){
actual = auxactual;
}
else {
contador -= rmt_v[nivel][auxactual];
//EN CASO DE ERROR REVISAR ESTE CASO (AL IR BAJANDO Y MOVERME A LA DERECHA ACA ESTOY ASUMIENDO QUE SIEMPRE  SE CUMPLE EL IF)
actual = auxactual-1;
}
	}

//REVISAR FORMULA EN CASO DE ERROR
//cout << " auxpos " << auxpos <<"  " <<contador <<"  " << obj <<"  " <<actual  << " " <<m_B[nivel].size()-1  <<endl;
int auxpos =  (actual-posinicial[nivel])*m_bs + (m_bs-1);
//cout << " auxpos " << auxpos <<"  " <<contador <<"  " << obj <<"  " <<actual  << " " <<m_B[nivel].size()-1  << " " <<  (actual-posinicial[nivel])*m_bs << endl;
int entro = 0;


if(auxpos > m_B[nivel].size()){

auxpos = m_B[nivel].size()-2;
//obj++;
}
//cout << " auxpos " << auxpos <<"  " <<contador <<"  " << obj <<"  " <<actual  << " " <<m_B[nivel].size()-1  <<endl;   


while(contador != obj && entro < m_bs){
int pos = selectwt(nivel,auxpos+1);
//cout <<"POS " << pos << endl;
if(m_B[0][pos] == 1)contador--;
else contador++;
auxpos--;
entro++;
}

//puts("RETORNARE FINAL");
return auxpos+1;
}








	size_type fopen2(size_type i, int nivel){
	int bloque = floor(i/m_bs);
	int top_nodes = pow(2, m_height[0])-1;
	int tope = bloque*m_bs + m_bs;
	int contador = 0;
	int resp = 0;
	for(int l =  i; l < tope; l++){
	if(m_B[0][l] == 1)contador++;
	else contador--;
	
	if(contador == 0){
		resp = l;
		break;
	}

	}
	if(contador == 0) return resp;

	int actual = top_nodes+bloque;
	while(actual != 0){
	if(actual %2 == 0)actual = (actual -2)/2;
	else{
  	if(contador + rmt_m[0][actual+1] <= 0){
		actual++;
		break;
		}
	contador   += rmt_v[0][actual+1];
	actual = (actual -1 )/2;
	}
//	if(contador + rmt_m[0][actual] <= 0)break;
//	contador += rmt_v[0][actual];
	}

	while(actual < posinicial[0]){
 	int auxactual = actual*2 +1;
 if(contador + rmt_m[0][auxactual]  <= 0)actual = actual*2 +1;
else {
contador += rmt_v[0][auxactual];
//EN CASO DE ERROR REVISAR ESTE CASO (AL IR BAJANDO Y MOVERME A LA DERECHA ACA ESTOY ASUMIENDO QUE SIEMPRE  SE CUMPLE EL IF)
actual = auxactual+1;
}
	}
//REVISAR FORMULA EN CASO DE ERROR
//cout << "TAM = "  << actual << "  " << posinicial[0] << "  " <<  m_num_leaves[0] << "  "  << m_num_internal[0]  << endl;
int auxpos =  (actual- posinicial[0])*m_bs;
int entro = 0;

while(contador > 0 && entro < m_bs){
if(m_B[0][auxpos] == 1)contador++;
else contador--;
auxpos++;
entro++;
}

//puts("RETORNARE AUXPOS -1");
return auxpos-1;
}



      size_type mate2(size_type i,int nivel, int tipoactual) {
          size_type pos_in_B = m_A_rank[nivel](i);
          if(m_A[nivel][i] == 1)pos_in_B++;
//cout << "tipo: " << tipoactual << endl;
	if(tipoactual != m_A[nivel][i])cout << "NUEVA MATE2 MALOOOO " << endl;
        if(m_A[nivel][i] == 1) {
//cout << "ENTRE " << endl;
          size_type match_in_B;
//          if(m_B[nivel][pos_in_B-1] == 1)puts("VALGO 1");
//cout << "valores: " << pos_in_B << " " << B[nivel].size();
//int pareo = B_select1[nivel](pos_in_B);
//cout << pos_in_B << " " << m_B[0].size() << " " <<m_B[nivel].size() << " "  << m_B[nivel][pos_in_B-1]  <<endl;
int realvalor;
if(pos_in_B>= m_B[nivel].size() )realvalor = 0;
else realvalor = Eat(pos_in_B,nivel);
 if(m_B[nivel][pos_in_B-1] != realvalor)cout << "EAT DISTINTO: " << endl; //<< Eat(pos_in_B,nivel)  << "  " <<  m_B[nivel][pos_in_B-1]  << " --  "  <<  endl; //pareo <<  "  "  << selectwt(nivel,pos_in_B)  <<
 // "  --  "     <<m_B[0][pareo]<< "  "  <<  nivel     <<  "  "   <<  m_B[nivel].size()   << "  "   <<  m_B[0].size()  <<"  " << pos_in_B <<endl;

//puts("ENTRARE IF");

if(realvalor == 1){
//puts("ENTREEEE");

            match_in_B = m_B_st[nivel].find_close(pos_in_B-1);
//puts("entre");
        int  match_in_B2 = fopen2(selectwt(nivel,pos_in_B),nivel);
        if(match_in_B != rankwt2(match_in_B2,nivel) )cout <<  "DISTINTOOOOOS: " << match_in_B << "  "  << match_in_B2 << "  "  <<   nivel <<"  "  << pos_in_B-1  <<"  "   <<  selectwt(nivel,pos_in_B) <<"  "  <<rankwt(selectwt(nivel,pos_in_B) ,nivel)  << "  " <<  m_B_st[0].find_open(match_in_B2)<< "  "  << rankwt2(match_in_B2,nivel) <<endl;
}
          else {
//	puts("ENTREEE A FOPEN");
                            match_in_B = m_B_st[nivel].find_open(pos_in_B-1);
//	int  match_in_B2 = fopen2(pos_in_B,nivel);
//	if(match_in_B2 != match_in_B)cout <<  "DISTINTOOOOOS: " << match_in_B << "  "  << match_in_B2 << endl;
          }
          return m_A_select1[nivel](match_in_B+1);
        }
        else
          {//		puts("ENTRE A STARS");
            size_type pos_in_B_star = i - pos_in_B+1;

            size_type match_in_B_star;
            if(m_B_star[nivel][pos_in_B_star-1] == 1){
            match_in_B_star = m_B_star_st[nivel].find_close(pos_in_B_star-1);
            return m_A_select0[nivel](match_in_B_star+1);
            }
            else{ 
            match_in_B_star = m_B_star_st[nivel].find_open(pos_in_B_star-1);
            return m_A_select0[nivel](match_in_B_star+1);
            }
          }
        return -1;
      }




/*

      pair<size_type,int> mate3(size_type i,int nivel, int tipoactual,int nuevaposB) {
          size_type pos_in_B = m_A_rank[nivel](i); //VERIFICAR
          if(tipoactual == 1){
pos_in_B++;
//nuevaposB++;
}
        if(tipoactual == 1) {
          size_type match_in_B;
int realvalor;
if(pos_in_B>= m_B[nivel].size() )realvalor = 0;
else realvalor = Eat(pos_in_B,nivel);
if(realvalor == 1){
 //           match_in_B = m_B_st[nivel].find_close(pos_in_B-1);
        int  match_in_B2 = fopen2(selectwt(nivel,pos_in_B),nivel);
        match_in_B = match_in_B2;
}
          else {
	int  match_in_B2 =  bsearch(pos_in_B-2,nivel,-1); //fopen2(selectwt(nivel,pos_in_B),nivel);
	if(match_in_B2+1 == m_B[nivel].size())match_in_B2 = 0;
match_in_B =	match_in_B2;
          }
          return make_pair(match_in_B,1);
        }
        else
          {
            size_type pos_in_B_star = i - pos_in_B+1;
            size_type match_in_B_star;
            if(/m_B_star[0][cselectwt(nivel,pos_in_B_star)] == 1){ 
int  match_in_B2_star = cfopen(selectwt(nivel,pos_in_B_star),nivel);
match_in_B_star = match_in_B2_star;
 return make_pair(match_in_B_star,0);
            }
            else{ 

int  match_in_B2_star = cfopen(selectwt(nivel,pos_in_B_star),nivel);
match_in_B_star = match_in_B2_star;
             return make_pair(match_in_B_star,0);
            }
          }
        return make_pair(-1,-1);
      }


*/


/*
      pair<size_type,int> mate3(size_type i,int nivel, int tipoactual,int nuevaposB) {
          size_type pos_in_B = m_A_rank[nivel](i);
          if(m_A[nivel][i] == 1){
pos_in_B++;
//nuevaposB++;
//cout << "ENTRE AL IF " << endl;
}
//else cout << "ENTRE AL ELSE " << endl;
//nuevaposB--;
//          if(nuevaposB != pos_in_B  && pos_in_B != rankwt(m_A_rank[0](nuevaposB),nivel) )cout << "POSB3 MALA " << nuevaposB << " " << pos_in_B << " " <<rankwt(m_A_rank[0](nuevaposB),nivel)  << endl; 
//cout << "tipo: " << tipoactual << endl;
	if(tipoactual != m_A[nivel][i])cout << "NUEVA MATE2 MALOOOO " << endl;
        if(m_A[nivel][i] == 1) {
//cout << "ENTRE " << endl;
          size_type match_in_B;
//          if(m_B[nivel][pos_in_B-1] == 1)puts("VALGO 1");
//cout << "valores: " << pos_in_B << " " << B[nivel].size();
//int pareo = B_select1[nivel](pos_in_B);
//cout << pos_in_B << " " << m_B[0].size() << " " <<m_B[nivel].size() << " "  << m_B[nivel][pos_in_B-1]  <<endl;
int realvalor;
if(pos_in_B>= m_B[nivel].size() )realvalor = 0;
else realvalor = Eat(pos_in_B,nivel);
 if(m_B[nivel][pos_in_B-1] != realvalor)cout << "EAT DISTINTO: " << endl; //<< Eat(pos_in_B,nivel)  << "  " <<  m_B[nivel][pos_in_B-1]  << " --  "  <<  endl; //pareo <<  "  "  << selectwt(nivel,pos_in_B)  <<
 // "  --  "     <<m_B[0][pareo]<< "  "  <<  nivel     <<  "  "   <<  m_B[nivel].size()   << "  "   <<  m_B[0].size()  <<"  " << pos_in_B <<endl;

//puts("ENTRARE IF");

if(realvalor == 1){
//puts("ENTREEEE");

            match_in_B = m_B_st[nivel].find_close(pos_in_B-1);
//puts("entre");
        int  match_in_B2 = fopen2(selectwt(nivel,pos_in_B),nivel);
        if(match_in_B != rankwt2(match_in_B2,nivel) )cout <<  "DISTINTOOOOOS 2: " << match_in_B << "  "  << match_in_B2 << "  "  <<   nivel <<"  "  << pos_in_B-1  <<"  "   <<  selectwt(nivel,pos_in_B) <<"  "  <<rankwt(selectwt(nivel,pos_in_B) ,nivel)  << "  " <<  m_B_st[0].find_open(match_in_B2)<< "  "  << rankwt2(match_in_B2,nivel) <<endl;
}
          else {
//	puts("ENTREEE A FOPEN");
                            match_in_B = m_B_st[nivel].find_open(pos_in_B-1);
	int  match_in_B2 =/* fopen2(pos_in_B,nivel)*/ /*fopen2(selectwt(nivel,pos_in_B-1),nivel)*/ /*  bsearch(pos_in_B-2,nivel,-1) ;
	if(match_in_B2+1 == m_B[nivel].size())match_in_B2 = 0;
	if(match_in_B2 != match_in_B)cout <<  "DISTINTOOOOOS: " << match_in_B << "  "    << bsearch(pos_in_B-2,nivel,-1)  << " "  
<< nivel  << " "  << m_B[nivel].size() << " " <<  pos_in_B-2 <<" "  <<  endl;
          }
//		cout << "PRUEBA: "<< m_A_select1[nivel](match_in_B+1)<< " " << selectwt(nivel,match_in_B+1); 
          return make_pair(match_in_B,1);
          //return m_A_select1[nivel](match_in_B+1);
        }
        else
          {//		puts("ENTRE A STARS");
            size_type pos_in_B_star = i - pos_in_B+1;

            size_type match_in_B_star;
	if(m_B_star[0][cselectwt(nivel,pos_in_B_star)] != m_B_star[nivel][pos_in_B_star-1])cout << "STAR MATE3 MALO" << endl;
            if(m_B_star[nivel][pos_in_B_star-1] == 1){
            match_in_B_star = m_B_star_st[nivel].find_close(pos_in_B_star-1);
int  match_in_B2_star = cfopen(cselectwt(nivel,pos_in_B_star),nivel);
if(match_in_B_star != crankwt(match_in_B2_star,nivel) ) cout << "STAR MALO 1 " <<   match_in_B_star << " " << crankwt(match_in_B2_star,nivel) <<" " << match_in_B2_star <<endl;  
//if(match_in_B_star != crankwt(match_in_B2_star,nivel) ) cout << "STAR MALO 1 " <<   match_in_B_star << " " << crankwt(match_in_B2_star,nivel) <<" " << match_in_B2_star <<endl;
 return make_pair(match_in_B_star,0);
//            return m_A_select0[nivel](match_in_B_star+1);
            }
            else{ 
            match_in_B_star = m_B_star_st[nivel].find_open(pos_in_B_star-1);

int  match_in_B2_star = cbsearch(pos_in_B_star-2,nivel,-1);//cfopen(cselectwt(nivel,pos_in_B_star),nivel);
if(match_in_B_star !=  match_in_B2_star ) cout << "STAR MALO 2" << match_in_B_star <<" " <<   match_in_B2_star  <<" "  << cbsearch(pos_in_B_star-1,nivel,-1)  <<endl ;

             return make_pair(match_in_B_star,0);
//            return m_A_select0[nivel](match_in_B_star+1);
            }
          }
        return make_pair(-1,-1);
      }

*/


      pair<size_type,int> mate3(size_type i,int nivel, int tipoactual,int nuevaposB) {
          size_type pos_in_B = m_A_rank[nivel](i); //VERIFICAR
          if(tipoactual == 1){
pos_in_B++;
//nuevaposB++;
}
        if(/*m_A[nivel][i]*/tipoactual == 1) {
          size_type match_in_B;
int realvalor;
if(pos_in_B>= m_B[nivel].size() )realvalor = 0;
else realvalor = Eat(pos_in_B,nivel);
if(realvalor == 1){
 //           match_in_B = m_B_st[nivel].find_close(pos_in_B-1);
        int  match_in_B2 = fopen2(selectwt(nivel,pos_in_B),nivel);
        match_in_B = rankwt2(match_in_B2,nivel);
}
          else {
	int  match_in_B2 = bsearch(pos_in_B-2,nivel,-1);   //  fopen2(pos_in_B,nivel);
if(match_in_B2+1 == m_B[nivel].size())match_in_B2 = 0;
match_in_B =	match_in_B2;
          }
          return make_pair(match_in_B,1);
        }
        else
          {
            size_type pos_in_B_star = i - pos_in_B+1;
            size_type match_in_B_star;
            if(m_B_star[0][cselectwt(nivel,pos_in_B_star)] == 1){ 
int  match_in_B2_star = cfopen(cselectwt(nivel,pos_in_B_star),nivel);
match_in_B_star = crankwt(match_in_B2_star,nivel);
 return make_pair(match_in_B_star,0);
            }
            else{ 

int  match_in_B2_star = cbsearch(pos_in_B_star-2,nivel,-1);
match_in_B_star = match_in_B2_star;
             return make_pair(match_in_B_star,0);
            }
          }
        return make_pair(-1,-1);
      }





      /* Assuming indices start with 0 */
      size_type next(size_type i) {
      //	puts("entre a next");
      	cout << i << endl;
	if(i > m_A[0].size())
	  return -1;
	
	if(m_A[0][i] == 0) {
	  return i+1;
	}
	else {
	//	puts("entre al else");
	  size_type pos_in_B = m_A_rank[0](i); // rank1
	 // cout << "pos in b "<< pos_in_B  <<endl;	
	  if(m_B[0][pos_in_B] == 1) {
	  	//cout << "pos in b "<< pos_in_B<<endl;	
	    return mate(i)+1;
	  }
	}
	return -1; //NUNCA DEBERIA RETORNAR ESTO, LA ESTRUCTURA ASEGURA QUE NO LLEGE A ESTE RETORNO
      }


size_type selectmBStartoA(int e, int nivel){

//Obtengo rank de m_B[nivel][pos]
//size_type posmB = rankwt(e,nivel); //inclusivo rank, creo que si?
//Paso desde m_B[nivel] a m_B[0]
size_type m_B0 = cselectwt(nivel,/*posmB+1*/ e+1);
//Obtengo rank de m_B[nivel][pos]
//size_type posmB = rankwt(m_B0,nivel);
//Paso de m_B[0] hasta m_A[0]
size_type posA = m_A_select0[0](m_B0+1);

//Obtengo cantidad de 0
size_type cantmBstar = posA - m_B0;
//cuento cuantos 0 pertenecen al nivel buscado
size_type cnivel = rankwt(cantmBstar,nivel);
// Respuesta es la suma de 1 perteneciente al nivel + suma de 0's
size_type retorno = cnivel +  e;
//cout << cnivel <<" " <<m_B_rank[nivel](e)   << " "  << posmB  <<  endl;
return retorno;
}



  /*
      size_type degree(size_type v) {
	if(v >= m_vertices)
	  return 0;
	
	size_type dg = 0;
	size_type nxt = first(v);
	while(nxt <= 2*m_edges+2 && nxt >= 0) {
		cout << nxt << endl;
	  nxt = next(nxt);
	  dg++;
	}
	return dg;
      }
*/

//Dado un parentesis me retorna la region correspondiente.

size_type vertex2(size_type e, int nivel){
      size_type pos_in_A = m_A_rank[nivel](e);
	  if(m_A[nivel][e] == 1)pos_in_A++;      	

      	if(m_B[nivel][pos_in_A-1] == 1 && m_A[nivel][e] == 1){
      		return m_B_st[nivel].rank(pos_in_A-1)-1;
      	}
      	else if(m_B[nivel][pos_in_A-1] == 0 && m_A[nivel][e] == 1){
      		size_type match_pos = m_B_st[nivel].find_open(pos_in_A-1);
      		return m_B_st[nivel].rank(match_pos)-1;
      	}
      	size_type pos_in_B = e- pos_in_A+1;

      	 if(m_A[nivel][e] == 0 && m_B_star[nivel][pos_in_B-1] == 1) {
     		size_type match_pos = mate(e,nivel);
     		size_type pos = m_A_rank[nivel](match_pos);
     		size_type pos2 = m_A_select1[nivel](pos+1); 
     		if(m_B[nivel][m_A_rank[nivel](pos2)] == 1){

     			pos2 = m_B_st[nivel].parent_t(m_A_rank[nivel](pos2));
     			size_type pos4 = m_A_select1[nivel](pos2+1);
				size_type pos3 =m_A_rank[nivel](pos4);
				if(m_A[nivel][pos4] == 1)pos3++;
     			return m_B_rank[nivel](pos3-1);
     		}
     		else {
     			pos2 = mate(pos2,nivel);
     			size_type pos3 =m_A_rank[nivel](pos2);
     			if(m_A[nivel][pos2] == 1)pos3++;
     			return m_B_rank[nivel](pos3-1);
     		}
      		
      	}
	else if(m_A[nivel][e] == 0 && m_B_star[nivel][pos_in_B-1] == 0){

     		size_type match_pos = mate(e,nivel);
     		size_type pos = m_A_rank[nivel](match_pos);
			if(m_A[nivel][match_pos] == 1)pos++;

     		size_type pos2 = m_A_select1[nivel](pos+1); 

     		if(m_B[nivel][m_A_rank[nivel](pos2)] == 1){

     			pos2 = m_B_st[nivel].parent_t(m_A_rank[nivel](pos2));

     			size_type pos4 = m_A_select1[nivel](pos2+1);
				size_type pos3 =m_A_rank[nivel](pos4);
				if(m_A[nivel][pos4] == 1)pos3++;
     			return m_B_rank[nivel](pos3-1);
     		}
     		else {

     			pos2 = mate(pos2,nivel);
     			size_type pos3 =m_A_rank[nivel](pos2);
     			if(m_A[nivel][pos2] == 1)pos3++;

     			return m_B_rank[nivel](pos3-1);
     		}
      	}
      	return -1;
      }


size_type selectmBtoAF(int e, int nivel){

//Obtengo rank de m_B[nivel][pos]
//size_type posmB = rankwt(e,nivel); //inclusivo rank, creo que si?
//Paso desde m_B[nivel] a m_B[0]
//size_type m_B0 = selectwt(nivel,/*posmB+1*/ e+1);
//Obtengo rank de m_B[nivel][pos]
//size_type posmB = rankwt(m_B0,nivel);
//Paso de m_B[0] hasta m_A[0]
size_type posA = m_A_select1[0](e+1);

//Obtengo cantidad de 0
size_type cantmBstar = posA - e;
//cuento cuantos 0 pertenecen al nivel buscado
size_type cnivel = crankwt(cantmBstar,nivel);
// Respuesta es la suma de 1 perteneciente al nivel + suma de 0's
size_type retorno = cnivel +  rankwt(e,nivel);
//cout << cnivel <<" " <<m_B_rank[nivel](e)   << " "  << posmB  <<  endl;
return retorno;
}


size_type mArankFromA(unsigned int x, unsigned int posx){
if(x == 0 )return m_A_rank[0](posx);
else return crankwtba(cselectwt(x,posx+1),x);
}


size_type vertex3(size_type e, int nivel, int inicial, int nposA){
//cout << "ENTRE A VERTEX 3 	" << endl;
      size_type pos_in_A = m_A_rank[nivel](e);
  if(nposA != pos_in_A) cout << "NPOSA MALA: " << nposA << " " << pos_in_A << endl;
//if( mArankFromA(nivel,e) != pos_in_A ) cout << "RANK MALO: " <<  mArankFromA(nivel,e) << " " << pos_in_A << " "  << nivel <<endl; 
	if(m_A[nivel][e] != inicial) cout << "PRUEBA CON INICIAL MALA" << endl;
	  if(m_A[nivel][e] == 1)pos_in_A++;      	
//cout << "PRUEBO EAT " <<" "  <<pos_in_A <<endl;
	  	if(m_B[nivel][pos_in_A-1]  !=  Eat(pos_in_A,nivel))cout << "VERTEX3 EAT MALO " << endl;
//cout << "tam: " << m_B[nivel].size() << " " << m_A[nivel].size() << " " << pos_in_A << " " << e << endl; 
      	if(m_B[nivel][pos_in_A-1] == 1 && m_A[nivel][e] == 1){
  //    		puts("RANK3");
	if(m_B_st[nivel].rank(pos_in_A-1) != rankmb(pos_in_A,nivel)) cout << "MALOO  " << m_B_st[nivel].rank(pos_in_A-1) << "  " <<  rankmb(pos_in_A,nivel)<< "  " << m_B_rank[nivel](pos_in_A) <<endl;
      		return m_B_st[nivel].rank(pos_in_A-1)-1;
      	}
      	else if(m_B[nivel][pos_in_A-1] == 0 && m_A[nivel][e] == 1){
//puts("ENTRE AL IF NUMERO 2");
      		size_type match_pos = m_B_st[nivel].find_open(pos_in_A-1);
//	puts("ENTRE AL IF NUMERO 2");
		if(match_pos != bsearch(pos_in_A,nivel,0)) cout << "MALOOO  " << match_pos  << "  " << bsearch(pos_in_A,nivel,0)+1 << endl; 
      		if(m_B_st[nivel].rank(match_pos)-1 != rankmB(match_pos,nivel)) cout << "IF1 MALO: " << rankmB(match_pos,nivel) << " " <<  m_B_st[nivel].rank(match_pos)-1<<endl;
      		return m_B_st[nivel].rank(match_pos)-1;
      	}
//puts("NO ENTRE A NADA");
      	size_type pos_in_B = e- pos_in_A+1;
//if(m_B_star[nivel][pos_in_B-1] != Eat2(pos_in_B-1,nivel)) cout << "EAT2 MALO\n";
if(m_B_star[nivel][pos_in_B-1] != Eat2(pos_in_B,nivel)) cout << "EAT22 MALO\n";
      	 if(m_A[nivel][e] == 0 && m_B_star[nivel][pos_in_B-1] == 1) {
     		size_type match_pos = mate2(e,nivel,inicial);
     	if(selectmBStartoA(mate3(e,nivel,inicial,nposA).first,nivel) != match_pos )cout << "comp malo: " << match_pos << " " <<selectmBStartoA(mate3(e,nivel,inicial,nposA).first,nivel)  <<endl; 
		size_type pos = m_A_rank[nivel](match_pos);
//    		if(pos != mate3(e,nivel,inicial,nposA).first) cout << "comp 2 malo: " << pos << " " <<mate3(e,nivel,inicial,nposA).first  << endl;

if(pos != cont1Ma(mate3(e,nivel,inicial,nposA).first,nivel)) cout << "comp 2 malo: " << pos << " " <<cont1Ma(mate3(e,nivel,inicial,nposA).first,nivel)  << endl;
     		size_type pos2 = m_A_select1[nivel](pos+1); 
int pos2a = m_A_rank[nivel](pos2);

//if(pos2 != selectmBtoA(rankwt(selectwt(nivel,pos+1),nivel),nivel))cout << "IMPRESIONES: " << pos +1 << " " << pos2<< " "  
//<< selectmBtoA(rankwt(selectwt(nivel,pos+1),nivel),nivel)<< " "  << selectwt(nivel,pos+1)<<" "<< rankwt(selectwt(nivel,pos+1),nivel) << " " <<m_B[nivel].size()   <<endl;
//if(m_B[nivel].size() == pos +1 ) cout << "IMPRESIONES 2 : " << pos << " " << pos2<< " "
//<< selectmBtoA(rankwt(selectwt(nivel,pos),nivel),nivel)<< " "  << ultimos[nivel]<<" "<< rankwt(selectwt(nivel,pos),nivel) << endl;
//cout << "pos's A: " << pos << " " << pos2a << endl;
//cout << "PREGUNTARE CON " << pos2 << "  " << nivel <<" " << pos2a+1 << " " << m_B[0].size()<<" " <<m_B[nivel].size()  << " " << m_B[nivel][m_A_rank[nivel](pos2)] <<endl;
int comparador,sientro = 0;
if(pos2a+1 >= m_B[nivel].size()){
comparador = 0;
sientro = 1;
}
else comparador = m_B[0][(selectwt(nivel,pos2a+1))];
 if(m_B[nivel][m_A_rank[nivel](pos2)]  !=  comparador)cout << "VERTEX3 EAT  2 MALO "<< pos2a <<" "  << comparador <<" "<< Eat(pos2a,nivel) <<" "<< (selectwt(nivel,pos2a+1))<<" "<<m_B[nivel].size() <<endl;
     		if(m_B[nivel][m_A_rank[nivel](pos2)] == 1){
//			    puts("ENTRE AL IF 2");
//			cout << pos2a << endl;
     			pos2 = m_B_st[nivel].parent_t(m_A_rank[nivel](pos2));
//			int posclose  = m_B_st[nivel].find_close(pos2a);
//			if(posclose != fopen2(selectwt(nivel,pos2a+1),nivel)) cout << "posclose: " << posclose << " " << fopen2(selectwt(nivel,pos2a+1),nivel) << endl;
//			cout << "ENTRE AL IF 2 "<< rankwt2(bsearch(selectwt(nivel,pos2a+1),nivel,-2),nivel) <<"  "  << pos2a << endl;
//if(pos2 != bsearch(pos2a,nivel,-2)) cout << "MALOO " <<pos2<< "  " << bsearch(pos2a,nivel,-2) <<"  "  << m_B[0][selectwt(nivel,pos2a+1)]<<"  " << nivel << " " << pos2a  <<"  "  <<posclose << endl;
     			size_type pos4 = m_A_select1[nivel](pos2+1);
//		cout << "pos4: " << pos4 << " " << selectmBtoA(pos2,nivel);
				size_type pos3 =m_A_rank[nivel](pos4);
				if(pos2 != pos3) cout <<"pos3 malo: " << pos2 << " "<< pos3 <<endl; 
				if(m_A[nivel][pos4] == 1)pos3++;
				else cout << "EN ALGUN MOMENTO NO SOY 1" << endl;
 if(m_B_rank[nivel](pos3-1) != rankmb(pos3-1,nivel) ) cout << "MALOO  21  " << m_B_rank[nivel](pos3-1) << "  " << rankmb(pos3-1,nivel)   <<endl;
     			return m_B_rank[nivel](pos3-1);
     		}
     		else {
//		    puts("ENTRE AL ELSE");
			size_type posaux = pos2;
                        if( pos != m_A_rank[nivel](pos2)) cout << "mate3pos:  " <<  pos << " " << m_A_rank[nivel](pos2) << endl;
     			pos2 = mate2(pos2,nivel,1);
if(pos2 !=  selectmBtoA(mate3(posaux,nivel,1,pos).first,nivel)) cout << "COMPMALO 2 : " << pos2 << " " <<selectmBtoA(mate3(posaux,nivel,1,pos).first,nivel) <<endl;
	
//     			if( pos != m_A_rank[nivel](pos2)) cout << "mate3pos:  " <<  pos << " " << m_A_rank[nivel](pos2) << endl;
     			size_type pos3 =m_A_rank[nivel](pos2);
//			if(pos3 !=  mate3(posaux,nivel,1,pos).first) cout << "COMPMALO: " << pos3 << " " <<mate3(posaux,nivel,1,pos).first <<endl;
     			if(m_A[nivel][pos2] == 1)pos3++;
			 else cout << "EN ALGUN MOMENTO NO SOY 1" << endl; 
			if(m_B_rank[nivel](pos3-1) != rankmb(pos3-1,nivel) ) cout << "MALOO  2  " << m_B_rank[nivel](pos3-1) << "  " << rankmb(pos3-1,nivel)   <<endl;
     			return m_B_rank[nivel](pos3-1);
     		}
      		
      	}
	else if(m_A[nivel][e] == 0 && m_B_star[nivel][pos_in_B-1] == 0){
	if( m_B_star[nivel][pos_in_B-1] != m_B_star[0][cselectwt(nivel,pos_in_B)]) cout << "MALOOOO B_STAR " <<m_B_star[0][cselectwt(nivel,pos_in_B)-1]  <<" "  <<  m_B_star[0][cselectwt(nivel,pos_in_B)] <<" " <<  nivel << " "<< crankwt(cselectwt(nivel,pos_in_B),nivel)  <<endl;
     		size_type match_pos = mate2(e,nivel,inicial);
     		if(selectmBStartoA(mate3(e,nivel,inicial,nposA).first,nivel) != match_pos) cout << "COMPFINAL1MALO: " << match_pos <<  " " << selectmBStartoA(mate3(e,nivel,inicial,nposA).first,nivel) << endl ;
     		size_type pos = m_A_rank[nivel](match_pos);
     		if(cont1Ma(mate3(e,nivel,inicial,nposA).first,nivel) != pos) cout << "COMPFINAL2MALO: " << pos << " " << cont1Ma(mate3(e,nivel,inicial,nposA).first,nivel) << endl;
			if(m_A[nivel][match_pos] == 1){
pos++;
cout << "EN ALGUN MOMENTO SOY 1" << endl;
}
			//else else cout << "EN ALGUN MOMENTO NO SOY 1" << endl; 
//		if(pos2 != selectmBtoA(pos,nivel))cout "COMPPOS2: " << pos2 << " " << selectmBtoA(pos,nivel) << endl;
     		size_type pos2 = m_A_select1[nivel](pos+1); 
 if(pos2 != selectmBtoA(pos,nivel))cout<< "COMPPOS2: " << pos2 << " " << selectmBtoA(pos,nivel) << endl;
if(Eat(pos+1,nivel) != m_B[nivel][m_A_rank[nivel](pos2)])cout << "EATFINAL MALO\n";
//		cout << Eat(pos,nivel) << " " << Eat(pos+1,nivel) << " " << m_B[nivel][m_A_rank[nivel](pos2)] << endl;
     		if(m_B[nivel][m_A_rank[nivel](pos2)] == 1){
//     puts("ENTRE AL IF 3");  
int pos2a = m_A_rank[nivel](pos2);
if(pos2a != pos)cout << "pos2a: " << pos2a << " " << pos<< endl;
     			pos2 = m_B_st[nivel].parent_t(m_A_rank[nivel](pos2));
if(pos2 != bsearch(pos2a,nivel,-2)) cout << "MALOO  3 " <<pos2<< "  " << bsearch(pos2a,nivel,-2) <<"  "  << m_B[0][selectwt(nivel,pos2a+1)]<<"  " << nivel << " " << pos2a  << endl;  
     			size_type pos4 = m_A_select1[nivel](pos2+1);
				size_type pos3 =m_A_rank[nivel](pos4);
if(pos4 != selectmBtoA(pos2,nivel))cout << "SIGUIENTES: " << pos4 << " " << pos3 << " " <<  selectmBtoA(pos2,nivel) <<endl;
if(pos2 != pos3)cout << "pos3FINAL: " << pos2 << " " << pos3 << endl;
				if(m_A[nivel][pos4] == 1)pos3++;
else cout << "EN ALGUN MOMENTO NO SOY 1" << endl;
 if(m_B_rank[nivel](pos3-1) != rankmb(pos3-1,nivel) ) cout << "MALOO  22  " << m_B_rank[nivel](pos3-1) << "  " << rankmb(pos3-1,nivel)   <<endl;
     			return m_B_rank[nivel](pos3-1);
     		}
     		else {
  //   puts("ENTRE AL ELSE 3");  
		size_type  posauxfinal = pos2;
                        if(m_A_rank[nivel](pos2) != pos) cout << "rankfinal: " << pos << " " << m_A_rank[nivel](pos2) << endl;
     			pos2 = mate2(pos2,nivel,1);
	//		if(m_A_rank[nivel](pos2) != pos) cout << "rankfinal: " << pos << " " << m_A_rank[nivel](pos2) << endl;
     			size_type pos3 =m_A_rank[nivel](pos2);
 if(selectmBtoA(mate3(posauxfinal,nivel,1,pos).first,nivel) != pos2)cout << "MATE3FINAL 1: " << selectmBtoA(mate3(posauxfinal,nivel,1,pos).first,nivel) << " " << pos2 << endl; 
			if(mate3(posauxfinal,nivel,1,pos).first != pos3)cout << "MATE3FINAL 2: " << mate3(posauxfinal,nivel,1,pos).first << " " << pos3 << endl;
     			if(m_A[nivel][pos2] == 1)pos3++;
else cout << "EN ALGUN MOMENTO NO SOY 1" << endl;
if(m_B_rank[nivel](pos3-1) != rankmb(pos3-1,nivel) ) cout << "MALOO  4  " << m_B_rank[nivel](pos3-1) << "  " << rankmb(pos3-1,nivel)   <<endl; 
     			return m_B_rank[nivel](pos3-1);
     		}
      	}
      	return -1;
      }






pair<int,int> MapeoA(int posini, int nivel,int nuevoB, int nuevoestado){
int tipo;
//cout << "INICIE CONTEO\n";
int contcerrados1 = m_A_rank[nivel](posini+1);
//if(posini == 1195569) cout << contcerrados1 << endl;
/*
if(m_A[nivel][posini] == 1){
int posinA = m_A_select1[0](nuevoB+1);
int poscorch = crankwt(posinA,nivel);
if(nuevoB != contcerrados1 || poscorch != corchetes) cout << "Mapeo Malo: " <<nuevoB <<" " <<contcerrados1 << " " << poscorch << " " << corchetes << endl; 
}
*/
/*
if(m_A[nivel][posini+1] == 1){

int pruebaA = selectwt(nivel,posini+1);
int pruebaB = rankwt(pruebaA,nivel);
if(m_A_rank[nivel](posini+1) != pruebaB) cout << "MALO: " << nivel << " " << m_A_rank[nivel](posini+1) << " " << pruebaA << " " << pruebaB << endl;

}
*/

int corchetes = posini   - contcerrados1+1;
//if(posini == 1195569) cout << corchetes <<" " << nuevoB<<endl;
int correcto;


if(/*nuevoestado >= 1 &&*/ nuevoB != 0) {
int posinA = m_A_select1[0](selectwt(nivel,nuevoB+1)+1);
int poscorch = crankwt(posinA-m_A_rank[0](posinA+1)+1,nivel);
if(nuevoB+1 != contcerrados1 || (poscorch != corchetes && corchetes != nuevoestado && abs(corchetes-nuevoestado) > 1)   ) cout << "Mapeo Malo 1 : " <<nuevoB+1 <<" " <<contcerrados1 << " " << poscorch << " " << corchetes << " " << nuevoestado<<endl;
} 

else if(nuevoestado == 0 && nuevoB != 0){
int posinA = m_A_select0[0](cselectwt(nivel,nuevoB+1)+1);
int poscorch = rankwt(m_A_rank[0](posinA+1)+1,nivel);
if(nuevoB+1 != corchetes || poscorch != contcerrados1) cout << "Mapeo Malo 2 : " <<nuevoB+1 <<" " <<contcerrados1 << " " << poscorch << " "<< corchetes <<" "<<cselectwt(nivel,nuevoB+1)<<" " <<posinA
<<" " <<crankwt(nuevoB+1,nivel) << " " << rankwt(m_A_rank[0](crankwt(nuevoB+1,nivel))+1,nivel) <<" " <<m_B_star[nivel].size()  <<endl;
}
if(corchetes >= m_B_star[nivel].size()-1 ){
correcto = m_A_select1[0](selectwt(nivel,contcerrados1+1)+1);
tipo = 1;
}
else if(contcerrados1 >= m_B[nivel].size()-1){
correcto = m_A_select0[0](cselectwt(nivel,corchetes+1)+1);
tipo = 0;
}
else{
correcto = (m_A_select1[0](selectwt(nivel,contcerrados1+1)+1) < m_A_select0[0](cselectwt(nivel,corchetes+1)+1) || m_A_select0[0](cselectwt(nivel,corchetes+1)+1)<
 m_A_select1[0](selectwt(nivel, m_A_rank[nivel](posini+1))+1)  )? m_A_select1[0](selectwt(nivel,contcerrados1+1)+1) : m_A_select0[0](cselectwt(nivel,corchetes+1)+1);
 tipo =  (m_A_select1[0](selectwt(nivel,contcerrados1+1)+1) < m_A_select0[0](cselectwt(nivel,corchetes+1)+1) || m_A_select0[0](cselectwt(nivel,corchetes+1)+1)<
 m_A_select1[0](selectwt(nivel, m_A_rank[nivel](posini+1))+1)  )? 1 :0;
} 
//cout << "RETORNO\n";
return make_pair(correcto,tipo);

}


      vector <int> lista_vecinos(size_type v, int nivel){
      	size_type nxt = first(v,nivel)+1,rankposini;
if(nxt != selectmBtoA(first2(v,nivel),nivel)+1)cout << "if1 malo "  <<  nxt << " " << selectmBtoA(first2(v,nivel),nivel)+1 << endl;
	pair <size_type,int > nxt2; 
      	pair<size_type,int> limite = mate3(first(v,nivel),nivel,1, first2(v,nivel)+1);
	int posini2 = first(v,nivel);
      	vector <int> retorno;
	pair<size_type,int> posini = make_pair(first2(v,nivel),1);
	int simb = 1;
int valordeprueba = 0;
int valorinsert = 0;
int tipoinsert = 1;
size_type nxtimp = nxt;
//tipoinsert = MapeoA(posini2,nivel,posini.first,posini.second).second;
/*if(tipoinsert == 1)*/valorinsert = rankwt(m_A_rank[0]( MapeoA(posini2,nivel,/*posini.first,posini.second*/0,0).first),nivel)-1;
//else valorinsert =crankwt(MapeoA(posini2,nivel,/*posini.first,posini.second*/0,0).first  - m_A_rank[0]( MapeoA(posini2,nivel,/*posini.first,posini.second*/0,0).first+1) -1,nivel);
//crankwt(MapeoA(posini2,nivel,posini.first,posini.second).first - m_A_rank[0](MapeoA(posini2,nivel,posini.first,posini.second).first+1),nivel);
nxt2 = make_pair(valorinsert , 1 );
if(MapeoA(posini2,nivel,/*posini.first,posini.second*/0,0).second == 0)
valordeprueba = crankwt(MapeoA(posini2,nivel,/*posini.first,posini.second*/0,0).first-m_A_rank[0]( MapeoA(posini2,nivel,/*posini.first,posini.second*/0,0).first),nivel);
bool inicial = true;
rankposini = posini.first+1;
size_type valorrank = first2(v,nivel), tipoA = selectmBtoA(first2(v,nivel)+1,nivel);
//cout << "INICIAL= " << valorrank<<endl;
if(  selectmBtoA(limite.first,nivel)  != mate2(first(v,nivel),nivel,1) ) cout << "LIMITE MALO: " << mate2(first(v,nivel),nivel,1) << " " <<   selectmBtoA(limite.first,nivel)  <<endl;
      	while(nxt < mate2(first(v,nivel),nivel,1)){
      		size_type x;
      		size_type pos_in_A = m_A_rank[nivel](nxt);
		if(valorrank+1 != pos_in_A)cout << "VALORRANK MALO " << valorrank+1 << " " << pos_in_A << endl;
if(inicial){
inicial = false;
//tipoA = m_A[0][MapeoA(posini2,nivel,posini.first,posini.second).first];
//if(m_A[nivel][nxt] != m_A[0][MapeoA(posini2,nivel,posini.first,posini.second).first])cout<< "FUNCION MALA" << endl;
}
else {
//cout << "CANT01 "  <<cantCeros(posini.first,nivel)  <<endl;
//if(mate2(nxt,nivel,1)+1 >= mate2(first(v,nivel),nivel,1))tipoA = -1;
//else{
//if(cantCeros(posini.first,nivel) == cantCeros(posini.first+1,nivel)) tipoA =1;
//else tipoA = 0;
//cout << "CANT02 "  <<cantCeros(posini.first+1,nivel)  <<endl;
//}
//if(m_A[nivel][nxt] == 1)cout << "ceros: " << cantCeros(posini.first,nivel) << " " << cantCeros(posini.first+1,nivel) << endl;
//tipoA = cantCeros(posini.first,nivel);
}
//if(tipoA != m_A[nivel][nxt]) cout << "TIPOA MALO:" << endl;
//if(MapeoA(posini2,nivel,posini.first,posini.second).second != m_A[nivel][nxt]) cout << "MAPEOMALO " << endl;
//if(m_A[nivel][nxt] != m_A[0][MapeoA(posini2,nivel,posini.first,posini.second).first])cout<< "FUNCION MALA" << endl;
if(tipoA == nxt &&  m_A[nivel][nxt] != 1)cout <<"IDEA MALA"<< endl;
      		if(m_A[nivel][nxt] == 1)pos_in_A++;
//		else cout << "NO SOY 1" << endl;
//          if(first2(v,nivel)+1 != pos_in_A)cout << "POSA MALA " << first2(v,nivel)+1 << " " << pos_in_A << endl;
if(m_B[nivel][pos_in_A-1]  !=  Eat(pos_in_A,nivel))cout << "V1 MALO " << endl;
      		if(m_A[nivel][nxt] == 1 && m_B[nivel][pos_in_A-1] == 1){
//	if()
      			x = vertex4(nxt,nivel,MapeoA(posini2,nivel,0,0).second,rankposini);
if(rankposini != m_A_rank[nivel](nxt)) cout << "RANK MALO0: " << posini.first+1 << " " << m_A_rank[nivel](nxt) << " " << rankposini	 <<endl;
posini = mate3(nxt,nivel,1,MapeoA(posini2,nivel,nxt2.first,valordeprueba/*0,0*/).first+1);
posini2 =mate2(nxt,nivel,1);
if(posini2 != selectmBtoA(posini.first,nivel)) cout<<"posini malo: " << posini2 << " " <<selectmBtoA(posini.first,nivel) << endl; 
simb = 1;
nxt = mate2(nxt,nivel,1)+1;
rankposini = posini.first+1;
tipoA = selectmBtoA(posini.first+1,nivel);
valorrank = posini.first;
if(nxt != selectmBtoA(posini.first,nivel)+1) cout <<"nxt malo: "   <<nxt << " " << selectmBtoA(posini.first,nivel)+1 << endl; 

valorinsert = rankwt(m_A_rank[0]( MapeoA(posini2,nivel,posini.first,posini.second/*0,0/*nxt2.first,valordeprueba*/).first),nivel)-1;
//if(MapeoA(posini2,nivel,/*posini.first,posini.second*/0,0).second == 0)
valordeprueba = crankwt(MapeoA(posini2,nivel,posini.first,posini.second).first-m_A_rank[0]( MapeoA(posini2,nivel,posini.first,posini.second).first),nivel);
nxt2 = make_pair(valorinsert , 1 );
      		}
      		else{
//cout << m_A[nivel][nxt] << " " << m_B[nivel][pos_in_A-1] << endl;
//posini = mate3(nxt,nivel,MapeoA(posini2,nivel,0,0).second,MapeoA(posini2,nivel,0,0).first+1); 
//posini2 =mate2(nxt,nivel,MapeoA(posini2,nivel,0,0).second);
//posini2 = nxt;
//cout << " Mapeo: " << MapeoA(posini2,nivel,0,0).first <<" "  << nxt << endl;
x = vertex4(nxt,nivel,MapeoA(posini2,nivel,0,0).second,rankposini);
if(rankposini != m_A_rank[nivel](nxt)) cout << "RANK MALO: " << posini.first+1 << " " << m_A_rank[nivel](nxt) <<" "  << rankposini  <<endl;
//posini = make_pair(nxt,MapeoA(posini2,nivel,/*posini.first,posini.second*/0,0).second);
posini = make_pair(nxt,MapeoA(posini2,nivel,nxt2.first,valordeprueba/*0,0*/).second);
if(rankposini != m_A_rank[nivel](nxt)) cout << "RANK MALO3: " << posini.first+1 << " " << m_A_rank[nivel](nxt) <<" "  << rankposini <<endl;
//posini = make_pair(nxt,MapeoA(posini2,nivel,0,0).second);
//tipoinsert = MapeoA(posini2,nivel,/*posini.first,posini.second*/0,0).second;
//if(tipoinsert == 1)valorinsert = rankwt(m_A_rank[0]( MapeoA(posini2,nivel,/*posini.first,posini.second*/0,0).first),nivel)+11;
//else valorinsert = crankwt(MapeoA(posini2,nivel,/*posini.first,posini.second*/0,0).first+1 - m_A_rank[0](MapeoA(posini2,nivel,/*posini.first,posini.second*/0,0).first),nivel)+11;
//nxt2 = make_pair(valorinsert , tipoinsert );
valorinsert = rankwt(m_A_rank[0]( MapeoA(posini2,nivel,/*posini.first,posini.second*/nxt2.first,valordeprueba).first),nivel)-1;
nxt2 = make_pair(valorinsert , 1 );
posini2 = nxt;
//		posini = make_pair(nxt,MapeoA(posini2,nivel,posini.first,posini.second).second);
      			nxt ++;
//rankposini++;
//valorinsert =rankwt(m_A_rank[0]( MapeoA(posini2,nivel,/*posini.first,posini.second*/0,0).first),nivel)-1;
//nxt2 = make_pair(valorinsert , 1 );
//cout << "V0: " << posini2 <<" "  << m_A[nivel].size() <<endl;      	
//if(nxt < mate2(first(v,nivel),nivel,1))//cout << "V1 :" << MapeoA(posini2,nivel,/*posini.first,posini.second*/0,0).first << endl;
//cout << "V2: " << m_A_rank[0]( MapeoA(posini2,nivel,/*posini.first,posini.second*/0,0).first) << endl;
if(nxt < mate2(first(v,nivel),nivel,1))if(MapeoA(posini2,nivel,/*posini.first,posini.second*/nxt2.first,valordeprueba).second == 0)
valordeprueba = crankwt(MapeoA(posini2,nivel,/*posini.first,posini.second*/0,0).first-m_A_rank[0]( MapeoA(posini2,nivel,/*posini.first,posini.second*/0,0).first),nivel);

	}
      		retorno.push_back(x);
      	}
int malo = 0;
      	if(v != 0){
      		nxt = m_B_st[nivel].parent_t(m_B_select1[nivel](v+1));
//if( m_B_select1[nivel](v+1) != selectmb(nivel,v+1))	cout << m_B_select1[nivel](v+1) << " " <<  selectmb(nivel,v+1)  <<  " "  << nivel  << endl;
//cout << m_B_select1[nivel](v) << endl;
//cout << m_B_select1[nivel](v+1) << endl;
//cout << bsearch(m_B_select1[nivel](v+1),nivel,-2) << endl;
//return retorno;
if(nxt != bsearch(m_B_select1[nivel](v+1),nivel,-2)) cout << "V2 MALO " <<" " << nxt << " "  << bsearch(m_B_select1[nivel](v+1),nivel,-2)   <<endl;
//return retorno;
      		int auxtemp = m_A_rank[nivel](m_A_select1[nivel](nxt+1));
if(auxtemp != nxt) cout<< "MALO: " << auxtemp << " " << nxt+1 <<" "  << m_A_select1[nivel](nxt+1) <<endl;
//return retorno;
      		size_type x = m_B_rank[nivel](auxtemp);
if(m_B_rank[nivel](auxtemp) != rankmb(auxtemp,nivel)) cout << "V3 MALO " << endl ;
if(m_B[nivel][auxtemp]  !=  Eat(auxtemp+1,nivel))cout << "V4 EAT MALO " << endl;
      		if(m_B[nivel][auxtemp] == 1)x++;
      		if(x-1 < limite.first)retorno.push_back(x-1);
      	}
      	return retorno;
      }


      //Funciones auxiliar para hacer coincidir las regiones consultadas en el baseline
	int getMapeo(size_type pos, unsigned int nivel){
    int retorno = m_B_select1[nivel](pos+1); // posicion en M_B donde esta el x-ecimo parentesis abierto
    retorno = B_select1[nivel](retorno+1); //posicion en B del parentesis abierto que representa X
	retorno = m_B_rank[0](retorno);
//return mapeo[m_vertices*nivel+retorno]; 
		return t.getposMapeo(m_vertices*nivel+retorno);
	}
	int getReverso(size_type pos, unsigned int nivel){
//return reverso[m_vertices*nivel+pos];
		return t.getposReverso(m_vertices*nivel+pos);
	}
//Retorna ID region correspondiente
        unsigned int rankwt(unsigned int posx, unsigned int x){
//	puts("ENTREEEE");
        int retorno = posx;
        tuple<int,int> mitad = wt.lex_smaller_count(posx,x);
//      retorno -= get<0>(mitad);
        retorno -= get<1>(mitad); 
        return retorno;
        }   



        unsigned int rankwt2(unsigned int posx, unsigned int x){
double time;
contadores+=1.0;

	int pos = 0,vactual = 0,vactual2 = 0,romper = 0,divisor = 2, total =pow(2,m_height[0])*m_bs;
   for(int i = 0; i < m_height[0]; i++){
	total = total >> 1;
                if(total + vactual2 >= posx){
                        pos = (pos<<1)+1;
                }
                else{
                        vactual+= resumen[x][2*pos+1];
                        vactual2+= total;
                        pos = (pos<<1)+2;
                }
//	divisor= divisor << 1;
        }

int actual = vactual;
int posactual =  (pos-posinicial2)*m_bs;
pos = posactual;
//pos += 142;

//cout << "LLEGUE CON "<< posactual << "  " << vactual2 << "  " <<B[0].size() <<"  " << posx << "  " <<pos <<endl;
for( /* int i = 0; i < m_bs &&*/; posactual < posx;/* i++,*/posactual++){
        if(plano[posactual] >= x){
                actual++;
                pos = posactual;
        }
}
 //       time = (double(t1-t0)/1000);
//unsigned int actual2 = rankwt(posx,x);
//if(actual != actual2)cout<<"ACTUAL = " << actual << "  " << actual2 << endl;


        return actual;
        }


	unsigned int selectwt(unsigned int x,unsigned int posx){
/*
double time;
unsigned t0, t1;
contadoreh += 1.0;
t0 = clock();
*/

/*
    unsigned int r = B[0].size() -1,l = 0,m;
    while(r > l){
		m = (l+r)/2;
		unsigned int mitad = rankwt(m+1,x);		if(mitad >= posx)r = m;
		else l = m+1;
	}

*/
int pos = 0,vactual = 0,romper = 0;
/*
unsigned t0, t1;
double time;
        t0=clock();
*/
//if(x == 1){
//puts("entre al if");
long long int mov;
	for(int i = 0; i < m_height[0]; i++){
		mov = pos<<1;
		if(resumen[x][mov+1]+vactual >= posx){
			pos = mov+1;
		}
		else{
			vactual+= resumen[x][mov+1];
			pos = mov+2;
		}
	}
/*  
      t1=clock();
        time = (double(t1-t0)/1000);
        tiempoh +=time;
*/
int actual = vactual;
//int actual = 0;
int posactual =  (pos-posinicial2)*m_bs;
//cout << actual << "  "  << posx  <<"  " << pos <<"  "  << posinicial<< endl;
//if(resumen[pos] > 0) actual = resumen[pos];
//else actual = resumen[pos];
pos = posactual;


    unsigned int r = posactual+ m_bs,l = posactual,m,totalentro = 0;
//	t0 = clock();
/*
while(r > l){
	
	
}
*/


    while(r > l){
                m = (l+r)/2;
                unsigned int mitad = rankwt(m+1,x);
                if(mitad >= posx)r = m;
                else l = m+1;
//		totalentro++;
        }


/*
        t1=clock();
        time = (double(t1-t0)/1000);
//	if(time > tiempohm)tiempohm = time;
	tiempohm+= time;
	if(time < tiempoh)tiempoh = time;
//        tiempos +=time;

//cout << "ENTRE = "<< totalentro<<endl;

*/
/*
int maximo = B[0].size();
for(int i = 0; i < m_bs && posactual < maximo; i++,posactual++){
	if(actual == posx)break;
	if(plano[posactual] >= x){
		actual++;
		pos = posactual;
	}
//	romper++;
//	if(romper == 4)break;
}
*/
//}


/*
        t1=clock();
        time = (double(t1-t0)/1000);
//      if(time > tiempohm)tiempohm = time;
        tiempohm+= time;
        if(time < tiempoh)tiempoh = time;
*/
/*
long long int lo = 0, hi = B[0].size(),pos,ranklo,rankhi,rankpos;
int retorno = -1;
        ranklo = rankwt(lo,x);
        rankhi = rankwt(hi,x);
  while (lo <= hi && posx >= ranklo && posx <= rankhi+1){
	if(ranklo == rankhi)pos = lo;	
	else pos = lo + (((double)(hi - lo) /
           (rankhi - ranklo)) * (posx - ranklo));
//	if(pos > hi) pos = hi;
//	cout << lo << "  " << pos << "  " << hi << "  " <<l << endl; 
	rankpos = rankwt(pos,x);
	if(rankpos == posx)return rankpos;
       if (rankpos < posx){
		lo = pos + 1;
		ranklo = rankpos;
		}
        else{ 
		hi = pos-1;
		rankhi = rankpos;
		}
    }
	retorno = hi;
        if(hi > 0){
            if (rankwt(hi-1,x) == posx){
                retorno = hi-1;
                }
        }
*/

//if(l != pos ) cout << "ERROR  " << l << "  " << pos << "  "  << posx <<endl; 
		return l;
	
//return pos;
	}



	unsigned int cselectwt(unsigned int x,unsigned int posx){
int pos = 0,vactual = 0,romper = 0;
long long int mov;
	for(int i = 0; i < cm_height[0]; i++){
		mov = pos<<1;
		if(cresumen[x][mov+1]+vactual >= posx){
			pos = mov+1;
		}
		else{
			vactual+= cresumen[x][mov+1];
			pos = mov+2;
		}
	}
int actual = vactual;
int posactual =  (pos-posinicial3)*m_bs;
//if(posx == 19)cout << "POSICION = " << pos << " " << posinicial3 << " " << posactual << endl;
pos = posactual;
    unsigned int r = posactual+ m_bs,l = posactual,m,totalentro = 0;
    while(r > l){
                m = (l+r)/2;
                unsigned int mitad = crankwt(m+1,x);
                if(mitad >= posx)r = m;
                else l = m+1;
        }
		return l;
	}


        unsigned int crankwt(unsigned int posx, unsigned int x){
//	puts("ENTREEEE");
        int retorno = posx;
        tuple<int,int> mitad = wt2.lex_smaller_count(posx,x);
//      retorno -= get<0>(mitad);
        retorno -= get<1>(mitad); 
        return retorno;
        }   





    unsigned int regionID(unsigned int x, unsigned int level){
	return m_B_rank[level](x);
	}
	//Retorna primer parentesis en nivel 0 contenido por la region entregada
    unsigned int goUp(unsigned int x, unsigned int level){
    //cout << x << " " << level << endl;
//	puts("entre");

/*
    int posicion = m_B_select1[level](x+1); // posicion en M_B donde esta el x-ecimo parentesis abierto
   posicion = selectwt(level,posicion+1);
*/

	int posicion = m_B_select1[level](x+1);	
	posicion = selectwt(level,posicion+1);



//   if(level > 1) posicion2 = B_select1[level](posicion2+1);
  // 	else posicion2 = BP_select1[level](posicion2+1); //posicion en B del parentesis abierto que representa X
    //puts("retorno");
//	if(posicion != posicion2)std::cout <<"posiciones = " << posicion << "	" << posicion2 << "  datos  " << x  <<"  "  << level <<   "\n";
//	puts("sali");
	return posicion;
	}




	//Retorna region que contiene a la region entregada, region entregada siempre pertenece al nivel 0
    unsigned int goDown(unsigned int x, unsigned int level){
    	int posB = m_B_select1[0](x+1);

  //    	if(level > 1){
      		if(plano[posB] >= level){
  //  		int retornar2 = B_rank[level](posB); //Obtengo posicion en m_B
      		int retornar = rankwt(posB,level);
//		if(retornar2 != retornar) cout << "ERROR1  " << retornar  << "  " << retornar2<< endl;
		  bool pab = (m_B[level][retornar] == 1) ? true: false; //Compruebo si parentesis es abierto o cerrado 
		  if(pab)return retornar;
		  else {
    			retornar = m_B_st[level].find_open(retornar); // obtengo el parentesis que abre
    			return retornar;
		   }
      	}
    	else{
    		int sumar2 = rankwt(posB,level);
    		int parentesiscercano =selectwt(level,sumar2+1); //Posicion siguiente parentesis en B
	//		int parentesiscercano2 =B_select1[level](B_rank[level](posB)+1);
//		if(parentesiscercano != parentesiscercano2) cout << "ERROR2  " << parentesiscercano  << "  " << parentesiscercano2<< endl;
    		sumar2 = rankwt(parentesiscercano,level);
//		if(sumar2 != rankwt(parentesiscercano,level))cout<<"ERROR RANKT1: " << sumar2 <<"  " << rankwt(parentesiscercano,level) <<endl;
    		int retornar = sumar2; // Posicion del parentesis siguiente en m_b
    		bool abierto = (m_B[level][retornar] == 1) ? true: false;  //Compruebo si es abierto o cerrado
    		int sumar3 = rankwt(posB,level);
    		int sumar4 = rankwt(B[0].size()-1,level);  		

       //     if(sumar3 != rankwt(posB,level))cout<<"ERROR RANKT2: " << sumar3 <<"  " << rankwt(posB,level) <<endl;
        //    if(sumar4 != rankwt(B[0].size()-1,level))cout<<"ERROR RANKT3: " << sumar4 <<"  " << rankwt(B[0].size()-1,level) <<endl;

  //            if(sumar3 != B_rank[level](posB)) cout << "ERROR3  " <<sumar3  << "  " << B_rank[level](posB)<< endl;

//           if(sumar4 != B_rank[level](B[level].size()-1)) cout << "ERROR4  " << sumar4  << "  " << B_rank[level](B[level].size()-1)<< endl;

    		if(abierto && sumar3 != sumar4 ){
    			retornar = m_B_st[level].parent_t(retornar); //en caso de ser abierto obtengo el padre en m_B
    			return retornar; //retorno directamente
    		}
    		else {
    		retornar = m_B_st[level].find_open(retornar); // Caso de ser cerrado obtengo parentesis que abre
    		return retornar;		
    		}


    	}
/*   }
		else {
		if(BP[level][posB] == 1){
    		int retornar = BP_rank[level](posB); //Obtengo posicion en m_B
		  	bool pab = (m_B[level][retornar] == 1) ? true: false; //Compruebo si parentesis es abierto o cerrado 
		   if(pab){
		//   	      		puts("entre 1");
		   	 return retornar;
		   }
		   
		   else {
    		retornar = m_B_st[level].find_open(retornar); // obtengo el parentesis que abre
    		return retornar;
		   }   
		    	
      	}  	
    	else { 
    		int parentesiscercano =BP_select1[level](BP_rank[level](posB)+1); //Posicion siguiente parentesis en B
    		int retornar = BP_rank[level](parentesiscercano); // Posicion del parentesis siguiente en m_b
    		bool abierto = (m_B[level][retornar] == 1) ? true: false;  //Compruebo si es abierto o cerrado

    		if(abierto && BP_rank[level](posB) != BP_rank[level](BP[level].size()-1)){
    			retornar = m_B_st[level].parent_t(retornar); //en caso de ser abierto obtengo el padre en m_B
    			return retornar; //retorno directamente
    		}
    		else {
    		retornar = m_B_st[level].find_open(retornar); // Caso de ser cerrado obtengo parentesis que abre
    		return retornar;			
    		}


    	}	
		}
*/
	}


    unsigned int goLevel(unsigned int x, unsigned int levels, unsigned int levelt){
    unsigned int regions = regionID(goUp(x,levels),0);
    unsigned int regiont = regionID(goDown(regions,levelt),levelt);
   //         cout << regions << endl;
   //			 return -1;
	return regiont;
	}

	bool Inside(unsigned int x, unsigned int levelx, unsigned int y, unsigned int levely){
		//puts("entre");
		unsigned int comprobar = goLevel(x,levelx,levely);
		//puts("sali");
		if(comprobar == y)return true;
		else return false;
	}

	bool Touches(unsigned int x, unsigned int levelx, unsigned int y, unsigned int levely){
//		puts("ENTRE AL TOUCH");
		if(x == y && levelx == levely)return true;
		if(levely < levelx){
			bool comprobar = Inside(y,levely,x,levelx);
			vector <int> vecinos = (lista_vecinos2(y,levely));
			for(int i = 0; i < vecinos.size();i++){
				unsigned int contenedor = goLevel(vecinos[i],levely,levelx);
				if(contenedor == x  && !comprobar)return true;
				if(contenedor != x  && comprobar)return true;
	
			}
		}
		else if(levely == levelx){
			vector <int> vecinos = (lista_vecinos2(y,levely));
			for(int i = 0; i < vecinos.size();i++){
				if(vecinos[i] == x )return true;
			}


		}
		else {
			bool comprobar = Inside(x,levelx,y,levely);
			vector <int> vecinos = (lista_vecinos2(x,levelx));
	//		cout << "tam de vecinos " << vecinos.size() << " " << x << " " << levelx << endl;
			for(int i = 0; i < vecinos.size();i++){
	//							cout << "vecino " << vecinos[i] <<endl;
				unsigned int contenedor = goLevel(vecinos[i],levelx,levely);
				if(contenedor == y  && !comprobar)return true;
				if(contenedor != y  && comprobar)return true;

			}
		}
		return false;
	} 

	vector <int> compvecinos(unsigned int x, int levelx){
		vector <int> vecinos = (lista_vecinos(x,levelx));
		return vecinos;

	}




	void frecacumulada(){
	for(int i = 3; i < 4; i++){
		for(int j = 1; j <= B[0].size();j++)printf("%d ",rankwt(j,i));
		}
	puts("");
	}

	double obtenertiempoh(){
		double tiempoaux = tiempoh;
		tiempoh = 1000.0;
		return tiempoaux;
	}
        double obtenertiempos(){
                double tiempoaux = tiempos;
                tiempos = 1000.0;
                return tiempoaux;
        }       



        double obtenertiempohm(){
                double tiempoaux = tiempobloque/contadores;
	contadoreh = 0.0;
                tiempoh = 0.0;
		contadoreh = 0.0;
                return tiempoaux;
        }
        double obtenertiemposm(){
                double tiempoaux = tiemposm/contadores;
		contadores =0.0;
                tiempos = 0.0;
                return tiempoaux;
        }


	vector <int> Contained(unsigned int x, int levelx, int levely){
		unsigned int a = goLevel(x,levelx,levely);
		size_type nxt = first(a,levely);
		size_type limit = mate(nxt,levely);
		size_type limit2 = mate(first(0,levely),levely);
		limit2 = m_B_rank[levely](m_A_rank[levely](limit2));
		vector <int> retorno;
		retorno.push_back(vertex2(nxt,levely));
		size_type pos_in_B = m_B_rank[levely](m_A_rank[levely](nxt))+2;
		if(pos_in_B > limit2)return retorno;
		nxt = m_A_select1[levely](m_B_select1[levely](pos_in_B)+1);
		int contador = 0;
		while(nxt < limit){
			size_type pos_in_A = m_B_rank[levely](m_A_rank[levely](nxt));
			size_type pos = m_B_select1[levely](pos_in_A+1);
			if(levely > 1)pos = B_select1[levely](pos+1);
			else pos = BP_select1[levely](pos+1);
			if(levelx > 1){
				if(B[levelx][pos] ==1){
				nxt = mate(nxt,levely);
				size_type pos_in_B = m_B_rank[levely](m_A_rank[levely](nxt))+1;
				if(pos_in_B > limit2)break;
		nxt = m_A_select1[levely](m_B_select1[levely](pos_in_B)+1);
			}
			else if(B[levelx][pos] == 0){
					//			cout << vertex2(nxt,levely) << endl;
				retorno.push_back(vertex2(nxt,levely));
				size_type pos_in_B = m_B_rank[levely](m_A_rank[levely](nxt))+2;
				if(pos_in_B > limit2)break;
				nxt = m_A_select1[levely](m_B_select1[levely](pos_in_B)+1);
			} 
		}
		else{
				if(BP[levelx][pos] ==1){
				nxt = mate(nxt,levely);
				size_type pos_in_B = m_B_rank[levely](m_A_rank[levely](nxt))+1;
				if(pos_in_B > limit2)break;
		nxt = m_A_select1[levely](m_B_select1[levely](pos_in_B)+1);
			}
			else if(BP[levelx][pos] == 0){
					//			cout << vertex2(nxt,levely) << endl;
				retorno.push_back(vertex2(nxt,levely));
				size_type pos_in_B = m_B_rank[levely](m_A_rank[levely](nxt))+2;
				if(pos_in_B > limit2)break;
				nxt = m_A_select1[levely](m_B_select1[levely](pos_in_B)+1);
			} 			
		} 

			contador++;

		}
		return retorno;

	}


        vector <int> Contained2(unsigned int x, int levelx, int levely){
                unsigned int a = goLevel(x,levelx,levely);
		int lugar = 0;
                size_type nxt = first(a,levely),nxt2 = first2(a,levely)+2,nxt3,nxt4,meter;
		
//		cout << nxt2  << " " << m_A_select1[levely](nxt2) <<" "  << nxt  <<endl;
	if(nxt != selectmBtoA(nxt2-2,levely))  cout << nxt2  << " " << m_A_select1[levely](nxt2-1) <<" "  << nxt  <<" "  << selectmBtoA(nxt2-2,levely) << endl;
                size_type limit = mate(nxt,levely);
		size_type limit3 = mate3(nxt,levely,1,nxt2-2).first;
		if(nxt2-2 != 	m_A_rank[levely](nxt))cout << "LIMITE: " << nxt2-2<< " " <<  m_A_rank[levely](nxt) << endl; 
		int sumar = 0;
                size_type limit2 = mate(first(0,levely),levely);
                limit2 = m_B_rank[levely](m_A_rank[levely](limit2));
                vector <int> retorno;
//		return retorno;
                retorno.push_back(vertex2(nxt,levely));
                size_type pos_in_B = m_B_rank[levely](m_A_rank[levely](nxt))+2;
//		cout <<"PRUEBA: " <<first2(a,levely) << " " << m_A_rank[levely](nxt) << " " << pos_in_B << " "  << m_B_rank[levely](first2(a,levely))+2   <<endl; 
		if(pos_in_B != rankmB(nxt2-2,levely)+2 )cout << "MALO: " << pos_in_B << " " << rankwtba( selectwt(levely,nxt2-2),levely)+2  <<" "  << m_A_rank[levely](nxt)   << " " <<nxt2-2   <<  
        " " <<   selectwt(levely,nxt2-1)    << " "  <<  selectwtba(levely,pos_in_B)  << " "  << rankwtba(selectwtba(levely,pos_in_B),levely)  << " " << levely <<endl;
                if(pos_in_B > limit2)return retorno;
//			return retorno;
  	              nxt = m_A_select1[levely](m_B_select1[levely](pos_in_B)+1);
//			return retorno;
			nxt4 = mBselect(pos_in_B,levely);
		    nxt3 =  selectmBtoA(mBselect(pos_in_B,levely),levely);
			meter = nxt3;
		if(m_B_select1[levely](pos_in_B) != /*rankwt(selectwtba(levely,pos_in_B),levely)*/ mBselect(pos_in_B,levely)) cout 
<< m_B_select1[levely](pos_in_B) << " " <</* rankwt(selectwtba(levely,pos_in_B),levely)+1 <<" "  <<*/  mBselect(pos_in_B,levely) << " " << levely <<endl;
			if(nxt != nxt3) cout << nxt <<" " << nxt3 << " " << m_B_select1[levely](pos_in_B) << " " << selectwtba(levely,nxt2) << " " << levely << " "   << 
   rankwt(selectwtba(levely,pos_in_B),levely)  <<endl;
//		cout << "Prueba: " << nxt<< " " <<   m_B_select1[levely](pos_in_B)+1  <<" " <<  pos_in_B <<" " << m_A_rank[levely](nxt) <<endl;
                int contador = 0;
	if(m_B_select1[levely](pos_in_B) != m_A_rank[levely](nxt)) cout << "SEGUNDO MALO\n";
//cout << "ENTRO AL WHILE " << endl;
//return retorno;
                while(true){

if(nxt >= limit){
//cout << "ROMPO IF: " << limit3  << " " <<  limit   <<" "  << nxt << " " << selectmBtoA(limit3,levely)  <<endl ;
break;
}
//cout << "ENTRO AL WHILE " << nxt  <<endl;
 if(m_B_select1[levely](pos_in_B) != m_A_rank[levely](nxt)) cout << "SEGUNDO MALO " <<  m_B_select1[levely](pos_in_B) << " "<<  m_A_rank[levely](nxt)<< endl;
                        size_type pos_in_A = m_B_rank[levely](m_A_rank[levely](nxt));
		//	if(pos_in_A != nxt4) cout << "nxt1 " << nxt << " " << nxt4 << " " << pos_in_A <<endl;
			if(pos_in_A != pos_in_B-1 ) cout <<"IF1 = " << pos_in_A << " " << m_A_rank[levely](nxt) << " "  << pos_in_B <<endl;
			lugar = 1;
                        size_type pos = m_B_select1[levely](pos_in_A+1);
			pos = selectwt(levely,pos+1);
//                        if(levely > 1)pos = B_select1[levely](pos+1);
//                        else pos = BP_select1[levely](pos+1);
               //         if(levelx > 1){
				if(plano[pos] >= levelx){
//                                if(B[levelx][pos] ==1){
//      	 	cout << "mate1" << endl;
//		nxt4 = nxt2;
//if(nxt != m_A_select1[levely](nxt2))  cout << "distintos"<<endl;//<< nxt2  << " " << m_A_select1[levely](nxt2) <<" "  << nxt  <<endl;
				sumar = 0;
				nxt2 = mate3(nxt,levely,1,nxt2).first+1;
				nxt3 = nxt2;
  //    	 	cout << "mate2" << endl;
                                nxt = mate2(nxt,levely,1);
if(nxt != m_A_select1[levely](nxt2)) cout<<"NXT= "   << nxt << " " << nxt2 << endl; 
//				nxt3 = mate3(nxt,levely,1,nxt4).first+1;
                                 pos_in_B = m_B_rank[levely](m_A_rank[levely](nxt))+1;
//				if (m_A_rank[levely](nxt) != nxt)cout << "RANK: " << nxt<< " " << m_A_rank[levely](nxt) << endl;
                                if(pos_in_B > limit2){
cout << "ENTRE A POS:IN;B" << endl;
break;
}
 if(pos_in_B != m_B_rank[levely](nxt2)+1) cout << "PRIMERO MALO 1 " << pos_in_B  << " " <<m_B_rank[levely](nxt2)+1   << endl;  
                nxt = m_A_select1[levely](m_B_select1[levely](pos_in_B)+1);
		 nxt3 =  selectmBtoA(mBselect(pos_in_B,levely),levely);
		   if(nxt != nxt3) cout << "nxt22 " << nxt << " " << nxt3 << endl;
		nxt4 =  mBselect(pos_in_B,levely);
		meter = selectmBtoA (mBselect(pos_in_B,levely),levely);
  // if(nxt != nxt4) cout << "nxt2 " << nxt << " " << nxt4 << endl;
		if(nxt !=  selectmBtoA (mBselect(pos_in_B,levely),levely) ) cout << nxt << " " << selectmBtoA (mBselect(pos_in_B,levely),levely) << " " << mBselect(pos_in_B,levely) 
		<< " " << m_B_select1[levely](pos_in_B) << endl;
	lugar = 1;
                        }
else if(plano[pos] < levelx){
  //                      else if(B[levelx][pos] == 0){
                                        //                      cout << vertex2(nxt,levely) << endl;
 if(nxt !=  meter ) cout <<"PUSH "  <<nxt << " " << meter << endl;
                                retorno.push_back(vertex2(nxt,levely));
                                 pos_in_B = m_B_rank[levely](m_A_rank[levely](nxt))+2;
//int sumar = 0;
//if(lugar == 1)sumar = 1;
// if(pos_in_B != m_B_rank[levely](nxt2)+sumar) cout << "PRIMERO MALO 2 " <<  pos_in_B  << " " <<m_B_rank[levely](nxt2+1)+1  << " " <<m_B_rank[levely](nxt2+1)   <<" " << lugar << " "  <<sumar    <<endl;  
                                if(pos_in_B > limit2)break;
				nxt2 = m_B_select1[levely](pos_in_B)+1;
                                nxt = m_A_select1[levely](m_B_select1[levely](pos_in_B)+1);
				nxt4 =  mBselect(pos_in_B,levely);
				meter = selectmBtoA (mBselect(pos_in_B,levely),levely);
  // if(nxt != nxt4) cout << "nxt3 " << nxt << " " << nxt4 << endl;
              if(nxt !=  selectmBtoA (mBselect(pos_in_B,levely),levely) ) cout << nxt << " " << selectmBtoA (mBselect(pos_in_B,levely),levely) << " " << mBselect(pos_in_B,levely) << endl;

//				nxt3 = ;
//if(nxt != m_A_select1[levely](nxt2))  cout <<"if malo: " << nxt2  << " " << m_A_select1[levely](nxt2) <<" "  << nxt  <<endl;
//lugar = 2;
                        }
              //  }
     /*           else{
                                if(BP[levelx][pos] ==1){
                                nxt = mate(nxt,levely);
                                size_type pos_in_B = m_B_rank[levely](m_A_rank[levely](nxt))+1;
                            //    if(pos_in_B > limit2)break;
                nxt = m_A_select1[levely](m_B_select1[levely](pos_in_B)+1);
                        }
                        else if(BP[levelx][pos] == 0){
                                        //                      cout << vertex2(nxt,levely) << endl;
                                retorno.push_back(vertex2(nxt,levely));
                                size_type pos_in_B = m_B_rank[levely](m_A_rank[levely](nxt))+2;
                                if(pos_in_B > limit2)break;
                                nxt = m_A_select1[levely](m_B_select1[levely](pos_in_B)+1);
                        }
                }*/

                        contador++;

                }
//		cout << "end" << endl;
                return retorno;

        }

size_type rankmB(size_type v, int nivel){

if(nivel == 0) return m_B_rank[0](v);
else return rankwtba(selectwt(nivel,v+1) ,nivel);
}


        unsigned int rankwtba(unsigned int posx, unsigned int x){
//	puts("ENTREEEE");
        int retorno = posx;
        tuple<int,int> mitad = wtba.lex_smaller_count(posx,x);
//      retorno -= get<0>(mitad);
        retorno -= get<1>(mitad); 
        return retorno;
        }   


        unsigned int crankwtba(unsigned int posx, unsigned int x){
//      puts("ENTREEEE");
        int retorno = posx;
        tuple<int,int> mitad = cwtba.lex_smaller_count(posx,x);
//      retorno -= get<0>(mitad);
        retorno -= get<1>(mitad);
        return retorno;
        }


	unsigned int selectwtba(unsigned int x,unsigned int posx){
int pos = 0,vactual = 0,romper = 0;
long long int mov;
	for(int i = 0; i < m_heightba[0]; i++){
		mov = pos<<1;
		if(resumen2[x][mov+1]+vactual >= posx){
			pos = mov+1;
		}
		else{
			vactual+= resumen2[x][mov+1];
			pos = mov+2;
		}
	}
int actual = vactual;
int posactual =  (pos-posinicial4)*m_bs;
pos = posactual;
    unsigned int r = posactual+ m_bs,l = posactual,m,totalentro = 0;
    while(r > l){
                m = (l+r)/2;
                unsigned int mitad = rankwtba(m+1,x);
                if(mitad >= posx)r = m;
                else l = m+1;
        }
		return l;
	}





        unsigned int cselectwtba(unsigned int x,unsigned int posx){
int pos = 0,vactual = 0,romper = 0;
long long int mov;
        for(int i = 0; i < cm_heightba[0]; i++){
                mov = pos<<1;
                if(cresumen2[x][mov+1]+vactual >= posx){
                        pos = mov+1;
                }
                else{
                        vactual+= cresumen2[x][mov+1];
                        pos = mov+2;
                }
        }
int actual = vactual;
int posactual =  (pos-posinicial5)*m_bs;
pos = posactual;
    unsigned int r = posactual+ m_bs,l = posactual,m,totalentro = 0;
    while(r > l){
                m = (l+r)/2;
                unsigned int mitad = crankwtba(m+1,x);
                if(mitad >= posx)r = m;
                else l = m+1;
        }
                return l;
        }







size_type selectmBtoA(int e, int nivel){

//Obtengo rank de m_B[nivel][pos]
//size_type posmB = rankwt(e,nivel); //inclusivo rank, creo que si?
//Paso desde m_B[nivel] a m_B[0]
size_type m_B0 = selectwt(nivel,/*posmB+1*/ e+1);
//Obtengo rank de m_B[nivel][pos]
size_type posmB = rankwt(m_B0,nivel);
//Paso de m_B[0] hasta m_A[0]
size_type posA = m_A_select1[0](m_B0+1);

//Obtengo cantidad de 0
size_type cantmBstar = posA - m_B0;
//cuento cuantos 0 pertenecen al nivel buscado
size_type cnivel = crankwt(cantmBstar,nivel);
// Respuesta es la suma de 1 perteneciente al nivel + suma de 0's
size_type retorno = cnivel +  e;
//cout << cnivel <<" " <<m_B_rank[nivel](e)   << " "  << posmB  <<  endl;
return retorno;
}




vector <int> Contained3(unsigned int x, int levelx, int levely){
        unsigned int a = goLevel(x,levelx,levely);
		int lugar = 0;
        size_type nxt ,nxt2 = first2(a,levely)+2,nxt3,nxt4,meter;
        nxt = selectmBtoA(nxt2-2,levely);
        size_type limit3 = mate3(nxt,levely,1,nxt2-2).first;
        size_type limit = selectmBtoA(limit3,levely) ;
		int sumar = 0;
        size_type limit2 = mate(first(0,levely),levely);
        limit2 = m_B_rank[levely](m_A_rank[levely](limit2));
        vector <int> retorno;
        retorno.push_back(vertex2(nxt,levely));
        size_type pos_in_B =rankmB(nxt2-2,levely)+2  ;
                if(pos_in_B > limit2)return retorno;
             nxt3 =  selectmBtoA(mBselect(pos_in_B,levely),levely);
  	              nxt = nxt3;
			nxt4 = mBselect(pos_in_B,levely);
			meter = nxt3;
                int contador = 0;
                while(true){
if(nxt >= limit){
break;
}
    size_type pos_in_A = pos_in_B-1;
	lugar = 1;
    size_type pos =  mBselect(pos_in_A+1,levely); 
	pos = selectwt(levely,pos+1);
	if(plano[pos] >= levelx){
			sumar = 0;
			nxt2 = mate3(nxt,levely,1,nxt2).first+1;
			nxt3 = nxt2;
            nxt =selectmBtoA(nxt2-1,levely); 
        pos_in_B = rankmB(nxt2,levely)+1; 
        if(pos_in_B > limit2){
		break;
		}
		 nxt3 =  selectmBtoA(mBselect(pos_in_B,levely),levely);
		 nxt = nxt3;
		nxt4 =  mBselect(pos_in_B,levely);
		meter = selectmBtoA (mBselect(pos_in_B,levely),levely);
	lugar = 1;
                        }
else if(plano[pos] < levelx){
	retorno.push_back(vertex2(meter,levely)); //CAMBIAR A VERTEX3
      pos_in_B = m_B_rank[levely](m_A_rank[levely](nxt))+2;
              if(pos_in_B > limit2)break;
				nxt2 =mBselect(pos_in_B,levely)+1;
                nxt =  selectmBtoA(mBselect(pos_in_B,levely),levely);
				nxt4 =  mBselect(pos_in_B,levely);
				meter = selectmBtoA (mBselect(pos_in_B,levely),levely);
                        }
                        contador++;
                }
                return retorno;

        }





size_type cantCeros(int e, int nivel){
size_type m_B0 = selectwt(nivel,/*posmB+1*/ e+1);
//Obtengo rank de m_B[nivel][pos]
//size_type posmB = rankwt(m_B0,nivel);
//Paso de m_B[0] hasta m_A[0]
size_type posA = m_A_select1[0](m_B0+1);

//Obtengo cantidad de 0
size_type cantmBstar = posA - m_B0;
size_type cnivel = crankwt(cantmBstar,nivel);
return cnivel;
}



size_type cont1Ma(int e, int nivel){

//Obtengo rank de m_B[nivel][pos]
//size_type posmB = rankwt(e,nivel); //inclusivo rank, creo que si?
//Paso desde m_B[nivel] a m_B[0]
size_type m_B0 = cselectwt(nivel,/*posmB+1*/ e+1);
//Obtengo rank de m_B[nivel][pos]
//size_type posmB = rankwt(m_B0,nivel);
//Paso de m_B[0] hasta m_A[0]
size_type posA = m_A_select0[0](m_B0+1);

//Obtengo cantidad de 0
size_type cantmBstar = posA - m_B0;
//cuento cuantos 0 pertenecen al nivel buscado
size_type cnivel = rankwt(cantmBstar,nivel);
return cnivel;
}





      vector <int> lista_vecinos2(size_type v, int nivel){
      	size_type nxt =  selectmBtoA(first2(v,nivel),nivel)+1;
	pair <size_type,int > nxt2; 
      	pair<size_type,int> limite = mate3(first(v,nivel),nivel,1, first2(v,nivel)+1);
	int posini2 = first(v,nivel);
      	vector <int> retorno;
	pair<size_type,int> posini = make_pair(first2(v,nivel),1);
	int simb = 1;
int valordeprueba = 0;
int valorinsert = 0;
int tipoinsert = 1;
size_type nxtimp = nxt, rankposini = posini.first+1;
valorinsert = rankwt(m_A_rank[0]( MapeoA(posini2,nivel,0,0).first),nivel)-1;
nxt2 = make_pair(valorinsert , 1 );
if(MapeoA(posini2,nivel,0,0).second == 0)
valordeprueba = crankwt(MapeoA(posini2,nivel,0,0).first-m_A_rank[0]( MapeoA(posini2,nivel,0,0).first),nivel);
bool inicial = true;
size_type valorrank = first2(v,nivel), tipoA = selectmBtoA(first2(v,nivel)+1,nivel);
      	while(nxt <  selectmBtoA(limite.first,nivel) ){
      		size_type x;
      		size_type pos_in_A =  valorrank+1;
      		if(tipoA == nxt)pos_in_A++;
      		if(tipoA == nxt &&  Eat(pos_in_A,nivel)){
      			x = vertex4(nxt,nivel,MapeoA(posini2,nivel,0,0).second,rankposini);
posini = mate3(nxt,nivel,1,MapeoA(posini2,nivel,nxt2.first,valordeprueba/*0,0*/).first+1);
posini2 = selectmBtoA(posini.first,nivel);
simb = 1;
nxt =  selectmBtoA(posini.first,nivel)+1;
rankposini = posini.first+1;
tipoA = selectmBtoA(posini.first+1,nivel);
valorrank = posini.first;
valorinsert = rankwt(m_A_rank[0]( MapeoA(posini2,nivel,posini.first,posini.second).first),nivel)-1;
valordeprueba = crankwt(MapeoA(posini2,nivel,posini.first,posini.second).first-m_A_rank[0]( MapeoA(posini2,nivel,posini.first,posini.second).first),nivel);
nxt2 = make_pair(valorinsert , 1 );
      		}
      		else{
      			x = vertex4(nxt,nivel,MapeoA(posini2,nivel,0,0).second,rankposini);    	
posini = make_pair(nxt,MapeoA(posini2,nivel,nxt2.first,valordeprueba/*0,0*/).second);
valorinsert = rankwt(m_A_rank[0]( MapeoA(posini2,nivel,nxt2.first,valordeprueba).first),nivel)-1;
nxt2 = make_pair(valorinsert , 1 );
posini2 = nxt;
      			nxt ++;
if(nxt < mate2(first(v,nivel),nivel,1))if(MapeoA(posini2,nivel,nxt2.first,valordeprueba).second == 0)
valordeprueba = crankwt(MapeoA(posini2,nivel,0,0).first-m_A_rank[0]( MapeoA(posini2,nivel,0,0).first),nivel);
	}
      		retorno.push_back(x);
      	}
int malo = 0;
      	if(v != 0){
      		nxt = bsearch(m_B_select1[nivel](v+1),nivel,-2);
      		int auxtemp = nxt;
      		size_type x = rankmb(auxtemp,nivel);
      		if(Eat(auxtemp+1,nivel) == 1)x++;
      		if(x-1 < limite.first)retorno.push_back(x-1);
      	}
      	return retorno;
      }





size_type vertex4(size_type e, int nivel, int inicial, int nposA){
      size_type pos_in_A = m_A_rank[nivel](e);  
	  if(inicial) pos_in_A++;      	
      	if(Eat(pos_in_A,nivel) == 1 && inicial == 1){
      return rankmb(pos_in_A,nivel)-1;
      	}
      	else if( m_A[nivel][e] == 1){
      		size_type match_pos = bsearch(pos_in_A,nivel,0);
      		return rankmB(match_pos,nivel);
      	}
      	size_type pos_in_B = e- pos_in_A+1;

      	if(m_B_star[nivel][pos_in_B-1] != Eat2(pos_in_B,nivel)) cout << "EAT2 MALO\n";
      	 if(/*m_A[nivel][e] == 0 &&*/ m_B_star[nivel][pos_in_B-1] == 1) { //
     		size_type match_pos = selectmBStartoA(mate3(e,nivel,inicial,nposA).first,nivel);
     		size_type pos = cont1Ma(mate3(e,nivel,inicial,nposA).first,nivel);
     		size_type pos2 =  selectmBtoA(rankwt(selectwt(nivel,pos+1),nivel),nivel); 
		if(pos+1 == m_B[nivel].size())pos2 = ultimos[nivel];
int pos2a = pos; 
int comparador,sientro = 0;
if(pos2a+1 >= m_B[nivel].size()){
comparador = 0;
sientro = 1;
}
else comparador = m_B[0][(selectwt(nivel,pos2a+1))];
     		if(comparador == 1){
     			pos2 = bsearch(pos2a,nivel,-2); 
     			size_type pos4 =  selectmBtoA(pos2,nivel); 
				size_type pos3 = pos2;  
				if(pos2 != pos3) cout <<"pos3 malo: " << pos2 << " "<< pos3 <<endl; 
				pos3++;
				return rankmb(pos3-1,nivel);
     		}
     		else {
     			size_type posaux = pos2;
     			pos2 =selectmBtoA(mate3(posaux,nivel,1,pos).first,nivel); 
     			size_type pos3 =mate3(posaux,nivel,1,pos).first; 
     			pos3++;
     			return rankmb(pos3-1,nivel);
     		}
      		
      	}
	else if( m_B_star[0][cselectwt(nivel,pos_in_B)]== 0){ 

     		size_type match_pos =selectmBStartoA(mate3(e,nivel,inicial,nposA).first,nivel);  /*mate2(e,nivel,inicial);*/ //ARREGLAR
     		size_type pos =cont1Ma(mate3(e,nivel,inicial,nposA).first,nivel); //m_A_rank[nivel](match_pos); //ARREGLAR (es retorno de mate3)
     		size_type pos2 =selectmBtoA(pos,nivel); 
     		if(Eat(pos+1,nivel) == 1){ 
				int pos2a = pos; 
     			pos2 =  bsearch(pos2a,nivel,-2); 
     			size_type pos4 =selectmBtoA(pos2,nivel);
				size_type pos3 = pos2;
				pos3++;
 				return rankmb(pos3-1,nivel);
     		}
     		else {
     			size_type  posauxfinal = pos2;
     			pos2 =selectmBtoA(mate3(posauxfinal,nivel,1,pos).first,nivel); 
     			size_type pos3 = mate3(posauxfinal,nivel,1,pos).first; 
     			pos3++;
     			return rankmb(pos3-1,nivel);
     		}
      	}
      	return -1;
      }




      pair<size_type,int> mate4(size_type i,int nivel, int tipoactual,int nuevaposB) {
          size_type pos_in_B = m_A_rank[nivel](i); 
          if(tipoactual == 1){
pos_in_B++;
//nuevaposB++;
}
        if(tipoactual == 1) {
          size_type match_in_B;
int realvalor;
if(pos_in_B>= m_B[nivel].size() )realvalor = 0;
else realvalor = Eat(pos_in_B,nivel);
if(realvalor == 1){
        int  match_in_B2 = fopen2(selectwt(nivel,pos_in_B),nivel);
        match_in_B = rankwt2(match_in_B2,nivel);
}
          else {
	int  match_in_B2 =  bsearch(pos_in_B-2,nivel,-1); 
	if(match_in_B2+1 == m_B[nivel].size())match_in_B2 = 0;
match_in_B =	match_in_B2;
          }
          return make_pair(match_in_B,1);
        }
        else
          {
            size_type pos_in_B_star = i - pos_in_B+1;
            size_type match_in_B_star;
            if(m_B_star[0][cselectwt(nivel,pos_in_B_star)] == 1){ 
int  match_in_B2_star = cfopen(selectwt(nivel,pos_in_B_star),nivel);
match_in_B_star =crankwt(match_in_B2_star,nivel); 
 return make_pair(match_in_B_star,0);
            }
            else{ 

int  match_in_B2_star = cfopen(selectwt(nivel,pos_in_B_star),nivel);
match_in_B_star =crankwt(match_in_B2_star,nivel);
             return make_pair(match_in_B_star,0);
            }
          }
        return make_pair(-1,-1);
      }





	size_type cbsearch(size_type i, int nivel, int obj){
	int bloque = floor(i/m_bs);
	int top_nodes = pow(2, cm_height[nivel])-1;
	int tope = bloque*m_bs;
	int contador = 0;
	int resp = 0;
	for(int l =  i; l >= tope; l--){
	int pos = cselectwt(nivel,l+1);
	if(m_B_star[0][pos] == 1)contador--;
	else contador++;
	if(contador == obj){
		resp = l;
		break;
	}

	}
	if(contador == obj) return resp;
	int actual = top_nodes+bloque;
//cout << "PASE BUSQUEDA LINEAL CON " << actual << " " << i <<endl;
	while(actual != 0 && contador - rmt_v2[nivel][actual-1] + rmt_m2[nivel][actual-1] > obj){
	if(actual %2 != 0)actual = (actual -1)/2;
	else{
/*
  	if(contador + rmt_m[0][actual-1] >= obj){
		actual--;
		break;
		}
*/
	contador   -= rmt_v2[nivel][actual-1];
	actual = (actual -2 )/2;
	}

	}


//	cout << "ACTUAL = " << actual << endl;
//AGREGAR CONDICION CON COMPARAR CON RAIZ EN CASO DE QUE VALOR NO SE ENCUENTRA
actual-=1;
	while(actual <  cposinicial[nivel]){
 	int auxactual = actual*2 +2;
 if(contador - rmt_v2[nivel][auxactual] + rmt_m2[nivel][auxactual]  <= obj){
actual = auxactual;
}
else {
contador -= rmt_v2[nivel][auxactual];
//EN CASO DE ERROR REVISAR ESTE CASO (AL IR BAJANDO Y MOVERME A LA DERECHA ACA ESTOY ASUMIENDO QUE SIEMPRE  SE CUMPLE EL IF)
actual = auxactual-1;
}
	}

//REVISAR FORMULA EN CASO DE ERROR
//cout << " auxpos " << auxpos <<"  " <<contador <<"  " << obj <<"  " <<actual  << " " <<m_B[nivel].size()-1  <<endl;
int auxpos =  (actual-cposinicial[nivel])*m_bs + (m_bs-1);
//cout << " auxpos " << auxpos <<"  " <<contador <<"  " << obj <<"  " <<actual  << " " <<m_B[nivel].size()-1  << " " <<  (actual-posinicial[nivel])*m_bs << endl;
int entro = 0;


if(auxpos > m_B_star[nivel].size()){

auxpos = m_B[nivel].size()-2;
//obj++;
}
//cout << " auxpos " << auxpos <<"  " <<contador <<"  " << obj <<"  " <<actual  << " " <<m_B[nivel].size()-1  <<endl;   


while(contador != obj && entro < m_bs){
int pos = cselectwt(nivel,auxpos+1);
//cout <<"POS " << pos << endl;
if(m_B_star[0][pos] == 1)contador--;
else contador++;
auxpos--;
entro++;
}

//puts("RETORNARE FINAL");
return auxpos+1;
}



 };




}// end namespace sdsl
#endif

