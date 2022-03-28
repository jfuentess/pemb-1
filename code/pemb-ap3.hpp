#ifndef INCLUDED_SDSL_PEMB
#define INCLUDED_SDSL_PEMB

#include <sdsl/sdsl_concepts.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/rank_support_v.hpp>
#include <sdsl/select_support_mcl.hpp>
#include <sdsl/wt_helper.hpp>
#include <sdsl/util.hpp>
#include <sdsl/bp_support_sada.hpp>
#include <sdsl/bp_support_gg.hpp>
#include <sdsl/bp_support_g.hpp>
#include <algorithm> // for std::swap
#include <stdexcept>
#include <vector>
#include <stack>
#include <utility>
#include "auxiliar.hpp"
#include "baseline.hpp"
#include <complementary/Vertex.hpp>
#include <complementary/Edge.hpp>
#include <complementary/Tree.hpp>
#include <vector>
#include <sdsl/sd_vector.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/enc_vector.hpp>
#include <sdsl/dac_vector.hpp>
#include <sdsl/vlc_vector.hpp>
#include <sdsl/coder.hpp>
#include <sdsl/inv_perm_support.hpp>



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
 *
 *  \par References
- *      - Leo Ferres, José Fuentes Sepúlveda, Travis Gagie, Meng He and Gonzalo
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
		class cinvector = dac_vector<>,
	   class t_succ_tree    = bp_support_sada<>,
	   class t_succ_tree2    = bp_support_g<>,
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
	typedef cinvector				cintvectort;
        typedef n_rank                               nrank_1_type;
        typedef n_select1                            nselect_1_type;
        typedef n_select0                            nselect_0_type;
        typedef t_succ_tree                          succ_tree;
	typedef t_succ_tree2                         succ_tree2; 

protected:
        size_type         m_vertices  = 0;
        size_type         m_edges  = 0;
        size_type		  m_levels = 1;
        bit_vector_type *m_A;
        bit_vector_type arbol;
        bit_vector_type pmarcado;
        nbit_vector_type *B;
        bit_vector_type *BP;
        nbit_vector_type *BF;
        nbit_vector_type *B_F;
        bit_vector_type *m_B;
        rank_1_type       pmarcado_rank;
	rank_1_type       arbol_rank;
	select_1_type	arbol_select1;
	select_1_type	pmarcado_select1;
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
        succ_tree          pmarcado_st;
        succ_tree          arbol_st;
        succ_tree          *m_B_st;
        succ_tree          *m_B_star_st;
        auxiliar			t;
	cintvectort comprimido;
//	int *plano;
	intv plano;
	intv plano2;
	intv plano3;
	intv samples,directoprueba;
inv_perm_support<> perm;
        wt_hutu_int<>      wt;
        intv  *resumen;
	int *nlimites;
	int *mapeo, *reverso;
	unsigned long long totaldistancia = 0;	
	uint32_t samples_dens =	16;
	double tiempoh =1000.0,tiempohm= 0.0,contadoreh = 0.0;
	double tiempos = 1000.0,tiemposm = 0.0,contadores =0.0, tiempobloque =0.0;
        int32_t           m_height  = 0; // height of the rmMt
        uint32_t          m_num_internal  = 0; // height of the rmMt
        uint32_t          m_num_leaves  = 0; // Leaves of the rmMt
        uint32_t          m_bs = 128;
        uint32_t          posinicial = 0;

        void copy(const pemb& p) {
            m_vertices          = p.m_vertices;
            m_edges         = p.m_edges;
            m_levels = p.m_levels;
            m_A       = p.m_A;
            B       = p.B;
            BP       = p.BP;
            BF       = p.BF;
            arbol = p.arbol;
		arbol_select1 = p.arbol_select1;
            pmarcado= p.pmarcado;
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
	    arbol_st = p.arbol_st;
	    pmarcado_st = p.pmarcado_st;
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




	vector <int> cantidad,acumulado;


	map <pair<int,int>,int >mapa;
	vector <int> limites2;
	FILE *fp = fopen(argv[1], "r");
	FILE *fp2 = fopen(argv[2], "r");
	fscanf(fp2,"%u",&niveles);
	fscanf(fp,"%u %u",&n,&m);
	for(int i = 0; i < niveles; i++){
		fscanf(fp2,"%d",&aux);
		cantidad.push_back(aux);
		acumulado.push_back(tam);
		limites2.push_back(aux);
		tam += aux;
	}
			acumulado.push_back(tam);
	baseline g3(tam,niveles);
	int enlaces = 0;	
	for(int i = 0; i < m*2; i++){
		fscanf(fp,"%u %u",&a,&b);
		g3.setVecino(a,b);
		//cout << "1 " << a << " " << b << endl;
		//if(mapa[make_pair(a,b)]==1)puts("MALOOO");
		mapa[make_pair(a,b)]++;
		enlaces++;
	}



	map< pair <int,int>,int > comprobar;

	for(int i = 0; i < n; i++){
		fscanf(fp2,"%d",&aux);
		for(int j = 1; j < niveles; j++){
			fscanf(fp2,"%u",&b);
			if(comprobar[make_pair(b+acumulado[j],aux)] == 0){
			g3.setPadre(aux,b+acumulado[j]);
		//			cout << "2 " << aux << " " << b+acumulado[j] << endl;
			g3.setHijo(b+acumulado[j],aux);
		//						cout << "3 " << b+acumulado[j] << " " << aux << endl;
			enlaces += 2;
			}
			comprobar[make_pair(b+acumulado[j],aux)]++;
			aux = b+acumulado[j];
		}
	}
	comprobar.clear();
	int pos = 0;
	vector<int>cenlaces(10,0);
	//	puts("VIVO");
	for(int i = 0; i < tam-cantidad[niveles-1]; i++){
		if(i == acumulado[pos])pos++;
		a = g3.getPadre(i);
		for(int j = 0; j < g3.getSizeV(i);j++){
			b = g3.getvecino(i,j);
			if(comprobar[make_pair(a,g3.getPadre(b))] == 0 && a!=g3.getPadre(b)){	
				cenlaces[pos]++;
				comprobar[make_pair(a,g3.getPadre(b))]++;
				enlaces ++;
				g3.setVecino(a,g3.getPadre(b));
//						cout << "1 " << a << " " << g.getPadre(b) << endl;
//				if(mapa[make_pair(a,g.getPadre(b))]==1)puts("MALOOO");
				mapa[make_pair(a,g3.getPadre(b))]++;
			}
		}
	}


	vector <bool> arbolito(tam*2,false),marcado(tam*2,false),mvisitado(tam,false);
	vector <int> directo2(tam,0);
	int auxpos = 0,otroaux= 0,numabierto = 0;
	for(int i = acumulado[niveles-1]; i < acumulado[niveles];i++){
//		cout << i  - acumulado[niveles-1]<< endl;
		stack<int> s;
		s.push(i);
		while(!s.empty()){
			int t = s.top();
			if(mvisitado[t]){
				s.pop();
				auxpos++;
			}
			else {
		otroaux++;
				mvisitado[t] = true;
				arbolito[auxpos] = true;
				int open = min(numabierto*2 +1,t*2);
				open == (t*2) ? marcado[t*2] = true : marcado[numabierto*2 +1] = true;
				directo2[t] = auxpos;
//				if(t == 42+acumulado[niveles-1])cout << "42 = " << t*2 << "  " << numabierto*2+1 << endl;
				auxpos++;
				numabierto++;
				if(t >= acumulado[1])for(int j = 0; j < g3.getSizeH(t); j++)s.push(g3.gethijo(t,j));
			}
		}
	}

double prom = 0.0;
vector<long long int> promedios;
vector <long long int> acumu(10,0);
promedios.push_back(directo2[0]);
for(int i = 1; i < tam; i++){
prom+= (double)abs(directo2[i] - directo2[i-1])/(double)tam;
promedios.push_back(directo2[i] - directo2[i-1]);
}

for(int i = 0; i < promedios.size();i++){
if(promedios[i] < 0)promedios[i] = (promedios[i]*-2) -1;
else promedios[i]*=2;
}
for(int i = 0;i < promedios.size();i++){


for(int j = 0; j < 10; j++)if(pow(10,j) > promedios[i]){
acumu[j]++;
break;
}

}
vector<int> directo(directo2.size());
for(int i = 0; i < directo2.size();i++)directo[i] = promedios[i]; 
fclose(fp2);

	FILE *fp11 = fopen(argv[2], "r");
	if (!fp11) {
    fprintf(stderr, "Error opening file \"%s\".\n", argv[2]);
    exit(EXIT_FAILURE);
	}

	n = g.vertices();
    fscanf(fp11,"%d",&niveles);
	unsigned int *limites = new unsigned int[niveles+1](); 
	unsigned int *totalnivel = new unsigned int[niveles](); 
	unsigned int *jerarquias = new unsigned int[niveles*n](); 
//	cantidad.clear();
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







	  //inicializo todo
	  m_vertices = g.vertices();
	  m_edges = g.edges();
	  m_levels = niveles;
     m_A = new bit_vector_type[niveles]();
    resumen = new intv[niveles];
     B = new nbit_vector_type[niveles]();
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
m_B_star_st = new succ_tree[niveles]();
nlimites = new int[niveles];

   mapeo = new int[6*m_vertices];
reverso = new int[6*m_vertices];


for(int i = 0; i < niveles; i++)nlimites[i] = limites[i+1]-limites[i];




	  t = g.dfs_spanning_tree_propio(inicial,jerarquias,limites,niveles,total);
	  vector <int> corch,parente;
	  vector <vector <char> >ParenteAux =t.getParentesis();
//	puts("aun vivo");
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
	  pos = 0;
	  	  	  bool *bits_aux = t.getBits();	
	  	  	  int *cierres_aux = t.getcierres();	



bit_vector_type auxA(total*2,0);
bit_vector_type auxB(total*2,0);
int compa = 0, compb = 0,compc = 0, compd = 0,compe = 0;
plano2.resize(total);
int mini = 10000000, maxim = 0;
double prome = 0.0;



long long int acul = 0;
long long sact = 0;
for(int ll = 4; ll < 5; ll++){
samples_dens = pow(2,ll);
samples.resize(directo2.size()/samples_dens);
sact  = 0;
for(int i = 0; i< directo2.size();i++){
if(i%samples_dens == 0){
samples[sact] = directo2[i];
sact++;
}
}
long long int diffact;
for(int i = 0; i < directo2.size();i++){
if(i%samples_dens == 0)diffact = samples[i/samples_dens];
directo[i] = directo2[i]-diffact;
}
for(int i = 0; i < directo.size();i++){
directo[i] >= 0? directo[i]*= 2:directo[i] = (directo[i]*-2)-1 ;
}
directoprueba.resize(directo2.size());
for(int i = 0; i < directo2.size();i++)directoprueba[i] = directo2[i];
cintvectort comp2(directo);


map<int,int> repetidito;
for(int i = 0; i < directo2.size();i++){
repetidito[directo2[i]]++;
}
comprimido.swap(comp2);
}


for(int i = 0; i < total*2; i++){
	if(arbolito[i]){auxA[i] = 1;compa++;}
}
vector<int> acumuladitos(9,0);
for(int i = 0; i < total; i++){plano2[i] = directo[i];
if(i != 0){
if(mini > abs(plano2[i]-plano2[i-1])) mini= abs(directo[i]-directo[i-1]);
if(maxim < abs(plano2[i]-plano2[i-1]))maxim = abs(directo[i]-directo[i-1]);
prome +=(double) ((double)abs(directo[i]-directo[i-1])/(double)(directo.size()-1));
int diff = abs(directo[i]-directo[i-1]);

for(int l = 0; l < 9; l++)if(pow(10,l)>diff){
acumuladitos[l]+=1;
break;
}
}
}
pmarcado.swap(auxB);
arbol.swap(auxA);
cintvectort comp2(directo);
comprimido.swap(comp2);


for(int i = 0; i < 18;i++){
long long int sumapre=0;
long long int act = comprimido[i];
if(i != 0)act%2 == 0?act = act >>1   : act = ( (act+1)>>1)*-1;
else act = act >> 1;

if( directo2[i] != (samples[i/samples_dens]+ act))cout<< i << " " <<comprimido[i]  << " " <<   directo2[i]<< "  " << samples[i/samples_dens]  <<"  " << sumapre  << " "  << directo[i]   <<endl;
}
int ant = -1;
for(int i = 0; i < total; i++){

if(ant == -1){
ant = plano2[i];
}
else if (ant > plano2[i]){
ant = plano2[i];
compc++;
}

else {
ant = -1;
compd++;
}

}

enc_vector<coder::elias_delta,4> ev(plano2);
int_vector<> v(plano2.size());
for(int kk = 0; kk < plano2.size();kk++)v[kk] = plano2[kk];
  util::bit_compress(v);
  util::bit_compress(plano2);

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
int_vector <> secuencia(vaux.size());
plano.resize(vaux.size()); //= new vint	[vaux.size()];
for(int i = 0; i < vaux.size();i++){
secuencia[i] = vaux[i];
plano[i] = vaux[i];
}
construct_im(wt,secuencia);
util::bit_compress(plano); 
          util::init_support(pmarcado_rank, &pmarcado);
        util::init_support(pmarcado_select1, &pmarcado);   
 util::init_support(arbol_select1, &arbol);
 util::init_support(arbol_rank, &arbol);
for(int i = 0; i < niveles; i++){
	  util::init_support(m_A_rank[i], &m_A[i]);
	  util::init_support(m_B_rank[i], &m_B[i]);
 	  util::init_support(B_rank[i], &B[i]);
 	  util::init_support(BP_rank[i], &BP[i]);
 	  util::init_support(BF_rank[i], &BF[i]);
 	  util::init_support(B_F_rank[i], &B_F[i]);
 	  util::init_support(m_A_select1[i], &m_A[i]);
  	  util::init_support(m_B_select1[i], &m_B[i]);
 	  util::init_support(B_select1[i], &B[i]); 
 	  util::init_support(BP_select1[i], &BP[i]); 
 	  util::init_support(B_F_select1[i], &B_F[i]); 
 	  util::init_support(BF_select1[i], &BF[i]); 
 	  util::init_support(BF_select0[i], &BF[i]); 
	  util::init_support(m_A_select0[i], &m_A[i]);
	  util::init_support(m_B_select0[i], &m_B[i]);
	  succ_tree B_local_st(&m_B[i]);
	  succ_tree B_star_local_st(&m_B_star[i]);

	  succ_tree local_arbol_st(&arbol);
	  arbol_st.swap(local_arbol_st);
	  arbol_st.set_vector(&arbol);
	  m_B_st[i].swap(B_local_st);
	  m_B_st[i].set_vector(&m_B[i]);
	  m_B_star_st[i].swap(B_star_local_st);
	  m_B_star_st[i].set_vector(&m_B_star[i]);
}


intv auxplano(directo2.size());

for(int i = 0; i < directo2.size();i++)auxplano[i] = arbol_st.rank(directo2[i]);
plano3.swap(auxplano);
inv_perm_support<> perm2(&plano3);
perm = perm2;
int suma= 0, anteri = 0;
for(int i = 0; i < 6; i++){
anteri = nlimites[i];
nlimites[i] = suma;
suma+= anteri;
}
util::bit_compress(directoprueba);
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
	  written_bytes += m_A_rank[i].serialize(out, child, "A_rank");
	  written_bytes += m_A_select1[i].serialize(out, child, "A_select1");
	  written_bytes += m_A_select0[i].serialize(out, child, "A_select0");
	  if(i > 1){
	  written_bytes += B[i].serialize(out, child, "B");
	  written_bytes += B_rank[i].serialize(out, child, "B_rank");
	  written_bytes += B_select1[i].serialize(out, child, "B_select1");
	  }
	  else if(i == 1){
	  written_bytes += BP[i].serialize(out, child, "B");
	  written_bytes += BP_rank[i].serialize(out, child, "B_rank");
	  written_bytes += BP_select1[i].serialize(out, child, "B_select1");	  	
	  }
	  written_bytes += m_B[i].serialize(out, child, "mB");
	  written_bytes += m_B_rank[i].serialize(out, child, "mB_rank");
	  written_bytes += m_B_select1[i].serialize(out, child, "mB_select1");
	  written_bytes += m_B_st[i].serialize(out, child, "mB_succ_tree");
	  written_bytes += m_B_select0[i].serialize(out, child, "B_select0");
	  
	  written_bytes += m_B_star[i].serialize(out, child, "mB_star");
	  written_bytes += m_B_star_st[i].serialize(out, child, "mB_star_succ_tree");
	  }
          written_bytes += arbol.serialize(out, child, "Arbol");    
          written_bytes += plano2.serialize(out, child, "PLANO");

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
       size_type first(size_type v, int nivel) {
	 if(v >= 0) {
	   size_type pos = m_B_select1[nivel](v+1);
	   size_type edge = 1;
	    edge = m_A_select1[nivel](pos+1);
	     return edge;
	 } else
	   return -1;
       }

      /* Assuming indices start with 0 */
       //dado una posicion con parentesis que abre me da la posicion donde cierra, es circular

       void iniresumen(int32_t bs){
	int valor;
	for(int l = 0; l < m_levels; l++)resumen[l].width(32);
       	m_bs = bs;
       	m_num_leaves = ceil((double)B[0].size()/m_bs);
        m_height = ceil(log(m_num_leaves)/log(2));
        m_num_internal =(pow(2,m_height+1)-1);
        for(int l = 0; l < m_levels; l++)resumen[l].resize(m_num_internal);
	for(int l = 0; l < m_levels; l++)util::set_to_value(resumen[l], 0);
	int comprobar = 0;
	for(int l = 0; l < m_levels; l++){
	valor = l;
        for(int i = 0; i < m_num_leaves; i++){
        	int pos = i*m_bs;
        	int cont = 0;
        	for(int j = 0; j <  m_bs && pos +j < B[0].size();j++){
        		if(wt[pos+j] >= valor)cont++;
        	}
        	int top_nodes = pow(2, m_height)-1;
		if(i == valor)posinicial = top_nodes;
        	resumen[l][top_nodes +i] = cont;
		comprobar += cont;
        }
        for(int lvl=m_height-1; lvl >= 0 ; lvl--){
	    int num_curr_nodes = pow(2, lvl);
	    int node = 0;
	    	for( ;node < num_curr_nodes  ;node++){
		      int pos = (pow(2,lvl)-1) + node;
		      resumen[l][pos] = (resumen[l][2*pos+1] + resumen[l][2*pos+2]);
	    	}
       	}
}


	for(int i = 0; i < m_levels; i++){
	}
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

      vector <int> lista_vecinos(size_type v, int nivel){
      	size_type nxt = first(v,nivel)+1, limite = mate(first(v,nivel),nivel);
      	vector <int> retorno;
      	while(nxt < limite){
      		size_type x;
      		size_type pos_in_A = m_A_rank[nivel](nxt);
      		if(m_A[nivel][nxt] == 1)pos_in_A++;
      		if(m_A[nivel][nxt] == 1 && m_B[nivel][pos_in_A-1] == 1){
      			x = vertex2(nxt,nivel);
      			nxt = mate(nxt,nivel)+1;
      		}
      		else{
      			x = vertex2(nxt,nivel);
      			nxt ++;
      		}
      		retorno.push_back(x);
      	}

      	if(v != 0){
      		nxt = m_B_st[nivel].parent_t(m_B_select1[nivel](v+1));
      		int auxtemp = m_A_rank[nivel](m_A_select1[nivel](nxt+1));
      		size_type x = m_B_rank[nivel](auxtemp);
      		if(m_B[nivel][auxtemp] == 1)x++;
      		if(x-1 < limite)retorno.push_back(x-1);
      	}
      	return retorno;
      }


      //Funciones auxiliar para hacer coincidir las regiones consultadas en el baseline
	int getMapeo(size_type pos, unsigned int nivel){
    int retorno = m_B_select1[nivel](pos+1); // posicion en M_B donde esta el x-ecimo parentesis abierto
    retorno = B_select1[nivel](retorno+1); //posicion en B del parentesis abierto que representa X
	retorno = m_B_rank[0](retorno);
		return t.getposMapeo(m_vertices*nivel+retorno);
	}
	int getReverso(size_type pos, unsigned int nivel){
		return t.getposReverso(m_vertices*nivel+pos);
	}
//Retorna ID region correspondiente

    unsigned int regionID(unsigned int x, unsigned int level){
	return m_B_rank[level](x);
	}
	//Retorna primer parentesis en nivel 0 contenido por la region entregada

        bool Inside(unsigned int x, unsigned int levelx, unsigned int y, unsigned int levely){

           	int posy,posx,posf,posx2;
posx = directoprueba[x];
posy = directoprueba[y];
                posf = /*posx-1000;*/arbol_st.find_close(posy);
                if(posf >=posx && posy <=posx ) return true;
                return false;  
              }   

	bool Touches(unsigned int x, unsigned int levelx, unsigned int y, unsigned int levely,unsigned int xx, unsigned int yy){
		if(x == y && levelx == levely)return true;
		if(levely < levelx){
	bool comprobar = Inside2(yy, levely, xx , levelx);
			vector <int> vecinos = (lista_vecinos(y,levely));
			for(int i = 0; i < vecinos.size();i++){
				bool contenido = Inside2(getMapeo(vecinos[i],levely)+nlimites[levely],levely,getMapeo(x,levelx)+nlimites[levelx],levelx);

				if(contenido  && !comprobar)return true;
				if(!contenido  && comprobar)return true;
	
			}
			return false;
		}
		else if(levely == levelx){
			vector <int> vecinos = (lista_vecinos(y,levely));
			for(int i = 0; i < vecinos.size();i++){
				if(vecinos[i] == x )return true;
			}


		}
		else {
			bool comprobar = Inside2(xx,levelx,yy,levely);

			vector <int> vecinos = (lista_vecinos(x,levelx));
			for(int i = 0; i < vecinos.size();i++){
bool contenido = Inside2(vecinos[i]+nlimites[levelx],levelx,y+nlimites[levely],levely);
                              if(contenido)return true;
				if(contenido && !comprobar)return true;
				if(!contenido  && comprobar)return true;

			}
		}
		return false;
	} 

	vector <int> compvecinos(unsigned int x, int levelx){
		vector <int> vecinos = (lista_vecinos(x,levelx));
		return vecinos;

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
          int posx = /*arbol_st.select(plano3[x])*/ arbol_select1(plano3[x]),diflevel = levelx - levely;
          int posini =  /*arbol_st.select(plano3[x]+diflevel)*/arbol_select1(plano3[x]+diflevel);
          int posfinal = arbol_st.find_close(posx),posaux = posini;
          vector <int> retorno;
	int limite2 = arbol_st.rank(posfinal);
                while(posaux < posfinal){
		int agregar= /*arbol_st.rank(posaux)*/arbol_rank(posaux);
                  retorno.push_back(perm[agregar]);
                  int posaux2 = arbol_st.find_close(posaux)+1;
		int arbolrank =/* arbol_st.rank(posaux2)*/arbol_rank(posaux2);
		if(arbolrank == limite2)break;
                  if(arbol[posaux2] == 0){
                      int totalabierto = arbolrank;//arbol_st.rank(posaux2);
                      int nextopen= /*arbol_st.select(totalabierto+1)*/arbol_select1(totalabierto+1);
                      int totalcerrado = nextopen - posaux2+1;
                      posaux = /*arbol_st.select(totalabierto + totalcerrado)*/arbol_select1(totalabierto + totalcerrado);
                  }
		else posaux= posaux2;
                }
                return retorno;
        }

  };




}// end namespace sdsl
#endif



8 porciento







if(pe->Touches(pe->getReverso(regionA[i],nivelA[i]),nivelA[i],pe->getReverso(regionB[i],nivelB[i]),nivelB[i],regionA[i]+limites[nivelA[i]],regionB[i]+limites[nivelB[i]]))
