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
#include "auxiliar.hpp"
#include <complementary/Vertex.hpp>
#include <complementary/Edge.hpp>
#include <complementary/Tree.hpp>
//#include "Graph.hpp"
#include "graph3.hpp"
#include <vector>
#include <sdsl/sd_vector.hpp>
#include <sdsl/rrr_vector.hpp>
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

  template<class t_bitvector   = bit_vector,
  		class n_bitvector = bit_vector,
//class nn_bitvector = rrr_vector <>,
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
        typedef t_rank                               rank_1_type;
        typedef t_select1                            select_1_type;
        typedef t_select0                            select_0_type;
        typedef n_bitvector                          nbit_vector_type;
//	typedef nn_bitvector                         nnbit_vector_type;
        typedef n_rank                               nrank_1_type;
        typedef n_select1                            nselect_1_type;
        typedef n_select0                            nselect_0_type;
        typedef t_succ_tree                          succ_tree;

protected:
	vector<int> histo;
	long long totalentro = 0, totalmaximo = 0;
        size_type         m_vertices  = 0;
        size_type         m_edges  = 0;
        size_type		  m_levels = 1;
        bit_vector_type *m_A;
        nbit_vector_type *B;
        bit_vector_type *BP;
        nbit_vector_type *BF;
        nbit_vector_type *B_F;
        bit_vector_type *m_B;
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
	int *mapeo,*reverso;
	unsigned int sep;

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
//sep default = 1;
     	pemb(Graph g,unsigned int inicial, unsigned int separador, char **argv) {
	int total, tam;
	unsigned int n,m,niveles,aux,mitotal = 0,tam = 0,a,b,c;
FILE *fp11 = fopen(argv[2], "r");
	if (!fp11) {
    fprintf(stderr, "Error opening file \"%s\".\n", argv[2]);
    exit(EXIT_FAILURE);
	}

	n = g2.vertices();
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
	total = mitotal;
	int tam = n;





	  //inicializo todo
		histo.resize(9,0);
     		sep = separador;
	  m_vertices = g.vertices();
	  m_edges = g.edges();
	  m_levels = niveles;
     m_A = new bit_vector_type[niveles]();
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
   mapeo = new int[6*m_vertices];
reverso = new int[6*m_vertices];


	  t = g.dfs_spanning_tree_propio(inicial,jerarquias,limites,niveles,total);
	  vector <int> corch,parente;
	  vector <vector <char> >ParenteAux =t.getParentesis();
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
	  	for(int j = 0; j < tam*2; j++){
	  		if(bits_aux[j + i*tam*2] == true){
	  			B_l[j] = 1;
	  		}
	  	}
	  	BP[i].swap(B_l);
	  }
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





for(int i = 0; i < niveles; i++){
	  util::init_support(m_A_rank[i], &m_A[i]);
	  util::init_support(m_B_rank[i], &m_B[i]);
 	  util::init_support(B_rank[i], &B[i]);
 	  util::init_support(BP_rank[i], &BP[i]);
 	  util::init_support(m_A_select1[i], &m_A[i]);
  	  util::init_support(m_B_select1[i], &m_B[i]);
 	  util::init_support(B_select1[i], &B[i]); 
 	  util::init_support(BP_select1[i], &BP[i]); 
	  util::init_support(m_A_select0[i], &m_A[i]);
	  util::init_support(m_B_select0[i], &m_B[i]);
	  succ_tree B_local_st(&m_B[i]);
	  succ_tree B_star_local_st(&m_B_star[i]);
	  m_B_st[i].swap(B_local_st);
	  m_B_st[i].set_vector(&m_B[i]);
	  m_B_star_st[i].swap(B_star_local_st);
	  m_B_star_st[i].set_vector(&m_B_star[i]);
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

	  for(int i = 0; i < m_levels;i++){
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
	  written_bytes += BP[i].serialize(out, child, "B2");
	  written_bytes += BP_rank[i].serialize(out, child, "B2_rank");
	  written_bytes += BP_select1[i].serialize(out, child, "B2_select1");	  	
	  }
	  written_bytes += m_B[i].serialize(out, child, "mB");
	  written_bytes += m_B_rank[i].serialize(out, child, "mB_rank");
	  written_bytes += m_B_select1[i].serialize(out, child, "mB_select1");
	  written_bytes += m_B_st[i].serialize(out, child, "mB_succ_tree");
	  written_bytes += m_B_select0[i].serialize(out, child, "B_select0");
	  
	  written_bytes += m_B_star[i].serialize(out, child, "mB_star");
	  written_bytes += m_B_star_st[i].serialize(out, child, "mB_star_succ_tree");
	  }
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
    unsigned int goUp(unsigned int x, unsigned int level){
    //cout << x << " " << level << endl;
    int posicion = m_B_select1[level](x+1); // posicion en M_B donde esta el x-ecimo parentesis abierto
   // puts("vivo");
   if(level > 1) posicion = B_select1[level](posicion+1);
   	else posicion = BP_select1[level](posicion+1); //posicion en B del parentesis abierto que representa X
    //puts("retorno");
	return posicion;
	}




	//Retorna region que contiene a la region entregada, region entregada siempre pertenece al nivel 0
    unsigned int goDown(unsigned int x, unsigned int level){
    	int posB = m_B_select1[0](x+1);

      	if(level > sep){if(B[level][posB] == 1){
    		int retornar = B_rank[level](posB); //Obtengo posicion en m_B
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
    		int parentesiscercano =B_select1[level](B_rank[level](posB)+1); //Posicion siguiente parentesis en B
    		int retornar = B_rank[level](parentesiscercano); // Posicion del parentesis siguiente en m_b
    		bool abierto = (m_B[level][retornar] == 1) ? true: false;  //Compruebo si es abierto o cerrado

    		if(abierto && B_rank[level](posB) != B_rank[level](B[level].size()-1)){
    			retornar = m_B_st[level].parent_t(retornar); //en caso de ser abierto obtengo el padre en m_B 	
		return retornar; //retorno directamente
    		}
    		else {
    		retornar = m_B_st[level].find_open(retornar); // Caso de ser cerrado obtengo parentesis que abre  	
	return retornar;		
    		}
    	}
    }
		else {
		if(BP[level][posB] == 1){
      	//	puts("entre");
    		int retornar = BP_rank[level](posB); //Obtengo posicion en m_B
		  	bool pab = (m_B[level][retornar] == 1) ? true: false; //Compruebo si parentesis es abierto o cerrado 
		   if(pab){
		//   	      		puts("entre 1");
		   	 return retornar;
		   }  
		   else {
//	retornar = 100;    
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
//retornar = 100;  
  			return retornar; //retorno directamente
    		}
    		else {
    		retornar = m_B_st[level].find_open(retornar); // Caso de ser cerrado obtengo parentesis que abre
//	retornar = 100;    	
	return retornar;			
    		}
    	}	
		}

	}


    unsigned int goLevel(unsigned int x, unsigned int levels, unsigned int levelt){
    unsigned int regions = regionID(goUp(x,levels),0);
    unsigned int regiont = regionID(goDown(regions,levelt),levelt);
	return regiont;
	}

	bool Inside(unsigned int x, unsigned int levelx, unsigned int y, unsigned int levely){
		unsigned int comprobar = goLevel(x,levelx,levely);
		if(comprobar == y)return true;
		else return false;
	}

	bool Touches(unsigned int x, unsigned int levelx, unsigned int y, unsigned int levely){
		if(x == y && levelx == levely)return true;
		if(levely < levelx){
			bool comprobar = Inside(y,levely,x,levelx);
			vector <int> vecinos = (lista_vecinos(y,levely));
			for(int i = 0; i < vecinos.size();i++){
				unsigned int contenedor = goLevel(vecinos[i],levely,levelx);
				if(contenedor == x  && !comprobar)return true;
				if(contenedor != x  && comprobar)return true;
	
			}
		}
		else if(levely == levelx){
			vector <int> vecinos = (lista_vecinos(y,levely));
			for(int i = 0; i < vecinos.size();i++){
				if(vecinos[i] == x )return true;
			}


		}
		else {
			bool comprobar = Inside(x,levelx,y,levely);
			vector <int> vecinos = (lista_vecinos(x,levelx));
			for(int i = 0; i < vecinos.size();i++){
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
		int contador = 0,totalaux=0;
		while(nxt < limit){
			size_type pos_in_A = m_B_rank[levely](m_A_rank[levely](nxt));
			size_type pos = m_B_select1[levely](pos_in_A+1);
			if(levely > 1)pos = B_select1[levely](pos+1);
			else pos = BP_select1[levely](pos+1);
			if(levelx > sep){
				if(B[levelx][pos] ==1){
			totalaux++;
				nxt = mate(nxt,levely);
				size_type pos_in_B = m_B_rank[levely](m_A_rank[levely](nxt))+1;
				if(pos_in_B > limit2)break;
		nxt = m_A_select1[levely](m_B_select1[levely](pos_in_B)+1);
			}
			else if(B[levelx][pos] == 0){
			//	totalaux++;
//				totalentro++;
				retorno.push_back(vertex2(nxt,levely));
				size_type pos_in_B = m_B_rank[levely](m_A_rank[levely](nxt))+2;
				if(pos_in_B > limit2)break;
				nxt = m_A_select1[levely](m_B_select1[levely](pos_in_B)+1);
			} 
		}
		else{
				if(BP[levelx][pos] ==1){
		totalaux++;
				nxt = mate(nxt,levely);
				size_type pos_in_B = m_B_rank[levely](m_A_rank[levely](nxt))+1;
				if(pos_in_B > limit2)break;
		nxt = m_A_select1[levely](m_B_select1[levely](pos_in_B)+1);
			}
			else if(BP[levelx][pos] == 0){
			//	totalaux++;
//				totalentro++;
				retorno.push_back(vertex2(nxt,levely));
				size_type pos_in_B = m_B_rank[levely](m_A_rank[levely](nxt))+2;
				if(pos_in_B > limit2)break;
				nxt = m_A_select1[levely](m_B_select1[levely](pos_in_B)+1);
			} 			
		} 

			contador++;

		}
		for(int i = 0 ; i< histo.size();i++){
		 if(totalaux <= pow(10,i)){
//			cout << totalaux <<endl;
			histo[i]++;
			break;
			}
		}
//		cout << totalaux << endl;
		if(totalaux > totalmaximo)totalmaximo = totalaux;
		return retorno;

	}

	void cuantoentro(){
//	puts("acumulado");
	for(int i = 0; i < histo.size();i++)cout << histo[i] << endl;
	cout << "ENTRE: " << totalentro << " " <<  totalmaximo <<endl;
	}
  };








}// end namespace sdsl
#endif

