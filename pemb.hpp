#ifndef INCLUDED_SDSL_PEMB
#define INCLUDED_SDSL_PEMB

#include "sdsl_concepts.hpp"
#include "int_vector.hpp"
#include "rank_support_v.hpp"
#include "select_support_mcl.hpp"
#include "wt_helper.hpp"
#include "util.hpp"
#include "bp_support_sada.hpp"
#include <algorithm> // for std::swap
#include <stdexcept>
#include <vector>
#include <stack>
#include <utility>
#include "auxiliar.hpp"
#include "../complementary/Vertex.hpp"
#include "../complementary/Edge.hpp"
#include "../complementary/Tree.hpp"
#include "../complementary/Graph.hpp"
#include <vector>
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
	   class t_succ_tree    = bp_support_sada<>,
	   class t_rank        = typename t_bitvector::rank_1_type,
	   class t_select1     = typename t_bitvector::select_1_type,
	   class t_select0     = typename t_bitvector::select_0_type>
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
        typedef t_succ_tree                          succ_tree;

protected:
        size_type         m_vertices  = 0;
        size_type         m_edges  = 0;
        size_type		  m_levels = 1;
        bit_vector_type *m_A;
        bit_vector_type *B;
        bit_vector_type *m_B;
        rank_1_type       *m_A_rank;
        rank_1_type       *m_B_rank;
        rank_1_type       *B_rank;
        select_1_type     *m_A_select1;
        select_1_type     *m_B_select1;
        select_1_type     *B_select1;
        select_0_type     *m_A_select0;
        select_0_type     *m_B_select0;
        bit_vector_type    *m_B_star;
        succ_tree          *m_B_st;
        succ_tree          *m_B_star_st;
        auxiliar			t;

        void copy(const pemb& p) {
            m_vertices          = p.m_vertices;
            m_edges         = p.m_edges;
            m_levels = p.m_levels;
            m_A       = p.m_A;
            B       = p.B;
            B_rank  = p.B_rank;
            B_select1     = p.B_select1;
            m_B_rank  = p.m_B_rank;
            m_B_select1     = p.m_B_select1;
            m_A_rank  = p.m_A_rank;
            m_A_select1     = p.m_A_select1;
            m_A_select0     = p.m_A_select0;
	    m_A_rank.set_vector(&m_A);
	    m_A_select1.set_vector(&m_A);
	    m_A_select0.set_vector(&m_A);
		m_B_select0.set_vector(&m_B);
	    m_B = p.m_B;
	    m_B_star = p.m_B_star;
	    m_B_st = p.m_B_st;
	    m_B_star_st = p.m_B_star_st;
	    m_B_st.set_vector(m_B);
	    m_B_star_st.set_vector(m_B_star);
        }


    public:

        //! Default constructor
        pemb() {};

     	pemb(Graph g,unsigned int inicial, unsigned int *jerarquias, unsigned int *limites, unsigned int niveles,unsigned int total, unsigned int tam) {
	  m_vertices = g.vertices();
	  m_edges = g.edges();
	  m_levels = niveles;
     m_A = new bit_vector_type[niveles]();
     B = new bit_vector_type[niveles]();
m_B = new bit_vector_type[niveles]();
m_B_star = new bit_vector_type[niveles]();
m_A_rank = new rank_1_type[niveles]();
m_A_select1 = new select_1_type[niveles]();
B_rank = new rank_1_type[niveles]();
B_select1 = new select_1_type[niveles]();
m_B_rank = new rank_1_type[niveles]();
m_B_select1 = new select_1_type[niveles]();
m_A_rank = new rank_1_type[niveles]();
m_A_select1 = new select_1_type[niveles]();
m_A_select0 = new select_0_type[niveles]();
m_B_select0 = new select_0_type[niveles]();
m_B_st = new succ_tree[niveles]();
m_B_star_st = new succ_tree[niveles]();
   //  bit_vector_type *m_A_aux = new bit_vector_type[niveles]();
	  t = g.dfs_spanning_tree_propio(inicial,jerarquias,limites,niveles,total);
	 //  puts("inicialice t");
	  vector <int> corch,parente;
	  vector <vector <char> >ParenteAux =t.getParentesis();
	  for(int i = 0; i < niveles; i++){
	  	int aux1 = 0, aux2= 0;
	  //	printf("ENTRE AL NIVEL %d\n",i);
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
	 // puts("llegue al parentesis");
	  	  	  bool *bits_aux = t.getBits();	

	  for(int i = 0; i < niveles; i++){
	  //puts("inicio for");
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
	  				  				  				//  				printf("%d ",j);	  		
	  			B_star_local[aux2] = 1;
	  			aux2++;	  			
	  		}
	  		if(ParenteAux[i][j] == ']' ){
	  				  				  				 // 			printf("%d ",j);
	  			aux2++;	  			
	  		}
	  	}
	  //	puts("");
	  //	puts("termine for");
	  	m_A[i].swap(A_local);
	    m_B[i].swap(B_local);
	    m_B_star[i].swap(B_star_local);
	  }
	  for(int i = 0; i <  niveles; i++){
	  	bit_vector_type B_l(tam,0);
	  	for(int j = 0; j < tam; j++){
	  		if(bits_aux[j + i*tam] == true){
	  			B_l[j] = 1;
	  		}
	  	}
	  	B[i].swap(B_l);
	  }
//	  m_A->swap(m_A_aux);
/*
	  	  for(int i = 0; i < 1; i++){
	  	for(int j = 0 ; j <40;j++)cout <<B[i][j] << "";
	  	cout << endl; 
	  }
	  */
	 /* 
	  puts("imprimo m_A");
	  for(int i = 0; i < niveles; i++){
	  	for(int j = 0 ; j < m_A[i].size();j++)cout <<m_A[i][j] << "";
	  	cout << endl; 
	  }
	  	  puts("imprimo B");
	  for(int i = 0; i < niveles; i++){
	  	for(int j = 0 ; j < B[i].size();j++)cout <<B[i][j] << "";
	  	cout << endl; 
	  }
	  	  	  puts("imprimo m_B");
	  	  for(int i = 0; i < niveles; i++){
	  	for(int j = 0 ; j < m_B[i].size();j++)cout <<m_B[i][j] << "";
	  	cout << endl; 
	  }
	  	  	  	  puts("imprimo m_B_star");
	  	  for(int i = 0; i < niveles; i++){
	  	for(int j = 0 ; j < m_B_star[i].size();j++)cout <<m_B_star[i][j] << "";
	  	cout << endl; 
	  }
	  */
for(int i = 0; i < niveles; i++){
	  util::init_support(m_A_rank[i], &m_A[i]);
	  util::init_support(m_B_rank[i], &m_B[i]);
 	  util::init_support(B_rank[i], &B[i]);
 	  util::init_support(m_A_select1[i], &m_A[i]);
  	  util::init_support(m_B_select1[i], &m_B[i]);
 	  util::init_support(B_select1[i], &B[i]); 
	  util::init_support(m_A_select0[i], &m_A[i]);
	  util::init_support(m_B_select0[i], &m_B[i]);
	  succ_tree B_local_st(&m_B[i]);
	  succ_tree B_star_local_st(&m_B_star[i]);
	  m_B_st[i].swap(B_local_st);
	  m_B_st[i].set_vector(&m_B[i]);
	  m_B_star_st[i].swap(B_star_local_st);
	  m_B_star_st[i].set_vector(&m_B_star[i]);
}

//	  succ_tree B_local_st(&m_B);
//	  succ_tree B_star_local_st(&m_B_star);
//[0]
//	  m_B_st.swap(B_local_st);
//	  m_B_st.set_vector(&m_B);
//	  m_B_star_st.swap(B_star_local_st);
//	  m_B_star_st.set_vector(&m_B_star);
	//  puts("TERMINEEEE");
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
		B_rank       = std::move(g.B_rank);
		B_select1       = std::move(g.B_select1);
		m_A_rank        = std::move(g.m_A_rank);
		m_A_select1     = std::move(g.m_A_select1);
		m_B_rank        = std::move(g.m_B_rank);
		m_B_select1     = std::move(g.m_B_select1);
		m_A_select0     = std::move(g.m_A_select0);
		m_B_select0     = std::move(g.m_B_select0);
		m_A_rank.set_vector(&m_A);
		m_A_select1.set_vector(&m_A);
		m_B_rank.set_vector(&m_B);
		m_B_select1.set_vector(&m_B);
		m_A_select0.set_vector(&m_A);
		m_B_select0.set_vector(&m_B);

		m_B             = std::move(g.m_B);
		m_B_star        = std::move(g.m_B_star);

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

	  written_bytes += B[i].serialize(out, child, "B");
	  written_bytes += B_rank[i].serialize(out, child, "B_rank");
	  written_bytes += B_select1[i].serialize(out, child, "B_select1");

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
	    B_rank.load(in, &m_A);
	    B_select1.load(in, &m_A);


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
	   //printf("%lu ",pos);
	   size_type edge = 1;
	    edge = m_A_select1[nivel](pos+1);
	     return edge;
	 } else
	   return -1;
       }

      /* Assuming indices start with 0 */
       //dado una posicion con parentesis que abre me da la posicion donde cierra, es circular
      size_type mate(size_type i,int nivel) {
	if(m_A[nivel][i] == 1) {
	  size_type pos_in_B = m_A_rank[nivel](i); // rank1

	  // Simulating the match operation
	  size_type match_in_B;
	 // cout << "pos in b = "<< pos_in_B << endl;
	  if(m_B[nivel][pos_in_B] == 1)
	    match_in_B = m_B_st[nivel].find_close(pos_in_B);
	  else {
	  	//	puts("buscare b");
	  		    match_in_B = m_B_st[nivel].find_open(pos_in_B);
	  	//	    	 cout << "b = "<<match_in_B << endl;
	  }
	  return m_A_select1[nivel](match_in_B+1);
	}
	else
	  {

	    size_type pos_in_B_star = i - m_A_rank[nivel](i); // rank0
	  //  cout  << " i = " << i << " pos_in_B_star " << pos_in_B_star << endl;
	    // Simulating the match operation
	    size_type match_in_B_star;
	    if(m_B_star[nivel][pos_in_B_star] == 1){
	    match_in_B_star = m_B_star_st[nivel].find_close(pos_in_B_star);
	    return m_A_select0[nivel](match_in_B_star+1);
	    }
	    else{
	    match_in_B_star = m_B_star_st[nivel].find_open(pos_in_B_star);
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
	return -1; //ACA SE MUERE SI ENCADENO LAS FUNCIONES, EN TEORIA NO LAS DEBERIA ENCADENAR(?)
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
size_type vertex2(size_type e, int nivel){
      	size_type pos_in_A = m_A_rank[nivel](e);
      	//cout <<" pos in a " << pos_in_A << endl;

      	if(m_B[nivel][pos_in_A] == 1 && m_A[nivel][e] == 1)return m_B_st[nivel].rank(pos_in_A)-1;

      	else if(m_B[nivel][pos_in_A] == 0 && m_A[nivel][e] == 1){
      		size_type match_pos = m_B_st[nivel].find_open(pos_in_A);
      		return m_B_st[nivel].rank(match_pos)-1;
  //    		return -1;
      	}

      	size_type pos_in_B = e- pos_in_A;
      	

      	//cout << "pos in b "<< pos_in_B<<endl; 
      	 if(m_A[nivel][e] == 0 && m_B_star[nivel][pos_in_B] == 1) {
      	 //	return -1;
      	 //	cout << "ENTRE AL IF DE VERTEX"<<endl;
     		size_type match_pos = mate(e,nivel);
     	//	cout << "match_pos "<< match_pos<<endl;
     		size_type pos = m_A_rank[nivel](match_pos);

     		size_type pos2 = m_A_select1[nivel](pos+1); 
     	//	cout << "pos2 "<< pos2<<endl;
     		if(m_B[nivel][m_A_rank[nivel](pos2)] == 1){
     		//	puts("entre al if");
     			pos2 = m_B_st[nivel].parent_t(m_A_rank[nivel](pos2));
     	//		cout << "PADRE "<<  pos2 << endl;
     			//cout << "POS IN A" << m_A_select1[0](pos2+1) ;

     			return m_B_rank[nivel](m_A_rank[nivel](m_A_select1[nivel](pos2+1)));
     		}
     		else {
     			pos2 = mate(pos2,nivel);

     			return m_B_rank[nivel](m_A_rank[nivel](pos2));
     		}
/*
     		size_type pos2 = m_B_select0[0](  (pos-m_B_rank[0](pos))+1); 
     		cout << "pos "<< pos<<" rank " << m_B_rank[0](pos)  <<" pos2 " <<pos2 <<endl;
     		pos_in_A = m_A_select1[0](pos2+1);
     		cout << "POS FINAL "<< pos_in_A << " MATE "<< mate(pos_in_A) << endl;
     		pos_in_A = m_A_rank[0](mate(pos_in_A));
*/
  //    		return m_B_st[0].rank(pos_in_A)-1;
//return m_B_st[0].rank(pos)-1;
      		
      	}
  	//cout << "POSICION = "<< pos_in_B << " VALOR " << m_B_star[0][pos_in_B]<<endl;
	else if(m_A[nivel][e] == 0 && m_B_star[nivel][pos_in_B] == 0){
	//	else {
     		size_type match_pos = mate(e,nivel);
     		size_type pos = m_A_rank[nivel](match_pos);

     	//	cout << "match_pos "<< match_pos<<endl;
     		size_type pos2 = m_A_select1[nivel](pos+1); 
     	//	cout << "pos2 "<< pos2<<endl;
     		if(m_B[nivel][m_A_rank[nivel](pos2)] == 1){
     			pos2 = m_B_st[nivel].parent_t(m_A_rank[nivel](pos2));
     			 	//		cout << "el retorno " << m_B_rank[nivel](m_A_rank[nivel](m_A_select1[nivel](pos2+1))) << endl;
     			return m_B_rank[nivel](m_A_rank[nivel](m_A_select1[nivel](pos2+1)));
     		}
     		else {
     			pos2 = mate(pos2,nivel);
     			return m_B_rank[nivel](m_A_rank[nivel](pos2));
     		}
/*
     		size_type pos2 = m_B_select0[0](  (pos-m_B_st[0].rank(pos))+1); 
     		cout << "pos "<< pos<<" rank " << m_B_st[0].rank(pos)  <<" pos2 " <<pos2 <<endl;
     		pos_in_A = m_A_select1[0](pos2+1);
     		//cout << "POS FINAL "<< pos_in_A << " MATE "<< mate(pos_in_A) << endl;
     		pos_in_A = m_A_rank[0](mate(pos_in_A));
      		return m_B_st[0].rank(pos_in_A)-1;


      	//	return m_B_st[0].rank(pos_in_A)-1;
     		cout << "POS " << pos << "  " << endl;
     	//	return m_B_st[0].rank(pos)-1;

     	*/
      	}
      	return -1;
      }
/*
      size_type vertex(size_type e) {
      	//	cout << "first = "<< first(e)<<endl;
      	size_type pos_in_A = m_A_rank[0](e+1); // rank1
      	//cout << " pos_in_A " << pos_in_A<<  " e " << e <<endl;
      	if(m_A[0][e] == 1) {
      	  if(m_B[0][pos_in_A] == 0) {
      	    size_type match_pos;
      	    if(m_B[0][pos_in_A] == 1)
      	      match_pos = m_B_st[0].find_close(pos_in_A);
      	    else 
      	      match_pos = m_B_st[0].find_open(pos_in_A);
	    //	cout << "match_pos = " << match_pos <<endl;
      	    return m_B_st[0].rank(match_pos);
      	  }
      	  else {
	    size_type par = m_B_st[0].parent_t(pos_in_A);
      	    return m_B_st[0].rank(par);
      	  }
      	}
      	else {
      	//	puts("entre al else");
      	  if(m_B[0][pos_in_A] == 1) {
      	    return m_B_st[0].rank(pos_in_A);
      	  }
      	  else {
      	    size_type match_pos;
      	    if(m_B[0][pos_in_A] == 1)
      	      match_pos = m_B_st[0].find_close(pos_in_A);
      	    else 
      	      match_pos = m_B_st[0].find_open(pos_in_A);
	   //  	cout << "match_pos = " << match_pos <<endl;
	    size_type par = m_B_st[0].parent_t(match_pos);
	     //	cout << "par = " << par <<endl;
      	    return m_B_st[0].rank(par);
      	  }
      	}
      }
*/
      // Traversal of the neighbors of vertex v. If necessary, return all values
      // of variable x to obtain the neighbors in a appropriate format
      // me paro en la primera ocurrencia de v, trabajo en general con posiciones de la secuencia
      /*
      void list_neighbors(size_type v) {
	if(v >= m_vertices)
	  return;
	size_type nxt = first(v);
//	cout << "limite = " << nxt << endl;
	while(nxt <= 2*m_edges+3) {
	  if(nxt <= 2*m_edges+3) {
	    size_type mt = mate(nxt);
	  //  puts("ejecutare vertex");
	    size_type x = vertex(mt); // Print neighbor
	    cout <<"mate = " << mt<< " x = "<<x<<endl;
	//    cout << x <<" ";
	  }
	  nxt = next(nxt);	
//	    cout <<"next = " << nxt <<endl;	  
	}
	puts("");
      }
**/
      vector <int> lista_vecinos(size_type v, int nivel){
      	size_type nxt = first(v,nivel)+1, limite = mate(first(v,nivel),nivel);
      //	cout << " limite = " << limite << endl;
      	vector <int> retorno;
      	while(nxt < limite){
      	//	cout << "NEXT "<< nxt <<endl;
      		size_type x;
      		size_type pos_in_A = m_A_rank[nivel](nxt);
      		if(m_A[nivel][nxt] == 1 && m_B[nivel][pos_in_A] == 1){
      	//			cout << " NEXT = " << nxt << endl;
      	//		puts("entre al if") ;
      			x = vertex2(nxt,nivel);
      			nxt = mate(nxt,nivel)+1;
      		}
      		else{
      	//		puts("entre al else") ;
      	//			cout << " NEXT = " << nxt << endl;
      			x = vertex2(nxt,nivel);
      			nxt ++;
      		}
      	//	cout<<x << " ";
      		retorno.push_back(x);
      	}
      	if(v != 0){
      		nxt = m_B_st[nivel].parent_t(m_B_select1[nivel](v+1));
      	//	cout << m_B_select1[nivel](v+1) << " " << nxt;
      		size_type x = m_B_rank[nivel](m_A_rank[nivel](m_A_select1[nivel](nxt+1)));
      		//cout << x;
      		if(x < limite)retorno.push_back(x);
      	}
      	//cout << endl;
      	return retorno;
      }

      // Traversal of the face where the edge e belongs. If necessary, return
      // all values of variable curr_vertex to obtain the vertices of the face
      // in a appropriate format
/*
      void face(size_type e) {
	if(e >= 2*m_edges)
	  return;
	char flag = 1;
	size_type nxt = e;
	size_type mt;
	size_type init_vertex = vertex(nxt);
	size_type curr_vertex = -1;

	while(curr_vertex != init_vertex || flag) {
	  if(nxt >= 2*m_edges) {
	    nxt = first(vertex(mt));
	  }
	  
	  flag = 0;
	  mt = mate(nxt);
	  curr_vertex = vertex(mt);
	  nxt = next(mt);
	}
      }
      */
//Las funciones que yo implemente van aca


	int getMapeo(size_type pos, unsigned int nivel){
		int retorno = B_select1[nivel](pos+1);
		return t.getposMapeo(m_vertices*nivel+retorno);
	}
	int getReverso(size_type pos, unsigned int nivel){
		return t.getposReverso(m_vertices*nivel+pos);
	}

    unsigned int regionID(unsigned int x, unsigned int level){
	//printf("rank = %ld\n",m_B_rank[level](x));
	return m_B_rank[level](x);
	}
    unsigned int goUp(unsigned int x, unsigned int level){
//	printf("select = %ld\n",B_select1[level](x+1));
	//printf("rank = %ld\n",m_B_rank[level](x));
//	printf("rank = %ld\n",m_B_rank[level](m_B_select1[0](B_select1[level](x+1)+1)));
	//m_B_select1[0](B_select1[level](x));
	return m_B_select1[0](B_select1[level](x+1)+1);
	}



/*
    unsigned int goDown(unsigned int x, unsigned int level){
    int posicion = m_B_select1[0](x+1);
    int q = m_B_select1[0](x+1),p = -1,i=x;
    int closeq = m_B_st[0].find_close(q);
    bool entre = false;
    int contador = 0;
    int actual = m_B_rank[0](posicion);
    int y;
    contador = 0;
    if(B[level][m_B_rank[0](posicion)] == 0)while(true){
    	entre = true;
    	y = B_select1[level](B_rank[level](actual));
    	actual = y;
    	posicion = m_B_select1[0](y+1);
    	p = m_B_select1[0](y+1);
    	if(B[level][y] == 1 && p < q && closeq < m_B_st[0].find_close(p) )break;
    	contador++;
    }
    printf("contador = %d\n",contador);
	unsigned int retorno;
	if(!entre )retorno = (B_rank[level](x));
	else retorno = B_rank[level](m_B_rank[0](posicion));
	return m_B_select1[level](retorno+1);
	}
*/



    unsigned int goDown(unsigned int x, unsigned int level){
	//printf("Down = %ld\n",B_rank[level](B_rank[level](m_B_select1[0](x+1))+1));
    int posicion = m_B_select1[0](x+1);
    bool entre = false;
    int auxiliar4;
    while(B[level][m_B_rank[0](posicion)] == 0){
    	entre = true;
    	posicion = m_B_st[0].parent_t(posicion);

    }
    auxiliar4 = B_rank[level](m_B_rank[0](posicion));
	unsigned int retorno;
	if(!entre )retorno = (B_rank[level](x));
	else retorno = B_rank[level](m_B_rank[0](posicion));
	return m_B_select1[level](retorno+1);
	}







    unsigned int goLevel(unsigned int x, unsigned int levels, unsigned int levelt){
    unsigned int regions = regionID(goUp(x,levels),0);
    unsigned int regiont = regionID(goDown(regions,levelt),levelt);
   // printf("%u\n",regiont);
	//printf("level = %u\n",regionID(goDown(regionID(goUp(x,levels),levels),levelt),levelt));
	return regiont;
	}
	bool Inside(unsigned int x, unsigned int levelx, unsigned int y, unsigned int levely){
		unsigned int comprobar = goLevel(x,levelx,levely);
		//printf("%u ",comprobar);
		//if(comprobar == y)printf("contenido\n");
		//else printf("no contenido\n");
		if(comprobar == y)return true;
		else return false;
	}
	void comprobarvecinos(){
		for(int i = 0; i < 2; i++){
	//		cout << "VECINOS NODO " << i << endl;
	//		lista_vecinos(i,2);
		}
	}	

	bool Touches(unsigned int x, unsigned int levelx, unsigned int y, unsigned int levely){
	//	puts("ENTRE");
		if(levely < levelx){
	//		puts("ENTRE ACA");
			bool comprobar = Inside(y,levely,x,levelx);
			vector <int> vecinos = (lista_vecinos(y,levely));
	//		puts("AUN VIVO");
			for(int i = 0; i < vecinos.size();i++){
				unsigned int contenedor = goLevel(vecinos[i],levely,levelx);
			//	cout << "VECINO " << vecinos[i] << " CONTENEDOR " << contenedor <<endl;
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
		//	puts("ENTRE AL ELSE");
			bool comprobar = Inside(x,levelx,y,levely);
			vector <int> vecinos = (lista_vecinos(x,levelx));
		//	puts("AUN VIVO");
		//	printf("tam = %d\n",vecinos.size());
			for(int i = 0; i < vecinos.size();i++){
		//		printf("vecino = %d\n",vecinos[i]);
				unsigned int contenedor = goLevel(vecinos[i],levelx,levely);
			//				cout << "VECINO " << vecinos[i] << " CONTENEDOR " << contenedor<<endl;
				if(contenedor == y  && !comprobar)return true;
				if(contenedor != y  && comprobar)return true;

			}
		}
		return false;
	} 

	vector <int> Contained(unsigned int x, int levelx, int levely){
		unsigned int a = goLevel(x,levelx,levely);
		size_type nxt = first(a,levely);
		size_type limit = mate(nxt,levely);
		size_type limit2 = mate(first(0,levely),levely);
		limit2 = m_B_rank[levely](m_A_rank[levely](limit2));
		//cout << limit2 << endl;
		vector <int> retorno;
		retorno.push_back(vertex2(nxt,levely));
		//cout << "PUSHEE " << vertex2(nxt,levely) << endl;
		size_type pos_in_B = m_B_rank[levely](m_A_rank[levely](nxt))+2;
		if(pos_in_B > limit2)return retorno;
	//	cout << m_B_select1[levely](pos_in_B+2) << endl;
		nxt = m_A_select1[levely](m_B_select1[levely](pos_in_B)+1);
	//	nxt = m_A_select1[levely](	 m_B_select1[levely](m_B_rank[levely](m_A_rank[levely](nxt))+1) +1);
	//	cout << m_A_rank[levely](nxt) << " " << m_B_select1[levely](m_A_rank[levely](nxt+1)+1) << endl;}
		int contador = 0;
		while(nxt < limit){
			size_type pos_in_A = m_B_rank[levely](m_A_rank[levely](nxt));
			size_type pos = B_select1[levely](pos_in_A+1);
		//	cout  << nxt<< "  " <<  pos_in_A <<"   "  <<pos << endl;
			if(B[levelx][pos] ==1){
				//sputs("ENTRE");
//
				nxt = mate(nxt,levely);
		//		cout << nxt <<endl;
				size_type pos_in_B = m_B_rank[levely](m_A_rank[levely](nxt))+1;
				//cout << pos_in_B << endl;
				if(pos_in_B > limit2)break;
		//		cout << m_A_rank[levely](nxt) << endl;
		//		cout << pos_in_B << endl;
		//		cout <<m_B_select1[levely](pos_in_B+1) << endl;
		//		cout << m_A_select1[levely](m_B_select1[levely](pos_in_B+1)+1) << endl;
		//		cout << m_B_select1[levely]( m_B_rank[levely](m_A_rank[levely](nxt+1))+1) <<  "  "<< m_A_rank[levely](nxt+1)<< endl;
		nxt = m_A_select1[levely](m_B_select1[levely](pos_in_B)+1);
			}
			else if(B[levelx][pos] == 0){
				retorno.push_back(vertex2(nxt,levely));
			//	cout << "PUSHEE " << vertex2(nxt,levely) <<endl;
				size_type pos_in_B = m_B_rank[levely](m_A_rank[levely](nxt))+2;
				if(pos_in_B > limit2)break;
				nxt = m_A_select1[levely](m_B_select1[levely](pos_in_B)+1);
			} 
			contador++;

		}
	//	sort(retorno.begin(),retorno.end());
		return retorno;

	}
	
	void comprobarfirst(){
		for(int i = 0; i < 2; i++)cout <<" "<<getMapeo(i,2);
		puts("");
	}	
	/*
	void comprobarmate(){
		for(int i = 0; i < 88; i++)if(m_A[0][i] == 0){
			int aux = mate(i);
			if(aux >= 0 )cout << " " << aux;
		}
			puts("");
			//mate(i);
	}	
	void comprobarvertex1(){
		for(int i = 1; i < 87; i++){
			int aux = vertex2(i);
			if(aux > -1)cout << " indice "<< i<<" aux " << aux << endl;
		}
			//puts("");
	}
	void comprobarnext(){
		for(int i = 0; i < 21; i++)cout << " "<<vertex(next(first(i)))<<endl;
			puts("");
	}	
	*/
 // Aca deberia terminar la seccion de funciones que yo implemente
  };





//first corregido, falta comprobar para demas niveles, la logica si deberia funcionar
//mate corregido igual para level 0


}// end namespace sdsl
#endif


//FUNCION FIRST RECIBE UN ENTERO QUE CORRESPONDE AL ID DEL NODO AL CUAL SE LE QUIERE SABER SU POSICION EN EL ARREGLO M_A
//