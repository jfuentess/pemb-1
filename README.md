Compact Representations of Spatial Hierarchical Structures with Support for Topological Queries
=========

Requirements
------------

The implementations provided requires:

* Installing the following version of SDSL: https://github.com/jfuentess/sdsl-lite

#### Fast and Compact Hierarchical Planar Embeddings

This project provides three succinct data structures to represent Hierarchical planar embeddings. 

The structures can be used as follows:

Structure based on approach 1

```cpp
#include "graph3.hpp"
#include "auxiliar.hpp"
#include "pemb-ap1.hpp"
#include <complementary/utils.hpp>



using namespace std;
using namespace sdsl;

pemb<bit_vector, sd_vector<> > *pe;

int main(int argc, char **argv){
  // argv[1] is the path to a file with the planar embedding.
  // argv[2] is the path to a file with the hierarchical description.
        int total;
        unsigned int n,niveles,aux,mitotal = 0;
        Graph g = read_graph_from_file(argv[1]);
pe = new pemb<bit_vector , sd_vector<> >(g,0,1,argv);


if(pe->Inside(2,0,7,8))cout << "Region 2 Inside in Region 8"\n;

if(pe->Touches(2,0,7,8))cout << "Region 2 Touches with Region 8"\n;


vector <int> regions = pe->Contained(8,2,1);
cout << "Region 8 contains the following regions at level 1\n";
four(int i = 0; i < regions.size();i++)cout << regions[i] << " ";
cout << endl;
return 0;
}
```

Structure based on approach 2 
```cpp
#include "graph3.hpp"
#include "auxiliar.hpp"
#include "pemb-ap2wt.hpp"
// to use one of the other variant based on approach two, change the include
//#include "pemb-ap2plain.hpp"
//#include "pemb-ap2rmmt.hpp"
#include <complementary/utils.hpp>



using namespace std;
using namespace sdsl;

pemb<> *pe;

int main(int argc, char **argv){
  // argv[1] is the path to a file with the planar embedding.
  // argv[2] is the path to a file with the hierarchical description.
        int total;
        unsigned int n,niveles,aux,mitotal = 0;
        Graph g = read_graph_from_file(argv[1]);
pe = new pemb<>(g,0,1,argv); 
pe->iniresumen(256); // Size of block
if(pe->Inside(2,0,7,8))cout << "Region 2 Inside in Region 8"\n;

if(pe->Touches(2,0,7,8))cout << "Region 2 Touches with Region 8"\n;


vector <int> regions = pe->Contained(8,2,1);
cout << "Region 8 contains the following regions at level 1\n";
four(int i = 0; i < regions.size();i++)cout << regions[i] << " ";
cout << endl;
return 0;
}
```

Structure based on approach 2 with reduced memory
```cpp
#include "graph4.hpp"
#include "pemb-ap2rmmt.hpp"
#include <complementary/utils.hpp>



using namespace std;
using namespace sdsl;

pemb<> *pe;

int main(int argc, char **argv){
  // argv[1] is the path to a file with the planar embedding.
  // argv[2] is the path to a file with the hierarchical description.
        int total;
        unsigned int n,niveles,aux,mitotal = 0;
        Graph g = read_graph_from_file(argv[1]);
pe = new pemb<>(g,0,argv); 
pe->iniresumen(256); // Size of block
if(pe->Inside(2,0,7,8))cout << "Region 2 Inside in Region 8"\n;

if(pe->Touches(2,0,7,8))cout << "Region 2 Touches with Region 8"\n;


vector <int> regions = pe->Contained(8,2,1);
cout << "Region 8 contains the following regions at level 1\n";
four(int i = 0; i < regions.size();i++)cout << regions[i] << " ";
cout << endl;
return 0;
}
```

Structure based on approach 3
```cpp
#include "graph3.hpp"
#include "auxiliar.hpp"
#include "pemb-ap3.hpp"
#include <complementary/utils.hpp>



using namespace std;
using namespace sdsl;

pemb<> *pe;

int main(int argc, char **argv){
  // argv[1] is the path to a file with the planar embedding.
  // argv[2] is the path to a file with the hierarchical description.
        int total;
        unsigned int n,niveles,aux,mitotal = 0;
        Graph g = read_graph_from_file(argv[1]);
pe = new pemb<>(g,0,argv);


if(pe->Inside(2,0,7,8))cout << "Region 2 Inside in Region 8"\n;

if(pe->Touches(2,0,7,8))cout << "Region 2 Touches with Region 8"\n;


vector <int> regions = pe->Contained(8,2,1);
cout << "Region 8 contains the following regions at level 1\n";
four(int i = 0; i < regions.size();i++)cout << regions[i] << " ";
cout << endl;
return 0;
}
```


to compile, just run:


```sh
g++ -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib program.cpp -o program -lsdsl
```

[STL]: http://www.sgi.com/tech/stl/ "Standard Template Library"
[pz]: http://pizzachili.di.unipi.it/ "Pizza&amp;Chli"
[d3js]: http://d3js.org "D3JS library"
[cmake]: http://www.cmake.org/ "CMake tool"
[MAKE]: http://www.gnu.org/software/make/ "GNU Make"
[gcc]: http://gcc.gnu.org/ "GNU Compiler Collection"
[DIVSUF]: https://github.com/y-256/libdivsufsort/ "libdivsufsort"
[LS]: http://www.sciencedirect.com/science/article/pii/S0304397507005257 "Larson &amp; Sadakane Algorithm"
[GTEST]: https://code.google.com/p/googletest/ "Google C++ Testing Framework"
[SDSLCS]: http://simongog.github.io/assets/data/sdsl-cheatsheet.pdf "SDSL Cheat Sheet"
[SDSLLIT]: https://github.com/simongog/sdsl-lite/wiki/Literature "Succinct Data Structure Literature"
[TUT]: http://simongog.github.io/assets/data/sdsl-slides/tutorial "Tutorial"
[QSUFIMPL]: http://www.larsson.dogma.net/qsufsort.c "Original Qsufsort Implementation"
[JESL]: http://www.itu.dk/people/jesl/ "Homepage of Jesper Larsson"
[CF]: https://github.com/simongog/sdsl-lite/blob/master/COPYING "Licence"
[SEAPAPER]: http://arxiv.org/pdf/1311.1249v1.pdf "SDSL paper"
[HB]: https://github.com/simongog/sdsl-lite/blob/hybrid_bitvector/include/sdsl/hybrid_vector.hpp "Hybrid bitevctor"
[DOXYGENDOCS]: http://algo2.iti.kit.edu/gog/docs/html/index.html "API Reference"
