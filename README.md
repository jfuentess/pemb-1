# Code and Datasets 

Code and datasets for the replication of the work "Compact epresentations of spatial hierarchical structures with support for topological queries"

## Pre-requirements

* The implementations are based on the SDSL library. Before executing our codes, install the SDSL library from https://github.com/jfuentess/sdsl-lite and follow the instructions.


## Datasets
All the datasets corresponds to subsets of the [TIGER dataset](https://www2.census.gov/geo/tiger/TIGER2019/), provided by the U.S. Census Bureau. The files are divided into two groups:

#### Planar embeddings
* tiger_8states.pg.tar.gz: Planar graph embedding representing the geographic faces of the neighboring states of Nevada, Utah, Arizona, Colorado, New Mexico, Kansas, Oklahoma and Texas, with a total of 4,761,354 vertices, 10,326,904 edges and 1 connected component.

* tiger_usa.pg.tar.gz: Planar graph embedding representing the geographic faces of the whole continental part of U.S., with 19,735,874 vertices and 43,837,150 edges.

* whole_usa.pg.tar.gz: Planar graph embedding representing the geographic faces of the whole continental part of U.S., disconnected subregions, and Alaska, Hawaii and overseas U.S. islands. The dataset consists of 20,037,199 vertices, 44,503,624 edges and 98 connected components.

#### Hierarchy data

* tiger_8states.hierarchy.tar.gz: Hierarchy information of the dataset `tiger_8states.pg`. It considers five aggregation levels: states, counties, census tracts, census block groups and census blocks.

* tiger_usa.hierarchy.tar.gz: Hierarchy information of the dataset `tiger_usa.pg`. It considers five aggregation levels: states, counties, census tracts, census block groups and census blocks.

* tiger_usa_plus.hierarchy.tar.gz:  Artificial hierarchy information for the dataset `tiger_usa.pg`. It considers five aggregation levels: states, counties, census tracts, census block groups and census blocks. The artificial hierarchy is conformed by grouping from 1 up to 10 contiguous regions into one, across the five aggregation levels.
* whole_usa.hierarchy.tar.gz: Hierarchy information of the dataset `whole_usa.pg`. It considers five aggregation levels: states, counties, census tracts, census block groups and census blocks.

## Implementations

The implementation files are available in the file `code.zip`, and are divided into two groups:

#### Representations of spatial hierarchical structures
* pemb-ap1.hpp: Implementation of the approach 1.
* pemb-ap2plain.hpp: Implementation of the approach 2 which uses a compact planar embedding on each granularity level and range min-max trees to store hierarchy information. 
* pemb-ap2rmmt.hpp: Implementation of the approach 2 which uses a compact planar embedding for the highest level of detail, and a combination of range min-max trees and wavelet trees to store hierarchy information.
* pemb-ap2wt.hpp: Implementation of the approach 2: Implementation of the approach 2 which uses a compact planar embedding on each granularity level and wavelet trees to store hierarchy information.
* pemb-ap3.hpp: Implementation of the approach 3.
* baseline.hpp: Implementation of a baseline which uses a compact planar embedding on each granularity level and a non-compact array to store hierarchy information. 

#### Auxiliary files
* graph3.hpp: 
* graph4.hpp: 
* Node.hpp: 
* auxiliar.hpp: 
* auxiliar2.hpp: 

## Usage

Now, we show some examples of how to use our implementations

#### Approach 1 (pemb-ap1.hpp)

```cpp
#include "graph3.hpp"
#include "auxiliar.hpp"
#include "pemb-ap1.hpp"
#include <complementary/utils.hpp>

using namespace std;
using namespace sdsl;

pemb<bit_vector, sd_vector<> > *pe;

// argv[1] is the path to a file with the planar embedding.
// argv[2] is the path to a file with the hierarchy information.
int main(int argc, char **argv) {
  // Load the input planar embedding
  Graph g = read_graph_from_file(argv[1]);
  
  // Build the compact representation with starting vertex 0 and threshold 1
  pe = new pemb<bit_vector , sd_vector<> >(g, 0, 1, argv);

  if(pe->Inside(2,0,7,8))
    cout << "Region 2 Inside in Region 8" << endl;

  if(pe->Touches(2,0,7,8))
  cout << "Region 2 Touches with Region 8" << endl;
  
  vector <int> regions = pe->Contained(8,2,1);
  cout << "Region 8 contains the following regions at level 1" << endl;
  for(int i = 0; i < regions.size();i++)
    cout << regions[i] << " ";
  cout << endl;

  return 0;
}
```

#### Approach 2 (pemb-ap2wt.hpp and pemb-ap2plain.hpp)
```cpp
#include "graph3.hpp"
#include "auxiliar.hpp"
// To use the variant pemb-ap2plain.hpp, change to include #include "pemb-ap2plain.hpp"
#include "pemb-ap2wt.hpp"
#include <complementary/utils.hpp>

using namespace std;
using namespace sdsl;

pemb<> *pe;

// argv[1] is the path to a file with the planar embedding.
// argv[2] is the path to a file with the hierarchy information.
int main(int argc, char **argv){
  // Load the input planar embedding
  Graph g = read_graph_from_file(argv[1]);
  
  // Build the compact representation with starting vertex 0 and threshold 1
  pe = new pemb<>(g,0,1,argv); 
  
  pe->iniresumen(256); // Size of block
  if(pe->Inside(2,0,7,8))
    cout << "Region 2 Inside in Region 8" << endl;

  if(pe->Touches(2,0,7,8))
    cout << "Region 2 Touches with Region 8" << endl;

  vector <int> regions = pe->Contained(8,2,1);
  cout << "Region 8 contains the following regions at level 1" << endl;
  for(int i = 0; i < regions.size();i++)
    cout << regions[i] << " ";
  cout << endl;

  return 0;
}
```

#### Approach 2 (pemb-ap2rmmt.hpp)

```cpp
#include "graph4.hpp"
#include "pemb-ap2rmmt.hpp"
#include <complementary/utils.hpp>

using namespace std;
using namespace sdsl;

pemb<> *pe;

// argv[1] is the path to a file with the planar embedding.
// argv[2] is the path to a file with the hierarchy information.
int main(int argc, char **argv){

  // Load the input planar embedding
  Graph g = read_graph_from_file(argv[1]);
  
  // Build the compact representation with starting vertex 0
  pe = new pemb<>(g,0,argv); 

  pe->iniresumen(256); // Size of block
  
  if(pe->Inside(2,0,7,8))
    cout << "Region 2 Inside in Region 8" << endl;

  if(pe->Touches(2,0,7,8))
    cout << "Region 2 Touches with Region 8" << endl;

  vector <int> regions = pe->Contained(8,2,1);
  cout << "Region 8 contains the following regions at level 1" << endl;
  for(int i = 0; i < regions.size();i++)
    cout << regions[i] << " ";
  cout << endl;

  return 0;
}
```

#### Approach 3 (pemb-ap3.hpp)
```cpp
#include "graph3.hpp"
#include "auxiliar.hpp"
#include "pemb-ap3.hpp"
#include <complementary/utils.hpp>

using namespace std;
using namespace sdsl;

pemb<> *pe;

// argv[1] is the path to a file with the planar embedding.
// argv[2] is the path to a file with the hierarchy information.
int main(int argc, char **argv){
  // Load the input planar embedding
  Graph g = read_graph_from_file(argv[1]);
  
  // Build the compact representation with starting vertex 0
  pe = new pemb<>(g,0,argv); 

  if(pe->Inside(2,0,7,8))
    cout << "Region 2 Inside in Region 8" << endl;

  if(pe->Touches(2,0,7,8))
    cout << "Region 2 Touches with Region 8" << endl;

  vector <int> regions = pe->Contained(8,2,1);
  cout << "Region 8 contains the following regions at level 1" << endl;
  for(int i = 0; i < regions.size();i++)
    cout << regions[i] << " ";
  cout << endl;
  
  return 0;
}
```

## Compilation
to compile, just run, where `~/include` and `~/lib` correspond to the installation folders of the SDSL library:

```sh
g++ -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib program.cpp -o program -lsdsl
```
