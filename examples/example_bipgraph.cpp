#include <dai/bipgraph.h>

using namespace std;
using namespace dai;

int main() {
    // Create a list of edges
    vector<BipartiteGraph::Edge> edges;
    edges.reserve( 5 );
    edges.push_back( BipartiteGraph::Edge(0, 0) );
    edges.push_back( BipartiteGraph::Edge(1, 0) );
    edges.push_back( BipartiteGraph::Edge(2, 0) );
    edges.push_back( BipartiteGraph::Edge(1, 1) );
    edges.push_back( BipartiteGraph::Edge(2, 1) );

    // Create a bipartite graph with 3 nodes of type 1,
    // 2 nodes of type 2 and edge list edges.
    BipartiteGraph G( 3, 2, edges.begin(), edges.end() );

    // Display some information about G
    cout << "G has " << G.nr1() << " nodes of type 1, " << G.nr2() << " nodes of type 2 and " << G.nrEdges() << " edges." << endl << endl;

    // Iterate over all nodes n1 of type 1
    for( size_t n1 = 0; n1 < G.nr1(); n1++ ) {
        cout << "Node " << n1 << " of type 1 has " << G.nb1(n1).size() << " neighbors:" << endl;
        // Iterate over all neighbors n2 of n1
        foreach( const BipartiteGraph::Neighbor &n2, G.nb1(n1) ) {
            // The n2.iter'th neighbor of n1 is n2:
            assert( G.nb1(n1)[n2.iter] == n2 );

            // The n2.dual'th neighbor of n2 is n1:
            assert( G.nb2(n2)[n2.dual] == n1 );

            // n2 can be used as an abbreviation of n2.node:
            assert( static_cast<size_t>(n2) == n2.node );

            cout << "  the " << n2.iter << "'th neighbor is node " << n2 << " of type 2" << endl;
        }
        cout << endl;
    }
}
