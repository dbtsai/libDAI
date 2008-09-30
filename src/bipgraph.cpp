/*  Copyright (C) 2006-2008  Joris Mooij  [joris dot mooij at tuebingen dot mpg dot de]
    Radboud University Nijmegen, The Netherlands /
    Max Planck Institute for Biological Cybernetics, Germany

    This file is part of libDAI.

    libDAI is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    libDAI is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with libDAI; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/


#include <dai/bipgraph.h>


namespace dai {


using namespace std;


void BipartiteGraph::erase1( size_t n1 ) {
    assert( n1 < nr1() );
    // Erase neighbor entry of node n1
    _nb1.erase( _nb1.begin() + n1 );
    // Adjust neighbor entries of nodes of type 2
    for( size_t n2 = 0; n2 < nr2(); n2++ ) {
        for( size_t iter = 0; iter < nb2(n2).size(); ) {
            Neighbor &m1 = nb2(n2, iter);
            if( m1.node == n1 ) {
                // delete this entry, because it points to the deleted node
                nb2(n2).erase( nb2(n2).begin() + iter );
            } else if( m1.node > n1 ) {
                // update this entry and the corresponding dual of the neighboring node of type 1
                m1.iter = iter;
                m1.node--;
                nb1( m1.node, m1.dual ).dual = iter;
                iter++;
            } else {
                // skip
                iter++;
            }
        }
    }
}


void BipartiteGraph::erase2( size_t n2 ) {
    assert( n2 < nr2() );
    // Erase neighbor entry of node n2
    _nb2.erase( _nb2.begin() + n2 );
    // Adjust neighbor entries of nodes of type 1
    for( size_t n1 = 0; n1 < nr1(); n1++ ) {
        for( size_t iter = 0; iter < nb1(n1).size(); ) {
            Neighbor &m2 = nb1(n1, iter);
            if( m2.node == n2 ) {
                // delete this entry, because it points to the deleted node
                nb1(n1).erase( nb1(n1).begin() + iter );
            } else if( m2.node > n2 ) {
                // update this entry and the corresponding dual of the neighboring node of type 2
                m2.iter = iter;
                m2.node--;
                nb2( m2.node, m2.dual ).dual = iter;
                iter++;
            } else {
                // skip
                iter++;
            }
        }
    }
}


std::vector<size_t> BipartiteGraph::delta1( size_t n1, bool include ) const {
    // get all second-order neighbors
    std::vector<size_t> result;
    foreach( const Neighbor &n2, nb1(n1) )
        foreach( const Neighbor &m1, nb2(n2) )
            if( include || (m1 != n1) )
                result.push_back( m1 );
    // remove duplicates
    std::vector<size_t>::iterator it = std::unique( result.begin(), result.end() );
    result.erase( it, result.end() );
    return result;
}


std::vector<size_t> BipartiteGraph::delta2( size_t n2, bool include ) const {
    // store all second-order neighbors
    std::vector<size_t> result;
    foreach( const Neighbor &n1, nb2(n2) )
        foreach( const Neighbor &m2, nb1(n1) )
            if( include || (m2 != n2) )
                result.push_back( m2 );
    // remove duplicates
    std::vector<size_t>::iterator it = std::unique( result.begin(), result.end() );
    result.erase( it, result.end() );
    return result;
}


bool BipartiteGraph::isConnected() const {
    if( nr1() == 0 ) {
        return true;
    } else {
        std::vector<bool> incomponent1( nr1(), false );
        std::vector<bool> incomponent2( nr2(), false );

        incomponent1[0] = true;
        bool found_new_nodes;
        do {
            found_new_nodes = false;

            // For all nodes of type 2, check if they are connected with the (growing) component
            for( size_t n2 = 0; n2 < nr2(); n2++ )
                if( !incomponent2[n2] ) {
                    foreach( const Neighbor &n1, nb2(n2) ) {
                        if( incomponent1[n1] ) {
                            found_new_nodes = true;
                            incomponent2[n2] = true;
                            break;
                        }
                    }
                }

            // For all nodes of type 1, check if they are connected with the (growing) component
            for( size_t n1 = 0; n1 < nr1(); n1++ )
                if( !incomponent1[n1] ) {
                    foreach( const Neighbor &n2, nb1(n1) ) {
                        if( incomponent2[n2] ) {
                            found_new_nodes = true;
                            incomponent1[n1] = true;
                            break;
                        }
                    }
                }
        } while( found_new_nodes );

        // Check if there are remaining nodes (not in the component)
        bool all_connected = true;
        for( size_t n1 = 0; (n1 < nr1()) && all_connected; n1++ )
            if( !incomponent1[n1] )
                all_connected = false;
        for( size_t n2 = 0; (n2 < nr2()) && all_connected; n2++ )
            if( !incomponent2[n2] )
                all_connected = false;

        return all_connected;
    }
}


bool BipartiteGraph::isTree() const {
    using namespace std;
    vector<levelType> levels;

    bool foundCycle = false;
    size_t nr_1 = 0;
    size_t nr_2 = 0;

    if( nr1() == 0 || nr2() == 0 )
        return true;
    else {
        levelType newLevel;
        do {
            newLevel.ind1.clear();
            newLevel.ind2.clear();
            if( levels.size() == 0 ) {
                size_t n1 = 0;
                // add n1 to ind1
                newLevel.ind1 = vector<size_t>( 1, n1 );
                // add all neighbors of n1 to ind2
                newLevel.ind2.reserve( nb1(n1).size() );
                foreach( const Neighbor &n2, nb1(n1) )
                    newLevel.ind2.push_back( n2 );
            } else {
                const levelType &prevLevel = levels.back();
                // build newLevel.ind1
                foreach( size_t n2, prevLevel.ind2 ) { // for all n2 in the previous level
                    foreach( const Neighbor &n1, nb2(n2) ) { // for all neighbors n1 of n2
                        if( find( prevLevel.ind1.begin(), prevLevel.ind1.end(), n1 ) == prevLevel.ind1.end() ) { // n1 not in previous level
                            if( find( newLevel.ind1.begin(), newLevel.ind1.end(), n1 ) != newLevel.ind1.end() )
                                foundCycle = true; // n1 already in new level: we found a cycle
                            else
                                newLevel.ind1.push_back( n1 ); // add n1 to new level
                        }
                        if( foundCycle )
                            break;
                    }
                    if( foundCycle )
                        break;
                }
                // build newLevel.ind2
                foreach( size_t n1, newLevel.ind1 ) { // for all n1 in this level
                    foreach( const Neighbor &n2, nb1(n1) ) { // for all neighbors n2 of n1
                        if( find( prevLevel.ind2.begin(), prevLevel.ind2.end(), n2 ) == prevLevel.ind2.end() ) { // n2 not in previous level
                            if( find( newLevel.ind2.begin(), newLevel.ind2.end(), n2 ) != newLevel.ind2.end() )
                                foundCycle = true; // n2 already in new level: we found a cycle
                            else
                                newLevel.ind2.push_back( n2 ); // add n2 to new level
                        }
                        if( foundCycle )
                            break;
                    }
                    if( foundCycle )
                        break;
                } 
            }
            levels.push_back( newLevel );
            nr_1 += newLevel.ind1.size();
            nr_2 += newLevel.ind2.size();
        } while( ((newLevel.ind1.size() != 0) || (newLevel.ind2.size() != 0)) && !foundCycle );
        if( nr_1 == nr1() && nr_2 == nr2() && !foundCycle )
            return true;
        else
            return false;
    }
}


void BipartiteGraph::printDot( std::ostream& os ) const {
    using namespace std;
    os << "graph G {" << endl;
    os << "node[shape=circle,width=0.4,fixedsize=true];" << endl;
    for( size_t n1 = 0; n1 < nr1(); n1++ )
        os << "\tx" << n1 << ";" << endl;
    os << "node[shape=box,width=0.3,height=0.3,fixedsize=true];" << endl;
    for( size_t n2 = 0; n2 < nr2(); n2++ )
        os << "\ty" << n2 << ";" << endl;
    for( size_t n1 = 0; n1 < nr1(); n1++ )
        foreach( const Neighbor &n2, nb1(n1) )
            os << "\tx" << n1 << " -- y" << n2 << ";" << endl;
    os << "}" << endl;
}


void BipartiteGraph::check() const {
    size_t N1 = nr1();
    size_t N2 = nr2();
    for( size_t n1 = 0; n1 < N1; n1++ ) {
        size_t iter = 0;
        foreach( const Neighbor &n2, nb1(n1) ) {
            assert( n2.iter == iter );
            assert( n2.node < N2 );
            assert( n2.dual < nb2(n2).size() );
            assert( nb2(n2, n2.dual) == n1 );
            iter++;
        }
    }
    for( size_t n2 = 0; n2 < N2; n2++ ) {
        size_t iter = 0;
        foreach( const Neighbor &n1, nb2(n2) ) {
            assert( n1.iter == iter );
            assert( n1.node < N1 );
            assert( n1.dual < nb1(n1).size() );
            assert( nb1(n1, n1.dual) == n2 );
            iter++;
        }
    }
}


} // end of namespace dai
