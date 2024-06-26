#ifndef TEMPLAT_LATTICE_NPARTITEITERATORS_NPARTITEITERATIONTESTER_H
#define TEMPLAT_LATTICE_NPARTITEITERATORS_NPARTITEITERATIONTESTER_H
 
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */ 
   
// File info: Main contributor(s): Adrien Florio,  Year: 2019

#include "TempLat/util/tdd/tdd.h"

namespace TempLat {


    /** \brief A class which tests the npartite iterations.
     *
     * 
     * Unit test: make test-npartiteiterationtester
     **/

    class NPartiteLooperTester {
    public:
        /* Put public methods here. These should change very little over time. */
        NPartiteLooperTester() {

        }


    private:
        /* Put all member variables and private methods here. These may change arbitrarily. */



    public:
#ifdef TEMPLATTEST
        static inline void Test(TDDAssertion& tdd);
#endif
    };



} /* TempLat */

#ifdef TEMPLATTEST
#include "TempLat/lattice/loopers/npartiteloopers/npartiteloopertester_test.h"
#endif


#endif
