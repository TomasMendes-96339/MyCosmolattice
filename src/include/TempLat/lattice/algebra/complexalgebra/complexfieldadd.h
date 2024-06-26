#ifndef COSMOINTERFACE_COMPLEXFIELDALGEBRA_COMPLEXFIELDADD_H
#define COSMOINTERFACE_COMPLEXFIELDALGEBRA_COMPLEXFIELDADD_H
 
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */ 
   
// File info: Main contributor(s): Adrien Florio,  Year: 2019

#include "TempLat/util/tdd/tdd.h"
#include "real.h"
#include "imag.h"
#include "TempLat/lattice/algebra/complexalgebra/complexwrapper.h"
#include "TempLat/lattice/algebra/helpers/getstring.h"
#include "TempLat/lattice/algebra/complexalgebra/complexfieldbinaryoperator.h"
#include "TempLat/lattice/algebra/helpers/hasstaticgetter.h"

namespace TempLat {


    /** \brief A class which implements complex field addition.

     *
     * 
     * Unit test: make test-complexfieldadd
     **/

    template < class R, class T >
    class ComplexFieldAddition : public ComplexFieldBinaryOperator<R,T>{
    public:

        using ComplexFieldBinaryOperator<R,T>::mR;
        using ComplexFieldBinaryOperator<R,T>::mT;

        /* Put public methods here. These should change very little over time. */
        ComplexFieldAddition(const R& pR, const T& pT) :
                ComplexFieldBinaryOperator<R,T>(pR,pT){

        }

        auto ComplexFieldGet(Tag<0> t)
        {
            return Real(mR) + Real(mT);
        }
        auto ComplexFieldGet(Tag<1> t)
        {
            return Imag(mR) + Imag(mT);
        }

        auto ComplexFieldGet(Tag<0> t, ptrdiff_t i)
        {
            return mR.ComplexFieldGet(0_c,i) + mT.ComplexFieldGet(0_c,i);
        }
        auto ComplexFieldGet(Tag<1> t, ptrdiff_t i)
        {
            return mR.ComplexFieldGet(1_c,i) + mT.ComplexFieldGet(1_c,i);
        }

        void eval(ptrdiff_t i)
        {
            DoEval::eval(mR, i);
            DoEval::eval(mT, i);
        }

        std::string operatorString() const {
            return "+";
        }

    private:
        /* Put all member variables and private methods here. These may change arbitrarily. */


    };

    struct ComplexFieldAddTester{
    public:
#ifdef TEMPLATTEST
        static inline void Test(TDDAssertion& tdd);
#endif
    };


    template <typename R, typename T>
    typename std::enable_if<HasComplexFieldGet<R>::value && HasComplexFieldGet<T>::value, ComplexFieldAddition<R,T>>::type
    operator+(const R& r, const T& t)
    {
        return {r,t};
    }

    template <typename R, typename T>
    typename std::enable_if<!HasStaticGetter<R>::value && HasComplexFieldGet<T>::value, ComplexFieldAddition<ComplexFieldWrapper<R,ZeroType>,T>>::type
    operator+(const R& r, const T& t)
    {
        return {Complexify(r,ZeroType()),t};
    }

    template <typename R, typename T>
    typename std::enable_if<!HasStaticGetter<T>::value && HasComplexFieldGet<R>::value, ComplexFieldAddition<R,ComplexFieldWrapper<T,ZeroType>>>::type
    operator+(const R& r, const T& t)
    {
        return {r,Complexify(t,ZeroType())};
    }



} /* TempLat */

#ifdef TEMPLATTEST
#include "TempLat/lattice/algebra/complexalgebra/complexfieldadd_test.h"
#endif


#endif
