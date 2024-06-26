#ifndef COSMOINTERFACE_ENERGIES_H
#define COSMOINTERFACE_ENERGIES_H
 
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */ 
   
// File info: Main contributor(s): Daniel G. Figueroa, Adrien Florio, Francisco Torrenti,  Year: 2020

#include "CosmoInterface/evolvers/evolvertype.h"
#include "TempLat/lattice/algebra/operators/operators.h"

namespace TempLat {


    /** \brief A class which computes the energy when working with our verlet-type algorithms. Assumes that the averages have already been computed.
     *
     * 
     **/

    class Energies {
    public:
        /* Put public methods here. These should change very little over time. */
        Energies() {

        }
        
        // The following functions compute the different volume-averaged contributions to the energy density.
        // The required field functionals are computed in fieldfunctionals.h, and then passed to these functions (via the "fldf" parameter) to rescale them by the appropriate factors.
        
        // EXAMPLE: The kinetic energy of a given field can be called as "Energies::kineticS(model,FieldFunctionals::pi2S(model,i))" (see energiesmeasurer.h for other examples)
        // In this case, FieldFunctionals::pi2S gives <pi^2>, with pi=a^(3-alpha)*f' the momenta, and f' the time-derivative. Then, Energies::kineticS multiplies this by (1/2)*a^(-6) to obtain the energy.

        template<class Model, class T>  
        static inline auto kineticS(Model& model, T fldf) {  // scalar singlet: kinetic energy
            return 0.5 * fldf* pow<-6>(model.aI);
        }

        template<class Model, class T>  
        static inline auto gradientS(Model& model, T fldf) {  // scalar singlet: gradient energy
            return 0.5 * fldf * pow<-2>(model.aI);
        }

        //------------------------ADDED BY TOMÁS-------------------------------------------------
        //------------------for energy densities
        template<class Model, class T>  
        static inline auto kineticS2(Model& model, T fldf) {  // scalar singlet: kinetic energy
            return 8 * model.xi * 3 *(2*fldf)* pow<-3>(model.aI) * (model.aDotI/model.aI)*pow<2>(model.cStar)/pow<2>(model.cutoff)* pow<2>(model.aDotI/model.aI); 
        }

        template<class Model, class T>
        static inline auto laplacian2S(Model& model, T fldf) {  // scalar singlet: laplacian energy
            return -8 * model.xi *fldf * pow<-2>(model.aI) * pow<2>(model.cStar)/pow<2>(model.cutoff)* pow<2>(model.aDotI/model.aI);
        }

        //--------------------for pressure
        template<class Model, class T>
        static inline auto PkineticS1(Model& model, T fldf) {  // scalar singlet: kinetic21
            return fldf* pow<-6>(model.aI) *0.5;
        }
        template<class Model, class T>
        static inline auto PkineticS2(Model& model, T fldf) {  // scalar singlet: kinetic22
            return fldf* pow<-6>(model.aI) *(-16*model.xi*pow<2>(model.aDotI/model.aI)*pow<2>(model.cStar/model.cutoff));
        }
        
        template<class Model, class T>  
        static inline auto PgradientS(Model& model, T fldf) {  // scalar singlet: gradient energy
            return -(1/6)* fldf * pow<-2>(model.aI);
        }

        template<class Model, class T>  
        static inline auto Pfld4S1(Model& model, T fldf) {  // scalar singlet: phi^4 (potential term)
            return - model.lambda * fldf * (0.25);
        }

        template<class Model, class T>  
        static inline auto Pfld2S1(Model& model, T fldf) {  // scalar singlet: m^2 * phi^2 (potential term)
            return - pow<2>(model.M/model.cStar) * fldf * (0.5);
        }

        template<class Model, class T>  
        static inline auto Pfld4S2(Model& model, T fldf) {  // scalar singlet: phi^4 (interaction)
            return model.lambda * fldf * (16*model.xi*pow<2>(model.aDotI/model.aI)*pow<2>(model.cStar/model.cutoff));
        }


        template<class Model, class T>
        static inline auto Plaplacian2S(Model& model, T fldf) {  // scalar singlet: \nabla^2\chi^2 (interaction)
            //return 16 *model.xi/3 *pow<2>(model.cStar)/pow<2>(model.cutoff)* (-0.5 * (1+0.333) * pow<2>(model.aDotI)/model.aI) *pow<-3>(model.aI) * fldf; //substituted addot with friedmann eqs. 
            return -8*model.xi/3*pow<2>(model.cStar/model.cutoff)*(1+3*0.333)*pow<2>(model.aDotI/model.aI)*pow<-2>(model.aI)*fldf;
        }

        template<class Model, class T>
        static inline auto PfldPiS(Model& model, T fldf) {  // scalar singlet: \chi\chi' (remember f' = pi/a^3) (interaction)
            //return -16 *model.xi *pow<2>(model.cStar)/pow<2>(model.cutoff) * (model.aDotI/model.aI) * (2*(-0.5 * (1+0.333) * pow<2>(model.aDotI)/model.aI)/model.aI-3*pow<2>(model.aDotI/model.aI)) * fldf * pow<-3>(model.aI);
            return 16*model.xi*pow<2>(model.cStar/model.cutoff)*(4+3*0.333)*pow<3>(model.aDotI/model.aI)*fldf;        
        }

        template<class Model, class T>
        static inline auto PfldLaplacianS(Model& model, T fldf) {  // scalar singlet: \chi\nabla^2\chi (interaction)
            return -16 * model.xi* pow<2>(model.cStar/model.cutoff)* pow<2>(model.aDotI/model.aI) * pow<-2>(model.aI) * fldf;
        }

        template<class Model, class T>
        static inline auto Pfld2S(Model& model, T fldf) {  // scalar singlet: \chi^2 (interaction)
            
            //return (16 *model.xi/3 * pow<2>(model.cStar)/pow<2>(model.cutoff)*(-72*pow<6>(model.aDotI/model.aI)*model.xi*(1+3*0.333)*pow<2>(model.cStar)/pow<2>(model.cutoff))+16*pow<2>(model.m)*2*model.xi*pow<2>(model.aDotI/model.aI)*pow<2>(model.cStar)/pow<2>(model.cutoff))*fldf;
            //return -(24*16) *pow<2>(12) *pow<2>(48*12) * pow<6>(model.aDotI/model.aI) * (1+3*0.333) * fldf;
            return -16*model.xi*pow<2>(model.cStar/model.cutoff)*pow<2>(model.aDotI/model.aI)*(24*model.xi*pow<2>(model.cStar/model.cutoff)*pow<4>(model.aDotI/model.aI)*(1+3*0.333)-pow<2>(model.M/model.cStar))*fldf;
     }  

        //---------------------------------------------------------------------------------------
        template<class Model, class T>  
        static inline auto kineticCS(Model& model, T fldf) {  // complex scalar: kinetic energy
            return  fldf * pow<-6>(model.aI);
        }

        template<class Model, class T>
        static inline auto gradientCS(Model& model, T fldf) {   // complex scalar: gradient energy
            return  fldf * pow<-2>(model.aI);
        }
		
        template<class Model, class T>
        static inline auto kineticSU2Doublet(Model& model, T fldf) {  // SU(2) doublet: kinetic energy
            return  fldf * pow<-6>(model.aI);
        }

        template<class Model, class T>
        static inline auto gradientSU2Doublet(Model& model, T fldf) {   // SU(2) doublet: gradient energy
            return fldf * pow<-2>(model.aI);
        }
		
        template<class Model, class T>    // U(1) gauge field: electric energy
        static inline auto electricU1(Model& model, T fldf) {  //This expects field already rescaled by the coupling constant. That ways, works like scalar, for different gauge fields.
            return 0.5 * pow<2>( model.omegaStar / model.fStar) * pow<-4>(model.aI) * fldf;
        }

        template<class Model, class T>
        static inline auto magneticU1(Model& model, T fldf) {   // U(1) gauge field: magnetic energy
            return 0.5 * pow<2>( model.omegaStar / model.fStar)  * pow<-4>(model.aI) * fldf;
        }
		
        template<class Model, class T>   // SU(2) gauge field: electric energy
        static inline auto electricSU2(Model& model, T fldf)
        {
            return  0.5 * pow<2>( model.omegaStar / model.fStar)  * pow<-4>(model.aI) * fldf;
        }

        template<class Model, class T>    // SU(2) gauge field: magnetic energy
        static inline auto magneticSU2(Model& model, T fldf)
        {
            return  0.5 * pow<2>( model.omegaStar / model.fStar)  * pow<-4>(model.aI) * fldf;
        }
        
        // The below functions are the same with fldf=model.pi2AvI,model.grad2AvI...
        template<class Model>
        static inline auto kineticS(Model& model) {
            return kineticS(model, model.pi2AvI);
        }

        template<class Model>
        static inline auto gradientS(Model& model) {
            return gradientS(model, model.grad2AvI);
        }

        //-----------------ADDED BY TOMÁS-------------------------------------------------------------------- 

        template<class Model>
        static inline auto kineticS2(Model& model) {
            return kineticS2(model, model.pi2S2AvI);
        }

        template<class Model>
        static inline auto laplacian2S(Model& model) {
            return laplacian2S(model, model.lap2AvI);
        }

        //----------------For Pressure-----------------------------------------------------------------------
        template<class Model>
        static inline auto PkineticS1(Model& model) {
            return PkineticS1(model, model.Ppi2AvI);
        }
        template<class Model>
        static inline auto PkineticS2(Model& model) {
            return PkineticS2(model, model.Ppi2AvI);
        }

        template<class Model>
        static inline auto PgradientS(Model& model) {
            return PgradientS(model, model.grad2AvI);
        }

        /*template<class Model>
        static inline auto Pfld4S1(Model& model) {
            return Pfld4S1(model, model.Pfld4AvI);
        }*/

        template<class Model>
        static inline auto Pfld4S2(Model& model) {
            return Pfld4S2(model, model.Pfld4AvI);
        }

        template<class Model>
        static inline auto Plaplacian2S(Model& model) {
            return Plaplacian2S(model, model.Plap2AvI);
        }

        template<class Model>
        static inline auto PfldPiS(Model& model) {
            return PfldPiS(model, model.PfldPiAvI);
        }

        template<class Model>
        static inline auto PfldLaplacianS(Model& model) {
            return PfldLaplacianS(model, model.PfldlapAvI);
        }

        template<class Model>
        static inline auto Pfld2S(Model& model) {
            return Pfld2S(model, model.Pfld2AvI);
        }

        /*template<class Model>
        static inline auto Pfld2S1(Model& model) {
            return Pfld2S1(model, model.Pfld2AvI);
        }*/
      //----------------------------------------------------------------------------------------------

         template<class Model>
        static inline auto kineticSSI(Model& model) {
            return kineticS(model, model.pi2AvSI);
        }

        template<class Model>
        static inline auto gradientSSI(Model& model) {
            return gradientS(model, model.grad2AvSI);
        }
		
        template<class Model>
        static inline auto kineticCS(Model& model) {
            return kineticCS(model, model.CSpi2AvI);
        }

        template<class Model>
        static inline auto gradientCS(Model& model) {
            return gradientCS(model, model.CSgrad2AvI);
        }
		
		 template<class Model>
        static inline auto kineticCSSI(Model& model) {
            return kineticCS(model, model.CSpi2AvSI);
        }

        template<class Model>
        static inline auto gradientCSSI(Model& model) {
            return gradientCS(model, model.CSgrad2AvSI);
        }
		
        template<class Model>
        static inline auto kineticSU2Doublet(Model& model) {
            return kineticSU2Doublet(model, model.SU2DblPi2AvI);
        }

        template<class Model>
        static inline auto gradientSU2Doublet(Model& model) {
            return gradientSU2Doublet(model, model.SU2DblGrad2AvI);
        }
		
		template<class Model>
        static inline auto kineticSU2DoubletSI(Model& model) {
            return kineticSU2Doublet(model, model.SU2DblPi2AvSI);
        }

        template<class Model>
        static inline auto gradientSU2DoubletSI(Model& model) {
            return gradientSU2Doublet(model, model.SU2DblGrad2AvSI);
        }
		
        template<class Model>
        static inline auto electricU1(Model& model) {
            return electricU1(model, model.U1pi2AvI);
        }

        template<class Model>
        static inline auto magneticU1(Model& model) {
            return magneticU1(model, model.U1Mag2AvI);
        }
        
        template<class Model>
        static inline auto electricU1SI(Model& model) {
            return electricU1(model, model.U1pi2AvSI);
        }

        template<class Model>
        static inline auto magneticU1SI(Model& model) {
            return magneticU1(model, model.U1Mag2AvSI);
        }

        template<class Model>
        static inline auto electricSU2(Model& model)
        {
            return electricSU2(model, model.SU2pi2AvI);
        }

        template<class Model>
        static inline auto magneticSU2(Model& model)
        {
            return   magneticSU2(model, model.SU2Mag2AvI);
        }
		
		template<class Model>
        static inline auto electricSU2SI(Model& model)
        {
            return electricSU2(model, model.SU2pi2AvSI);
        }

        template<class Model>
        static inline auto magneticSU2SI(Model& model)
        {
            return   magneticSU2(model, model.SU2Mag2AvSI);
        }
		
        template <class Model>
        static inline auto rho(Model& model)  // Total energy density (sum of all contributions)
        {
            auto Eks = (model.Ns > 0 ?  kineticS(model) : 0);
            auto Eks2 = (model.Ns > 0 ?  kineticS2(model) : 0); //ADDED BY TOMÁS
            auto Ekcs =  (model.NCs > 0 ? kineticCS(model) : 0);
            auto EkSU2Dbl =  (model.NSU2Doublet > 0 ? kineticSU2Doublet(model) : 0);
            auto Egs =  (model.Ns > 0 ? gradientS(model) : 0 );
            auto Els =  (model.Ns > 0 ? laplacian2S(model) : 0 ); //ADDED BY TOMÁS
            auto Egcs =  (model.NCs > 0 ? gradientCS(model) : 0 );
            auto EgSU2Dbl =  (model.NSU2Doublet > 0 ? gradientSU2Doublet(model) : 0 );
            auto EelU1 =  (model.NU1 > 0 ? electricU1(model) : 0);
            auto EmagU1 =  (model.NU1 > 0 ? magneticU1(model) : 0);
            auto EelSU2 =  (model.NSU2 > 0 ? electricSU2(model) : 0);
            auto EmagSU2 =  (model.NSU2 > 0 ? magneticSU2(model) : 0);

            return  (Eks + Eks2 + Ekcs + EkSU2Dbl
                    + Egs + Els + Egcs + EgSU2Dbl
                    + EelU1 + EmagU1
                    + EelSU2 + EmagSU2
                    + model.potAvI);   // we also include the potential energy to the sum
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
#include "CosmoInterface/definitions/energies_test.h"
#endif


#endif
