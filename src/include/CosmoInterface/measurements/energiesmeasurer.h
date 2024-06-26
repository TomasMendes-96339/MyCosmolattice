#ifndef COSMOINTERFACE_MEASUREMENTS_ENERGIESMEASURER_H
#define COSMOINTERFACE_MEASUREMENTS_ENERGIESMEASURER_H
 
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */ 
   
// File info: Main contributor(s): Daniel G. Figueroa, Adrien Florio, Francisco Torrenti,  Year: 2020

#include "CosmoInterface/runparameters.h"
#include "CosmoInterface/measurements/meansmeasurer.h"
#include "CosmoInterface/measurements/measurementsIO/filesmanager.h"
#include "TempLat/util/templatvector.h"
#include "TempLat/util/rangeiteration/sum_in_range.h"
#include "TempLat/util/rangeiteration/tagliteral.h"
#include "TempLat/lattice/algebra/helpers/getngrid.h"
#include "CosmoInterface/definitions/energies.h"
#include "CosmoInterface/definitions/hubbleconstraint.h"
#include "CosmoInterface/definitions/gaugederivatives.h"
#include "CosmoInterface/definitions/fieldfunctionals.h"


namespace TempLat {

    /** \brief A class which contains measurements of energies and scale factor.
     *
     **/

    template <typename T>
    class EnergiesMeasurer {
    public:
        /* Put public methods here. These should change very little over time. */
        template <typename Model>
        EnergiesMeasurer(Model& model, FilesManager& filesManager, const RunParameters<T>& par, bool append) :
               amIRoot( model.getToolBox()->amIRoot()),
                expansion(par.expansion),
                fixedBackground(par.fixedBackground), // boolean: if true, expansion is given by fixed background
                Etot0(0),  // Initial total energy
                energies(filesManager, "energies", amIRoot, append, getEnergyHeaders(model)),  // Output file for volume-average energies.
                energyCons(filesManager,  "energy_conservation", amIRoot, append, getEnergyConsHeaders(), fixedBackground)   // Output file for checking energy conservation.
                {
        }

        template <class Model>
        void measure(Model& model, T t, bool saveEtot)
        {
            energies.addAverage(t);  // add to file
            T Etot = 0;  // stores total energy
            T Egrad = 0;  // auxiliary variable, stores grad energy
            T Ekin = 0;  // auxiliary variable, stores kinetic energy
            T Elap = 0; //ADDED BY TOMÁS 
            T Ekin2 = 0; //ADDED BY TOMÁS 
            //----------implementing pressure 
            T Ptot = 0;
            T Pkin1 = 0;
            T Pkin12 = 0;
            T Pgrad = 0;
            //T Pfld41 = 0;
            T Pfld42 = 0;
            T Plap = 0;
            T Pkin2 = 0;
            T Plap2 = 0;
            T Pfld2 = 0;
            //T Pfld21 = 0;
        
            // The "energies" functions contain the appropriate scale factor rescaling. Here we compute the energy species by species

            // Scalar singlets
            ForLoop(i,0, Model::Ns-1,
                    Ekin = average(Energies::kineticS(model,FieldFunctionals::pi2S(model,i)));
                    Egrad =  average(Energies::gradientS(model,FieldFunctionals::grad2S(model,i)));
                    Ekin2 =  average(Energies::kineticS2(model,FieldFunctionals::pi2S2(model,i)));
                    Elap = average(Energies::laplacian2S(model,FieldFunctionals::lap2S(model,i)));

                    Etot += Ekin + Egrad + Ekin2 + Elap;  // add to total energy
                    energies.addAverage(Ekin);
                    energies.addAverage(Egrad);
                    energies.addAverage(Ekin2); //ADDED BY TOMÁS
                    energies.addAverage(Elap); //ADDED BY TOMÁS

                //------------ Pressure ------------------------ ALL ADDED BY TOMÁS

                    Pkin1 = average(Energies::PkineticS1(model,FieldFunctionals::pi2S(model,i)));
                    Pkin12 = average(Energies::PkineticS2(model,FieldFunctionals::pi2S(model,i)));
                    Pgrad =  average(Energies::PgradientS(model,FieldFunctionals::grad2S(model,i)));
                    //Pfld21 =  average(Energies::Pfld2S1(model,FieldFunctionals::fld2S(model,i)));
                    //Pfld41 =  average(Energies::Pfld4S1(model,FieldFunctionals::fld4S(model,i)));
                    Pfld42 =  average(Energies::Pfld4S2(model,FieldFunctionals::fld4S(model,i)));
                    Plap = average(Energies::Plaplacian2S(model,FieldFunctionals::lap2S(model,i)));
                    Pkin2 =  average(Energies::PfldPiS(model,FieldFunctionals::pi2S2(model,i)));
                    Plap2 =  average(Energies::PfldLaplacianS(model,FieldFunctionals::fldLapS(model,i)));
                    Pfld2 =  average(Energies::Pfld2S(model,FieldFunctionals::fld2S(model,i)));

                    Ptot += Pkin1 + Pkin12 + (-Egrad/3) + /*Pfld41 +*/ Pfld42 + Plap + Pkin2 + Plap2 + Pfld2 /*+ Pfld21*/;

                    energies.addAverage(Pkin1);
                    energies.addAverage(-Egrad/3);
                    //energies.addAverage(Pfld41); //quartic potential 
                    //energies.addAverage(Pfld21); //quadratic potential
                    //energies.addAverage(Plap);
                    //energies.addAverage(Pkin2);
                    //energies.addAverage(Plap2);
                    //energies.addAverage(Pfld2);
                    //energies.addAverage(Pkin12);
                    //energies.addAverage(Pfld42);
                    energies.addAverage(Plap+Pkin2+Plap2+Pfld2+Pkin12+Pfld42); //Interaction

            );
            
            // Complex scalars
            ForLoop(i,0, Model::NCs-1,
                    Ekin = average(Energies::kineticCS(model,FieldFunctionals::pi2CS(model,i)));
                    Egrad =  average(Energies::gradientCS(model,FieldFunctionals::grad2CS(model,i)));
                    Etot += Ekin + Egrad;   // add to total energy
                    energies.addAverage(Ekin);
                    energies.addAverage(Egrad);
            );
            
            // SU2 Doublets
            ForLoop(i,0, Model::NSU2Doublet-1,
                    Ekin = average(Energies::kineticSU2Doublet(model, FieldFunctionals::pi2SU2Doublet(model,i)));
                    Egrad =  average(Energies::gradientSU2Doublet(model,FieldFunctionals::grad2SU2Doublet(model,i)));
                    Etot += Ekin + Egrad;  // add to total energy
                    energies.addAverage(Ekin);
                    energies.addAverage(Egrad);
            );

            // U1 gauge fields
            ForLoop(i,0, Model::NU1-1,
                    Ekin = average(Energies::electricU1(model,FieldFunctionals::pi2U1(model,i)));
                    Egrad =  average(Energies::magneticU1(model,FieldFunctionals::B2U1(model,i)));
                    Etot += Ekin + Egrad;  // add to total energy
                    energies.addAverage(Ekin);
                    energies.addAverage(Egrad);
            );

			// SU2 gauge fields
            ForLoop(i,0, Model::NSU2-1,
                    Ekin = average(Energies::electricSU2(model,FieldFunctionals::pi2SU2(model,i)));
                    Egrad = average(Energies::magneticSU2(model,FieldFunctionals::B2SU2(model,i)));
                    Etot += Ekin + Egrad;  // add to total energy
                    energies.addAverage(Ekin);
                    energies.addAverage(Egrad);
            );
            
            //CHANGED BY TOMÁS 
            // Potential
            T potTerm = 0;
            ForLoop(i, 0, Model::NPotTerms-2, //Chaged here (Only want (minimaly coupled) potential terms)
                    potTerm = average(model.potentialTerms(i));
                    energies.addAverage(potTerm);
                    Etot += potTerm;
                    Ptot -= potTerm;
            );

            energies.addAverage(Etot);
            energies.addAverage(Ptot);
            energies.save();

            if(!fixedBackground) {  // Energy cannot be checked if expansion is fixed
            
                // We now check energy conservation:
                energyCons.addAverage(t);
                if (saveEtot) Etot0 = Etot;  // Saves the initial total energy before the first iteration

                if (expansion) {   // If self-consistent expansion, energy conservation is checked via the first Friedmann equation
                    auto hubbleLaw = HubbleConstraint::get(model);
                    energyCons.addAverage(hubbleLaw[0]);
                    energyCons.addAverage(hubbleLaw[1]);
                    energyCons.addAverage(hubbleLaw[2]);
                } 
                else {  // If no expansion, energy must be approximately constant during the evolution
                    energyCons.addAverage(abs(1.0 - Etot / Etot0));
                }

                energyCons.save();
            }
        }

    private:
    
     	// Returns string with the header of the energies file.
        template <typename Model>
        std::vector<std::string> getEnergyHeaders(Model& model) const 
        {
            std::vector<std::string> ret;
            ret.emplace_back("t");
            ForLoop(i,0, Model::Ns-1,
                    ret.emplace_back("E^kin_scal" + std::to_string(i));
                            ret.emplace_back("E^grad_scal" + std::to_string(i));
                                    ret.emplace_back("E^lap_scal" + std::to_string(i)); //ADDED BY TOMÁS
            );
            ForLoop(i,0, Model::NCs-1,
                    ret.emplace_back("E^kin_cmplxscal" + std::to_string(i));
                            ret.emplace_back("E^grad_cmplxscal" + std::to_string(i));
            );
            ForLoop(i,0, Model::NSU2Doublet-1,
                    ret.emplace_back("E^kin_SU2matter");
                            ret.emplace_back("E^grad_SU2matter");
            );
            ForLoop(i,0, Model::NU1-1,
                    ret.emplace_back("E^kin_U1" + std::to_string(i));
                            ret.emplace_back("E^grad_U1" + std::to_string(i));
            );
            ForLoop(i,0, Model::NSU2-1,
                    ret.emplace_back("E^kin_SU2");
                            ret.emplace_back("E^grad_SU2");
            );
            ForLoop(i,0, Model::NPotTerms-1,
                    ret.emplace_back("Vpot_term_" + std::to_string(i) );
            );

            ret.emplace_back("E_tot");

            return ret;
        }

		// Returns header for energy conservation file.
        std::vector<std::string> getEnergyConsHeaders() const   
        {
            std::vector<std::string> ret;
            ret.emplace_back("t");
            if(expansion){
                ret.emplace_back("rel_diff_friedmann");
                ret.emplace_back("LHS_friedmann");
                ret.emplace_back("RHS_friedmann");
            }else{
                ret.emplace_back("energy_cons");
            }

            return ret;
        }

        /* Put all member variables and private methods here. These may change arbitrarily. */
        const bool amIRoot;
        const bool expansion, fixedBackground;
        T Etot0;

        MeasurementsSaver<T> energies;
        //MeasurementsSaver<T> pressures;
        MeasurementsSaver<T> energyCons;

    };

    class EnergiesMeasurerTester{
    public:
#ifdef TEMPLATTEST
        static inline void Test(TDDAssertion& tdd);
#endif
    };



} /* TempLat */

#ifdef TEMPLATTEST
#include "CosmoInterface/measurements/energiesmeasurer_test.h"
#endif


#endif
