/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Copyright (C) 2022-2025 VTT Technical Research Centre of Finland Ltd
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::fv::HeterogeneousShrinkingCoreReaction

Description
    Heterogeneous reactions between Eulerian phases using shrinking grain core
    model.
    A few notes:
    1. Assumes all grains in a cell have the same core radius (conversion).
    2. The model starts to behave in a weird way when volume fraction of the
    solid phase is smaller than 1e-5 (see 0d validation). The reaction rate
    is then decreased by alpha/1e-5, which is not an intended behaviour. This is
    happening probably due to floating-point precision issues.

Usage
    Example usage:
    \verbatim
    Abad_pre_CH4
    {
        type            reactionControlled;
        reaction        "4 Fe2O3 + CH4 = 8 FeO + CO2 + 2 H2O";
        
        rg      1.25e-6;    // [m] grain radius

        A       9.8;        // [mol^(1-n)m^(3n-2)s^-1] chemical kinetic constant
        Ea      135.2e3;    // [J/mol] activation energy
        n       1;          // [-] reaction order
    }
    \endverbatim

SourceFiles
    HeterogeneousShrinkingCoreReaction.C

\*---------------------------------------------------------------------------*/

#ifndef HeterogeneousShrinkingCoreReaction_H
#define HeterogeneousShrinkingCoreReaction_H

#include "fvMesh.H"
#include "fluidMulticomponentThermo.H"
#include "specieCoeffs.H"
#include "fvMatrices.H"
#include "constants.H"

#ifndef namespaceFoam
#define namespaceFoam
    using namespace Foam;
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
           Class HeterogeneousShrinkingCoreReaction Declaration
\*---------------------------------------------------------------------------*/

class HeterogeneousShrinkingCoreReaction
{
    // Protected Data
protected:

    //- Name of the reaction
    const string name_;

    //- Type of the reaction
    const word type_;

    //- Name of the gaseous phase
    const word gasPhaseName_;

    //- Name of the solid phase
    const word solidPhaseName_;

    //- Name of the inerts component in solids 
    const word inertsName_;

    //- String enconding the reaction
    const string reactionString_;

    //- Reference to the thermo of the gas phase
    const fluidMulticomponentThermo& gasThermo_;

    //- Reference to the thermo of the solid phase
    const fluidMulticomponentThermo& solidThermo_;

    //- Molar densities [kg/kmol]
    dimensionedScalar MW_SR;
    dimensionedScalar MW_SP;
    dimensionedScalar MW_GR;
    PtrList<dimensionedScalar> MW_GPs;

    //- Total density of unreacted grains [kg/m^3]
    dimensionedScalar rhou_;

    //- Boolean indicating if the reacton is reduction (true) or 
    // oxidation (false)
    bool reduction_;

    // Gaseous reactant mass source terms [kg/m^3/s]
    volScalarField SGR_;
    
    // Solid reactant mass source terms [kg/m^3/s]
    volScalarField SSR_;
    
    // Solid product mass source terms [kg/m^3/s]
    volScalarField SSP_;

    // List of gaseous products mass source terms [kg/m^3/s]
    PtrList<volScalarField> SGPs_;

    // Net reaction rate according to reaction equation [kmol/m^3/s]
    volScalarField NRR_;
    // volScalarField dNRRdYSR_; // derivative of the NRR with respect to Y_SR

    specieCoeffs scGR_;

    specieCoeffs scSR_;

    specieCoeffs scSP_;

    List<specieCoeffs> scGPs_;

    // Switch for TASIF correction
    bool correctTASIF_;

    // Correction factor for net reaction rates (included in NRR_, so to
    // extract raw NRR, divide NRR_ by this correction factor) [-]
    volScalarField NRRAlphaCorr_;

    // Correction coefficients for NRRAlphaCorr_ [-]
    scalar k1_, k2_;

    //- Field names to apply source to
    wordList fieldNames_;
    wordList reactantSpeciesNames_; // Separate list for limiting reaction rates

public:

    // Constructors

        //- Construct from components
        HeterogeneousShrinkingCoreReaction
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict,
            const word gasPhaseName,
            const word solidPhaseName,
            const word inertsName
        );

    //- Destructor
    virtual ~HeterogeneousShrinkingCoreReaction()
    {}


    // Member Functions

        inline const string& name() const{return name_;}

        wordList getFieldNames() const {
            return fieldNames_;
        }

        wordList getReactantSpeciesNames() const {
            return reactantSpeciesNames_;
        }

        inline bool isReduction() const {
            return reduction_;
        }

        inline specieCoeffs getGaseousReactantCoeffs() const {
            return scGR_;
        }
        inline specieCoeffs getSolidReactantCoeffs() const {
            return scSR_;
        }
        inline specieCoeffs getSolidProductCoeffs() const {
            return scSP_;
        }
        inline List<specieCoeffs> getGaseousProductsCoeffsList() const {
            return scGPs_;
        }

        const volScalarField fieldSource(const word& fieldName) const;

        virtual void calculateNetReactionRate
        (
            const volScalarField& rFraction
        ) = 0;

        // Correct

            //- Update the source terms
            void correct(const volScalarField& conversion);

            void downscale(const volScalarField& downscalingFactor);

            void addSup
            (
                const volScalarField& rho,
                const volScalarField& field,
                fvMatrix<scalar>& eqn
            ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const HeterogeneousShrinkingCoreReaction&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
