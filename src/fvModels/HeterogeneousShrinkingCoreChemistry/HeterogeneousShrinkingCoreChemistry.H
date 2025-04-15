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
    Foam::fv::HeterogeneousShrinkingCoreChemistry

Description
    Heterogeneous reactions between Eulerian phases using shrinking grain core
    model.
    A few notes:
    1. Assumes all grains in a cell have the same core radius (conversion).
    2. The model starts to behave in a weird way when volume fraction of the
    solid phase is smaller than 1e-5 (see 0d validation). The reaction rate
    is then decreased by alpha/1e-5, which is not an intended behaviour. This is
    happening probably due to floating-point precision issues.

    TODOs/what is not yet here:
    1. Energy equations source terms
    2. More than one gas products/reactants

Usage
    Example usage:
    \verbatim
    HeterogeneousShrinkingCoreChemistryAbad
    {
        type            HeterogeneousShrinkingCoreChemistry;
        libs            ("libHeterogeneousChemistryCLC.so");
        Ro              0.033;      // [-] oxygen carrying capacity
        gasPhaseName    gas;
        solidPhaseName  solids;
        inertsName      inerts;
        
        reactions
        {
            ...
        }
    }
    \endverbatim

SourceFiles
    HeterogeneousShrinkingCoreChemistry.C

\*---------------------------------------------------------------------------*/

#ifndef HeterogeneousShrinkingCoreChemistry_H
#define HeterogeneousShrinkingCoreChemistry_H

#include "fvModel.H"
#include "fluidMulticomponentThermo.H"
#include "specieCoeffs.H"
#include "HeterogeneousShrinkingCoreReaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
           Class HeterogeneousShrinkingCoreChemistry Declaration
\*---------------------------------------------------------------------------*/

class HeterogeneousShrinkingCoreChemistry
:
    public fvModel
{
    // Private Data

    //- Name of the gaseous phase
    const word gasPhaseName_;

    //- Name of the solid phase
    const word solidPhaseName_;
    
    //- Name of the inerts component in solids 
    const word inertsName_;

    //- Name of the reduced component in solids 
    const word reducedName_;

    //- Name of the oxidized component in solids 
    const word oxidizedName_;

    //- Oxygen carrying capacity [-]
    const scalar Ro_;

    //- Reference to the thermo of the gas phase
    const fluidMulticomponentThermo& gasThermo_;

    //- Reference to the thermo of the solid phase
    const fluidMulticomponentThermo& solidThermo_;

    //- Conversion (reaction coordinate). 0 - fully oxidized, 1 - fully reduced
    volScalarField conversion_;

    //- Total density in oxidized state [kg/m^3]
    dimensionedScalar rhoo_;

    //- Index of the inerts component in solids 
    label inertsIndex_;
    
    //- Index of the reduced component in solids 
    label oxidizedIndex_;

    //- Index of the oxidized component in solids 
    label reducedIndex_;

    //- Density ratio between reduced and oxidized state
    scalar rho_ratio_;

    //- Field names to apply source to
    wordList fieldNames_;

    //- Reactions
    PtrList<HeterogeneousShrinkingCoreReaction> reactions_;

    //- Downscaling for reaction rates to ensure mass conservation
    PtrList<volScalarField::Internal> downscalingFactors_;


    void update() const;

public:

    //- Runtime type information
    TypeName("HeterogeneousShrinkingCoreChemistry");


    // Constructors

        //- Construct from components
        HeterogeneousShrinkingCoreChemistry
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


        //- Disallow default bitwise copy construction
        HeterogeneousShrinkingCoreChemistry
        (
            const HeterogeneousShrinkingCoreChemistry&
        ) = delete;


    //- Destructor
    virtual ~HeterogeneousShrinkingCoreChemistry()
    {}


    // Member Functions

        // Checks

            //- Return the list of fields for which the option adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Correct

            //- Update the source terms
            virtual void correct() override;

            virtual void addSup
            (
                const volScalarField& field,
                fvMatrix<scalar>& eqn
            ) const override
            {
                NotImplemented;
            }

            virtual void addSup
            (
                const volScalarField& rho,
                const volScalarField& field,
                fvMatrix<scalar>& eqn
            ) const override;

            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                const volScalarField& field,
                fvMatrix<scalar>& eqn
            ) const override
            {
                this->addSup(rho, field, eqn);
            }

        // Mesh changes

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);

            //- Update for mesh motion
            virtual bool movePoints();

        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const HeterogeneousShrinkingCoreChemistry&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
