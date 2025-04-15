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
    Foam::fv::HSCRAbadReactionRateAndDiffusionControlled

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
    Abad_O2
    {
        active      true;
        reaction    "4 FeO + O2 = 2 Fe2O3";
        type        reactionRateAndDiffusionControlled;

        rg      0.48e-6;    // [m] grain radius
        Ea      11.8e3;     // [J/mol] activation energy
        n       1;          // [-] reaction order
        Edif    77.4e3;     // [J/mol] activation energy 
        rp      160e-6;     // [m] particle radius, used by diffusion controlled 
        // reactions
        Xcrit   0.25;       // [-] Critical conversion at which mode switches 
        // from reaction rate controlled to diffusion controlled 

        A       5.27350965e-5; // [m/s] chemical kinetic constant (=ks0)
        De0     9.03088528e-6; // [m^2/s] pre-exponential factor
        b       4;             // [mol/mol] legacy stoichiometric coeff
    }
    \endverbatim

SourceFiles
    HSCRAbadReactionRateAndDiffusionControlled.C

\*---------------------------------------------------------------------------*/

#ifndef HSCRAbadReactionRateAndDiffusionControlled_H
#define HSCRAbadReactionRateAndDiffusionControlled_H

#include "specieCoeffs.H"
#include "HeterogeneousShrinkingCoreReaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
           Class HSCRAbadReactionRateAndDiffusionControlled Declaration
\*---------------------------------------------------------------------------*/

class HSCRAbadReactionRateAndDiffusionControlled 
    : public HeterogeneousShrinkingCoreReaction
{
    // Private Data

    // //- String enconding the reaction
    // const string reactionString_;

    //- Kinetic parameter: grain radius [m]
    const dimensionedScalar rg_;

    //- Kinetic parameter: chemical kinetic constant [mol^(1-n)m^(3n-2)s^-1]
    const dimensionedScalar A_;

    //- Activation energy [J/mol]
    const dimensionedScalar Ea_;

    //- Activation temperature [K]
    const dimensionedScalar Ta_;

    // //- Kinetic parameter: average stoichiometric coefficient [mol/mol]
    // const scalar b_;

    //- Kinetic parameter: reaction order [-]
    const scalar n_;

    //- Diffusion pre-exponential factor [m^2/s]
    const dimensionedScalar De0_;

    //- Diffusion activation energy [J/mol]
    const dimensionedScalar Edif_;

    //- Diffusion activation temperature [K]
    const dimensionedScalar Tdif_;

    //- Kinetic parameter: particle radius [m]
    const dimensionedScalar rp_;

    //- Critical conversion at which mode switches from reaction rate controlled
    // to diffusion controlled 
    const scalar Xcrit_;

    //- Volume of a single grain [m^3]
    const dimensionedScalar Vgtot_;

    //- Volume of a particle [m^3]
    const dimensionedScalar Vptot_;

public:

    // Constructors

        //- Construct from components
        HSCRAbadReactionRateAndDiffusionControlled
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
    virtual ~HSCRAbadReactionRateAndDiffusionControlled()
    {}


    // Member Functions

        virtual void calculateNetReactionRate(const volScalarField& rFraction);

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const HSCRAbadReactionRateAndDiffusionControlled&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
