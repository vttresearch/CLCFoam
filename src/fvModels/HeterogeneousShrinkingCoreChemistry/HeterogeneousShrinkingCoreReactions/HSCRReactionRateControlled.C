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

\*---------------------------------------------------------------------------*/

#include "HSCRReactionRateControlled.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "thermodynamicConstants.H"
#include "HeterogeneousShrinkingCoreReaction.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

HSCRReactionRateControlled::HSCRReactionRateControlled
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const word gasPhaseName,
    const word solidPhaseName,
    const word inertsName
)
:
    HeterogeneousShrinkingCoreReaction(
        name, 
        modelType, 
        mesh, 
        dict, 
        gasPhaseName, 
        solidPhaseName, 
        inertsName
    ),
    rg_(dimLength,dict.lookup<scalar>("rg")),
    A_(dimLength/dimTime,dict.lookup<scalar>("A")),
    Ea_(dimEnergy/dimMoles,dict.lookup<scalar>("Ea")*1000), 
    // TODO: better dimensions handling (also NRR_ = has 1000* now)
    Ta_(Ea_/Foam::constant::physicoChemical::RR),
    // b_(dict.lookup<scalar>("b")),
    n_(dict.lookup<scalar>("n")),
    Vgtot_(dimVolume,4./3.*constant::mathematical::pi*rg_.value()*rg_.value()*rg_.value())
{
    
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void HSCRReactionRateControlled::calculateNetReactionRate(const volScalarField& rFraction)
{
    const volScalarField& p = gasThermo_.p(); // Pressure [Pa]
    // Note that the gas temperature is used.
    const volScalarField& Tg = gasThermo_.T(); // Gas phase temperature [K]
    const volScalarField& alphas = p.mesh().lookupObject<volScalarField>
        (
            IOobject::groupName("alpha", solidPhaseName_)
        ); // Solid phase volume fraction [-]
    const volScalarField MW_GMIX(gasThermo_.W()); //- Mean molar weight [kg/kmol]

    //- Grain radius calculation [m]
    const volScalarField r(rg_*rFraction);

    // Concentration of gaseous reactant (assuming ideal gas) [kmol/m^3]
    const volScalarField cGR
    (
        gasThermo_.rho()/MW_GR*max(gasThermo_.Y(scGR_.index),
        dimensionedScalar("0", dimless, 0.0))
    );

    // Arrhenius reaction rate
    const volScalarField k(A_*exp(-Ta_/Tg));
    NRR_ = (4*constant::mathematical::pi*r*r)*k*pow(cGR*1000, n_)
         / 1000*pow(dimensionedScalar(dimMoles/dimVolume, 1),1.-n_)
         * alphas/Vgtot_;
    
    // dNRRdYSR_ = 
    //     (
    //         b_ * (4 * constant::mathematical::pi * r * r) * k 
    //         * pow(gasThermo_.rho() / MW_GR, n_) 
    //         * pow
    //         (
    //             max
    //             (
    //                 gasThermo_.Y(scGR_.index), 
    //                 dimensionedScalar(dimless, SMALL)
    //             ), 
    //             n_ - 1
    //         ) 
    //         * n_ 
    //         * pow(dimensionedScalar(dimMoles / dimVolume, 1), 1. - n_) 
    //         * alphas / Vgtot_
    //     );
}

// ************************************************************************* //