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

#include "HSCRAbadReactionRateAndDiffusionControlled.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "thermodynamicConstants.H"
#include "HeterogeneousShrinkingCoreReaction.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

HSCRAbadReactionRateAndDiffusionControlled::HSCRAbadReactionRateAndDiffusionControlled
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
    Ta_(Ea_/Foam::constant::physicoChemical::RR),
    // b_(dict.lookup<scalar>("b")),
    n_(dict.lookup<scalar>("n")),
    De0_(dimLength*dimLength/dimTime,dict.lookup<scalar>("De0")),
    Edif_(dimEnergy/dimMoles,dict.lookup<scalar>("Edif")*1000),
    Tdif_(Edif_/Foam::constant::physicoChemical::RR),
    rp_(dimLength,dict.lookup<scalar>("rp")),
    Xcrit_(dict.lookup<scalar>("Xcrit")),
    Vgtot_(
        dimVolume,
        4./3.*constant::mathematical::pi*rg_.value()*rg_.value()*rg_.value()
    ),
    Vptot_(
        dimVolume,
        4./3.*constant::mathematical::pi*rp_.value()*rp_.value()*rp_.value()
    )
{
    
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void HSCRAbadReactionRateAndDiffusionControlled::calculateNetReactionRate
(
    const volScalarField& rFraction
)
{
    const volScalarField& p = gasThermo_.p(); // Pressure [Pa]
    const volScalarField& Tg = gasThermo_.T(); // Gas phase temperature [K]
    const volScalarField& alphas = p.mesh().lookupObject<volScalarField>
        (
            IOobject::groupName("alpha", solidPhaseName_)
        ); // Solid phase volume fraction [-]
    const volScalarField MW_GMIX(gasThermo_.W()); //- Mean molar weight [kg/mol]

    //- Conversion
    const volScalarField conversion(1-pow(rFraction, 3));
     // TODO utilize saved field

    //- Grain core radius calculation [kg/mol]
    const volScalarField rgc(rg_*rFraction);

    //- Particle core radius calculation [kg/mol]
    const volScalarField rpc(rp_*rFraction);

    // Concentration of gaseous reactant (assuming ideal gas) [mol/m^3]
    const volScalarField cGR(
        gasThermo_.rho()/MW_GR*max(gasThermo_.Y(scGR_.index),
        dimensionedScalar("0", dimless, 0.0))
    );

    // Arrhenius reaction rate
    const volScalarField kR(A_*exp(-Ta_/Tg));

    // Diffusion rate
    const volScalarField kD(De0_*exp(-Tdif_/Tg));

    forAll(NRR_, celli) {
        if (conversion[celli] < Xcrit_) {
            NRR_[celli] = 4*constant::mathematical::pi*kR[celli]
                        * Foam::pow(cGR[celli], n_)*alphas[celli]
                        / Vgtot_.value()* (rgc[celli]*rgc[celli]);
            // dNRRdYSR_[celli] = 4 * b_ * constant::mathematical::pi * kR[celli]
            //     * Foam::pow(gasThermo_.rho()()[celli] / MW_GR.value(), n_)
            //     * Foam::pow(gasThermo_.Y(scGR_.index)[celli], n_ - 1) * n_
            //     * alphas[celli] / Vgtot_.value()
            //     * (rgc[celli] * rgc[celli]);
        } else {
            NRR_[celli] = 4*constant::mathematical::pi*kD[celli]*cGR[celli]
                        * alphas[celli]/Vptot_.value() 
                        * (rpc[celli]*rp_.value()/(rp_.value()-rpc[celli]));
            // dNRRdYSR_[celli] = 4 * b_ * constant::mathematical::pi * kD[celli]
            //     * gasThermo_.rho()()[celli] / MW_GR.value()
            //     * alphas[celli] / Vptot_.value()
            //     * (rpc[celli] * rp_.value() / (rp_.value() - rpc[celli]));
        }
    }

    // return NRR;
}

// ************************************************************************* //
