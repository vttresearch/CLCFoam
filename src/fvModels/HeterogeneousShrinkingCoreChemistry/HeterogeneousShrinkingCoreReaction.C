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

#include "HeterogeneousShrinkingCoreReaction.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "thermodynamicConstants.H"
// #include "fvmSup.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

HeterogeneousShrinkingCoreReaction::HeterogeneousShrinkingCoreReaction
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
    name_(name),
    type_(dict.lookup<word>("type")),
    gasPhaseName_(gasPhaseName),
    solidPhaseName_(solidPhaseName),
    inertsName_(inertsName),
    reactionString_(dict.lookup<string>("reaction")),
    gasThermo_
    (
        mesh.lookupObject<fluidMulticomponentThermo>
        (
            IOobject::groupName
            (
               physicalProperties::typeName,
               gasPhaseName_
            )
        )
    ),
    solidThermo_
    (
        mesh.lookupObject<fluidMulticomponentThermo>
        (
            IOobject::groupName
            (
                physicalProperties::typeName,
                solidPhaseName_
            )
        )
    ),
    MW_SR(dimMass*inv(dimMoles), 0),
    MW_SP(dimMass*inv(dimMoles), 0),
    MW_GR(dimMass*inv(dimMoles), 0),
    rhou_(dimDensity, 0),
    SGR_
    (
        IOobject
        (
            "SGR_" + name,
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimMass/dimTime/dimVolume, 0)
    ),
    SSR_
    (
        IOobject
        (
            "SSR_" + name,
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimMass/dimTime/dimVolume, 0)
    ),
    SSP_
    (
        IOobject
        (
            "SSP_" + name,
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimMass/dimTime/dimVolume, 0)
    ),
    NRR_
    (
        IOobject
        (
            "NRR_" + name,
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimMoles/dimTime/dimVolume, 0)
    ),
    correctTASIF_(dict.lookupOrDefault<bool>("correctTASIF", false)),
    NRRAlphaCorr_
    (
        IOobject
        (
            "NRRAlphaCorr_" + name,
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 0)
    )
    // dNRRdYSR_
    // (
    //     IOobject
    //     (
    //         "NRR_" + name,
    //         mesh.time().timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh,
    //     dimensionedScalar(dimMoles/dimTime/dimVolume, 0)
    // )
{
    if (correctTASIF_) {
        k1_ = (dict.lookup<scalar>("k1"));
        k2_ = (dict.lookup<scalar>("k2"));
    }

    const speciesTable& gasSpecies = gasThermo_.species();
    const speciesTable& solidSpecies = solidThermo_.species();
    label inertsIndex = -1;
    label isolids = gasSpecies.size();
    speciesTable species(gasSpecies);
    forAll(solidSpecies, i) {
        if (solidSpecies[i] == inertsName_) {
            // Info << tab << "Index of inert species is " << i << endl;
            inertsIndex = i;
        }
        species.append(solidSpecies[i]);
    }

    IStringStream reaction_stream(reactionString_);
    List<specieCoeffs> lhs;
    List<specieCoeffs> rhs;
    specieCoeffs::setLRhs(reaction_stream, species, lhs, rhs);

    if (lhs.size() != 2 || rhs.size() < 1) {
        FatalErrorInFunction << "LHS of the equation must contain one solid and"
            << " one gaseous component, RHS should contain one solid component"
            << " and any number of gaseous" << abort(FatalError);
    }
    // ---=== LHS processing ===---
    if (lhs[0].index < isolids) {
        if (lhs[1].index >= isolids) {
            scGR_ = lhs[0]; // gaseous reactant
            scSR_ = lhs[1]; // solid reactant
        } else {
            FatalErrorInFunction << "No solid reactant found on the LHS of "
                << "reaction " << name_ 
                << abort(FatalError);
        }
    } else {
        if (lhs[1].index < isolids) {
            scGR_ = lhs[1]; // gaseous reactant
            scSR_ = lhs[0]; // solid reactant
        } else {
            FatalErrorInFunction << "No gaseous reactant found on the LHS of "
                << "reaction " << name_ 
                << abort(FatalError);
        }
    }
    scSR_.index -= isolids;

    // ---=== RHS processing ===---
    bool found_solid = false;
    forAll(rhs, i) {
        if (rhs[i].index >= isolids) { // if is solid
            if (!found_solid) {
                found_solid = true;
                scSP_ = rhs[i];
                scSP_.index -= isolids;
            } else {
                FatalErrorInFunction << "RHS can have only one solid in " 
                     << "reaction " << name_ << abort(FatalError);
            }
        } else {
            scGPs_.append(rhs[i]);
        }
    }

    Info << tab << "Gas side: " << scGR_.stoichCoeff  
         << gasSpecies[scGR_.index] << " -> ";
    forAll(scGPs_, i){
        Info << scGPs_[i].stoichCoeff  << gasSpecies[scGPs_[i].index];
        if (i < scGPs_.size()-1) {
            Info << " + ";
        }
    }
    Info << endl;

    Info << tab << "Solid side: " << scSR_.stoichCoeff 
         << solidSpecies[scSR_.index]
         << " -> " << scSP_.stoichCoeff  << solidSpecies[scSP_.index] << endl;

    SGPs_.setSize(scGPs_.size());
    forAll(scGPs_, i){
        SGPs_.set
        (
            i,
            new volScalarField(
                IOobject
                (
                    "SGP_" + gasSpecies[scGPs_[i].index] + "_" + name,
                    mesh.time().name(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar(dimMass/dimTime/dimVolume, 0)
            )
        );
    }

    fieldNames_.setSize(7+scGPs_.size());

    fieldNames_[0] = IOobject::groupName("rho", gasPhaseName_);
    fieldNames_[1] = IOobject::groupName("rho", solidPhaseName_);
    fieldNames_[2] = gasThermo_.he().name();
    fieldNames_[3] = solidThermo_.he().name();

    fieldNames_[4] = IOobject::groupName
    (
        solidSpecies[scSR_.index], 
        solidPhaseName_
    );
    fieldNames_[5] = IOobject::groupName
    (
        solidSpecies[scSP_.index],
        solidPhaseName_
    );

    fieldNames_[6] = IOobject::groupName
    (
        gasSpecies[scGR_.index],
        gasPhaseName_
    );
    forAll(scGPs_, i){
        fieldNames_[7+i] = IOobject::groupName
        (
            gasSpecies[scGPs_[i].index],
            gasPhaseName_
        );
    }

    reactantSpeciesNames_.setSize(2);
    reactantSpeciesNames_[0] = IOobject::groupName
    (
        solidSpecies[scSR_.index],
        solidPhaseName_
    );
    reactantSpeciesNames_[1] = IOobject::groupName
    (
        gasSpecies[scGR_.index],
        gasPhaseName_
    );

    // see line 182 in basicSpecieMixture.H
    const scalar rhoi  = solidThermo_.rhoi(inertsIndex, 101325, 300);
    const scalar rhoSR = solidThermo_.rhoi(scSR_.index, 101325, 300);
    const scalar rhoSP = solidThermo_.rhoi(scSP_.index, 101325, 300);
    rhou_ = dimensionedScalar(dimDensity, max(rhoSR,rhoSP));

    if (rhoi != max(rhoSR,rhoSP)) {
        FatalErrorInFunction << "According to the model assumptions,"
            << " density of inerts ("  << rhoi 
            << ") should be equal to the max density of a particle (" 
            << max(rhoSR,rhoSP) << ")" << abort(FatalError);
    }

    reduction_ = (rhoSR > rhoSP);
    if (reduction_) {
        Info << tab << "Reduction reaction" << endl;
    } else {
        Info << tab << "Oxidation reaction" << endl;
    }
    if (rhoSR == rhoSP){
        FatalErrorInFunction << "Densities of reactants cannot be equal"
            << abort(FatalError);
    }

    //- Molar weights [kg/mol]
    MW_SR = solidThermo_.Wi(scSR_.index);
    MW_SP = solidThermo_.Wi(scSP_.index);
    MW_GR = gasThermo_.Wi(scGR_.index);
    MW_GPs.setSize(scGPs_.size());
    forAll(MW_GPs, i){
        MW_GPs.set
        (
            i,
            new dimensioned<scalar>
            (
                "MW_GP_" + gasSpecies[scGPs_[i].index], // TODO name
                // MW_GR.dimensions(),
                gasThermo_.Wi(scGPs_[i].index)///1000
            )
        );
    }

    if ( 
        mag( rhoSR / rhoSP - (scSR_.stoichCoeff * MW_SR.value() ) 
        / (scSP_.stoichCoeff * MW_SP.value())) > 1e-6
    ){
        FatalErrorInFunction << "According to the model assumptions,"
            << " ratio between densities ("  << rhoSR/rhoSP 
            << ") must be equal to the molar mass ratio (" 
            << (scSR_.stoichCoeff*MW_SP)/(scSR_.stoichCoeff*MW_SP)
            << ")" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void HeterogeneousShrinkingCoreReaction::correct
(
    const volScalarField& conversion
)
{
    volScalarField rFraction(pow(conversion, 1./3.)); // r = rFraction*rg_
    if (! reduction_) {
        rFraction = pow(1-conversion, 1./3.);
    }
    
    calculateNetReactionRate(rFraction);

    // Info << "[" << name_ << "] min=" << gMin(NRRSR) 
    //     << ", max=" << gMax(NRRSR) << endl;
    // Warning: limiting for this particular reaction
    // const dimensionedScalar dt(dimTime, NRR_.mesh().time().deltaTValue());
    // const volScalarField maxNRRSR =
    //     max( (gasThermo_.Y(scGR_.index)*gasThermo_.rho()/dt)
    //        / (scGR_.stoichCoeff/scSR_.stoichCoeff)/MW_GR,
    //        dimensionedScalar(dimMoles/dimVolume/dimTime,0) );
    // // Info << NRRSR << endl;
    // const volScalarField NRRSR_corr = min( NRR_, maxNRRSR );

    if (correctTASIF_) {
        const volScalarField& alphas =
            rFraction.mesh().lookupObject<volScalarField>
         (IOobject::groupName("alpha", solidPhaseName_));

        NRRAlphaCorr_ = 
        (
          1/(
            1
            + neg(1e-4-alphas) * (alphas-1e-4)
            * k1_*max(((alphas-1e-4)/k2_)*(1-(alphas-1e-4)/k2_), scalar(0))
            )
        );

        NRR_ *= NRRAlphaCorr_;
    }


    // Source terms [kg/s/m^3]
    SSR_ = scSR_.stoichCoeff*MW_SR*NRR_;
    SSP_ = scSP_.stoichCoeff*MW_SP*NRR_;
    SGR_ = scGR_.stoichCoeff*MW_GR*NRR_;
    forAll(scGPs_, i){
        SGPs_[i] = scGPs_[i].stoichCoeff*MW_GPs[i]*NRR_;
    }
    
    // SSR_ = scSR_.stoichCoeff*MW_SR*NRRSR_corr;
    // SSP_ = scSP_.stoichCoeff*MW_SP*NRRSR_corr;
    // SGR_ = scGR_.stoichCoeff*MW_GR*NRRSR_corr;
    // forAll(scGPs_, i){
    //     SGPs_[i] = scGPs_[i].stoichCoeff*MW_GPs[i]*NRRSR_corr;
    // }
}

void HeterogeneousShrinkingCoreReaction::downscale
(
    const volScalarField& downscalingFactor
){
    SSR_ *= downscalingFactor;
    SSP_ *= downscalingFactor;
    SGR_ *= downscalingFactor;
    forAll(scGPs_, i){
        SGPs_[i] *= downscalingFactor;
    }
    NRR_ *= downscalingFactor;
}


void HeterogeneousShrinkingCoreReaction::addSup
(
    const volScalarField& rho,
    const volScalarField& field,
    fvMatrix<scalar>& eqn
) const
{
    const speciesTable& gasSpecies = gasThermo_.species();
    const speciesTable& solidSpecies = solidThermo_.species();
    const word fieldName = field.name();
    if (fieldName == IOobject::groupName("rho", gasPhaseName_)) {
        eqn -= SGR_;
        forAll(scGPs_, i){
            eqn += SGPs_[i];
        }
    } else if (fieldName == IOobject::groupName("rho", solidPhaseName_)) {
        eqn += -SSR_ + SSP_;
    } else if (fieldName == gasThermo_.he().name()) {
        // Reaction heat is negative for exothermic reaction.
        eqn += SGR_*gasThermo_.hfi(scGR_.index);
        forAll(scGPs_, i){
            eqn -= SGPs_[i]*gasThermo_.hfi(scGPs_[i].index);
        }
    } else if (fieldName == solidThermo_.he().name()) {
        eqn -= - SSR_*solidThermo_.hfi(scSR_.index)
               + SSP_*solidThermo_.hfi(scSP_.index);
    } else if 
    (
        fieldName 
     == IOobject::groupName(solidSpecies[scSR_.index], solidPhaseName_)
    ) {
        eqn -= SSR_;
    } else if 
    (
        fieldName 
     == IOobject::groupName(solidSpecies[scSP_.index], solidPhaseName_)
    ) {
        eqn += SSP_;
    } else if 
    (
        fieldName 
     == IOobject::groupName(gasSpecies[scGR_.index], gasPhaseName_)
    ) {
        eqn -= SGR_;
        // Trying to linearize...
        // const volScalarField& rhoGR  = gasThermo_.rhoi(scGR_.index, p, Tg);
        // const volScalarField& YGR    = gasThermo_.Y(scGR_.index);
        // eqn -= fvm::Sp(dNRRdYSR_, YGR);
    } else {
        forAll(scGPs_, i){
            if 
            (
                fieldName 
             == IOobject::groupName(gasSpecies[scGPs_[i].index], gasPhaseName_)
            ) {
                eqn += SGPs_[i];
            }
        }
    }
}


const volScalarField HeterogeneousShrinkingCoreReaction::fieldSource
(
    const word& fieldName
) const
{
    const speciesTable& gasSpecies = gasThermo_.species();
    const speciesTable& solidSpecies = solidThermo_.species();
    
    if (fieldName == IOobject::groupName("rho", gasPhaseName_)) {
        volScalarField result(-SGR_);
        forAll(scGPs_, i){
            result += SGPs_[i];
        }
        return result;
    } else if (fieldName == IOobject::groupName("rho", solidPhaseName_)) {
        return -SSR_ + SSP_;
    } else if (fieldName == gasThermo_.he().name()) {
        volScalarField result(SGR_*gasThermo_.hfi(scGR_.index));
        forAll(scGPs_, i){
            result -= SGPs_[i]*gasThermo_.hfi(scGPs_[i].index);
        }
        return result;
    } else if (fieldName == solidThermo_.he().name()) {
        return   SSR_*solidThermo_.hfi(scSR_.index)
               - SSP_*solidThermo_.hfi(scSP_.index);
    } else if 
    (
        fieldName 
     == IOobject::groupName(solidSpecies[scSR_.index], solidPhaseName_)
    ) {
        return -SSR_;
    } else if 
    (
        fieldName 
     == IOobject::groupName(solidSpecies[scSP_.index], solidPhaseName_)
    ) {
        return SSP_;
    } else if
    (
        fieldName == IOobject::groupName(gasSpecies[scGR_.index], gasPhaseName_)
    ) {
        return -SGR_;
    } else {
        forAll(scGPs_, i){
            if 
            (
                fieldName 
             == IOobject::groupName(gasSpecies[scGPs_[i].index], gasPhaseName_)
            ) {
               return SGPs_[i];
            }
        }
    }
    return volScalarField
        (
            IOobject
            (
                "zeroFieldSource",
                SSR_.mesh().time().name(),
                SSR_.mesh()
            ),
            SSR_.mesh(),
            dimensionedScalar(dimMass/dimTime/dimVolume, 0.0)
        );
}

// ************************************************************************* //
