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

#include "HeterogeneousShrinkingCoreChemistry.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "thermodynamicConstants.H"
#include "HeterogeneousShrinkingCoreReaction.H"
#include "HSCRReactionRateControlled.H"
#include "HSCRAbadReactionRateAndDiffusionControlled.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(HeterogeneousShrinkingCoreChemistry, 0);

    addToRunTimeSelectionTable
    (
        fvModel,
        HeterogeneousShrinkingCoreChemistry,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::HeterogeneousShrinkingCoreChemistry
::HeterogeneousShrinkingCoreChemistry
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict),
    gasPhaseName_(dict.lookup<word>("gasPhaseName")),
    solidPhaseName_(dict.lookup<word>("solidPhaseName")),
    inertsName_(dict.lookup<word>("inertsName")),
    reducedName_(dict.lookup<word>("reducedName")),
    oxidizedName_(dict.lookup<word>("oxidizedName")),
    Ro_(dict.lookup<scalar>("Ro")),
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
    conversion_
    (
       IOobject
       (
           "conversion",
           mesh.time().timeName(),
           mesh,
           IOobject::NO_READ,
           IOobject::AUTO_WRITE
       ),
       mesh,
       dimensionedScalar(dimless, 0)
    ),
    rhoo_(dimDensity, 0)
{
    const speciesTable& solidSpecies = solidThermo_.species();
    inertsIndex_   = -1;
    oxidizedIndex_ = -1;
    reducedIndex_  = -1;
    forAll(solidSpecies, i) {
        if (solidSpecies[i] == inertsName_) inertsIndex_ = i;
        else if (solidSpecies[i] == oxidizedName_) oxidizedIndex_ = i;
        else if (solidSpecies[i] == reducedName_) reducedIndex_ = i;
    }

    const dimensionedScalar& Wo = solidThermo_.Wi(oxidizedIndex_);
    const dimensionedScalar& Wr = solidThermo_.Wi(reducedIndex_);

    // Assuming density does not change with T and p
    const scalar rhoi = solidThermo_.rhoi(inertsIndex_,   101325, 300);
    const scalar rhoo = solidThermo_.rhoi(oxidizedIndex_, 101325, 300);
    const scalar rhor = solidThermo_.rhoi(reducedIndex_,  101325, 300);

    scalar Romax = (rhoo - rhor)/rhoo;
    Info << tab << "Oxygen carrying capacity = " << Ro_ 
         << " (max = " << Romax << ")" << endl;
    if (Ro_ > Romax) {
        FatalErrorInFunction 
            << "Oxygen carrying capacity is larger than the maximum: "
           << Ro_ << " > " << Romax << abort(FatalError);
    }

    if (rhoi != rhoo) {
        FatalErrorInFunction << "According to the model assumptions,"
            << " density of inerts ("  << rhoi 
            << ") should be equal to the max density of a particle (" 
            << rhoo << ")" << abort(FatalError);
    }
    if (rhoo <= rhor){
        FatalErrorInFunction << "Density of the reduced reactant must be "
            << "smaller than the density of the oxidized one"
            << abort(FatalError);
    }
    rhoo_ = dimensionedScalar(dimDensity, rhoo);

    const dimensionedScalar vi = 1 - Ro_/Romax;
    Info << tab << "Volume fraction of reacting part:       " 
         << vi.value() << endl;
    Info << tab << "Volume fraction of inerts in unreacted: " 
         << 1 - vi.value() << endl;

    const dimensionedScalar rhomo = (1-vi)*rhoo/Wo;
    const dimensionedScalar rhomr = (1-vi)*rhor/Wr;
    Info << tab << "Molar density of oxidized: " 
         << rhomo.value() << " mol/m^3" << endl;
    Info << tab << "Molar density of reduced:  "
         << rhomr.value() << " mol/m^3" << endl;


    // Create reactions
    const dictionary& reactionsDict(dict.subDict("reactions"));

    label count = 0;
    forAllConstIter(dictionary, reactionsDict, iter)
    {
        if (iter().isDict())
        {
            count++;
        }
    }

    label i = 0;
    forAllConstIter(dictionary, reactionsDict, iter)
    {
        // reactions_.setSize(count);
        if (iter().isDict() 
            && iter().dict().lookupOrDefault<bool>("active", true))
        {
            reactions_.setSize(i+1);
            downscalingFactors_.setSize(i+1);

            const word& name = iter().keyword();
            const dictionary& reactionDict = iter().dict();

            const word   reactionType(reactionDict.lookup("type"));

            Info << name << endl;
            if (reactionType == "reactionRateControlled") {
                reactions_.set(
                    i,
                    new HSCRReactionRateControlled(
                        name,
                        reactionType,
                        mesh,
                        reactionDict,
                        gasPhaseName_,
                        solidPhaseName_,
                        inertsName_
                    )
                );
            } else if (reactionType == "reactionRateAndDiffusionControlled") {
                reactions_.set(
                    i,
                    new HSCRAbadReactionRateAndDiffusionControlled(
                        name,
                        reactionType,
                        mesh,
                        reactionDict,
                        gasPhaseName_,
                        solidPhaseName_,
                        inertsName_
                    )
                );
            } else {
                FatalErrorInFunction << "Unknown reaction type "
                    << reactionType << abort(FatalError);
            }
            

            downscalingFactors_.set(
                i,
                new volScalarField::Internal
                (
                    IOobject
                    (
                        "downscalingFactors_" + name,
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensionedScalar(dimless, 1.0)
                )
            );

            i++;
        }
    }
    reactions_.setSize(i);
    
    // Check if relation between ox and red solids is always the same
    // and calculate density ratio between ox and red
    // This is checked once more in constructor of every reaction
    rho_ratio_ = rhoo/rhor;
    forAll(reactions_, i) {
        scalar MW_SR = 
            solidThermo_.Wi(reactions_[i].getSolidReactantCoeffs().index).value();
        scalar MW_SP = 
            solidThermo_.Wi(reactions_[i].getSolidProductCoeffs().index).value();
        scalar stoich_SR = reactions_[i].getSolidReactantCoeffs().stoichCoeff;
        scalar stoich_SP = reactions_[i].getSolidProductCoeffs().stoichCoeff;
    
        scalar rho_ratioi = reactions_[i].isReduction() ? 
            (stoich_SR*MW_SR) / (stoich_SP*MW_SP)
            :
            (stoich_SP*MW_SP) / (stoich_SR*MW_SR);
        if (mag(rho_ratioi - rho_ratio_) > 1e-6) {
            FatalErrorInFunction << "Relations between solid species' "
                << "densities are not consistent between different reactions. "
                << "Physical properties give " << rho_ratio_ << " whereas "
                << " molar weights and stoichiometric coefficients of "
                << "reaction " << reactions_[i].name() << " give " 
                << rho_ratioi
                << abort(FatalError);
        }
    }

    label nfields = 0;
    fieldNames_.setSize(0);
    forAll(reactions_, i) {
        const wordList fieldNamesReactioni_ = reactions_[i].getFieldNames();
        forAll(fieldNamesReactioni_, j) {
            bool found = false;
            forAll(fieldNames_, k) {
                if (fieldNamesReactioni_[j] == fieldNames_[k]) {
                    found = true;
                    break;
                }
            }
            if (! found) {
                fieldNames_.setSize(nfields+1);
                fieldNames_[nfields] = fieldNamesReactioni_[j];
                nfields++;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::HeterogeneousShrinkingCoreChemistry::addSupFields() const
{
    return fieldNames_;
}


void Foam::fv::HeterogeneousShrinkingCoreChemistry::correct()
{
    const volScalarField& Yox  = solidThermo_.Y(oxidizedIndex_);
    const volScalarField& Yred = solidThermo_.Y(reducedIndex_);
    conversion_ = 
        min(
            max(
                Yox/(Yox+rho_ratio_*Yred),
                dimensionedScalar(dimless,0)
            ), 
            dimensionedScalar(dimless,1)
        );

    forAll(reactions_, i) {
        reactions_[i].correct(conversion_);
    }

    // Enforce correction for reactions to ensure mass does not go negative
    // PtrList<volScalarField>& downscalingFactors; // make class attribute
    const dimensionedScalar dt(dimTime, conversion_.mesh().time().deltaTValue());
    
    forAll(downscalingFactors_, ri) {
        forAll(downscalingFactors_[ri], celli) {
            downscalingFactors_[ri][celli] = 1;
        }
    }
    forAll(fieldNames_, fieldi) {
        // Quite ugly solution to get proper dimensions
        volScalarField::Internal totalSourceOfFieldI(
            reactions_[0].fieldSource(fieldNames_[fieldi])
        );

        forAll(reactions_, ri) {
            if (ri>0) {
                //  generation of species also affects total source
                // (but not int his mechanism)
                const wordList& reactantSpeciesNames = 
                    reactions_[ri].getReactantSpeciesNames();
                label rfi = findIndex(reactantSpeciesNames, fieldNames_[fieldi]);
                if (rfi != -1) {
                    totalSourceOfFieldI += 
                        reactions_[ri].fieldSource(
                            fieldNames_[fieldi]).internalField();
                }
            }
        }
        
        volScalarField maxSourceOfFieldI =
            mesh().lookupObject<volScalarField>(fieldNames_[fieldi])
          * (
                fieldNames_[fieldi].find(gasPhaseName_) != string::npos ? 
                gasThermo_.rho() : solidThermo_.rho()
            )
          * (
                fieldNames_[fieldi].find(gasPhaseName_) != string::npos ? 
                mesh().lookupObject<volScalarField>
                (
                    IOobject::groupName("alpha", gasPhaseName_)
                )
              : mesh().lookupObject<volScalarField>
                (
                IOobject::groupName("alpha", solidPhaseName_)
                )
            )
          / dt;

        forAll(downscalingFactors_, ri) {
            const label reactantIndex = findIndex(
                reactions_[ri].getReactantSpeciesNames(), fieldNames_[fieldi]
            );
            if (reactantIndex != -1) {
                forAll(downscalingFactors_[ri], celli) {
                    scalar newFactor = min
                    (
                        maxSourceOfFieldI[celli] 
                      / max(-totalSourceOfFieldI[celli],SMALL),
                        1
                    );
                    if (newFactor < downscalingFactors_[ri][celli]) {
                        downscalingFactors_[ri][celli] = newFactor;
                    }
                }
            }
        }
        
    }

    forAll(reactions_, i) {
        volScalarField minDownscalingFactor
        (
            IOobject
            (
                "minDownscalingFactor" + i,
                mesh()
            ),
            mesh(),
            dimensionedScalar(dimless, 1)
        );
        forAll(reactions_[i].getReactantSpeciesNames(), rfieldi) {
            forAll(minDownscalingFactor, celli) {
                minDownscalingFactor[celli] = min
                (
                    minDownscalingFactor[celli], 
                    downscalingFactors_[i][celli]
                );
            }
        }
        reactions_[i].downscale(minDownscalingFactor);
    }
}


void Foam::fv::HeterogeneousShrinkingCoreChemistry::addSup
(
    const volScalarField& rho,
    const volScalarField& field,
    fvMatrix<scalar>& eqn
) const
{
    forAll(reactions_, i) {
        reactions_[i].addSup(rho, field, eqn);
    }
}


bool Foam::fv::HeterogeneousShrinkingCoreChemistry::movePoints()
{
    return true;
}


void Foam::fv::HeterogeneousShrinkingCoreChemistry::topoChange(const polyTopoChangeMap&)
{}


void Foam::fv::HeterogeneousShrinkingCoreChemistry::mapMesh(const polyMeshMap&)
{}


void Foam::fv::HeterogeneousShrinkingCoreChemistry::distribute(const polyDistributionMap&)
{}


bool Foam::fv::HeterogeneousShrinkingCoreChemistry::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
