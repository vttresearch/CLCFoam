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

#include "unityLewisFourierNoPatchMassDiffusion.H"
#include "fvmLaplacian.H"
#include "fvcSnGrad.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarThermophysicalTransportModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicThermophysicalTransportModel>
unityLewisFourierNoPatchMassDiffusion<BasicThermophysicalTransportModel>::unityLewisFourierNoPatchMassDiffusion
(
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo
)
:
    unityLewisFourierNoPatchMassDiffusion
    (
        typeName,
        momentumTransport,
        thermo
    )
{}


template<class BasicThermophysicalTransportModel>
unityLewisFourierNoPatchMassDiffusion<BasicThermophysicalTransportModel>::unityLewisFourierNoPatchMassDiffusion
(
    const word& type,
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo
)
:
    laminarThermophysicalTransportModel<BasicThermophysicalTransportModel>
    (
        type,
        momentumTransport,
        thermo
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicThermophysicalTransportModel>
const dictionary&
unityLewisFourierNoPatchMassDiffusion<BasicThermophysicalTransportModel>::coeffDict() const
{
    return dictionary::null;
}


template<class BasicThermophysicalTransportModel>
bool unityLewisFourierNoPatchMassDiffusion<BasicThermophysicalTransportModel>::read()
{
    return true;
}


template<class BasicThermophysicalTransportModel>
tmp<surfaceScalarField>
unityLewisFourierNoPatchMassDiffusion<BasicThermophysicalTransportModel>::q() const
{
    // tmp<surfaceScalarField> tqsf = surfaceScalarField::New
    // (
    //     IOobject::groupName
    //     (
    //         "q",
    //         this->momentumTransport().alphaRhoPhi().group()
    //     ),
    //    -fvc::interpolate
    //     (
    //         this->alpha()*this->thermo().kappa()/this->thermo().Cpv()
    //     )
    //    *fvc::snGrad(this->thermo().he())
    // );

    // const fvPatchList& patches = this->thermo().mesh().boundary();
    // surfaceScalarField::Boundary& qsfBf = tqsf.ref().boundaryFieldRef();
    // forAll(patches, patchi)
    // {
    //     if (qsfBf[patchi].coupled())
    //     {
    //         continue;
    //     }
    //     if(qsfBf[patchi].patch().type() != "wall")
    //     {
    //         qsfBf[patchi] *= 0;
    //     }
    // }
    // return tqsf;
    return surfaceScalarField::New
    (
        IOobject::groupName
        (
            "q",
            this->momentumTransport().alphaRhoPhi().group()
        ),
       -fvc::interpolate
        (
            this->alpha()*this->thermo().kappa()/this->thermo().Cpv()
        )
       *fvc::snGrad(this->thermo().he())
    );
}


template<class BasicThermophysicalTransportModel>
tmp<fvScalarMatrix>
unityLewisFourierNoPatchMassDiffusion<BasicThermophysicalTransportModel>::
divq(volScalarField& he) const
{
    // volScalarField alphahe
    // (
    //     volScalarField::New
    //     (
    //         "alphahe",
    //         this->thermo().kappa()/this->thermo().Cpv()
    //     )
    // );

    // volScalarField alpha(this->alpha()*alphahe);

    // // volScalarField::Boundary& kappaBf = kappa.ref().boundaryFieldRef();
    // volScalarField::Boundary& alphaBf = alpha.boundaryFieldRef();
    // const fvPatchList& patches = this->thermo().mesh().boundary();
    // forAll(patches, patchi)
    // {
    //     if (alphaBf[patchi].coupled())
    //     {
    //         continue;
    //     }
    //     if(alphaBf[patchi].patch().type() != "wall")
    //     {
    //         alphaBf[patchi] *= 0;
    //     }
    // }

    // return -fvm::laplacian(alpha, he);
    volScalarField alphahe
    (
        volScalarField::New
        (
            "alphahe",
            this->thermo().kappa()/this->thermo().Cpv()
        )
    );

    return -fvm::laplacian(this->alpha()*alphahe, he);
}


template<class BasicThermophysicalTransportModel>
tmp<surfaceScalarField>unityLewisFourierNoPatchMassDiffusion<BasicThermophysicalTransportModel>::j
(
    const volScalarField& Yi
) const
{
    const fvPatchList& patches = Yi.mesh().boundary();
    tmp<volScalarField> dEff = this->DEff(Yi);
    volScalarField::Boundary& dEffBf = dEff.ref().boundaryFieldRef();
    forAll(patches, patchi)
    {
        if (dEffBf[patchi].coupled())
        {
            continue;
        }
        if(dEffBf[patchi].patch().type() != "wall")
        {
            dEffBf[patchi] *= 0;
        }
    }
    return surfaceScalarField::New
    (
        IOobject::groupName
        (
            "j(" + Yi.name() + ')',
            this->momentumTransport().alphaRhoPhi().group()
        ),
       -fvc::interpolate(this->alpha()*dEff)
       *fvc::snGrad(Yi)
    );
}


template<class BasicThermophysicalTransportModel>
tmp<fvScalarMatrix>
unityLewisFourierNoPatchMassDiffusion<BasicThermophysicalTransportModel>::
divj(volScalarField& Yi) const
{
    const fvPatchList& patches = Yi.mesh().boundary();
    tmp<volScalarField> dEff = this->DEff(Yi);
    volScalarField::Boundary& dEffBf = dEff.ref().boundaryFieldRef();
    forAll(patches, patchi)
    {
        if (dEffBf[patchi].coupled())
        {
            continue;
        }
        if(dEffBf[patchi].patch().type() != "wall")
        {
            dEffBf[patchi] *= 0;
        }
    }
    return -fvm::laplacian(this->alpha()*dEff, Yi);
}


template<class BasicThermophysicalTransportModel>
void unityLewisFourierNoPatchMassDiffusion<BasicThermophysicalTransportModel>::predict()
{
    laminarThermophysicalTransportModel
    <
        BasicThermophysicalTransportModel
    >::predict();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace laminarThermophysicalTransportModels
} // End namespace Foam

// ************************************************************************* //
