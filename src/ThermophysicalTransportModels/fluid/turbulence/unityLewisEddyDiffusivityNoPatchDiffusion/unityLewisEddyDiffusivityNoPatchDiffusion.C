/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "unityLewisEddyDiffusivityNoPatchDiffusion.H"
#include "fvcSnGrad.H"
#include "fvmLaplacian.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulenceThermophysicalTransportModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// template<class TurbulenceThermophysicalTransportModel>
// void unityLewisEddyDiffusivityNoPatchDiffusion<TurbulenceThermophysicalTransportModel>::
// correctAlphat()
// {
//     alphat_ =
//         this->momentumTransport().rho()
//        *this->momentumTransport().nut()/Prt_;
//     alphat_.correctBoundaryConditions();
// }


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class TurbulenceThermophysicalTransportModel>
unityLewisEddyDiffusivityNoPatchDiffusion<TurbulenceThermophysicalTransportModel>::
unityLewisEddyDiffusivityNoPatchDiffusion
(
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo,
    const bool allowDefaultPrt
)
:
    unityLewisEddyDiffusivityNoPatchDiffusion
    (
        typeName,
        momentumTransport,
        thermo,
        allowDefaultPrt
    )
{
}


template<class TurbulenceThermophysicalTransportModel>
unityLewisEddyDiffusivityNoPatchDiffusion<TurbulenceThermophysicalTransportModel>::
unityLewisEddyDiffusivityNoPatchDiffusion
(
    const word& type,
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo,
    const bool allowDefaultPrt
)
:   
    unityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>
    (
        type,
        momentumTransport,
        thermo,
        allowDefaultPrt
    )
    // TurbulenceThermophysicalTransportModel
    // (
    //     type,
    //     momentumTransport,
    //     thermo
    // ),

    // Prt_
    // (
    //     allowDefaultPrt
    //   ? dimensioned<scalar>::lookupOrAddToDict
    //     (
    //         "Prt",
    //         this->coeffDict_,
    //         1
    //     )
    //   : dimensioned<scalar>
    //     (
    //         "Prt",
    //         dimless,
    //         this->coeffDict_
    //     )
    // ),

    // alphat_
    // (
    //     IOobject
    //     (
    //         IOobject::groupName
    //         (
    //             "alphat",
    //             this->momentumTransport().alphaRhoPhi().group()
    //         ),
    //         momentumTransport.time().name(),
    //         momentumTransport.mesh(),
    //         IOobject::MUST_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     momentumTransport.mesh()
    // )
{
    Info << "CONSTRUCTING!!!!!!!!!!!!!!!!!!" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// template<class TurbulenceThermophysicalTransportModel>
// bool unityLewisEddyDiffusivityNoPatchDiffusion<TurbulenceThermophysicalTransportModel>::read()
// {
//     if (TurbulenceThermophysicalTransportModel::read())
//     {
//         Prt_.readIfPresent(this->coeffDict());

//         return true;
//     }
//     else
//     {
//         return false;
//     }
// }


template<class TurbulenceThermophysicalTransportModel>
tmp<surfaceScalarField>
unityLewisEddyDiffusivityNoPatchDiffusion<TurbulenceThermophysicalTransportModel>::q() const
{
    tmp<surfaceScalarField> tqsf =  surfaceScalarField::New
    (
        IOobject::groupName
        (
            "q",
            this->momentumTransport().alphaRhoPhi().group()
        ),
       -fvc::interpolate(this->alphaEff()*this->alpha())
       *fvc::snGrad(this->thermo().he())
    );
    const fvPatchList& patches = this->thermo().mesh().boundary();
    surfaceScalarField::Boundary& qsfBf = tqsf.ref().boundaryFieldRef();
    forAll(patches, patchi)
    {
        Info << "    q:Patch " << qsfBf[patchi].patch().name() << " of type " << qsfBf[patchi].patch().type() << endl;

        if (qsfBf[patchi].coupled())
        {
            continue;
        }
        if(qsfBf[patchi].patch().type() != "wall")
        {
            Info << "Disabling diffusion for a patch " << qsfBf[patchi].patch().name() << " of type " << qsfBf[patchi].patch().type() << endl;
            qsfBf[patchi] *= 0;
        }
    }
    return tqsf;
}


template<class TurbulenceThermophysicalTransportModel>
tmp<scalarField>
unityLewisEddyDiffusivityNoPatchDiffusion<TurbulenceThermophysicalTransportModel>::q
(
    const label patchi
) const
{
    const fvPatch& patch = this->thermo().mesh().boundary()[patchi];

    if (patch.coupled() || patch.type() == "wall")
    {
        return
          - (
                this->alpha().boundaryField()[patchi]
               *this->alphaEff(patchi)
               *this->thermo().he().boundaryField()[patchi].snGrad()
            );
    }

    return tmp<scalarField>
    (
        new scalarField(patch.size(), scalar(0))
    );
}


template<class TurbulenceThermophysicalTransportModel>
tmp<fvScalarMatrix>
unityLewisEddyDiffusivityNoPatchDiffusion<TurbulenceThermophysicalTransportModel>::divq
(
    volScalarField& he
) const
{
    volScalarField alphaEffhe
    (
        volScalarField::New
        (
            "alphaEffhe",
            this->alpha()*this->alphaEff()
        )
    );
    volScalarField::Boundary& alphaEffheBf = alphaEffhe.boundaryFieldRef();
    const fvPatchList& patches = this->thermo().mesh().boundary();
    forAll(patches, patchi)
    {
        Info << "    divq:Patch " << alphaEffheBf[patchi].patch().name() << " of type " << alphaEffheBf[patchi].patch().type() << endl;
        if (alphaEffheBf[patchi].coupled())
        {
            continue;
        }
        if(alphaEffheBf[patchi].patch().type() != "wall")
        {
            Info << "Disabling diffusion for a patch " << alphaEffheBf[patchi].patch().name() << " of type " << alphaEffheBf[patchi].patch().type() << endl;
            alphaEffheBf[patchi] *= 0;
        }
    }
    return -fvm::laplacian(alphaEffhe, he);
}


template<class TurbulenceThermophysicalTransportModel>
tmp<surfaceScalarField>
unityLewisEddyDiffusivityNoPatchDiffusion<TurbulenceThermophysicalTransportModel>::j
(
    const volScalarField& Yi
) const
{
    const fvPatchList& patches = Yi.mesh().boundary();
    tmp<volScalarField> dEff = this->DEff(Yi);
    volScalarField::Boundary& dEffBf = dEff.ref().boundaryFieldRef();
    forAll(patches, patchi)
    {
        Info << "    j:Patch " << dEffBf[patchi].patch().name() << " of type " << dEffBf[patchi].patch().type() << endl;

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


template<class TurbulenceThermophysicalTransportModel>
tmp<scalarField>
unityLewisEddyDiffusivityNoPatchDiffusion<TurbulenceThermophysicalTransportModel>::j
(
    const volScalarField& Yi,
    const label patchi
) const
{
    const fvPatch& patch = Yi.mesh().boundary()[patchi];

    if (patch.coupled() || patch.type() == "wall")
    {
        return
          - (
                this->alpha().boundaryField()[patchi]
               *this->DEff(Yi, patchi)
               *Yi.boundaryField()[patchi].snGrad()
            );
    }

    return tmp<scalarField>
    (
        new scalarField(patch.size(), scalar(0))
    );
}


template<class TurbulenceThermophysicalTransportModel>
tmp<fvScalarMatrix>
unityLewisEddyDiffusivityNoPatchDiffusion<TurbulenceThermophysicalTransportModel>::divj
(
    volScalarField& Yi
) const
{
    const fvPatchList& patches = Yi.mesh().boundary();
    tmp<volScalarField> dEff = this->DEff(Yi);
    volScalarField::Boundary& dEffBf = dEff.ref().boundaryFieldRef();
    forAll(patches, patchi)
    {
        Info << "    divj:Patch " << patches[patchi].name() << " of type " << patches[patchi].type() << endl;
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


// template<class TurbulenceThermophysicalTransportModel>
// void unityLewisEddyDiffusivityNoPatchDiffusion<TurbulenceThermophysicalTransportModel>::
// predict()
// {
//     TurbulenceThermophysicalTransportModel::predict();
//     correctAlphat();
// }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbulenceThermophysicalTransportModels
} // End namespace Foam

// ************************************************************************* //
