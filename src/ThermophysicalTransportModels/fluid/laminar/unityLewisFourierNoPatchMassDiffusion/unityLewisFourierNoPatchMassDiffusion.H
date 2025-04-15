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
    Foam::laminarThermophysicalTransportModels::unityLewisFourierNoPatchMassDiffusion

Description
    unityLewisFourierNoPatchMassDiffusion's energy gradient heat flux model for laminar flow.
    Specie fluxes are computed assuming a unity turbulent Lewis number.

SourceFiles
    unityLewisFourierNoPatchMassDiffusion.C

\*---------------------------------------------------------------------------*/

#ifndef unityLewisFourierNoPatchMassDiffusion_H
#define unityLewisFourierNoPatchMassDiffusion_H

#include "laminarThermophysicalTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarThermophysicalTransportModels
{

/*---------------------------------------------------------------------------*\
                           Class unityLewisFourierNoPatchMassDiffusion Declaration
\*---------------------------------------------------------------------------*/

template<class BasicThermophysicalTransportModel>
class unityLewisFourierNoPatchMassDiffusion
:
    public laminarThermophysicalTransportModel
    <
        BasicThermophysicalTransportModel
    >
{

public:

    typedef typename BasicThermophysicalTransportModel::alphaField
        alphaField;

    typedef typename BasicThermophysicalTransportModel::momentumTransportModel
        momentumTransportModel;

    typedef typename BasicThermophysicalTransportModel::thermoModel
        thermoModel;


    //- Runtime type information
    TypeName("unityLewisFourierNoPatchMassDiffusion");


    // Constructors

        //- Construct from a momentum transport model and a thermo model
        unityLewisFourierNoPatchMassDiffusion
        (
            const momentumTransportModel& momentumTransport,
            const thermoModel& thermo
        );

        //- Construct from a type name, a momentum transport model and a thermo
        //  model
        unityLewisFourierNoPatchMassDiffusion
        (
            const word& type,
            const momentumTransportModel& momentumTransport,
            const thermoModel& thermo
        );


    //- Destructor
    virtual ~unityLewisFourierNoPatchMassDiffusion()
    {}


    // Member Functions

        //- Const access to the coefficients dictionary
        virtual const dictionary& coeffDict() const;

        //- Read thermophysicalTransport dictionary
        virtual bool read();

        //- Effective mass diffusion coefficient
        //  for a given specie mass-fraction [kg/m/s]
        virtual tmp<volScalarField> DEff(const volScalarField& Yi) const
        {
            return volScalarField::New
            (
                "DEff",
                this->thermo().kappa()/this->thermo().Cp()
            );
        }

        //- Effective mass diffusion coefficient
        //  for a given specie mass-fraction for patch [kg/m/s]
        virtual tmp<scalarField> DEff
        (
            const volScalarField& Yi,
            const label patchi
        ) const
        {
            return
                this->thermo().kappa().boundaryField()[patchi]
               /this->thermo().Cp().boundaryField()[patchi];
        }

        //- Return the heat flux [W/m^2]
        virtual tmp<surfaceScalarField> q() const;

        //- Return the source term for the energy equation
        virtual tmp<fvScalarMatrix> divq(volScalarField& he) const;

        //- Return the specie flux for the given specie mass-fraction [kg/m^2/s]
        virtual tmp<surfaceScalarField> j(const volScalarField& Yi) const;

        //- Return the source term for the given specie mass-fraction equation
        virtual tmp<fvScalarMatrix> divj(volScalarField& Yi) const;

        //- Correct the unityLewisFourierNoPatchMassDiffusion viscosity
        virtual void predict();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace laminarThermophysicalTransportModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "unityLewisFourierNoPatchMassDiffusion.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
