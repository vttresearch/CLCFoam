/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

Application
    testThemo

Description
    Small code for manual editing and testing thermo properties

    
    
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"
#include "fluidMulticomponentThermo.H" 
#include "Function1.H"
#include "IOmanip.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;


int main(int argc, char *argv[])
{
    argList::noParallel();
    

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    argList::noParallel();
    argList::validArgs.append("phase");
    argList::validArgs.append("specie");
    argList::addOption
    (
        "Tmin",
        "value",
        "Temperature to start from"
    );
    argList::addOption
    (
        "Tend",
        "value",
        "Temperature to end"
    );
    argList::addOption
    (
        "Tstep",
        "value",
        "Temperature step"
    );
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    scalar Tmin=300, Tend=1500, Tstep=50;
    args.optionReadIfPresent("Tmin", Tmin);
    args.optionReadIfPresent("Tend", Tend);
    args.optionReadIfPresent("Tstep", Tstep);
    scalar p = 100000;

    

    const string phase = args.argRead<string>(1);
    const string specie = args.argRead<string>(2);

    autoPtr<fluidMulticomponentThermo> thermo_(fluidMulticomponentThermo::New(mesh, phase)); 
    fluidMulticomponentThermo& thermo = thermo_();
    label specieID(thermo.species()[specie]);

    scalar MW = thermo.Wi(specieID).value()/1000;

    unsigned int w = IOstream::defaultPrecision() + 7;
    OFstream fileStream(specie + ".dat");
    fileStream
        << setw(w) << "T" 
        << setw(w) << "Cp" 
        << setw(w) << "HE"
        << endl;
    for (scalar T = Tmin; T <= Tend; T+=Tstep) {
        fileStream << setw(w) << T;
        fileStream << setw(w) << thermo.Cpi(specieID,p,T);
        fileStream << setw(w) << thermo.hei(specieID,p,T);
        fileStream << endl;

    }
    Info << "t   = " << thermo.T().time().userTimeValue() << endl;
    Info << "Cp  = " << thermo.Cp()[0] << endl;
    Info << "Cv  = " << thermo.Cv()[0] << endl;
    Info << "rho = " << thermo.rho()()[0] << endl;
    Info << "kappa = " << thermo.kappa()[0] << " W/m/K" << endl;
    Info << "Hf  = " << thermo.hfi(specieID)/1000 << " kJ/kg" << endl;
    Info << "Hf  = " << thermo.hfi(specieID)*MW/1000 << " kJ/mol" << endl;
    return 0;
}


// ************************************************************************* //
