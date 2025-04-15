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

Description
    Some useful functions to set up a case with ilmenite oxygen carrier.


\*---------------------------------------------------------------------------*/

#ifndef flame_cpp_funcs_H
#define flame_cpp_funcs_H

namespace Foam
{

void _calcIlmeniteVolumeFractions(
    const scalar& rhoFe2O3, 
    const scalar& rhoFeO, 
    const scalar& rhoinerts,
    const scalar& Ro,
    const scalar& conversion,
    scalar& vFe2O3,
    scalar& vFeO,
    scalar& vi
) {
    // Calculation done under the assumption that the 
    // particle size does not change (see integral analysis)
    scalar Ro_max = 1 - rhoFeO / rhoFe2O3;
    vi = 1.0 - Ro / Ro_max;
    vFe2O3 = conversion*(1-vi);
    vFeO = 1-vi-vFe2O3;
    if (fabs(vFeO) < SMALL) {
        // Info << "[calcFunctions.H:_calcIlmeniteVolumeFractions] Setting vFeO from " << vFeO << " to 0" << endl;
        vi += vFeO;
        vFeO = 0;
    }
}

scalar calcIlmeniteDensity(
    const scalar& rhoFe2O3, 
    const scalar& rhoFeO, 
    const scalar& rhoinerts,
    const scalar& Ro,
    const scalar& conversion
) {
    scalar vFe2O3, vFeO, vi;
    _calcIlmeniteVolumeFractions(rhoFe2O3, rhoFeO, rhoinerts, Ro, conversion, vFe2O3, vFeO, vi);
    const scalar rhos = vFe2O3*rhoFe2O3 + vFeO*rhoFeO + vi*rhoinerts;
    return rhos;
}

List<scalar> calcIlmeniteMassFractions(
    const wordList& species,
    const scalar& rhoFe2O3, 
    const scalar& rhoFeO, 
    const scalar& rhoinerts,
    const scalar& Ro,
    const scalar& conversion
) {
    scalar vFe2O3, vFeO, vi;
    _calcIlmeniteVolumeFractions(rhoFe2O3, rhoFeO, rhoinerts, Ro, conversion, vFe2O3, vFeO, vi);
    const scalar rhos = calcIlmeniteDensity(rhoFe2O3, rhoFeO, rhoinerts, Ro, conversion);

    List<scalar> result(species.size(), 0.0);
    result[findIndex(species, "Fe2O3")]  = vFe2O3*rhoFe2O3/rhos;
    result[findIndex(species, "FeO")]    = vFeO*rhoFeO/rhos;
    result[findIndex(species, "inerts")] = vi*rhoinerts/rhos;

    scalar sum = 0;
    forAll(result, i) {
        sum+=result[i];
        if (result[i] < 0) {
            FatalErrorInFunction
                << "Some of the mass fractions of ilmenite are negative: " << result[i]
                << exit(FatalError);
        }
    }
    if (abs(1.0 - sum) > SMALL) {
        FatalErrorInFunction
                << "Sum of specie's mass fractions for ilmenite "
                << "is not 1.0" << exit(FatalError);
    }

    return result;
}

List<scalar> dictToMassFractionsList(
    const wordList& species, 
    const dictionary& composition
) {
    // const word compositionName = "composition_H2O";
    // const wordList species = dict.lookup<wordList>("species_gas");
    // const dictionary& composition = dict.subDict(compositionName);
    scalar Ysum = 0;
    const wordList& speciesInMassFractionList(composition.toc());
    List<scalar> result(species.size(), 0.0);
    forAll(speciesInMassFractionList, i) {
        const word& specie_name = speciesInMassFractionList[i];
        const scalar& massFraction = composition.lookup<scalar>(composition.toc()[i]);
        const label speciei = findIndex(species, specie_name);
        if (speciei == -1) {
            FatalErrorInFunction
                << "Specie " << specie_name 
                << " not found in species list "
                << exit(FatalError);
        }
        result[speciei] = massFraction;
        Ysum += massFraction;
    }
    
    // createWithValues<List<scalar>>(species.size(), 0, indices, setValues)

    if (abs(1.0 - Ysum) > SMALL) { // TODO: remove hardcoding
        FatalErrorInFunction
                << "Sum of specie's mass fractions for  " << composition.name()
                << "is not 1.0" << exit(FatalError);
    }
    
    return result;
}

List<scalar> moleFractionsToMassFractions(
    const List<scalar>& moleFractions,
    const List<scalar>& molarWeights
) {
    scalar Wavg = 0;
    forAll(moleFractions, i) {
        Wavg += moleFractions[i]*molarWeights[i];
    }
    List<scalar> massFractions(moleFractions.size(), 0.0);
    forAll(massFractions, i) {
        massFractions[i] = molarWeights[i]/Wavg*moleFractions[i];
    }
    return massFractions;
}

// scalar LnMin_to_m3s(
//     const List<scalar>& massFractions,
//     const scalar& Vn
// ) {

// }

// scalar getSpecieValue(
//     const wordList& species,
//     const scalarList& values,
//     const word& specie
// ) {
//     return species[findIndex(species, specie)];
// }



} // end of namespace Foam

#endif
