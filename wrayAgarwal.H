/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

Class
    Foam::RASModels::wrayAgarwal

Group
    grpRASTurbulence

See also
    Foam::wrayAgarwal

SourceFiles
    wrayAgarwal.C

\*---------------------------------------------------------------------------*/

#ifndef wrayAgarwal_H
#define wrayAgarwal_H

// #include "kOmegaSSTBase.H"
#include "kOmegaSST.H"
// #include "RASModel.H"
// #include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                          Class wrayAgarwal Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class wrayAgarwal
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{
    // Disallow default bitwise copy construct and assignment
        wrayAgarwal(const wrayAgarwal&);
        void operator=(const wrayAgarwal&);

protected:
    dimensionedScalar C1kw_;
    dimensionedScalar C1ke_;

    dimensionedScalar sigmakw_;
    dimensionedScalar sigmake_;

    dimensionedScalar C2kw_;
    dimensionedScalar C2ke_;

    dimensionedScalar kappa_;
    dimensionedScalar Cw_;

    const volScalarField& y_;

    volScalarField R_;

    tmp<volScalarField> fmu();
    tmp<volScalarField> f1(const volScalarField& S);
    tmp<volScalarField> C1(const volScalarField& f1);
    tmp<volScalarField> blend(const volScalarField& f1, const dimensionedScalar var1, const dimensionedScalar var2) const;

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;

    //- Runtime type information
    TypeName("wrayAgarwal");

    // Constructors

        //- Construct from components
        wrayAgarwal
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~wrayAgarwal()
    {}

    virtual bool read();

    virtual tmp<volScalarField> k() const;

    virtual tmp<volScalarField> epsilon() const;

    tmp<volScalarField> DREff(const volScalarField& S) const;
    
    virtual void correctNut();
    
    virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "wrayAgarwal.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
