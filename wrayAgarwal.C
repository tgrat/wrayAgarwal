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

\*---------------------------------------------------------------------------*/

#include "wrayAgarwal.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
wrayAgarwal<BasicTurbulenceModel>::wrayAgarwal
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),
    C1kw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1kw",
            this->coeffDict_,
            0.0829
        )
    ),
    C1ke_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1ke",
            this->coeffDict_,
            0.1127
        )
    ),
   sigmakw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmakw",
            this->coeffDict_,
            0.72
        )
    ),
    sigmake_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmake",
            this->coeffDict_,
            1.0
        )
    ),
    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            0.41
        )
    ),
    Cw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw",
            this->coeffDict_,
            8.54
        )
    ),
    y_(wallDist::New(this->mesh_).y()),
    R_
    (
        IOobject
        (
            IOobject::groupName("R", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{
   C2kw_ = C1kw_/sqr(kappa_) + sigmakw_;
   C2ke_ = C1ke_/sqr(kappa_) + sigmake_;

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
wrayAgarwal<BasicTurbulenceModel>::fmu()
{
    volScalarField chi3(pow3(R_/this->nu()));

    return chi3/(chi3 + pow3(Cw_));
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
wrayAgarwal<BasicTurbulenceModel>::f1(const volScalarField& S)
{
    volScalarField arg1
    (
        (1.0 + y_*sqrt(R_*S)/this->nu())/
        (1.0 + sqr(max(y_*sqrt(R_*S), 1.5*R_)/(20.0*this->nu())))
    );

    return min(tanh(pow4(arg1)), scalar(0.9));    
    // return tmp<volScalarField>
    // (
    //     new volScalarField("f1", dimensionSet(0, 0, 0, 0, 0), min(tanh(pow4(arg1)), 0.9))
    // );
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
wrayAgarwal<BasicTurbulenceModel>::blend
(
    const volScalarField& f1,
    const dimensionedScalar var1,
    const dimensionedScalar var2
) const
{
    return f1*(var1 - var2) + var2;
}

template<class BasicTurbulenceModel>
bool
wrayAgarwal<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        C1kw_.readIfPresent(this->coeffDict());
        C1ke_.readIfPresent(this->coeffDict());
        sigmakw_.readIfPresent(this->coeffDict());
        sigmake_.readIfPresent(this->coeffDict());
        kappa_.readIfPresent(this->coeffDict());
        Cw_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
wrayAgarwal<BasicTurbulenceModel>::k() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "k",
                this->runTime_.timeName(),
                this->mesh_
            ),
            this->mesh_,
            dimensionedScalar("0", dimensionSet(0, 2, -2, 0, 0), 0)
        )
    );
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
wrayAgarwal<BasicTurbulenceModel>::epsilon() const
{
    WarningInFunction
        << "Turbulence kinetic energy dissipation rate not defined for "
        << "Wray Agarwal model. Returning zero field"
        << endl;

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "epsilon",
                this->runTime_.timeName(),
                this->mesh_
            ),
            this->mesh_,
            dimensionedScalar("0", dimensionSet(0, 2, -3, 0, 0), 0)
        )
    );
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
wrayAgarwal<BasicTurbulenceModel>::DREff(const volScalarField& tf1) const
{
    return blend(tf1, sigmakw_, sigmake_)*R_ + this->nu();
}

template<class BasicTurbulenceModel>
void
wrayAgarwal<BasicTurbulenceModel>::correct()
{
   if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();
   
    volScalarField S(sqrt(2*magSqr(symm(fvc::grad(this->U_)))));
    S = max(S, dimensionedScalar("SMALL", S.dimensions(), 1.0e-16));
   
    const volVectorField gradS(fvc::grad(S));    

    const volScalarField tf1(f1(S));
    const volScalarField C1(blend(tf1, C1kw_, C1ke_));

    tmp<fvScalarMatrix> REqn
    (
        fvm::ddt(alpha, rho, R_)
      + fvm::div(alphaRhoPhi, R_)
      - fvm::laplacian(alpha*rho*DREff(tf1), R_)
     ==
        alpha*rho*C1*R_*S
      + alpha*rho*tf1*C2kw_*R_*(fvc::grad(R_) & gradS)/S
      - fvm::Sp(alpha*rho*(1.0 - tf1)*C2ke_*R_*(gradS & gradS)/sqr(S), R_)
      + fvOptions(alpha, rho, R_)
    );

    REqn.ref().relax();
    fvOptions.constrain(REqn.ref());
    solve(REqn);
    fvOptions.correct(R_);
    bound(R_, dimensionedScalar("0", R_.dimensions(), 0.0));
    R_.correctBoundaryConditions();

    correctNut();
}

template<class BasicTurbulenceModel>
void
wrayAgarwal<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = fmu()*R_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
