Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);





dimensionedVector g
(
    "g",
    dimAcceleration,
    transportProperties
);

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


#include "createPhi.H"


Info<< "Reading field alpha\n" << endl;

// Read first phase volume fraction



const wordList phases(transportProperties.get<wordList>("phases")); 
PtrList<volScalarField> alphaList(phases.size()-1);
PtrList<dimensionedScalar> rhoList(phases.size());
PtrList<dimensionedScalar> nuList(phases.size());

forAll (phases, i)
{
    rhoList.set(
        i,
        new dimensionedScalar
        (
            "rho"+std::to_string(i+1),
            dimDensity,
            transportProperties
        )
    );
    nuList.set(
        i,
        new dimensionedScalar
        (
            "nu"+std::to_string(i+1),
            dimViscosity,
            transportProperties
        )
    );
    if (i < phases.size() - 1)
        {
            alphaList.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "alpha" + std::to_string(i+1),
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );
        };
};

volScalarField alphaSum = 0 * alphaList[0];
volScalarField rhoMix =  0 * alphaList[0] * rhoList[0];
volScalarField muMix = 0 * alphaList[0] * rhoList[0] * nuList[0];
for (int i = 0; i < alphaList.size(); i++)
{
    alphaSum += alphaList[i];
    rhoMix += alphaList[i] * rhoList[i];
    muMix += alphaList[i] * nuList[i] * rhoList[i];
};
rhoMix += (1 - alphaSum) * rhoList.last();
muMix += (1 - alphaSum) * nuList.last() * rhoList.last();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rhoMix
);

// Create mixture dynamic viscosity
volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    muMix
);


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());
