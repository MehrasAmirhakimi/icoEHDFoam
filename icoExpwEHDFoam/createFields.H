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

dimensionedScalar nu
(
    "nu",
    dimViscosity,
    transportProperties.lookup("nu")
);
/*
dimensionedScalar epsilonr
(
    transportProperties.lookup("epsilonr")
);
*/
dimensionedScalar epsilon
(
    transportProperties.lookup("epsilon")
);


dimensionedScalar bpos
(
    transportProperties.lookup("bpos")
);

dimensionedScalar bneg
(
    transportProperties.lookup("bneg")
);
/*
dimensionedScalar sig
(
    transportProperties.lookup("sig")
);
*/
dimensionedScalar neq
(
    transportProperties.lookup("neq")
);

dimensionedScalar D
(
    transportProperties.lookup("D")
);
dimensionedScalar rho0
(
    transportProperties.lookup("rho0")
);

dimensionedScalar kB
(
    transportProperties.lookup("kB")
);
dimensionedScalar electron
(
    transportProperties.lookup("electron")
);

dimensionedScalar T0
(
    transportProperties.lookup("T0")
);

double pinum = constant::mathematical::pi;
// dimensionedScalar neq = sig / (bpos - bneg);
dimensionedScalar kr = (bpos - bneg) / epsilon;
dimensionedScalar kdN = kr * neq * neq;


Info << "Reading field p\n" << endl;
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


Info << "Reading field U\n" << endl;
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


Info << "Reading field phii\n" << endl;
volScalarField phii
(
    IOobject
    (
        "phii",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info << "Reading field pos\n" << endl;
volScalarField pos
(
    IOobject
    (
        "pos",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading field neg\n" << endl;
volScalarField neg
(
    IOobject
    (
        "neg",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);




Info << "Calculating field chargeDensity\n" << endl;
volScalarField chargeDensity
(
    IOobject
    (
        "chargeDensity",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    pos - neg
);

Info << "Calculating field source\n" << endl;
volScalarField source
(
    IOobject
    (
        "source",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    kdN - kr * pos * neg
);

Info << "Calculating field posFlux\n" << endl;
surfaceScalarField posFlux
(
    IOobject
    (
        "posFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    -1 * bpos * mesh.magSf() * fvc::snGrad(phii)
);

Info << "Calculating field negFlux\n" << endl;
surfaceScalarField negFlux
(
    IOobject
    (
        "negFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    -1 * bneg * mesh.magSf() * fvc::snGrad(phii)
);

Info << "Calculating field force\n" << endl;
volVectorField force
(
    IOobject
    (
        "force",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -chargeDensity * fvc::grad(phii) / rho0
);
Info << "Calculating field elec\n" << endl;
volVectorField elec
(
    IOobject
    (
        "elec",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(phii)
);

Info << "Calculating field omega\n" << endl;
volScalarField omega
(
    IOobject
    (
        "omega",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    sqrt(pow(electron, 3) / (4 * pinum * epsilon * pow(kB * T0, 2)) * mag(elec))
);
Info << "Calculating field fw\n" << endl;
volScalarField fw
(
    IOobject
    (
        "fw",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("fw", dimless, scalar(1))
    //       (2*sqrt(pow(electron,3)*mag(elec)/epsilon0)/kbT)
);



#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

