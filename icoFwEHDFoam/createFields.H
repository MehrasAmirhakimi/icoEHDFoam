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


double p1 = 0.5;
double p2 = 0.87890594;
double p3 = 0.51498869;
double p4 = 0.15084934;
double p5 = 0.2658733e-1;
double p6 = 0.301532e-2;
double p7 = 0.32411e-3;
double q1 = 0.39894228;
double q2 = 0.3988024e-1;
double q3 = -0.362018e-2;
double q4 = 0.163801e-2;
double q5 = -0.1031555e-1;
double q6 = 0.2282967e-1;
double q7 = -0.2895312e-1;
double q8 = 0.1787654e-1;
double q9 = -0.420059e-2;
double yy = 0.0;
double ax = 0.0;
double bessi1 = 0.0;
double twoOmega = 0.0;


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

