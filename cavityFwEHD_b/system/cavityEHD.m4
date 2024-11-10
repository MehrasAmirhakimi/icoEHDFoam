/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  6
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//Run using:
//m4 -P cavityEHD.m4 > blockMeshDict

//m4 definitions:
m4_changecom(//)m4_changequote([,])
m4_define(calc, [m4_esyscmd(perl -e 'use Math::Trig; printf ($1)')])
m4_define(VCOUNT, 0)
m4_define(vlabel, [[// ]Vertex $1 = VCOUNT m4_define($1, VCOUNT)m4_define([VCOUNT], m4_incr(VCOUNT))])


//*******************************************************************************

//Geometry

m4_define(L, 3.5)

//Number of cells:
m4_define(m, 40)
m4_define(n, 20)


//x grading:
m4_define(xGrading, 1)

//y grading:
m4_define(yGrading, 1)

//*******************************************************************************



convertToMeters 0.01;

vertices
(

(0 0 0) vlabel(B0)
(calc(2*L/5) 0 0) vlabel(B1)
(calc(3*L/5) 0 0) vlabel(B2)
(L 0 0) vlabel(B3)
(0 calc(2*L/5) 0) vlabel(B4)
(calc(2*L/5) calc(2*L/5) 0) vlabel(B5)
(calc(3*L/5) calc(2*L/5) 0) vlabel(B6)
(L calc(2*L/5) 0) vlabel(B7)
(0 L 0) vlabel(B8)
(calc(2*L/5) L 0) vlabel(B9)
(calc(3*L/5) L 0) vlabel(B10)
(L L 0) vlabel(B11)

(0 0 calc(L/10)) vlabel(B12)
(calc(2*L/5) 0 calc(L/10)) vlabel(B13)
(calc(3*L/5) 0 calc(L/10)) vlabel(B14)
(L 0 calc(L/10)) vlabel(B15)
(0 calc(2*L/5) calc(L/10)) vlabel(B16)
(calc(2*L/5) calc(2*L/5) calc(L/10)) vlabel(B17)
(calc(3*L/5) calc(2*L/5) calc(L/10)) vlabel(B18)
(L calc(2*L/5) calc(L/10)) vlabel(B19)
(0 L calc(L/10)) vlabel(B20)
(calc(2*L/5) L calc(L/10)) vlabel(B21)
(calc(3*L/5) L calc(L/10)) vlabel(B22)
(L L calc(L/10)) vlabel(B23)

);


// Defining blocks:
blocks
(

	//Blocks :

        // block0
    hex (B0 B1 B5 B4 B12 B13 B17 B16) (m m 1) simpleGrading (xGrading yGrading 1) 
    hex (B1 B2 B6 B5 B13 B14 B18 B17) (n m 1) simpleGrading (xGrading yGrading 1)
    hex (B2 B3 B7 B6 B14 B15 B19 B18) (m m 1) simpleGrading (xGrading yGrading 1)
    hex (B4 B5 B9 B8 B16 B17 B21 B20) (m calc(m+n) 1) simpleGrading (xGrading yGrading 1)
    hex (B6 B7 B11 B10 B18 B19 B23 B22) (m calc(m+n) 1) simpleGrading (xGrading yGrading 1)
 
);

edges
();

// Defining patches:
boundary
(
    insulatedWall
    {
        type wall;
        faces
        (
            (B0 B12 B16 B4)
            (B4 B16 B20 B8)
            (B8 B20 B21 B9)
            (B10 B22 B23 B11)
            (B11 B23 B19 B7)
            (B7 B19 B15 B3)
        );
    }
    electrodeHV
    {
        type wall;
        faces
        (
            (B5 B9 B21 B17)
            (B5 B17 B18 B6)
            (B10 B6 B18 B22)
        );
    }
    electrodeG
    {
        type wall;
        faces
        (
            (B0 B1 B13 B12)
            (B1 B2 B14 B13)
            (B2 B3 B15 B14)
        );
    }

    frontAndBack
    {
        type empty;
        faces
        (
            (B12 B13 B17 B16)
            (B13 B14 B18 B17)
            (B14 B15 B19 B18)
            (B16 B17 B21 B20)
            (B18 B19 B23 B22)

            (B0 B4 B5 B1)
            (B1 B5 B6 B2)
            (B2 B6 B7 B3)
            (B4 B8 B9 B5)
            (B6 B10 B11 B7)
        );
    }
);

mergePatchPairs 
(
);

// ************************************************************************* //
