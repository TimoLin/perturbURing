/*---------------------------------------------------------------------------*\
Date: August 24, 2006
Author: Eugene de Villiers
Source: http://www.cfd-online.com/Forums/openfoam-solving/58905-les-turbulent-pipe-flow.html#post192079
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Application
    perturbUCylinder

Description
    initialise channel velocity with superimposed streamwise streaks.
    To be used to force channelOodles to transition and reach a fully
    developed flow sooner.

    Reads in perturbUDict.

    EdV from paper:
        Schoppa, W. and Hussain, F.
        "Coherent structure dynamics in near wall turbulence",
        Fluid Dynamics Research, Vol 26, pp119-139, 2000.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "Random.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"

    const vectorField centers(mesh.C());

    Info<< "Time = " << runTime.value() << endl;


    IOobject Uheader
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    // Check U exists
    if (Uheader.headerOk())
    {
        Info<< "    Reading U" << endl;
        volVectorField U(Uheader, mesh);
        const scalar Retau = 300;

        IOdictionary transportProperties
        (
            IOobject
            (
                "transportProperties",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        dimensionedScalar nu
        (
            transportProperties.lookup("nu")
        );
        dimensionedVector Ubar
        (
            transportProperties.lookup("Ubar")
        );
        dimensionedScalar diameterOuter
        (
            transportProperties.lookup("dOuter")
        );
        dimensionedScalar diameterInner
        (
            transportProperties.lookup("dInner")
        );
        
        // Circle ring half height
        const scalar dOuter = diameterOuter.value()/1000;
        const scalar dInner = diameterInner.value()/1000;
        const scalar d = (dOuter-dInner)/2;
        const scalar utau = Retau*nu.value()/d;
        //wall normal circulation
        const scalar duplus = Ubar.value()[0]*0.5/utau;
        //spanwise wavenumber: spacing z+ = 200
        const scalar betaPlus = 2.0*constant::mathematical::pi*(1.0/200.0);
        //const scalar sigma = 0.00055;
        const scalar sigma = 0.0002;
        //streamwise wave number: spacing x+ = 500
        const scalar alphaPlus = 2.0*constant::mathematical::pi*(1.0/800.0);
        const scalar epsilon = Ubar.value()[0]/20.0;
        //Random perturbation(1234567);
        
        forAll(centers, celli)
        {
            //scalar deviation=1.0 + 0.2*perturbation.GaussNormal();
            scalar& Ux(U[celli].x());
            vector cCenter = centers[celli];
            scalar r = ::sqrt(::sqr(cCenter.y()) + ::sqr(cCenter.z()));
            //scalar r=rr;//-innerd/2;
            //scalar ringr=mag(r-(outerd+innerd)/4);
            scalar rHalf = mag(r-(dOuter+dInner)/4);

            // Laminar parabolic profile
            Ux = 2*mag(Ubar.value())*(1-::sqr(rHalf/d));
            // Turbulent 1/7 power law profile. 
            // Turbulent profile is  not recommend according to Dr. Eugene de Villiers
            //Ux = 1.218*mag(Ubar.value())*std::pow(max((1-2*rHalf/1.01/d),SMALL),0.142857);

            // Wall distance
            scalar y = mag(min(dOuter/2-r,r-dInner/2));

            r = r*Retau/d; //zPlus->theta*r
            y = y*Retau/d; //yPlus
 
            scalar theta = ::atan(cCenter.y()/cCenter.z());
            scalar x = cCenter.x()*Retau/d;
 
            Ux = Ux + (utau*duplus/2.0)
                    *::cos(betaPlus*theta*r) *(y/30)
                    *::exp(-sigma*::sqr(y) + 0.5);
            //Ux = Ux*deviation; //White noise
            scalar utheta = epsilon*::sin(alphaPlus*x)*y
                            *::exp(-sigma*::sqr(y));
            vector tangential
            (
                 0, cCenter.y(), cCenter.z()
            );
            tangential = tangential ^ vector(1,0,0);
            tangential = tangential/mag(tangential);
 
            U[celli] = U[celli] + utheta*tangential;
 
        }
 
        U.write();
    }
    else
    {
         Info<< "    No U" << endl;
    }
 
    Info<< endl;
 
 
    return(0);
}


// ************************************************************************* //
