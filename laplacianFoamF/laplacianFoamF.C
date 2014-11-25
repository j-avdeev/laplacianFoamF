/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    interDyMFoam

Description
    Solver for 2 incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"
#include "subCycle.H"
#include "interfaceProperties.H"
#include "simpleControl.H"
#include "fvIOoptionList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"

    simpleControl simple(mesh);

    #include "createFields.H"
    #include "readTimeControls.H"
//     #include "createPrghCorrTypes.H"

//     volScalarField rAU
//     (
//         IOobject
//         (
//             "rAU",
//             runTime.timeName(),
//             mesh,
//             IOobject::READ_IF_PRESENT,
//             IOobject::AUTO_WRITE
//         ),
//         mesh,
//         dimensionedScalar("rAUf", dimTime/rho.dimensions(), 1.0)
//     );

//     #include "correctPhi.H"
//     #include "createUf.H"
//     #include "CourantNo.H"
//     #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readControls.H"
//         #include "alphaCourantNo.H"
//         #include "CourantNo.H"

//         #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (simple.loop())
        {
//             if (moveMeshOuterCorrectors)
//             {
                scalar timeBeforeMeshUpdate = runTime.elapsedCpuTime();

                mesh.update();
		
		Info << "mesh updated" << endl;
                if (mesh.changing())
                {
                    Info<< "Execution time for mesh.update() = "
                        << runTime.elapsedCpuTime() - timeBeforeMeshUpdate
                        << " s" << endl;

//                     gh = g & mesh.C();
//                     ghf = g & mesh.Cf();
                }


                if (mesh.changing() && checkMeshCourantNo)
                {
                    #include "meshCourantNo.H"
                }
//             }

//             #include "alphaControls.H"

//             if (pimple.firstIter() || alphaOuterCorrectors)
//             {
//                 twoPhaseProperties.correct();
// 
//                 #include "alphaEqnSubCycle.H"
//                 interface.correct();
//             }

//             #include "UEqn.H"
	    fvScalarMatrix TEqn
	    (
	      fvm::ddt(T) - fvm::laplacian(DT, T)
	    );
	    
	    while (simple.correctNonOrthogonal())
	    {
		TEqn.solve();
	    }
	    
	    

	
	List<scalar> nonDiag; // sums of non-diagonal elements

	nonDiag.resize(T.size());
	forAll(nonDiag, i)
	{
	    nonDiag[i] = 0.0;
	}
	
	const lduAddressing& addr = TEqn.lduAddr();
	const labelList& lowerAddr = addr.lowerAddr();
	const labelList& upperAddr = addr.upperAddr();
	
	forAll(lowerAddr, i)
	{
// 	    A[lowerAddr[i]][upperAddr[i]] = pEqn.upper()[i];
// 	    A[upperAddr[i]][lowerAddr[i]] = pEqn.lower()[i];
	    nonDiag[lowerAddr[i]] += TEqn.upper()[i];
	    nonDiag[upperAddr[i]] += TEqn.lower()[i];
	}
	
	forAll(T.boundaryField(),I)
	{
	    const fvPatch &ptch=T.boundaryField()[I].patch();
	    forAll(ptch,J)
	    {
		int w=ptch.faceCells()[J];
// 		nonDiag[w]+=UEqn.internalCoeffs()[I][J];
// // // 		Info << "nonDiag" << "[" << w << "]" << " = " << nonDiag[w]  << nl << endl;
	    }
	
	}
	
// F = alpha1;
// 	forAll(p_rghEqn.diag(),i)
	forAll(TEqn.diag(),i)
	    {

// 	      F[i] = alpha1Eqn.diag()[i]+nonDiag[i];
F[i] = TEqn.diag()[i];
// 	      F[i] = abs(nonDiag[i]);
// 	      F[i] = abs(p_rghEqn.diag()[i])+abs(nonDiag[i]);
// 	      Info << "p_rghEqn.diag()[i]=" << p_rghEqn.diag()[i] << nl << endl;
// 	      Info << "UEqn.diag()=" << UEqn.diag() << nl << endl; 
	    }
       
        Info << "nonDiag.size()=" << nonDiag.size() << nl << endl;
	Info << "U.size()=" << T.size() << nl << endl;
	Info << "F.size()=" << F.size() << nl << endl;
// 	Info << "p_rghEqn.size()=" << p_rghEqn.size() << nl << endl;	
            
            
        // Indicators for refinement. Note: before runTime++
        // only for postprocessing reasons.
//         tmp<volScalarField> tmagGradP = mag(fvc::grad(p));
// 	//! not good to recreate p_rghEqn here as in pEqn.H
// 	fvScalarMatrix p_rghEqn
// 		(
// 		    fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
// 		);

// std::cout << typeid(a).name() << '\n';
// Info << "tmagGradP = " << typeid(tmagGradP) << endl;
// Info << "mag(UEqn.diag()) = " << UEqn.diag().type() << endl;
// Info << "mag(fvc::grad(p)) = " << fvc::grad(p).type() << endl;
// 	tmp<volScalarField> tmagGradP = mag(UEqn.diag());
        volScalarField normalisedGradP
        (
            "normalisedGradP",
            F/max(F)
        );
        normalisedGradP.writeOpt() = IOobject::AUTO_WRITE;
//         tmagGradP.clear();

	F = normalisedGradP;            
        
	    #include "write.H"
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
