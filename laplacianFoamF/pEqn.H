{
    rAU = 1.0/UEqn.A();

    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn.H();

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (fvc::interpolate(HbyA) & mesh.Sf())
      + fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, Uf)
    );

    if (p_rgh.needReference())
    {
        fvc::makeRelative(phiHbyA, U);
        adjustPhi(phiHbyA, U, p_rgh);
        fvc::makeAbsolute(phiHbyA, U);
    }

    surfaceScalarField phig
    (
        (
            fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the fixedFluxPressure BCs to ensure flux consistency
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p_rgh.boundaryField(),
        (
            phiHbyA.boundaryField()
          - (mesh.Sf().boundaryField() & U.boundaryField())
        )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - p_rghEqn.flux();

            p_rgh.relax();

            U = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
        

// // // 	List<scalar> nonDiag; // sums of non-diagonal elements
// // // 
// // // 	nonDiag.resize(p.size());
// // // 	forAll(nonDiag, i)
// // // 	{
// // // 	    nonDiag[i] = 0.0;
// // // 	}
// // // 	
// // // 	const lduAddressing& addr = p_rghEqn.lduAddr();
// // // 	const labelList& lowerAddr = addr.lowerAddr();
// // // 	const labelList& upperAddr = addr.upperAddr();
// // // 	
// // // 	forAll(lowerAddr, i)
// // // 	{
// // // // 	    A[lowerAddr[i]][upperAddr[i]] = pEqn.upper()[i];
// // // // 	    A[upperAddr[i]][lowerAddr[i]] = pEqn.lower()[i];
// // // 	    nonDiag[lowerAddr[i]] += p_rghEqn.upper()[i];
// // // 	    nonDiag[upperAddr[i]] += p_rghEqn.lower()[i];
// // // 	}
// // // 	
// // // 	forAll(p.boundaryField(),I)
// // // 	{
// // // 	    const fvPatch &ptch=p.boundaryField()[I].patch();
// // // 	    forAll(ptch,J)
// // // 	    {
// // // 		int w=ptch.faceCells()[J];
// // // // 		nonDiag[w]+=UEqn.internalCoeffs()[I][J];
// // // // // // 		Info << "nonDiag" << "[" << w << "]" << " = " << nonDiag[w]  << nl << endl;
// // // 	    }
// // // 	
// // // 	}
// // // 	
// // // 
// // // 
// // // 
// // // // F = alpha1;
// // // // 	forAll(p_rghEqn.diag(),i)
// // // 	forAll(nonDiag,i)
// // // 	    {
// // // 
// // // 	      F[i] = -p_rghEqn.diag()[i]+nonDiag[i];
// // // // 	      F[i] = abs(nonDiag[i]);
// // // // 	      F[i] = abs(p_rghEqn.diag()[i])+abs(nonDiag[i]);
// // // // 	      Info << "p_rghEqn.diag()[i]=" << p_rghEqn.diag()[i] << nl << endl;
// // // // 	      Info << "nonDiag[i]=" << nonDiag[i] << nl << endl; 
// // // 	    }
// // //        
// // //         Info << "nonDiag.size()=" << nonDiag.size() << nl << endl;
// // // 	Info << "alpha1.size()=" << alpha1.size() << nl << endl;
// // // 	Info << "F.size()=" << F.size() << nl << endl;
// // // // 	Info << "p_rghEqn.size()=" << p_rghEqn.size() << nl << endl;
    }

    #include "continuityErrs.H"

    {
        Uf = fvc::interpolate(U);
        surfaceVectorField n(mesh.Sf()/mesh.magSf());
        Uf += n*(phi/mesh.magSf() - (n & Uf));
    }

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);

    p == p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }
    
    
}
