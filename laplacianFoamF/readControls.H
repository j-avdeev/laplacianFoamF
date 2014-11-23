#include "readTimeControls.H"

//     const dictionary& potentialFlow =
//         mesh.solutionDict().subDict("potentialFlow");

//     const int nNonOrthCorr =
//         potentialFlow.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);

bool correctPhi
(
    simple.dict().lookupOrDefault<Switch>("correctPhi", true)
);

bool checkMeshCourantNo
(
    simple.dict().lookupOrDefault<Switch>("checkMeshCourantNo", false)
);

bool moveMeshOuterCorrectors
(
    simple.dict().lookupOrDefault<Switch>("moveMeshOuterCorrectors", false)
);
