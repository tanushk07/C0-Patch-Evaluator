#include "Evaluator.h"

int main()
{
    /*
    // Test triangle on the unit sphere
    vector<Vector3> verts = { Vector3(1,0,0), Vector3(0,1,0), Vector3(0,0,1) };
    vector<Vector3> norms = verts;  

    coefficients c = MakeCoefficients(verts, norms);

    int passed = 0, total = 4;

    if (Test1_CheckOrthogonality(verts, norms))  passed++;
    if (Test2_VertexRecovery(c, verts))    passed++;
    if (Test3_NormalRecovery(c, norms))    passed++;
    if (Test4_SphereAccuracy())            passed++;

    cout << "\n Results: " << passed << "/" << total << " passed\n";
    
    Mesh two=MakeSphereMesh(5,50,50);
    WriteFlatObj(two, "Flat_Patch.obj");
    WriteNagataObj(two, "Nagata_Patch.obj", 10);
    
    return (passed == total) ? 0 : 1;
    */
    Evaluator evaluator;
    evaluator.RunSimplificationExperiment();
}
