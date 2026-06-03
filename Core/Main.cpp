#include "Evaluator.h"
#include "Tests.h"
#include "MakeSphereMesh.h"
#include "MakeTorusMesh.h"
#include "Vector3.h"
#include <iostream>

using namespace std;

int main()
{
    // Test triangle on the unit sphere
    vector<Vector3> verts = { Vector3(1,0,0), Vector3(0,1,0), Vector3(0,0,1) };
    vector<Vector3> norms = verts;  

    coefficients c = Evaluator::MakeCoefficients(verts, norms);

    int passed = 0, total = 4;

    Mesh two;
    MakeSphereMesh(5,50,50,two);
    
    if (Tests::Test1_CheckOrthogonality(verts, norms))  passed++;
    if (Tests::Test2_VertexRecovery(c, verts))    passed++;
    if (Tests::Test3_NormalRecovery(c, norms))    passed++;
    if (Tests::Test4_Accuracy(two,5))       passed++;
    cout << "\n Results: " << passed << "/" << total << " passed\n";
    
    Mesh TorusMesh;
    MakeTorusMesh(5,3,20,10,TorusMesh);
    
    Evaluator::WriteFlatObj(TorusMesh, "FlatTorus_Patch.obj");
    Evaluator::WriteNagataObj(TorusMesh, "NagataTorus_Patch.obj", 10);
    
    Evaluator::RunSimplificationExperiment(); 
    
    return (passed == total) ? 0 : 1;
}
