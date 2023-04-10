#include <iostream>
#include <string>
#include "src/OpenMesh_Typedef.h"
#include "src/OpenMesh_Simplification.h"
#include "src/OpenMesh_Denoising.h"

using namespace std;

int main() {
    string filename = "../data/armadillo.obj";
    iGame_OpenMesh_TriMesh tri_mesh;
    tri_mesh.request_face_normals();
    tri_mesh.request_vertex_normals();
    if (!OpenMesh::IO::read_mesh(tri_mesh, filename)) {
        cout << "Can not read obj" << endl;
        return -1;
    }
    if (!tri_mesh.is_trimesh()) {
        cout << "Support trimesh only!" << endl;
        return -1;
    }
    int idx_dot = filename.find_last_of('.');
    string filename_result = filename.substr(0, idx_dot) + "_result.obj";

    iGame_OpenMesh_Simplification app_simplification(tri_mesh);
    app_simplification.simplification_qem(tri_mesh.n_vertices() * 0.5);

//    iGame_OpenMesh_Denoising app_denoisng(tri_mesh);
//    app_denoisng.denoising_bnf();

    OpenMesh::IO::write_mesh(tri_mesh, filename_result);

    return 0;
}
