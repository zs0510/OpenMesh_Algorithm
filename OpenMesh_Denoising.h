#pragma once
#include <Eigen/Core>
#include <Kernel/iGame_OpenMesh.h>
#include <Surface_Algorithm/AABB_Tree.h>
#include <omp.h>

class iGame_OpenMesh_Denoising {
public:
	iGame_OpenMesh_Denoising(OM_TriMesh& _mesh, double _sigma_normal = 0.5f, size_t _iter_normal = 5, size_t _iter_vex = 15) : 
		mesh(_mesh), sigma_normal(_sigma_normal), iter_normal(_iter_normal), iter_vex(_iter_vex) {
		initMeshData();
	}
	void denoising_bnf();
private:

	OM_TriMesh& mesh;
	OpenMesh::FPropHandleT<OM_Pos> centroids;
	OpenMesh::FPropHandleT<double> areas;
	OpenMesh::FPropHandleT<OM_Vec> normals;
	OpenMesh::FPropHandleT<OM_Vec> new_normals;
	double sigma_normal = 0.5f, sigma_centroid = 0.f;
	size_t iter_normal = 15, iter_vex = 5;
	void initMeshData();
	void BilateralNormalFiltering();
	void VertexUpdating();

};

