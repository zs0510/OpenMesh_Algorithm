#pragma once
#include <Kernel/iGame_OpenMesh.h>
#include <omp.h>

class iGame_OpenMesh_Smoother {
private:
	iGame_OpenMesh_TriMesh& mesh;
public:
	iGame_OpenMesh_Smoother(iGame_OpenMesh_TriMesh& _mesh): mesh(_mesh) {
	
	}
	void smoothing(size_t iter_times);// average smoothing
};