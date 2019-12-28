#pragma once
#include "voxelization.h"

class VoxelConeTracing
{
public:
    static void init();
    static void terminate();

    static Voxelizer* voxelizer() { return m_voxelizer; }
    static Voxelizer* virtualVoxelizer() { return m_virtualVoxelizer; }

private:
    static Voxelizer* m_voxelizer;
    static Voxelizer* m_virtualVoxelizer;
};
