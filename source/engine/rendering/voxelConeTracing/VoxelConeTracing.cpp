#include "VoxelConeTracing.h"
#include "Downsampler.h"
#include "CopyAlpha.h"
#include "engine/rendering/util/ImageCleaner.h"
#include "Globals.h"

Voxelizer* VoxelConeTracing::m_voxelizer = nullptr;
Voxelizer* VoxelConeTracing::m_virtualVoxelizer = nullptr;

void VoxelConeTracing::init()
{
    ImageCleaner::init();
    Downsampler::init();
    CopyAlpha::init();

    m_voxelizer = new Voxelizer(VOXEL_RESOLUTION);
    m_virtualVoxelizer = new Voxelizer(VIRTUAL_VOXEL_RESOLUTION);
}

void VoxelConeTracing::terminate()
{
    if (m_voxelizer)
        delete m_voxelizer;
    if (m_virtualVoxelizer)
        delete m_virtualVoxelizer;
}
