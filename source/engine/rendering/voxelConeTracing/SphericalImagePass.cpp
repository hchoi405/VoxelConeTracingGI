#include "SphericalImagePass.h"
#include <engine/rendering/Texture3D.h>
#include <engine/rendering/architecture/RenderPipeline.h>
#include <engine/rendering/util/GLUtil.h>
#include <engine/rendering/Screen.h>
#include "Globals.h"
#include <engine/geometry/Transform.h>
#include <engine/resource/ResourceManager.h>
#include <engine/rendering/renderer/MeshRenderers.h>
#include "settings/VoxelConeTracingSettings.h"
#include "VoxelRegion.h"
#include "engine/util/ECSUtil/ECSUtil.h"
#include <SOIL2.h>
#include <fstream>

SphericalImagePass::SphericalImagePass()
    : RenderPass("SphericalImagePass")
{
    // Input position map
    GLfloat *posData = nullptr, *normalData = nullptr;
    posData = new GLfloat[m_height * m_width * 3];
    normalData = new GLfloat[m_height * m_width * 3];
    readData("../assets/cglab/matterport_pos_dasan613.txt", posData);
    readData("../assets/cglab/matterport_normal_dasan613.txt", normalData);

    m_positionMap = std::make_shared<Texture2D>();
    m_positionMap->create(m_width, m_height, GL_RGB32F, GL_RGB, GL_FLOAT, Texture2DSettings::S_T_REPEAT_MIN_MAG_LINEAR, posData);

    m_normalMap = std::make_shared<Texture2D>();
    m_normalMap->create(m_width, m_height, GL_RGB32F, GL_RGB, GL_FLOAT, Texture2DSettings::S_T_REPEAT_MIN_MAG_LINEAR, normalData);

    // Set up framebuffer
    m_framebuffer = std::make_unique<Framebuffer>(m_width, m_height, false);
    m_framebuffer->begin();

    // Output Spherical Map
    m_sphericalTexture = std::make_shared<Texture2D>();
    m_sphericalTexture->create(m_width, m_height, GL_RGB8, GL_RGB, GL_UNSIGNED_BYTE, Texture2DSettings::S_T_REPEAT_MIN_MAG_LINEAR);
    m_framebuffer->attachRenderTexture2D(m_sphericalTexture, GL_COLOR_ATTACHMENT0);

    m_framebuffer->end();

    m_shader = ResourceManager::getShader("shaders/voxelConeTracing/sphericalImagePass.vert", "shaders/voxelConeTracing/sphericalImagePass.frag");
}

void SphericalImagePass::readData(std::string filename, GLfloat *data)
{
    std::ifstream in(filename);
    if (!in.is_open())
    {
        std::cout << "Failed to open: " << std::endl;
        exit(-1);
    }
    in >> m_width >> m_height;
    uint64_t idx = 0;
    while (!in.eof())
    {
        in >>
            data[idx + 0] >>
            data[idx + 1] >>
            data[idx + 2];
        data[idx + 2] = -data[idx + 2];

        idx += 3;
    }
}

void SphericalImagePass::render() const
{
    m_framebuffer->begin();
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if (RENDERING_SETTINGS.cullBackFaces)
    {
        glFrontFace(GL_CW);
        glEnable(GL_CULL_FACE);
    }
    else
    {
        glDisable(GL_CULL_FACE);
    }

    m_shader->bind();

    Texture3D *voxelRadiance = m_renderPipeline->fetchPtr<Texture3D>("VoxelRadiance");
    Texture3D *voxelOpacity = m_renderPipeline->fetchPtr<Texture3D>("VoxelOpacity");
    auto clipRegions = m_renderPipeline->fetchPtr<std::vector<VoxelRegion>>("ClipRegions");

    GLint textureUnit = 100;
    m_shader->bindTexture2D(*m_positionMap, "u_positionMap", textureUnit++);
    m_shader->bindTexture2D(*m_normalMap, "u_normalMap", textureUnit++);
    m_shader->bindTexture3D(*voxelRadiance, "u_voxelRadiance", textureUnit++);
    m_shader->bindTexture3D(*voxelOpacity, "u_voxelOpacity", textureUnit++);

    // VCT
    m_shader->setFloat("u_voxelSizeL0", clipRegions->at(0).voxelSize);
    m_shader->setUnsignedInt("u_volumeDimension", VOXEL_RESOLUTION);
    m_shader->setVector("u_volumeCenterL0", clipRegions->at(0).getCenterPosWorld());
    m_shader->setFloat("u_stepFactor", GI_SETTINGS.stepFactor);
    m_shader->setFloat("u_occlusionDecay", GI_SETTINGS.occlusionDecay);
    m_shader->setFloat("u_traceStartOffset", GI_SETTINGS.traceStartOffset);

    // Spherical image
    m_shader->setVector("u_sphericalCenter", glm::vec3(1.2987847328186035, 1.3614389896392822, -1.2158539295196533));
    m_shader->setFloat("u_viewAperture", DEBUG_SETTINGS.viewAperture);
    m_shader->setFloat("u_hitpointOffset", DEBUG_SETTINGS.hitpointOffset);

    ECSUtil::renderEntities(m_shader.get());

    glDisable(GL_CULL_FACE);
    m_framebuffer->end();

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void SphericalImagePass::update()
{
    // Store to the sphericalTexture
    render();

    m_framebuffer->bind();
    glReadBuffer(GL_COLOR_ATTACHMENT0);
    GLsizei channels = 3;
    GLubyte *data = new GLubyte[m_width * m_height * channels];
    glReadPixels(0, 0, m_width, m_height, GL_RGB, GL_UNSIGNED_BYTE, data);

    for (int i = 0; i * 2 < m_height; ++i)
    {
        int idx0 = i * m_width * channels;
        int idx1 = (m_height - 1 - i) * m_width * channels;
        for (int j = m_width * channels; j > 0; --j)
            std::swap(data[idx0++], data[idx1++]);
    }

    if (!SOIL_save_image(
            "test.bmp",
            SOIL_SAVE_TYPE_BMP,
            m_width, m_height, channels,
            data))
    {
        std::cout << "Failed to save image." << std::endl;
    }
    delete data;
}
