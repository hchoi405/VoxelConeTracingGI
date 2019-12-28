#pragma once
#include <engine/rendering/architecture/RenderPass.h>
#include <engine/rendering/shader/Shader.h>
#include <engine/rendering/Framebuffer.h>
#include <memory>

class SphericalImagePass : public RenderPass
{
public:
    SphericalImagePass();
    void render() const;

    void update() override;

    void readData(std::string filename, GLfloat* data);

private:
    GLsizei m_width = 4096;
    GLsizei m_height = 2048;
    std::shared_ptr<Shader> m_shader;
    std::unique_ptr<Framebuffer> m_framebuffer;

    std::shared_ptr<Texture2D> m_sphericalTexture;
    std::shared_ptr<Texture2D> m_positionMap;
    std::shared_ptr<Texture2D> m_normalMap;
    
};
