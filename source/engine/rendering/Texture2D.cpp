#include "Texture2D.h"
#include <SOIL2.h>
#include <engine/util/Logger.h>
#include <engine/util/Timer.h>
#include "util/GLUtil.h"
#include "engine/util/convert.h"
#include <memory>
#include <fstream>

Texture2D::Texture2D(const std::string& path, Texture2DSettings settings) { load(path, settings); }

Texture2D::~Texture2D()
{
    if (isValid())
    {
        glDeleteTextures(1, &m_glId);
        m_glId = 0;
    }
}

void Texture2D::create(GLsizei width, GLsizei height, GLint internalFormat, GLenum format, GLenum pixelType, Texture2DSettings settings, const void* pixels)
{
    if (isValid())
    {
        glDeleteTextures(1, &m_glId);
        m_glId = 0;
    }

    m_target = GL_TEXTURE_2D;
    m_format = format;
    m_internalFormat = internalFormat;
    m_pixelType = pixelType;
    m_width = int(width);
    m_height = int(height);

    switch (format)
    {
    case GL_LUMINANCE:
    case GL_DEPTH_COMPONENT:
        m_channels = 1;
        break;
    case GL_LUMINANCE_ALPHA:
        m_channels = 2;
        break;
    case GL_RGB:
        m_channels = 3;
        break;
    case GL_RGBA:
        m_channels = 4;
        break;
    default: break;
    }

    glGenTextures(1, &m_glId);

    glBindTexture(m_target, m_glId);
    glTexImage2D(m_target, 0, m_internalFormat, m_width, m_height, 0, m_format, m_pixelType, pixels);

    applySettings(settings);
}

void Texture2D::createMultisampled(uint8_t numSamples, GLsizei width, GLsizei height, GLint internalFormat)
{
    assert(numSamples > 1);

    if (isValid())
    {
        glDeleteTextures(1, &m_glId);
        m_glId = 0;
    }

    m_target = GL_TEXTURE_2D_MULTISAMPLE;
    m_internalFormat = internalFormat;
    m_pixelType = 0;
    m_width = int(width);
    m_height = int(height);
    m_numSamples = numSamples;

    glGenTextures(1, &m_glId);

    glBindTexture(m_target, m_glId);
    glTexImage2DMultisample(m_target, GLsizei(numSamples), m_internalFormat, m_width, m_height, false);
}

void Texture2D::save(const std::string& path, bool isFBO, GLuint fboId, Texture2DSettings settings)
{
    if (!isValid()) {
        LOG("Texture is not valid. Failed to save: " + path);
        return;
    }

    unsigned int pixelTypeBytes = 0;
    switch(m_pixelType) {
    case GL_UNSIGNED_BYTE:
        pixelTypeBytes = 1;
        break;
    case GL_FLOAT:
        pixelTypeBytes = 4;
        break;
    default:
        LOG("Unknown pixelType for saving texture");
        return;
    }

    unsigned int bufSize = m_width * m_height * m_channels;
    // std::cout << "bufSize: " << bufSize << std::endl;
    void *data = new GLubyte[bufSize * pixelTypeBytes];

    if (isFBO)
    {
        glBindFramebuffer(GL_FRAMEBUFFER, fboId);
        glReadBuffer(GL_COLOR_ATTACHMENT5);
        glReadPixels(0, 0, m_width, m_height, m_format, m_pixelType, data);
    }
    else
    {
        // requires OpenGL 4.5
        glGetTextureImage(m_glId, 0, m_format, m_pixelType, bufSize, data);
    }
    GL_ERROR_CHECK();

    // invert y
    switch(m_pixelType) {
    case GL_UNSIGNED_BYTE:{
        auto dataPtr = (GLubyte*) data;
        for (int i = 0; i * 2 < m_height; ++i)
        {
            int idx0 = i * m_width * m_channels;
            int idx1 = (m_height - 1 - i) * m_width * m_channels;
            for (int j = m_width * m_channels; j > 0; --j)
                std::swap(dataPtr[idx0++], dataPtr[idx1++]);
        }
        if (!SOIL_save_image((path + ".bmp").c_str(), SOIL_SAVE_TYPE_BMP, m_width, m_height, m_channels, (GLubyte*)data))
        {
            LOG("Failed to save: " + path);
        }
        break;
    }
    case GL_FLOAT: {
        auto dataPtr = (GLfloat*) data;
        for (int i = 0; i * 2 < m_height; ++i)
        {
            int idx0 = i * m_width * m_channels;
            int idx1 = (m_height - 1 - i) * m_width * m_channels;
            for (int j = m_width * m_channels; j > 0; --j)
                std::swap(dataPtr[idx0++], dataPtr[idx1++]);
        }

        std::ofstream out("position.xyz");
        GLsizei rowSize = m_width * m_channels;
        out << m_width * m_height << std::endl;
        for (GLsizei y = 0; y < m_height; ++y) 
        {
            for (GLsizei x = 0; x < rowSize; x += m_channels)
            {
                for (GLsizei c = 0; c < m_channels; ++c)
                {
                    out << dataPtr[y * rowSize + x + c] << " ";
                }
                out << "\n";
            }
        }
        out.close();
        break;
    }
    default:
        LOG("Unknown pixelType for saving texture");
    }
    delete (GLubyte*)data;
}

void Texture2D::load(const std::string& path, Texture2DSettings settings)
{
    m_target = GL_TEXTURE_2D;

    if (isValid())
    {
        glDeleteTextures(1, &m_glId);
        m_glId = 0;
    }

    auto imgData = SOIL_load_image(path.c_str(), &m_width, &m_height, &m_channels, SOIL_LOAD_AUTO);
    GL_ERROR_CHECK();

    if (!imgData)
    {
        LOG("Failed to load: " + path);
        return;
    }

    // Using SOIL_create_OGL_texture() results in OpenGL errors due to the use of deprecated functionaly (GL_LUMINANCE)
    glGenTextures(1, &m_glId);
    glBindTexture(m_target, m_glId);

    switch (m_channels)
    {
    case 1:
    {
        m_format = GL_RED;
        m_internalFormat = GL_R16;
        GLint swizzleMask[] = {GL_RED, GL_RED, GL_RED, GL_ONE};
        glTexParameteriv(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_RGBA, swizzleMask);
        break;
    }
    case 2:
    {
        m_format = GL_RG;
        m_internalFormat = GL_RG16;
        GLint swizzleMask[] = {GL_RED, GL_RED, GL_RED, GL_GREEN};
        glTexParameteriv(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_RGBA, swizzleMask);
        break;
    }
    case 3:
        m_format = GL_RGB;
        m_internalFormat = GL_RGB8;
        break;
    case 4:
        m_format = GL_RGBA;
        m_internalFormat = GL_RGBA8;
        break;
    default: break;
    }

    // Invert_Y
    for (int i = 0; i * 2 < m_height; ++i)
    {
        int idx0 = i * m_width * m_channels;
        int idx1 = (m_height - 1 - i) * m_width * m_channels;
        for (int j = m_width * m_channels; j > 0; --j)
            std::swap(imgData[idx0++], imgData[idx1++]);
    }

    m_pixelType = GL_UNSIGNED_BYTE;
    glTexImage2D(m_target, 0, m_internalFormat, m_width, m_height, 0, m_format, m_pixelType, imgData);

    applySettings(settings);

    SOIL_free_image_data(imgData);
    GL_ERROR_CHECK();
}

void Texture2D::setParameteri(GLenum name, GLint value) const
{
    assert(isValid());
    assert(GL::isTextureBound(m_target, m_glId));

    glTexParameteri(m_target, name, value);
}

void Texture2D::resize(GLsizei width, GLsizei height)
{
    m_width = width;
    m_height = height;

    glBindTexture(m_target, m_glId);
    if (m_target == GL_TEXTURE_2D)
        glTexImage2D(m_target, 0, m_internalFormat, m_width, m_height, 0, m_format, m_pixelType, nullptr);
    else
        glTexImage2DMultisample(m_target, GLsizei(m_numSamples), m_internalFormat, m_width, m_height, false);

    glBindTexture(m_target, 0);
}

void Texture2D::applySettings(Texture2DSettings settings)
{
    switch (settings)
    {
    case Texture2DSettings::S_T_REPEAT_MIN_MAG_LINEAR:
        glTexParameteri(m_target, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(m_target, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(m_target, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(m_target, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        break;
    case Texture2DSettings::S_T_REPEAT_MIN_MIPMAP_LINEAR_MAG_LINEAR:
        glGenerateMipmap(GL_TEXTURE_2D);
        glTexParameteri(m_target, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(m_target, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(m_target, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(m_target, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        break;
    case Texture2DSettings::S_T_CLAMP_TO_BORDER_MIN_MIPMAP_LINEAR_MAG_LINEAR:
        glGenerateMipmap(GL_TEXTURE_2D);
        glTexParameteri(m_target, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
        glTexParameteri(m_target, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
        glTexParameteri(m_target, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(m_target, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        break;
    case Texture2DSettings::S_T_REPEAT_MIN_MAG_NEAREST:
        glTexParameteri(m_target, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(m_target, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(m_target, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(m_target, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        break;
    case Texture2DSettings::S_T_REPEAT_ANISOTROPIC:
        glGenerateMipmap(GL_TEXTURE_2D);
        glTexParameteri(m_target, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(m_target, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(m_target, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(m_target, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

        GLfloat fLargest;
        glGetFloatv(GL_MAX_TEXTURE_MAX_ANISOTROPY_EXT, &fLargest);
        glTexParameterf(m_target, GL_TEXTURE_MAX_ANISOTROPY_EXT, fLargest);
        break;
    case Texture2DSettings::S_T_CLAMP_TO_BORDER_MIN_MAX_NEAREST:
        glTexParameteri(m_target, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
        glTexParameteri(m_target, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
        glTexParameteri(m_target, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(m_target, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        break;
    case Texture2DSettings::S_T_CLAMP_TO_BORDER_MIN_MAX_LINEAR:
        glTexParameteri(m_target, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
        glTexParameteri(m_target, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
        glTexParameteri(m_target, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(m_target, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        break;
    case Texture2DSettings::Custom:
        break;
    default:
        break;
    }
}
