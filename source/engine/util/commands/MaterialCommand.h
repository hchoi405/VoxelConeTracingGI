#pragma once
#include <glm/glm.hpp>
#include <engine/geometry/Transform.h>
#include "engine/rendering/renderer/MeshRenderers.h"
#include <glm/ext.hpp>
#include <engine/util/math.h>
#include "Command.h"

/**
* Change material property from a start value to a destination value within a given time constraint.
*/
struct MaterialCommand : Command
{
    MaterialCommand() {}

    MaterialCommand(const ComponentPtr<MeshRenderer>& renderer, std::string property, const glm::vec4& startValue, const glm::vec4& destinationValue, float duration = 1.0f)
        : renderer(renderer), property(property), startValue(startValue), destinationValue(destinationValue), duration(duration)
    {
        if (this->duration <= 0.0f)
            this->duration = math::EPSILON;
    }

    void operator()(float deltaTime) override
    {
        if (t >= 1.0f || !renderer)
            return;

        t += deltaTime / duration;
        glm::vec4 curValue = glm::lerp(startValue, destinationValue, math::clamp(t, 0.0f, 1.0f));
        auto material = renderer->getMaterial(0);
        if (property == "u_shininess") {
            material->setFloat("u_shininess", curValue[0]);
        }
        else if (property == "u_color") {
            material->setColor("u_color", curValue);
        }
        else if (property == "u_emissionColor") {
            material->setColor("u_emissionColor", glm::vec3(curValue));
        }
        else if (property == "u_specularColor") {
            material->setColor("u_specularColor", glm::vec3(curValue));
        }
    }

    bool done() const override { return t >= 1.0f || !renderer; }

    void reset() override { t = 0.0f; }

    ComponentPtr<MeshRenderer> renderer;
    std::string property;
    glm::vec4 startValue;
    glm::vec4 destinationValue;
    float t{ 0.0f };
    float duration{ 1.0f };
};
