#pragma once
#include <glm/glm.hpp>
#include <engine/geometry/Transform.h>
#include <glm/ext.hpp>
#include <engine/util/math.h>
#include "Command.h"

/**
* Moves a transform from a start position to a destination position within a given time constraint.
*/
struct TransformCommand : Command
{
    TransformCommand() {}

    TransformCommand(const ComponentPtr<Transform>& transform, const glm::vec3& startPos, const glm::vec3& dstPos, const glm::quat startRot, const glm::quat dstRot, float duration = 1.0f)
        : transform(transform), startPos(startPos), dstPos(dstPos), startRot(startRot), dstRot(dstRot), duration(duration)
    {
        if (this->duration <= 0.0f)
            this->duration = math::EPSILON;
    }

    void operator()(float deltaTime) override
    {
        if (t >= 1.0f || !transform)
            return;

        t += deltaTime / duration;
        glm::vec3 curPos = glm::lerp(startPos, dstPos, math::clamp(t, 0.0f, 1.0f));
        transform->setPosition(curPos);
        glm::quat curRotation = glm::slerp(startRot, dstRot, math::clamp(t, 0.0f, 1.0f));
        transform->setRotation(curRotation);
    }

    bool done() const override { return t >= 1.0f || !transform; }

    void reset() override { t = 0.0f; }

    ComponentPtr<Transform> transform;
    glm::vec3 startPos, dstPos;
    glm::quat startRot, dstRot;
    float t{0.0f};
    float duration{1.0f};
};
