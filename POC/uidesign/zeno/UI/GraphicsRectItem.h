#pragma once


#include <zeno/UI/GraphicsWidget.h>


struct GraphicsRectItem : GraphicsWidget {
    void paint() const override {
        if (selected || pressed[0]) {
            glColor3f(0.75f, 0.5f, 0.375f);
        } else if (hovered) {
            glColor3f(0.375f, 0.5f, 1.0f);
        } else {
            glColor3f(0.375f, 0.375f, 0.375f);
        }
        glRectf(bbox.x0, bbox.y0, bbox.x0 + bbox.nx, bbox.y0 + bbox.ny);
    }
};
