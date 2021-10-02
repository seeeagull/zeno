#pragma once


#include <zeno/UI/Label.h>


struct Button : Label {
    SignalSlot on_clicked;

    void on_event(Event_Mouse e) override;
    void paint() const override;
};
