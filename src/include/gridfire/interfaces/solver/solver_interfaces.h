#pragma once

#include "fourdst/plugin/plugin.h"

class SolverPluginInterface : public fourdst::plugin::PluginBase {
public:
    using PluginBase::PluginBase;

    ~SolverPluginInterface() override = default;
    virtual void log_time(double t, double dt) = 0;
};