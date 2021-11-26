#pragma once

#include <deque>

#include <mitsuba/core/fwd.h>
#include <mitsuba/render/fwd.h>

NAMESPACE_BEGIN(mitsuba)

/**
 * \brief Spectrum container with additional time information. Whereas standard
 * integrators return a `std::pair<Spectrum, Mask>`, transient integrators
 * return multiple separated RadianceSamples.
 *
 * \remark Time is stored as optical path length (i.e. distance * medium_ior),
 * because speed of light is big
 */
template<typename Float, typename Spectrum>
struct MTS_EXPORT_RENDER RadianceSample {
    MTS_IMPORT_CORE_TYPES()

    Float opl;
    Spectrum radiance;
    Mask mask;

    RadianceSample(Float opl, Spectrum radiance, Mask mask)
        : opl(opl), radiance(radiance), mask(mask) {}

    void add_opl(Float opl) { opl += opl; }
    void add_opl(Float distance, Float medium_eta) {
        opl += distance * medium_eta;
    }

    decltype(auto) time() const { return opl / math::SpeedOfLight<Float>; }
    decltype(auto) time() { return opl / math::SpeedOfLight<Float>; }
};

/**
 * \brief Contains arbitrary time-tagged color information (i.e. a `std::deque`
 * of values instead of `Spectrum`)
 * 
 * \remark Time is stored as optical path length (i.e. distance * medium_ior),
 * because speed of light is big
 */
template<typename Float>
struct MTS_EXPORT_RENDER FloatSample {
    MTS_IMPORT_CORE_TYPES()

    Float opl;
    std::deque<Float> values;
    Mask mask;

    FloatSample(Float opl, Mask mask) : opl(opl), mask(mask) {}

    void push_front(Float value) { values.push_front(value); }

    decltype(auto) time() const { return opl / math::SpeedOfLight<Float>; }
    decltype(auto) time() { return opl / math::SpeedOfLight<Float>; }
};

NAMESPACE_END(mitsuba)