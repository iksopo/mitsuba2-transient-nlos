#pragma once

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/logger.h>
#include <mitsuba/core/object.h>
#include <mitsuba/core/rfilter.h>
#include <mitsuba/core/vector.h>
#include <mitsuba/render/film.h>
#include <mitsuba/render/sampler.h>
#include <mitsuba/render/fwd.h>

NAMESPACE_BEGIN(mitsuba)

/**
 * TODO: add documentation
 * */
template <typename Float, typename Spectrum>
class MTS_EXPORT_RENDER StreakFilm : public Film<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(Film, m_size, m_crop_size, m_crop_offset,
                    m_filter, m_high_quality_edges, bitmap)
    MTS_IMPORT_TYPES(ImageBlock, StreakImageBlock, ReconstructionFilter)

    /// Merge an image block into the film. This methods should be thread-safe.
    virtual void put(const StreakImageBlock *block) = 0;

    /// Return the bitmap object storing the developed contents of the film corresponding to the i-th slice x-t
    virtual ref<Bitmap> bitmap(int slice, bool raw) = 0;

    // =============================================================
    //! @{ \name Accessor functions
    // =============================================================

    size_t num_bins() const { return m_num_bins; }

    float bin_width_opl() const { return m_bin_width_opl; }

    float start_opl() const { return m_start_opl; }

    Float end_opl() const { return start_opl() + num_bins() * bin_width_opl(); }

    const ReconstructionFilter *time_reconstruction_filter() const {
        return m_time_filter.get();
    }

    //! @}
    // =============================================================

    virtual std::string to_string() const override;

    virtual ref<StreakImageBlock> getStreakImageBlock() const = 0;

    MTS_DECLARE_CLASS()
protected:

    /// Create a film
    StreakFilm(const Properties &props);

    /// Virtual destructor
    virtual ~StreakFilm();

protected:
    uint32_t m_num_bins;
    float m_bin_width_opl;
    float m_start_opl;
    ref<ReconstructionFilter> m_time_filter;
};

MTS_EXTERN_CLASS_RENDER(StreakFilm)
NAMESPACE_END(mitsuba)
