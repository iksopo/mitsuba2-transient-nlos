#include <mitsuba/core/logger.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/progress.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/radiancesample.h>
#include <mitsuba/render/sampler.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/sensor.h>
#include <mitsuba/render/streakfilm.h>
#include <mitsuba/render/transientintegrator.h>
#include <mutex>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <thread>

NAMESPACE_BEGIN(mitsuba)

MTS_VARIANT StreakFilm<Float, Spectrum>::StreakFilm(const Properties &props) : Base(props) {
    // NOTE(diego): old names are included for compatibility with old files
    m_num_bins = props.int_("num_bins", props.int_("time", 1000));
    m_bin_width_opl =
        props.float_("bin_width_opl", props.float_("exposure_time", -1.f));
    m_start_opl = props.float_("start_opl", props.float_("time_offset", -1.f));
    m_auto_detect_bins = props.bool_("auto_detect_bins", false);
    Assert(!m_auto_detect_bins || m_start_opl < 0.f);

    // Use the provided reconstruction filter, if any.
    for (auto &[name, obj] : props.objects(false)) {
        if(name != "tfilter") {
           continue;
        }
        auto *tfilter = dynamic_cast<ReconstructionFilter *>(obj.get());
        if (tfilter) {
            if (m_time_filter)
                Throw("A film can only have one reconstruction filter.");
            m_time_filter = tfilter;
            props.mark_queried(name);
        }
    }

    if (!m_time_filter) {
        // No reconstruction filter has been selected. Load a Gaussian filter by default
        m_time_filter =
            PluginManager::instance()->create_object<ReconstructionFilter>(Properties("gaussian"));
    }
}

MTS_VARIANT StreakFilm<Float, Spectrum>::~StreakFilm() {}

MTS_VARIANT void StreakFilm<Float, Spectrum>::auto_detect_bins(Scene *scene,
                                                               Sensor *sensor) {
    ThreadEnvironment env;
    ref<ProgressReporter> progress =
        new ProgressReporter("Auto-detecting StreakFilm histogram limits");
    std::mutex mutex;
    ref<TransientSamplingIntegrator> integrator =
        dynamic_cast<TransientSamplingIntegrator *>(scene->integrator());

    // Total number of samples
    size_t total_samples = 500000000, samples_done = 0;
    Float global_min_opl          = std::numeric_limits<Float>::max(),
          global_max_opl          = std::numeric_limits<Float>::min();
    ScalarFloat diff_scale_factor = rsqrt((ScalarFloat) total_samples);

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, total_samples, 1),
        [&](const tbb::blocked_range<size_t> &range) {
            ScopedSetThreadEnvironment set_env(env);
            ref<Sampler> sampler = sensor->sampler()->clone();
            for (auto i = range.begin(); i != range.end(); ++i) {
                Vector2f position_sample = sampler->next_2d();

                Point2f aperture_sample(.5f);
                Float max_opl           = math::Infinity<Float>;
                Float wavelength_sample = sampler->next_1d();

                auto [ray, ray_weight] = sensor->sample_ray_differential(
                    0.f, wavelength_sample, position_sample, aperture_sample);

                ray.scale_differential(diff_scale_factor);

                const Medium *medium = sensor->medium();
                std::vector<FloatTimeSample<Float, Mask>> aovs_record;
                std::vector<RadianceSample<Float, Spectrum, Mask>>
                    timed_samples_record;
                integrator->sample(scene, sampler, ray, medium, aovs_record,
                                   timed_samples_record, max_opl);

                size_t n_timed_samples = timed_samples_record.size();

                if (n_timed_samples == 0) {
                    continue;
                }

                Float local_min_opl = std::numeric_limits<Float>::max(),
                      local_max_opl = std::numeric_limits<Float>::min();
                for (size_t j = 0; j < n_timed_samples; j++) {
                    // pick samples with >0 radiance
                    Mask good = any_inner(
                        depolarize<Spectrum>(timed_samples_record[j].radiance) >
                        math::Epsilon<Float>);
                    good &= timed_samples_record[j].mask;
                    if (any_or<true>(good)) {
                        local_min_opl =
                            hmin(select(good, timed_samples_record[j].opl,
                                        std::numeric_limits<Float>::max()));
                        break;
                    }
                }
                // NOTE(diego): j is an unsigned value, cant check for >= 0
                for (size_t j = n_timed_samples - 1; j < n_timed_samples; j--) {
                    // pick samples with >0 radiance
                    Mask good = any_inner(
                        depolarize<Spectrum>(timed_samples_record[j].radiance) >
                        math::Epsilon<Float>);
                    good &= timed_samples_record[j].mask;
                    if (any_or<true>(good)) {
                        local_max_opl =
                            hmax(select(good, timed_samples_record[j].opl,
                                        std::numeric_limits<Float>::min()));
                        break;
                    }
                }

                /* Critical section: update min/max OPL */ {
                    global_max_opl = std::max(global_max_opl, local_max_opl);
                    global_min_opl = std::min(global_min_opl, local_min_opl);
                }
            }

            /* Critical section: update progress bar */ {
                std::lock_guard<std::mutex> lock(mutex);
                samples_done += range.end() - range.begin();
                progress->update(samples_done / (ScalarFloat) total_samples);
            }
        });

    Float initial_bin_width    = (global_max_opl - global_min_opl) / m_num_bins;
    size_t edge_padding_before = m_num_bins / 20; // bins
    size_t edge_padding_after  = m_num_bins / 5;  // bins
    // add some left-right padding just in case
    m_start_opl     = global_min_opl - edge_padding_before * initial_bin_width;
    Float m_end_opl = global_max_opl + edge_padding_after * initial_bin_width;
    if (m_bin_width_opl < 0.f) {
        m_bin_width_opl = (m_end_opl - m_start_opl) / m_num_bins;
    } else {
        Log(Info, "Using user-specified StreakFilm bin width: %i",
            m_bin_width_opl);
        m_end_opl = m_start_opl + m_bin_width_opl * m_num_bins;
    }

    Log(Info,
        "Auto-detected StreakFilm histogram OPL limits: [%i, %i] with bin "
        "width %i",
        m_start_opl, m_end_opl, m_bin_width_opl);
}

MTS_VARIANT std::string StreakFilm<Float, Spectrum>::to_string() const {
    std::ostringstream oss;
    oss << "StreakFilm[" << std::endl
        << "  size = " << m_size << "," << std::endl
        << "  crop_size = " << m_crop_size << "," << std::endl
        << "  crop_offset = " << m_crop_offset << "," << std::endl
        << "  num_bins = " << m_num_bins << "," << std::endl
        << "  bin_width_opl = " << m_bin_width_opl << "," << std::endl
        << "  start_opl = " << m_start_opl << "," << std::endl
        << "  auto_detect_bins = " << m_auto_detect_bins << "," << std::endl
        << "  high_quality_edges = " << m_high_quality_edges << "," << std::endl
        << "  filter = " << m_filter << std::endl
        << "  time_filter = " << m_time_filter << std::endl
        << "]";
    return oss.str();
}


MTS_IMPLEMENT_CLASS_VARIANT(StreakFilm, Film)
MTS_INSTANTIATE_CLASS(StreakFilm)
NAMESPACE_END(mitsuba)
