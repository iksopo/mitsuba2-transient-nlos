#include <enoki/stl.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/ray.h>
#include <mitsuba/render/ior.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/emitter.h>
#include <mitsuba/render/records.h>
#include <mitsuba/render/transientintegrator.h>
#include <random>

NAMESPACE_BEGIN(mitsuba)

template <typename Float, typename Spectrum>
class TransientPathIntegrator : public TransientMonteCarloIntegrator<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(TransientMonteCarloIntegrator, m_max_depth, m_rr_depth)
    MTS_IMPORT_TYPES(Scene, Sampler, Medium, Emitter, EmitterPtr, BSDF, BSDFPtr)

    TransientPathIntegrator(const Properties &props) : Base(props) {
        m_filter_depth = props.int_("filter_depth", -1);
        // avoid the case filter_depth >= max_depth
        Assert(m_filter_depth < m_max_depth);
        m_discard_direct_paths = props.bool_("discard_direct_paths", false);
        // avoid the case m_discard_direct_paths && m_filter_depth > 0
        Assert(!m_discard_direct_paths || m_filter_depth <= 0);
        m_nlos_emitter_sampling = props.bool_("nlos_emitter_sampling", false);
        m_nlos_hidden_geometry_sampling =
            props.bool_("nlos_hidden_geometry_sampling", false);
        m_nlos_hidden_geometry_sampling_do_mis =
            props.bool_("nlos_hidden_geometry_sampling_do_mis", false) &&
            m_nlos_hidden_geometry_sampling;
    }

    void emitter_sample(
        const Scene *scene, Sampler *sampler, BSDFContext &ctx,
        SurfaceInteraction3f &si, Mask &active_e,
        std::vector<RadianceSample<Float, Spectrum>> &timed_samples_record,
        const BSDFPtr &bsdf, const Spectrum &throughput, const Float &path_opl,
        const Float &current_ior, const int &depth) const {
        auto [ds, emitter_val] = scene->sample_emitter_direction(
            si, sampler->next_2d(active_e), true, active_e);
        active_e &= neq(ds.pdf, 0.f);

        // Query the BSDF for that emitter-sampled direction
        Vector3f wo       = si.to_local(ds.d);
        Spectrum bsdf_val = select(
            neq(bsdf, nullptr),
            si.to_world_mueller(bsdf->eval(ctx, si, wo, active_e), -wo, si.wi),
            0.0f);

        // Determine density of sampling that same direction
        // using BSDF sampling
        Float bsdf_pdf =
            select(neq(bsdf, nullptr), bsdf->pdf(ctx, si, wo, active_e), 0.0f);

        Float mis = select(ds.delta, 1.f, mis_weight(ds.pdf, bsdf_pdf));

        Spectrum radiance(0.f);
        radiance[active_e] += mis * throughput * bsdf_val * emitter_val;

        if (depth == m_filter_depth || (!m_discard_direct_paths || depth >= 2))
            timed_samples_record.emplace_back(path_opl + ds.dist * current_ior,
                                              radiance, active_e);
    }

    void complete_path(
        const Scene *scene, Sampler *sampler, RayDifferential3f &ray,
        SurfaceInteraction3f &si, Mask &active,
        std::vector<RadianceSample<Float, Spectrum>> &timed_samples_record,
        const Float max_path_opl, const Point3f nlos_laser_target,
        const int initial_depth, const Spectrum &_throughput, Float path_opl,
        Float current_ior, bool do_mis, Float pdf_hg) const {

        EmitterPtr emitter = si.emitter(scene);
        Spectrum throughput(_throughput);
        // MIS weight for intersected emitters (set by prev. iteration)
        Float emission_weight(1.f);

        for (int depth = initial_depth;; ++depth) {

            path_opl += si.distance(ray) * current_ior;

            // Transient samples are not going to get stored anyway, discard
            if (path_opl > max_path_opl)
                break;

            // ---------------- Intersection with emitters ----------------

            if (any_or<true>(neq(emitter, nullptr))) {
                Spectrum radiance(0.f);
                radiance[active] +=
                    emission_weight * throughput * emitter->eval(si, active);
                Mask path_finished = active;
                if (depth == m_filter_depth ||
                    (!m_discard_direct_paths || depth >= 2))
                    timed_samples_record.emplace_back(path_opl, radiance,
                                                      path_finished);
            }

            active &= si.is_valid();

            // Stop if we've exceeded the number of requested bounces, or
            // if there are no more active lanes. Only do this latter check
            // in GPU mode when the number of requested bounces is infinite
            // since it causes a costly synchronization.
            if ((uint32_t) depth >= (uint32_t) m_max_depth ||
                ((!is_cuda_array_v<Float> || m_max_depth < 0) && none(active)))
                break;

            // --------------------- Emitter sampling ---------------------

            BSDFContext ctx;
            BSDFPtr bsdf = si.bsdf(ray);
            Mask active_e =
                active && has_flag(bsdf->flags(), BSDFFlags::Smooth);

            if (likely(any_or<true>(active_e))) {
                if (!m_nlos_emitter_sampling) {
                    // same emitter sampling as before
                    emitter_sample(scene, sampler, ctx, si, active_e,
                                   timed_samples_record, bsdf, throughput,
                                   path_opl, current_ior, depth);
                } else {
                    // nlos scenes only have one laser emitter - standard
                    // emitter sampling techniques do not apply as most
                    // directions do not emit any radiance, it needs to be very
                    // lucky to bsdf sample the exact point that the laser is
                    // illuminating
                    //
                    // this modifies the emitter sampling so instead of directly
                    // sampling the laser we sample (1) the point that the laser
                    // is illuminating and then (2) the laser

                    // 1. Obtain direction to NLOS illuminated point
                    //    and test visibility with ray_test
                    Vector3f d = nlos_laser_target - si.p;
                    Float dist = norm(d);
                    d /= dist;
                    Ray3f ray_bsdf(si.p, d,
                                   math::RayEpsilon<Float> *
                                       (1.f + hmax(abs(si.p))),
                                   dist * (1.f - math::ShadowEpsilon<Float>),
                                   si.time, si.wavelengths);
                    // active rays are those that did NOT intersect
                    active_e &= !scene->ray_test(ray_bsdf, active_e);

                    // 2. Evaluate BSDF to desired direction
                    if (any_or<true>(active_e)) {
                        Vector3f wo       = si.to_local(d);
                        Spectrum bsdf_val = bsdf->eval(ctx, si, wo, active_e);
                        bsdf_val = si.to_world_mueller(bsdf_val, -wo, si.wi);

                        ray_bsdf.maxt = math::Infinity<Float>;
                        SurfaceInteraction3f si_bsdf =
                            scene->ray_intersect(ray_bsdf, active_e);
                        active_e &= si_bsdf.is_valid();
                        active_e &= any_inner(depolarize<Spectrum>(bsdf_val) >
                                              math::Epsilon<Float>);
                        Vector3f wl = si_bsdf.to_local(-d);
                        active_e &= Frame3f::cos_theta(wl) > 0.f;

                        // NOTE(diego): as points are not randomly chosen,
                        // we need to account for d^2 and cos term because of
                        // the solid angle projection of si.p to
                        // nlos_laser_target.
                        // This is like a point light, but the extra cos term
                        // as it is not a point light that emits in all
                        // directions :^)
                        // The incident cos term at
                        // nlos_laser_target will be taken into account by
                        // f_emitter_sample's bsdf
                        bsdf_val *= sqr(rcp(dist)) * Frame3f::cos_theta(wl);

                        if (any(active_e)) {
                            BSDFPtr bsdf_next = si_bsdf.bsdf(ray_bsdf);

                            // 3. Combine nlos + emitter sampling
                            emitter_sample(scene, sampler, ctx, si_bsdf,
                                           active_e, timed_samples_record,
                                           bsdf_next, throughput * bsdf_val,
                                           path_opl + dist * current_ior,
                                           current_ior, depth + 1);
                        }
                    }
                }
            }

            // ----------------------- BSDF sampling ----------------------

            // Sample BSDF * cos(theta)
            auto [bs, bsdf_val] =
                bsdf->sample(ctx, si, sampler->next_1d(active),
                             sampler->next_2d(active), active);
            bsdf_val = si.to_world_mueller(bsdf_val, -bs.wo, si.wi);

            throughput = throughput * bsdf_val;
            if (do_mis && depth == initial_depth) {
                // MIS weight for BSDF samples
                Float mis = mis_weight(bs.pdf, pdf_hg);
                throughput *= mis;
            }
            active &= any(neq(depolarize(throughput), 0.f));
            if (none_or<false>(active))
                break;

            current_ior *= bs.eta;

            // Intersect the BSDF ray against the scene geometry
            ray                          = si.spawn_ray(si.to_world(bs.wo));
            SurfaceInteraction3f si_bsdf = scene->ray_intersect(ray, active);

            /* Determine probability of having sampled that same
               direction using emitter sampling. */
            // NOTE(diego): NLOS scenes (should) only have point light emitters,
            //              so they cannot be reached with BSDF sampling
            emitter = si_bsdf.emitter(scene, active);
            DirectionSample3f ds(si_bsdf, si);
            ds.object = emitter;

            if (any_or<true>(neq(emitter, nullptr))) {
                Float emitter_pdf =
                    select(neq(emitter, nullptr) &&
                               !has_flag(bs.sampled_type, BSDFFlags::Delta),
                           scene->pdf_emitter_direction(si, ds), 0.f);

                emission_weight = mis_weight(bs.pdf, emitter_pdf);
            }

            si = std::move(si_bsdf);
        }
    }

    void
    sample(const Scene *scene, Sampler *sampler, const RayDifferential3f &ray_,
           const Medium * /* medium */,
           std::vector<FloatSample<Float>> & /* aovs_record */,
           std::vector<RadianceSample<Float, Spectrum>> &timed_samples_record,
           Float max_path_opl, Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::SamplingIntegratorSample, active);

        RayDifferential3f ray = ray_;

        // Index of refraction of the current medium (used for OPL
        // calculation)
        Float current_ior(lookup_ior("air"));

        Spectrum throughput(1.f);

        // Time associated to a path, measured in optical path length
        // NOTE(diego): this assumes that the ray's `time` variable is
        // measured in OPL, we'll just have to believe now :-)
        Float path_opl(ray.time);

        // ------------------ Laser sampling preparation ------------------

        Point3f nlos_laser_target;
        if (m_nlos_emitter_sampling) {
            auto emitters = scene->emitters();
            if (unlikely(emitters.size() != 1)) {
                Throw("NLOS emitter sampling is not implemented for scenes "
                      "with more than one emitter.");
            }
            Transform4f trafo =
                emitters[0]->world_transform()->eval(ray.time, true);
            Vector3f laser_dir = trafo.transform_affine(Vector3f(0, 0, 1));
            Ray3f ray_laser(trafo.translation(), laser_dir, ray.time);
            nlos_laser_target = scene->ray_intersect(ray_laser, true).p;
        }

        // ---------------------- First intersection ----------------------

        SurfaceInteraction3f si = scene->ray_intersect(ray, active);
        active &= si.is_valid();

        if (unlikely((!is_cuda_array_v<Float> || m_max_depth < 0) &&
                     none(active))) {
            return;
        }

        // NLOS geometry sampling preparation
        int initial_depth = 1;
        BSDFContext ctx;
        BSDFPtr bsdf   = si.bsdf(ray);
        Mask active_hg = active && has_flag(bsdf->flags(), BSDFFlags::Smooth);
        Float hg_pdf   = 1.f;
        if (m_nlos_hidden_geometry_sampling) {
            PositionSample3f ps_hg = scene->sample_hidden_geometry_position(
                si, sampler->next_2d(active_hg), active_hg);
            active_hg &= neq(ps_hg.pdf, 0.f);

            Vector3f d = ps_hg.p - si.p;
            Float dist = norm(d);
            d /= dist;
            RayDifferential3f ray_hg(
                si.p, d, math::RayEpsilon<Float> * (1.f + hmax(abs(si.p))),
                math::Infinity<Float>, si.time, si.wavelengths);
            SurfaceInteraction3f si_hg =
                scene->ray_intersect(ray_hg, active_hg);
            active_hg &= si_hg.is_valid();

            Vector3f wo       = si.to_local(d);
            Spectrum bsdf_val = bsdf->eval(ctx, si, wo, active_hg);
            bsdf_val          = si.to_world_mueller(bsdf_val, -wo, si.wi);

            active_hg &= any_inner(depolarize<Spectrum>(bsdf_val) >
                                   math::Epsilon<Float>);
            Vector3f wg = si_hg.to_local(-d);
            bsdf_val *= sqr(rcp(dist)) * Frame3f::cos_theta(wg);

            if (any(active_hg)) {
                Float mis = 1.f;

                if (m_nlos_hidden_geometry_sampling_do_mis) {
                    // Determine density of sampling that same direction
                    // using BSDF sampling
                    Float bsdf_pdf = bsdf->pdf(ctx, si, wo, active_hg);
                    hg_pdf         = ps_hg.pdf;

                    mis = select(ps_hg.delta, 1.f,
                                 mis_weight(ps_hg.pdf, bsdf_pdf));
                }

                complete_path(scene, sampler, ray_hg, si_hg, active_hg,
                              timed_samples_record, max_path_opl,
                              nlos_laser_target, initial_depth + 1,
                              mis * bsdf_val * throughput / hg_pdf,
                              path_opl + si.distance(ray) * current_ior,
                              current_ior, false, hg_pdf);
            }
        }

        if (!m_nlos_hidden_geometry_sampling ||
            m_nlos_hidden_geometry_sampling_do_mis) {
            complete_path(scene, sampler, ray, si, active, timed_samples_record,
                          max_path_opl, nlos_laser_target, initial_depth,
                          throughput, path_opl, current_ior,
                          m_nlos_hidden_geometry_sampling_do_mis, hg_pdf);
        }
    }

    //! @}
    // =============================================================

    std::string to_string() const override {
        return tfm::format("TransientPathIntegrator[\n"
                           "  max_depth = %i,\n"
                           "  rr_depth = %i\n"
                           "]",
                           m_max_depth, m_rr_depth);
    }

    Float mis_weight(Float pdf_a, Float pdf_b) const {
        pdf_a *= pdf_a;
        pdf_b *= pdf_b;
        return select(pdf_a > 0.f, pdf_a / (pdf_a + pdf_b), 0.f);
    }

    MTS_DECLARE_CLASS()
private:
    int m_filter_depth;
    bool m_discard_direct_paths;
    bool m_nlos_emitter_sampling,
        m_nlos_hidden_geometry_sampling,
        m_nlos_hidden_geometry_sampling_do_mis;
};

MTS_IMPLEMENT_CLASS_VARIANT(TransientPathIntegrator, TransientMonteCarloIntegrator)
MTS_EXPORT_PLUGIN(TransientPathIntegrator, "Transient Path Tracer integrator");
NAMESPACE_END(mitsuba)
