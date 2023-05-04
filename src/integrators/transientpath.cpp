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
        Assert(m_filter_depth == -1 || m_max_depth == -1 ||
               m_filter_depth < m_max_depth);
        m_discard_direct_paths = props.bool_("discard_direct_paths", false);
        // avoid the case m_discard_direct_paths && m_filter_depth == 1
        Assert(!m_discard_direct_paths || m_filter_depth != 1);
        m_nlos_laser_sampling = props.bool_("nlos_laser_sampling", false);
        m_nlos_hidden_geometry_sampling =
            props.bool_("nlos_hidden_geometry_sampling", false);
        m_nlos_hidden_geometry_sampling_do_rroulette =
            props.bool_("nlos_hidden_geometry_sampling_do_rroulette", false) &&
            m_nlos_hidden_geometry_sampling;
        m_nlos_hidden_geometry_sampling_includes_relay_wall =
            props.bool_("nlos_hidden_geometry_sampling_includes_relay_wall", true) &&
            m_nlos_hidden_geometry_sampling;
    }

    bool should_store_sample(const int depth) const {
        if (m_discard_direct_paths && depth == 1)
            return false;
        if (m_filter_depth != -1)
            return depth == m_filter_depth;
        return true;
    }

    void emitter_nee_sample(
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

        if (should_store_sample(depth))
            timed_samples_record.emplace_back(path_opl + ds.dist * current_ior,
                                              radiance, active_e);
    }

    void emitter_laser_sample(
        const Scene *scene, Sampler *sampler, BSDFContext &ctx,
        SurfaceInteraction3f &si, Mask &active_e,
        std::vector<RadianceSample<Float, Spectrum>> &timed_samples_record,
        const BSDFPtr &bsdf, const Spectrum &throughput, const Float &path_opl,
        const Float &current_ior, const int &depth) const {
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
        Vector3f d = m_nlos_laser_target - si.p;
        Float dist = norm(d);
        d /= dist;
        Ray3f ray_bsdf(
            si.p, d, math::RayEpsilon<Float> * (1.f + hmax(abs(si.p))),
            dist * (1.f - math::ShadowEpsilon<Float>), si.time, si.wavelengths);
        // active rays are those that did NOT intersect
        active_e &= !scene->ray_test(ray_bsdf, active_e);

        // 2. Evaluate BSDF to desired direction
        if (any_or<true>(active_e)) {
            Vector3f wo       = si.to_local(d);
            Spectrum bsdf_val = bsdf->eval(ctx, si, wo, active_e);
            bsdf_val          = si.to_world_mueller(bsdf_val, -wo, si.wi);

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
            // exists as it is not a point light that emits in all
            // directions :^)
            // The incident cos term at
            // nlos_laser_target will be taken into account by
            // emitter_nee_sample's bsdf
            bsdf_val *= sqr(rcp(dist)) * Frame3f::cos_theta(wl);

            if (any(active_e)) {
                BSDFPtr bsdf_next = si_bsdf.bsdf(ray_bsdf);

                // 3. Combine nlos + emitter sampling
                emitter_nee_sample(
                    scene, sampler, ctx, si_bsdf, active_e,
                    timed_samples_record, bsdf_next, throughput * bsdf_val,
                    path_opl + dist * current_ior, current_ior, depth + 1);
            }
        }
    }

    std::pair<BSDFSample3f, Spectrum>
    hidden_geometry_sample(const Scene *scene, const BSDFPtr &bsdf,
                           const BSDFContext &ctx,
                           const SurfaceInteraction3f &si, Float /* sample1 */,
                           const Point2f &sample2, Mask active) const {

        PositionSample3f ps_hg =
            scene->sample_hidden_geometry_position(si, sample2, active);
        active &= neq(ps_hg.pdf, 0.f);

        Vector3f d = ps_hg.p - si.p;
        Float dist = norm(d);
        d /= dist;
        RayDifferential3f ray_hg(
            si.p, d, math::RayEpsilon<Float> * (1.f + hmax(abs(si.p))),
            math::Infinity<Float>, si.time, si.wavelengths);
        SurfaceInteraction3f si_hg = scene->ray_intersect(ray_hg, active);
        active &= si_hg.is_valid();

        SurfaceInteraction3f si_test;
        Point3f p_test         = si_hg.p;
        Mask active_test       = active;
        uint num_intersections = 0;
        while (num_intersections == 0 || any(active_test)) {
            num_intersections++;
            Ray3f ray_test(p_test, d, si.time);
            ray_test.mint = math::RayEpsilon<Float> * (1.f + hmax(abs(si.p)));
            si_test       = scene->ray_intersect(ray_test, active_test);
            active_test &= si_test.is_valid();
            if (unlikely(num_intersections > 100)) {
                // Some rays get stuck in an infinite loop, creating a cycle
                // of p_test points. Just ignore those cases.
                active_test = false;
                break;
            }
            p_test = si_test.p;
        }

        Vector3f wo       = si.to_local(d);
        Spectrum bsdf_val = bsdf->eval(ctx, si, wo, active);
        bsdf_val          = si.to_world_mueller(bsdf_val, -wo, si.wi);

        Vector3f wg = si_hg.to_local(-d);
        Float travel_dist = norm(si_hg.p - si.p);
        bsdf_val *= sqr(rcp(travel_dist)) * Frame3f::cos_theta(wg);
        // discard low travel dist, produces high variance
        // the integrator will use bsdf sampling instead
        active &= travel_dist > 1.f;  
        active &=
            !any_inner(depolarize<Spectrum>(bsdf_val) < math::Epsilon<Float>);
            
        BSDFSample3f bs      = zero<BSDFSample3f>();
        bs.wo                = wo;
        bs.pdf               = ps_hg.pdf * num_intersections;
        bs.eta               = 1.f;
        bs.sampled_type      = +BSDFFlags::DeltaReflection;
        bs.sampled_component = 0;

        return { bs, select(active && bs.pdf > math::Epsilon<Float>, bsdf_val, 0.f) };
    }

    void
    sample(const Scene *scene, Sampler *sampler, const RayDifferential3f &ray_,
           const Medium * /* medium */,
           std::vector<FloatSample<Float>> & /* aovs_record */,
           std::vector<RadianceSample<Float, Spectrum>> &timed_samples_record,
           Float max_path_opl, Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::SamplingIntegratorSample, active);

        RayDifferential3f ray = ray_;

        // Index of refraction of the current medium (used for OPL calculation)
        Float current_ior(lookup_ior("air"));

        // MIS weight for intersected emitters (set by prev. iteration)
        Float emission_weight(1.f);

        Spectrum throughput(1.f);

        // Time associated to a path, measured in optical path length
        // NOTE(diego): this assumes that the ray's `time` variable is
        // measured in OPL, we'll just have to believe now :-)
        Float path_opl(ray.time);

        // ---------------------- First intersection ----------------------

        SurfaceInteraction3f si = scene->ray_intersect(ray, active);
        EmitterPtr emitter      = si.emitter(scene);

        if (unlikely((!is_cuda_array_v<Float> || m_max_depth < 0) &&
                     none(active))) {
            return;
        }

        for (int depth = 1;; ++depth) {

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
                if (should_store_sample(depth))
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
                if (m_nlos_laser_sampling) {
                    // laser sampling (nee modified for nlos scenes)
                    emitter_laser_sample(scene, sampler, ctx, si, active_e,
                                         timed_samples_record, bsdf, throughput,
                                         path_opl, current_ior, depth);
                } else {
                    // standard nee
                    emitter_nee_sample(scene, sampler, ctx, si, active_e,
                                       timed_samples_record, bsdf, throughput,
                                       path_opl, current_ior, depth);
                }
            }

            // ------------- Geometry sampling (BSDF or area) -------------

            bool do_bsdf_sampling;
            Float pdf_bsdf_method;
            if (m_nlos_hidden_geometry_sampling) {
                if (m_nlos_hidden_geometry_sampling_do_rroulette) {
                    Float BSDF_SAMPLING_Q = 0.5f;
                    do_bsdf_sampling =
                        sampler->next_1d(active) < BSDF_SAMPLING_Q;
                    pdf_bsdf_method = do_bsdf_sampling
                                          ? BSDF_SAMPLING_Q
                                          : (1.0f - BSDF_SAMPLING_Q);
                } else {
                    do_bsdf_sampling = false;
                    pdf_bsdf_method  = 1.f;
                }
            } else {
                do_bsdf_sampling = true;
                pdf_bsdf_method  = 1.f;
            }

            BSDFSample3f bs;
            Spectrum bsdf_val;
            if (!do_bsdf_sampling) {
                std::tie(bs, bsdf_val) = hidden_geometry_sample(
                    scene, bsdf, ctx, si, sampler->next_1d(active),
                    sampler->next_2d(active), active);
            }
            if (do_bsdf_sampling || all_inner(depolarize<Spectrum>(bsdf_val) <
                                              math::Epsilon<Float>)) {
                std::tie(bs, bsdf_val) =
                    bsdf->sample(ctx, si, sampler->next_1d(active),
                                 sampler->next_2d(active), active);
                bsdf_val = si.to_world_mueller(bsdf_val, -bs.wo, si.wi);
            }

            throughput = throughput * bsdf_val / pdf_bsdf_method;
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

    //! @}
    // =============================================================

    std::string to_string() const override {
        return tfm::format("TransientPathIntegrator[\n"
                           "  max_depth = %i,\n"
                           "  rr_depth = %i,\n"
                           "  filter_depth = %i,\n"
                           "  discard_direct_paths = %s,\n"
                           "  nlos_laser_sampling = %s,\n"
                           "  nlos_hidden_goemetry_sampling = %s,\n"
                           "  nlos_hidden_goemetry_sampling_do_rroulette = %s\n"
                           "]",
                           m_max_depth, m_rr_depth, m_filter_depth,
                           m_discard_direct_paths, m_nlos_laser_sampling,
                           m_nlos_hidden_geometry_sampling,
                           m_nlos_hidden_geometry_sampling_do_rroulette);
    }

    void prepare_integrator(const Scene* scene) override {
        if (m_nlos_laser_sampling) {
            auto emitters = scene->emitters();
            Assert(emitters[0]->world_transform()->size() == 1);
            if (unlikely(emitters.size() != 1)) {
                Throw("NLOS laser sampling is not implemented for scenes "
                      "with more than one emitter.");
            }
            Transform4f trafo =
                emitters[0]->world_transform()->eval(0.0f, true);
            Vector3f laser_dir = trafo.transform_affine(Vector3f(0, 0, 1));
            Ray3f ray_laser(trafo.translation(), laser_dir, 0.0f);
            m_nlos_laser_target = scene->ray_intersect(ray_laser, true).p;
        }
    }

public:
    bool hidden_geometry_sampling_includes_relay_wall() const override {
        return m_nlos_hidden_geometry_sampling_includes_relay_wall;
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
    bool m_nlos_laser_sampling,
        m_nlos_hidden_geometry_sampling,
        m_nlos_hidden_geometry_sampling_do_rroulette,
        m_nlos_hidden_geometry_sampling_includes_relay_wall;
    Point3f m_nlos_laser_target;
};

MTS_IMPLEMENT_CLASS_VARIANT(TransientPathIntegrator, TransientMonteCarloIntegrator)
MTS_EXPORT_PLUGIN(TransientPathIntegrator, "Transient Path Tracer integrator");
NAMESPACE_END(mitsuba)
