#include <mitsuba/core/properties.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>

NAMESPACE_BEGIN(mitsuba)

/**!

.. _bsdf-diffuse:

Smooth diffuse material (:monosp:`diffuse`)
-------------------------------------------

.. pluginparameters::

 * - reflectance
   - |spectrum| or |texture|
   - Specifies the diffuse albedo of the material (Default: 0.5)

The smooth diffuse material (also referred to as *Lambertian*)
represents an ideally diffuse material with a user-specified amount of
reflectance. Any received illumination is scattered so that the surface
looks the same independently of the direction of observation.

.. subfigstart::
.. subfigure:: ../../resources/data/docs/images/render/bsdf_diffuse_plain.jpg
   :caption: Homogeneous reflectance
.. subfigure:: ../../resources/data/docs/images/render/bsdf_diffuse_textured.jpg
   :caption: Textured reflectance
.. subfigend::
   :label: fig-diffuse

Apart from a homogeneous reflectance value, the plugin can also accept
a nested or referenced texture map to be used as the source of reflectance
information, which is then mapped onto the shape based on its UV
parameterization. When no parameters are specified, the model uses the default
of 50% reflectance.

Note that this material is one-sided---that is, observed from the
back side, it will be completely black. If this is undesirable,
consider using the :ref:`twosided <bsdf-twosided>` BRDF adapter plugin.
The following XML snippet describes a diffuse material,
whose reflectance is specified as an sRGB color:

.. code-block:: xml
    :name: diffuse-srgb

    <bsdf type="diffuse">
        <rgb name="reflectance" value="0.2, 0.25, 0.7"/>
    </bsdf>

Alternatively, the reflectance can be textured:

.. code-block:: xml
    :name: diffuse-texture

    <bsdf type="diffuse">
        <texture type="bitmap" name="reflectance">
            <string name="filename" value="wood.jpg"/>
        </texture>
    </bsdf>

*/
template <typename Float, typename Spectrum>
class SmoothDiffuseFakePolarized final : public BSDF<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(BSDF, m_flags, m_components)
    MTS_IMPORT_TYPES(Texture)

    SmoothDiffuseFakePolarized(const Properties &props) : Base(props) {
        m_reflectance = props.texture<Texture>("reflectance", .5f);
        m_flags = BSDFFlags::DiffuseReflection | BSDFFlags::FrontSide;
        m_components.push_back(m_flags);
    }

    std::pair<BSDFSample3f, Spectrum> sample(const BSDFContext &ctx,
                                             const SurfaceInteraction3f &si,
                                             Float /* sample1 */,
                                             const Point2f &sample2,
                                             Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);

        Float cos_theta_i = Frame3f::cos_theta(si.wi);
        BSDFSample3f bs = zero<BSDFSample3f>();

        active &= cos_theta_i > 0.f;
        if (unlikely(none_or<false>(active) ||
                     !ctx.is_enabled(BSDFFlags::DiffuseReflection)))
            return { bs, 0.f };

        bs.wo = warp::square_to_cosine_hemisphere(sample2);
        bs.pdf = warp::square_to_cosine_hemisphere_pdf(bs.wo);
        bs.eta = 1.f;
        bs.sampled_type = +BSDFFlags::DiffuseReflection;
        bs.sampled_component = 0;

        Normal3f m = si.n;

        UnpolarizedSpectrum value = m_reflectance->eval(si, active);

        if constexpr (is_polarized_v<Spectrum>) {
            /* Due to the coordinate system rotations for polarization-aware
               pBSDFs below we need to know the propagation direction of light.
               In the following, light arrives along `-wo_hat` and leaves along
               `+wi_hat`. */
            Vector3f wo_hat = ctx.mode == TransportMode::Radiance ? bs.wo : si.wi,
                     wi_hat = ctx.mode == TransportMode::Radiance ? si.wi : bs.wo;

            // Mueller matrix for specular reflection.
            Spectrum F = mueller::preserver<Float>();

            /* The Stokes reference frame vector of this matrix lies perpendicular
               to the plane of reflection. */
            Vector3f s_axis_in = normalize(cross(m, -wo_hat)),
                     s_axis_out = normalize(cross(m, wi_hat));

            /* Rotate in/out reference vector of F s.t. it aligns with the implicit
               Stokes bases of -wo_hat & wi_hat. */
            F = mueller::rotate_mueller_basis(F,
                                              -wo_hat, s_axis_in, mueller::stokes_basis(-wo_hat),
                                               wi_hat, s_axis_out, mueller::stokes_basis(wi_hat));
            return { bs, select(active && bs.pdf > 0.f, F * value, 0.f) };
        } else {
            return { bs, select(active && bs.pdf > 0.f, unpolarized<Spectrum>(value), 0.f) };
        }
        
    }

    Spectrum eval(const BSDFContext &ctx, const SurfaceInteraction3f &si,
                  const Vector3f &wo, Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        if (!ctx.is_enabled(BSDFFlags::DiffuseReflection))
            return 0.f;

        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);

        active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

        UnpolarizedSpectrum value =
            m_reflectance->eval(si, active) * math::InvPi<Float> * cos_theta_o;

        Normal3f m = si.n;

        if constexpr (is_polarized_v<Spectrum>) {
            /* Due to the coordinate system rotations for polarization-aware
               pBSDFs below we need to know the propagation direction of light.
               In the following, light arrives along `-wo_hat` and leaves along
               `+wi_hat`. */
            Vector3f wo_hat = ctx.mode == TransportMode::Radiance ? wo : si.wi,
                     wi_hat = ctx.mode == TransportMode::Radiance ? si.wi : wo;

            // Mueller matrix for specular reflection.
            Spectrum F = mueller::preserver<Float>();

            /* The Stokes reference frame vector of this matrix lies perpendicular
               to the plane of reflection. */
            Vector3f s_axis_in = normalize(cross(m, -wo_hat)),
                     s_axis_out = normalize(cross(m, wi_hat));

            /* Rotate in/out reference vector of F s.t. it aligns with the implicit
               Stokes bases of -wo_hat & wi_hat. */
            F = mueller::rotate_mueller_basis(F,
                                              -wo_hat, s_axis_in, mueller::stokes_basis(-wo_hat),
                                               wi_hat, s_axis_out, mueller::stokes_basis(wi_hat));
            return (F * value) & active;
        } else {
            return select(active, unpolarized<Spectrum>(value), 0.f);
        }
    }

    Float pdf(const BSDFContext &ctx, const SurfaceInteraction3f &si,
              const Vector3f &wo, Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        if (!ctx.is_enabled(BSDFFlags::DiffuseReflection))
            return 0.f;

        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);

        Float pdf = warp::square_to_cosine_hemisphere_pdf(wo);

        return select(cos_theta_i > 0.f && cos_theta_o > 0.f, pdf, 0.f);
    }

    void traverse(TraversalCallback *callback) override {
        callback->put_object("reflectance", m_reflectance.get());
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "SmoothDiffuseFakePolarized[" << std::endl
            << "  reflectance = " << string::indent(m_reflectance) << std::endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    ref<Texture> m_reflectance;
};

MTS_IMPLEMENT_CLASS_VARIANT(SmoothDiffuseFakePolarized, BSDF)
MTS_EXPORT_PLUGIN(SmoothDiffuseFakePolarized, "Smooth diffuse fake-polarized material")
NAMESPACE_END(mitsuba)
