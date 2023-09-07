#include <mitsuba/render/sensor.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/transform.h>
#include <mitsuba/core/bbox.h>

NAMESPACE_BEGIN(mitsuba)

// This code was fixed thanks to the following PR:
// https://github.com/mitsuba-renderer/mitsuba2/pull/279/files
// Thanks!

/**!

.. _sensor-perspective:

Ortographic pinhole camera (:monosp:`perspective`)
--------------------------------------------------

.. pluginparameters::

 * - to_world
   - |transform|
   - Specifies an optional camera-to-world transformation.
     (Default: none (i.e. camera space = world space))
 * - fov
   - |float|
   - Denotes the camera's field of view in degrees---must be between 0 and 180,
     excluding the extremes. Alternatively, it is also possible to specify a
     field of view using the :monosp:`focal_length` parameter.
 * - focal_length
   - |string|
   - Denotes the camera's focal length specified using *35mm* film
     equivalent units. Alternatively, it is also possible to specify a field of
     view using the :monosp:`fov` parameter. See the main description for further
     details. (Default: :monosp:`50mm`)
 * - fov_axis
   - |string|
   - When the parameter :monosp:`fov` is given (and only then), this parameter further specifies
     the image axis, to which it applies.

     1. :monosp:`x`: :monosp:`fov` maps to the :monosp:`x`-axis in screen space.
     2. :monosp:`y`: :monosp:`fov` maps to the :monosp:`y`-axis in screen space.
     3. :monosp:`diagonal`: :monosp:`fov` maps to the screen diagonal.
     4. :monosp:`smaller`: :monosp:`fov` maps to the smaller dimension
        (e.g. :monosp:`x` when :monosp:`width` < :monosp:`height`)
     5. :monosp:`larger`: :monosp:`fov` maps to the larger dimension
        (e.g. :monosp:`y` when :monosp:`width` < :monosp:`height`)

     The default is :monosp:`x`.
 * - near_clip, far_clip
   - |float|
   - Distance to the near/far clip planes. (Default: :monosp:`near_clip=1e-2` (i.e. :monosp:`0.01`)
     and :monosp:`far_clip=1e4` (i.e. :monosp:`10000`))
 * - principal_point_offset_x, principal_point_offset_y
   - |float|
   - Specifies the position of the camera's principal point relative to the center of the film.

.. subfigstart::
.. subfigure:: ../../resources/data/docs/images/render/sensor_perspective.jpg
   :caption: The material test ball viewed through a perspective pinhole camera. (:monosp:`fov=28`)
.. subfigure:: ../../resources/data/docs/images/render/sensor_perspective_large_fov.jpg
   :caption: The material test ball viewed through a perspective pinhole camera. (:monosp:`fov=40`)
.. subfigend::
   :label: fig-perspective

This plugin implements a simple idealizied perspective camera model, which
has an infinitely small aperture. This creates an infinite depth of field,
i.e. no optical blurring occurs.

By default, the camera's field of view is specified using a 35mm film
equivalent focal length, which is first converted into a diagonal field
of view and subsequently applied to the camera. This assumes that
the film's aspect ratio matches that of 35mm film (1.5:1), though the
parameter still behaves intuitively when this is not the case.
Alternatively, it is also possible to specify a field of view in degrees
along a given axis (see the :monosp:`fov` and :monosp:`fov_axis` parameters).

The exact camera position and orientation is most easily expressed using the
:monosp:`lookat` tag, i.e.:

.. code-block:: xml

    <sensor type="ortographic">
        <transform name="to_world">
            <!-- Move and rotate the camera so that looks from (1, 1, 1) to (1, 2, 1)
                and the direction (0, 0, 1) points "up" in the output image -->
            <lookat origin="1, 1, 1" target="1, 2, 1" up="0, 0, 1"/>
        </transform>
    </sensor>

 */

template <typename Float, typename Spectrum>
class OrtographicCamera final : public ProjectiveCamera<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(ProjectiveCamera, m_world_transform, m_needs_sample_3,
                    m_film, m_sampler, m_resolution, m_shutter_open,
                    m_shutter_open_time, m_near_clip, m_far_clip)
    MTS_IMPORT_TYPES()

    // =============================================================
    //! @{ \name Constructors
    // =============================================================


    Transform<Point<Float, 4>>
    ortographic_projection(const Vector<int, 2> &film_size,
                       const Vector<int, 2> &crop_size,
                       const Vector<int, 2> &crop_offset,
                       Float near_clip, Float far_clip) {

        using Vector2f = Vector<Float, 2>;
        using Vector3f = Vector<Float, 3>;
        using Transform4f = Transform<Point<Float, 4>>;

        Vector2f film_size_f = Vector2f(film_size),
                rel_size    = Vector2f(crop_size) / film_size_f,
                rel_offset  = Vector2f(crop_offset) / film_size_f;

        Float aspect = film_size_f.x() / film_size_f.y();

        /**
         * These do the following (in reverse order):
         *
         * 1. Create transform from camera space to [-1,1]x[-1,1]x[0,1] clip
         *    coordinates (not taking account of the aspect ratio yet)
         *
         * 2+3. Translate and scale to shift the clip coordinates into the
         *    range from zero to one, and take the aspect ratio into account.
         *
         * 4+5. Translate and scale the coordinates once more to account
         *     for a cropping window (if there is any)
         */
        return Transform4f::scale(
                Vector3f(1.f / rel_size.x(), 1.f / rel_size.y(), 1.f)) *
            Transform4f::translate(
                Vector3f(-rel_offset.x(), -rel_offset.y(), 0.f)) *
            Transform4f::scale(Vector3f(-0.5f, -0.5f * aspect, 1.f)) *
            Transform4f::translate(Vector3f(-1.f, -1.f / aspect, 0.f)) *
            Transform4f::orthographic(near_clip, far_clip);
    }

    OrtographicCamera(const Properties &props) : Base(props) {
        update_camera_transforms();
    }

    void update_camera_transforms() {
        m_camera_to_sample = ortographic_projection(
            m_film->size(), m_film->crop_size(), m_film->crop_offset(),
            m_near_clip, m_far_clip);

        m_sample_to_camera = m_camera_to_sample.inverse();

        // Position differentials on the near plane
        m_dx = m_sample_to_camera * ScalarPoint3f(1.f / m_resolution.x(), 0.f, 0.f) -
               m_sample_to_camera * ScalarPoint3f(0.f);
        m_dy = m_sample_to_camera * ScalarPoint3f(0.f, 1.f / m_resolution.y(), 0.f)
             - m_sample_to_camera * ScalarPoint3f(0.f);
        
        std::cout << "dx " << m_dx << " dy " << m_dy << std::endl;

        m_normalization = 1.f / m_image_rect.volume();
    }

    //! @}
    // =============================================================

    // =============================================================
    //! @{ \name Sampling methods (Sensor interface)
    // =============================================================

    std::pair<Ray3f, Spectrum> sample_ray(Float time, Float wavelength_sample,
                                          const Point2f &position_sample,
                                          const Point2f & /*aperture_sample*/,
                                          Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::EndpointSampleRay, active);

        auto [wavelengths, wav_weight] = sample_wavelength<Float, Spectrum>(wavelength_sample);
        Ray3f ray;
        ray.time = time;
        ray.wavelengths = wavelengths;

        // Compute the sample position on the near plane (local camera space).
        Point3f near_p = m_sample_to_camera *
                         Point3f(position_sample.x(), position_sample.y(), 0.f);

        ray.mint = m_near_clip;
        ray.maxt = m_far_clip;

        auto trafo = m_world_transform->eval(ray.time, active);
        ray.o = trafo.transform_affine(Point3f(near_p.x(), near_p.y(), 0.f));
        ray.d = normalize(trafo * Vector3f(0, 0, 1));
        ray.update();

        return std::make_pair(ray, wav_weight);
    }

    std::pair<RayDifferential3f, Spectrum>
    sample_ray_differential(Float time, Float wavelength_sample, const Point2f &position_sample,
                            const Point2f & /*aperture_sample*/, Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::EndpointSampleRay, active);

        auto [wavelengths, wav_weight] = sample_wavelength<Float, Spectrum>(wavelength_sample);
        RayDifferential3f ray;
        ray.time = time;
        ray.wavelengths = wavelengths;

        // Compute the sample position on the near plane (local camera space).
        Point3f near_p = m_sample_to_camera *
                         Point3f(position_sample.x(), position_sample.y(), 0.f);

        ray.mint = m_near_clip;
        ray.maxt = m_far_clip;

        auto trafo = m_world_transform->eval(ray.time, active);
        ray.o = trafo.transform_affine(near_p);
        ray.d = normalize(trafo * Vector3f(0, 0, 1));
        ray.update();

        ray.o_x = ray.o_y = ray.o;

        ray.d_x = trafo * normalize(Vector3f(near_p) + m_dx);
        ray.d_y = trafo * normalize(Vector3f(near_p) + m_dy);
        ray.has_differentials = true;

        return std::make_pair(ray, wav_weight);
    }

    ScalarBoundingBox3f bbox() const override {
        return m_world_transform->translation_bounds();
    }

    // =============================================================

    void traverse(TraversalCallback *callback) override {
        Base::traverse(callback);
        // TODO x_fov
    }

    void parameters_changed(const std::vector<std::string> &keys) override {
        Base::parameters_changed(keys);
        // TODO
    }

    std::string to_string() const override {
        using string::indent;

        std::ostringstream oss;
        oss << "OrtographicCamera[" << std::endl
            << "  near_clip = " << m_near_clip << "," << std::endl
            << "  far_clip = " << m_far_clip << "," << std::endl
            << "  film = " << indent(m_film) << "," << std::endl
            << "  sampler = " << indent(m_sampler) << "," << std::endl
            << "  resolution = " << m_resolution << "," << std::endl
            << "  shutter_open = " << m_shutter_open << "," << std::endl
            << "  shutter_open_time = " << m_shutter_open_time << "," << std::endl
            << "  world_transform = " << indent(m_world_transform) << std::endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    ScalarTransform4f m_camera_to_sample;
    ScalarTransform4f m_sample_to_camera;
    ScalarBoundingBox2f m_image_rect;
    ScalarFloat m_normalization;
    ScalarVector3f m_dx, m_dy;
};

MTS_IMPLEMENT_CLASS_VARIANT(OrtographicCamera, ProjectiveCamera)
MTS_EXPORT_PLUGIN(OrtographicCamera, "Ortographic Camera");
NAMESPACE_END(mitsuba)
