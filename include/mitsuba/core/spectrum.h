#pragma once

#include <utility>

#include <enoki/matrix.h>
#include <mitsuba/core/fwd.h>
#include <mitsuba/core/logger.h>
#include <mitsuba/core/simd.h>
#include <mitsuba/core/traits.h>
#include <mitsuba/core/math.h>
#include <mitsuba/render/mueller.h>

NAMESPACE_BEGIN(mitsuba)

#if !defined(MTS_WAVELENGTH_MIN)
#  define MTS_WAVELENGTH_MIN 360.f
#endif

#if !defined(MTS_WAVELENGTH_MAX)
#  define MTS_WAVELENGTH_MAX 830.f
#endif

// =======================================================================
//! @{ \name Data types for RGB data
// =======================================================================

template <typename Value_, size_t Size_ = 3>
struct Color : enoki::StaticArrayImpl<Value_, Size_, false, Color<Value_, Size_>> {
    using Base = enoki::StaticArrayImpl<Value_, Size_, false, Color<Value_, Size_>>;

    /// Helper alias used to implement type promotion rules
    template <typename T> using ReplaceValue = Color<T, Size_>;

    using ArrayType = Color;
    using MaskType = enoki::Mask<Value_, Size_>;

    decltype(auto) r() const { return Base::x(); }
    decltype(auto) r() { return Base::x(); }

    decltype(auto) g() const { return Base::y(); }
    decltype(auto) g() { return Base::y(); }

    decltype(auto) b() const { return Base::z(); }
    decltype(auto) b() { return Base::z(); }

    decltype(auto) a() const { return Base::w(); }
    decltype(auto) a() { return Base::w(); }

    ENOKI_ARRAY_IMPORT(Base, Color)
};

//! @}
// =======================================================================

// =======================================================================
//! @{ \name Data types for RGB data WITH additional time information
// =======================================================================

#define MTS_SPEED_OF_LIGHT 299792458.0f

template <typename Value_, size_t Size_ = 3>
struct TimedColor : enoki::StaticArrayImpl<Value_, Size_ + 1, false, TimedColor<Value_, Size_>> {
    using Base = enoki::StaticArrayImpl<Value_, Size_ + 1, false, TimedColor<Value_, Size_>>;

    /// Helper alias used to implement type promotion rules
    template <typename T> using ReplaceValue = TimedColor<T, Size_>;

    using ArrayType = TimedColor;
    using MaskType  = enoki::Mask<Value_, Size_ + 1>;

    TimedColor() = default;

    template <typename Value2_ = Value_, size_t Size2_ = Size_>
    TimedColor(const Color<Value2_, Size2_> &color, const Value2_ &time = 0)
        : Base(color.x(), color.y(), color.z(), time) {}

    decltype(auto) r() const { return Base::x(); }
    decltype(auto) r() { return Base::x(); }

    decltype(auto) g() const { return Base::y(); }
    decltype(auto) g() { return Base::y(); }

    decltype(auto) b() const { return Base::z(); }
    decltype(auto) b() { return Base::z(); }

    template <typename Value2_ = Value_, size_t Size2_ = Size_>
    operator Color<Value2_, Size2_>() const { return this->color(); }

    template <typename Value2_ = Value_, size_t Size2_ = Size_>
    Color<Value2_, Size2_> color() const { return enoki::head<Size2_>(*this); }

    decltype(auto) opl() const { return Base::operator[](Size_); }
    decltype(auto) opl() { return Base::operator[](Size_); }

    decltype(auto) time() const { return this->opl() / MTS_SPEED_OF_LIGHT; }
    decltype(auto) time() { return this->opl() / MTS_SPEED_OF_LIGHT; }

    void add_opl(Value_ opl) { this->opl() += opl; }
    void add_opl(Value_ distance, Value_ medium_eta) { this->opl() += distance * medium_eta; }

    ENOKI_ARRAY_IMPORT(Base, TimedColor)
};

//! @}
// =======================================================================

// =======================================================================
//! @{ \name Data types for spectral quantities with sampled wavelengths
// =======================================================================

template <typename Value_, size_t Size_ = 4>
struct Spectrum : enoki::StaticArrayImpl<Value_, Size_, false, Spectrum<Value_, Size_>> {
    using Base = enoki::StaticArrayImpl<Value_, Size_, false, Spectrum<Value_, Size_>>;

    // Never allow matrix-vector products involving this type (important for polarized rendering)
    static constexpr bool IsVector = false;

    /// Helper alias used to implement type promotion rules
    template <typename T> using ReplaceValue = Spectrum<T, Size_>;

    using ArrayType = Spectrum;
    using MaskType = enoki::Mask<Value_, Size_>;

    ENOKI_ARRAY_IMPORT(Base, Spectrum)
};

/// Return the (1,1) entry of a Mueller matrix. Identity function for all other-types.
template <typename T> depolarize_t<T> depolarize(const T& spectrum) {
    if constexpr (is_polarized_v<T>) {
        // First entry of the Mueller matrix is the unpolarized spectrum
        return spectrum(0, 0);
    } else {
        return spectrum;
    }
}

/**
 * Turn a spectrum into a Mueller matrix representation that only has a non-zero
 * (1,1) entry. For all non-polarized modes, this is the identity function.
 */
template <typename T> auto unpolarized(const T &spectrum) {
    if constexpr (is_polarized_v<T>) {
        T result = zero<T>();
        result(0, 0) = spectrum(0, 0);
        return result;
    } else {
        return spectrum;
    }
}

//! @}
// =======================================================================

// =======================================================================
//! @{ \name Masking support for color and spectrum data types
// =======================================================================

template <typename Value_, size_t Size_>
struct Color<enoki::detail::MaskedArray<Value_>, Size_>
    : enoki::detail::MaskedArray<Color<Value_, Size_>> {
    using Base = enoki::detail::MaskedArray<Color<Value_, Size_>>;
    using Base::Base;
    using Base::operator=;
    Color(const Base &b) : Base(b) { }
};

template <typename Value_, size_t Size_>
struct TimedColor<enoki::detail::MaskedArray<Value_>, Size_>
    : enoki::detail::MaskedArray<TimedColor<Value_, Size_>> {
    using Base = enoki::detail::MaskedArray<TimedColor<Value_, Size_>>;
    using Base::Base;
    using Base::operator=;
    TimedColor(const Base &b) : Base(b) { }
};

template <typename Value_, size_t Size_>
struct Spectrum<enoki::detail::MaskedArray<Value_>, Size_>
    : enoki::detail::MaskedArray<Spectrum<Value_, Size_>> {
    using Base = enoki::detail::MaskedArray<Spectrum<Value_, Size_>>;
    using Base::Base;
    using Base::operator=;
    Spectrum(const Base &b) : Base(b) { }
};

//! @}
// =======================================================================

#define MTS_CIE_MIN           360.f
#define MTS_CIE_MAX           830.f
#define MTS_CIE_SAMPLES       95

/* Scaling the CIE curves by the following constant ensures that
   a unit-valued spectrum integrates to a luminance of 1.0 */
#define MTS_CIE_Y_NORMALIZATION float(1.0 / 106.7502593994140625)

/// Table with fits for \ref cie1931_xyz and \ref cie1931_y
extern MTS_EXPORT_CORE const float *cie1931_x_data;
extern MTS_EXPORT_CORE const float *cie1931_y_data;
extern MTS_EXPORT_CORE const float *cie1931_z_data;

/// Allocate GPU memory for the CIE 1931 tables
extern MTS_EXPORT_CORE void cie_alloc();

/**
 * \brief Evaluate the CIE 1931 XYZ color matching functions given a wavelength
 * in nanometers
 */
template <typename Float, typename Result = Color<Float, 3>>
Result cie1931_xyz(Float wavelength, mask_t<Float> active = true) {
    using Int32       = int32_array_t<Float>;
    using Float32     = float32_array_t<Float>;
    using ScalarFloat = scalar_t<Float>;

    Float t = (wavelength - (ScalarFloat) MTS_CIE_MIN) *
              ((MTS_CIE_SAMPLES - 1) /
               ((ScalarFloat) MTS_CIE_MAX - (ScalarFloat) MTS_CIE_MIN));

    active &= wavelength >= (ScalarFloat) MTS_CIE_MIN &&
              wavelength <= (ScalarFloat) MTS_CIE_MAX;

    Int32 i0 = clamp(Int32(t), zero<Int32>(), Int32(MTS_CIE_SAMPLES - 2)),
          i1 = i0 + 1;

    Float v0_x = (Float) gather<Float32>(cie1931_x_data, i0, active),
          v1_x = (Float) gather<Float32>(cie1931_x_data, i1, active),
          v0_y = (Float) gather<Float32>(cie1931_y_data, i0, active),
          v1_y = (Float) gather<Float32>(cie1931_y_data, i1, active),
          v0_z = (Float) gather<Float32>(cie1931_z_data, i0, active),
          v1_z = (Float) gather<Float32>(cie1931_z_data, i1, active);

    Float w1 = t - Float(i0),
          w0 = (ScalarFloat) 1.f - w1;

    return Result(fmadd(w0, v0_x, w1 * v1_x),
                  fmadd(w0, v0_y, w1 * v1_y),
                  fmadd(w0, v0_z, w1 * v1_z)) & mask_t<Result>(active);
}

/**
 * \brief Evaluate the CIE 1931 Y color matching function given a wavelength in
 * nanometers
 */

template <typename Float>
Float cie1931_y(Float wavelength, mask_t<Float> active = true) {
    using Int32       = int32_array_t<Float>;
    using Float32     = float32_array_t<Float>;
    using ScalarFloat = scalar_t<Float>;

    Float t = (wavelength - (ScalarFloat) MTS_CIE_MIN) *
              ((MTS_CIE_SAMPLES - 1) /
               ((ScalarFloat) MTS_CIE_MAX - (ScalarFloat) MTS_CIE_MIN));

    active &= wavelength >= (ScalarFloat) MTS_CIE_MIN &&
              wavelength <= (ScalarFloat) MTS_CIE_MAX;

    Int32 i0 = clamp(Int32(t), zero<Int32>(), Int32(MTS_CIE_SAMPLES - 2)),
          i1 = i0 + 1;

    Float v0 = (Float) gather<Float32>(cie1931_y_data, i0, active),
          v1 = (Float) gather<Float32>(cie1931_y_data, i1, active);

    Float w1 = t - Float(i0),
          w0 = (ScalarFloat) 1.f - w1;

    return select(active, fmadd(w0, v0, w1 * v1), 0.f);
}

/// Spectral responses to XYZ.
template <typename Float, size_t Size>
Color<Float, 3> spectrum_to_xyz(const Spectrum<Float, Size> &value,
                                const Spectrum<Float, Size> &wavelengths,
                                mask_t<Float> active = true) {
    Array<Spectrum<Float, Size>, 3> XYZ = cie1931_xyz(wavelengths, active);
    return { hmean(XYZ.x() * value),
             hmean(XYZ.y() * value),
             hmean(XYZ.z() * value) };
}

////////////////////////////////////////////////////////////////////////////////
// sRGB to XYZ conversion and vice-versa
//
// NOTE(diego): depending on what we are compiling, we might want to define
//   (a) two overloaded conversion functions (for TimedColor and Color) or
//   (b) only one, as enoki::vectorize cannot assume the function type when it
//       is overloaded
// Case (b) is required when creating python bindings, and in this case the
// MTS_VARIANT_TRANSIENT flag should be set as we are compiling for a specific
// Spectrum variant
// The other case (a) is required for librender, which (I think) is
// Spectrum-independent, does not set MTS_VARIANT_TRANSIENT, does not vectorize
// and thus the function can be overloaded without problems

template <typename Float>
MTS_INLINE Color<Float, 3> _srgb_to_xyz(const Color<Float, 3> &rgb) {
    using ScalarMatrix3f = enoki::Matrix<scalar_t<Float>, 3>;
    const ScalarMatrix3f M(0.412453f, 0.357580f, 0.180423f,
                           0.212671f, 0.715160f, 0.072169f,
                           0.019334f, 0.119193f, 0.950227f);
    return M * rgb;
}

#if !defined(MTS_VARIANT_TRANSIENT) || MTS_VARIANT_TRANSIENT

/// Convert ITU-R Rec. BT.709 linear RGB to XYZ tristimulus values (Timed Remix)
template <typename Float>
TimedColor<Float, 3> srgb_to_xyz(const TimedColor<Float, 3> &timed, mask_t<Float> /*active*/ = true) {
    const Color<Float, 3> xyz = _srgb_to_xyz(timed.color());
    return TimedColor<Float, 3>(xyz, timed.time());
}

#endif

#if !defined(MTS_VARIANT_TRANSIENT) || !MTS_VARIANT_TRANSIENT

/// Convert ITU-R Rec. BT.709 linear RGB to XYZ tristimulus values
template <typename Float>
Color<Float, 3> srgb_to_xyz(const Color<Float, 3> &rgb, mask_t<Float> /*active*/ = true) {
    return _srgb_to_xyz(rgb);
}

#endif

template <typename Float>
MTS_INLINE Color<Float, 3> _xyz_to_srgb(const Color<Float, 3> &rgb) {
    using ScalarMatrix3f = enoki::Matrix<scalar_t<Float>, 3>;
    const ScalarMatrix3f M(3.240479f, -1.537150f, -0.498535f,
                          -0.969256f,  1.875991f,  0.041556f,
                           0.055648f, -0.204043f,  1.057311f);
    return M * rgb;
}

#if !defined(MTS_VARIANT_TRANSIENT) || MTS_VARIANT_TRANSIENT

/// Convert ITU-R Rec. BT.709 linear RGB to XYZ tristimulus values (Timed Remix)
template <typename Float>
TimedColor<Float, 3> xyz_to_srgb(const TimedColor<Float, 3> &timed, mask_t<Float> /*active*/ = true) {
    const Color<Float, 3> srgb = _xyz_to_srgb(timed.color());
    return TimedColor<Float, 3>(srgb, timed.time());
}

#endif

#if !defined(MTS_VARIANT_TRANSIENT) || !MTS_VARIANT_TRANSIENT

/// Convert XYZ tristimulus values to ITU-R Rec. BT.709 linear RGB
template <typename Float>
Color<Float, 3> xyz_to_srgb(const Color<Float, 3> &rgb, mask_t<Float> /*active*/ = true) {
    return _xyz_to_srgb(rgb);
}

#endif
////////////////////////////////////////////////////////////////////////////////

template <typename Float, size_t Size>
Float luminance(const Spectrum<Float, Size> &value,
                const Spectrum<Float, Size> &wavelengths,
                mask_t<Float> active = true) {
    return hmean(cie1931_y(wavelengths, active) * value);
}

template <typename Float> Float luminance(const Color<Float, 3> &c) {
    return c[0] * 0.212671f + c[1] * 0.715160f + c[2] * 0.072169f;
}

template <typename Value>
std::pair<Value, Value> sample_uniform_spectrum(const Value &sample) {
    return { sample * (MTS_CIE_MAX - MTS_CIE_MIN) + MTS_CIE_MIN,
             MTS_CIE_MAX - MTS_CIE_MIN };
}

template <typename Value>
Value pdf_uniform_spectrum(const Value & /* wavelength */) {
    return Value(1.f / (MTS_WAVELENGTH_MAX - MTS_WAVELENGTH_MIN));
}

/**
 * Importance sample a "importance spectrum" that concentrates the computation
 * on wavelengths that are relevant for rendering of RGB data
 *
 * Based on "An Improved Technique for Full Spectral Rendering"
 * by Radziszewski, Boryczko, and Alda
 *
 * Returns a tuple with the sampled wavelength and inverse PDF
 */
template <typename Value>
std::pair<Value, Value> sample_rgb_spectrum(const Value &sample) {
    if constexpr (MTS_WAVELENGTH_MIN == 360.f && MTS_WAVELENGTH_MAX == 830.f) {
        Value wavelengths =
            538.f - atanh(0.8569106254698279f -
                          1.8275019724092267f * sample) * 138.88888888888889f;

        Value tmp    = cosh(0.0072f * (wavelengths - 538.f));
        Value weight = 253.82f * tmp * tmp;

        return { wavelengths, weight };
    } else {
        // Fall back to uniform sampling for other wavelength ranges
        return sample_uniform_spectrum(sample);
    }
}

/**
 * PDF for the \ref sample_rgb_spectrum strategy.
 * It is valid to call this function for a single wavelength (Float), a set
 * of wavelengths (Spectrumf), a packet of wavelengths (SpectrumfP), etc. In all
 * cases, the PDF is returned per wavelength.
 */
template <typename Value> Value pdf_rgb_spectrum(const Value &wavelengths) {
    if constexpr (MTS_WAVELENGTH_MIN == 360.f && MTS_WAVELENGTH_MAX == 830.f) {
        Value tmp = sech(0.0072f * (wavelengths - 538.f));
        return select(wavelengths >= MTS_WAVELENGTH_MIN && wavelengths <= MTS_WAVELENGTH_MAX,
                      0.003939804229326285f * tmp * tmp, zero<Value>());
    } else {
        return pdf_uniform_spectrum(wavelengths);
    }
}

/// Helper function to sample a wavelength (and a weight) given a random number
template <typename Float, typename Spectrum>
std::pair<wavelength_t<Spectrum>, Spectrum> sample_wavelength(Float sample) {
    if constexpr (!is_spectral_v<Spectrum>) {
        ENOKI_MARK_USED(sample);
        // Note: wavelengths should not be used when rendering in RGB mode.
        return { std::numeric_limits<scalar_t<Float>>::quiet_NaN(), 1.f };
    } else {
        auto wav_sample = math::sample_shifted<wavelength_t<Spectrum>>(sample);
        return sample_rgb_spectrum(wav_sample);
    }
}

template <typename Scalar>
MTS_EXPORT_CORE void spectrum_from_file(const std::string &filename,
                                        std::vector<Scalar> &wavelengths,
                                        std::vector<Scalar> &values);

template <typename Scalar>
MTS_EXPORT_CORE Color<Scalar, 3> spectrum_to_rgb(const std::vector<Scalar> &wavelengths,
                                                 const std::vector<Scalar> &values,
                                                 bool bounded = true);

NAMESPACE_END(mitsuba)
