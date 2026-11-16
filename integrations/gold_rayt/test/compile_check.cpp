// Check of MeasuredHpBRDF.hpp against the GoLD_rayt headers: compiles for the RGB build and for the spectral build
// (/DRAYT_SPECTRAL_MODE), links with GoLD_rayt's GGX.cpp, loads the table and evaluates a few configurations.
// Run from integrations/gold_rayt:  test\check_rgb.exe [TABLE.hpbpca]   (default ../../out/pca/Aluminum.hpbpca)
#include "Materials/MeasuredHpBRDF.hpp"
#include <chrono>
#include <cstdio>

using namespace rayt;

static void show(const char* what, const Spectrum& s) {
    std::printf("  %-34s", what);
    for (int i = 0; i < pol::kN; ++i) std::printf(" %9.4f", double(s[i]));
    std::printf("\n");
}

int main(int argc, char** argv) {
    std::string err;
    auto t0 = std::chrono::steady_clock::now();
    auto mat = MeasuredHpBRDF::load(argc > 1 ? argv[1] : "../../out/pca/Aluminum.hpbpca", 0.1, &err);
    if (!mat) { std::printf("load failed: %s\n", err.c_str()); return 1; }
    const double loadS = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("loaded in %.2f s (%s)\n", loadS, kRenderModeName);
    const Vector3 a = mat->aovAlbedo(SurfaceInteraction{});
    std::printf("  directional albedo (RGB AOV)         %.3f %.3f %.3f\n", a.x, a.y, a.z);

    SurfaceInteraction rec; rec.n = rec.gn = Vector3(0, 0, 1); rec.p = Point3(0);
    SampledWavelengths lambda;   // 550 nm x4 by default
#if defined(RAYT_SPECTRAL_MODE)
    lambda.lambda = { 450, 550, 650, 800 };
#endif
    const double th = 45 * constants::PI / 180;
    const Vector3 wo(std::sin(th), 0, std::cos(th));            // to the viewer
    const Vector3 wiSpec(-std::sin(th), 0, std::cos(th));       // mirror direction
    const Vector3 wiOff = glm::normalize(Vector3(-0.3, 0.5, 0.8));
    show("f_r at the specular direction", mat->eval(rec, wo, wiSpec, lambda));
    show("f_r off specular", mat->eval(rec, wo, wiOff, lambda));
    show("f_r reciprocity (wo <-> wi)", mat->eval(rec, wiOff, wo, lambda));
    std::printf("  pdf(spec) %.4f  pdf(off) %.4f\n", mat->pdf(rec, wo, wiSpec, lambda), mat->pdf(rec, wo, wiOff, lambda));

    // polarization: the Mueller matrix (M/m00) at the specular configuration, component 0
    pol::ScatterEvent ev;
    if (mat->scatterMueller(rec, wo, wiSpec, lambda, ev)) {
        std::printf("  scatterMueller (spec, 45 deg): s = (%.3f %.3f %.3f)\n", ev.s.x, ev.s.y, ev.s.z);
        for (int r = 0; r < 4; ++r) {
            std::printf("    ");
            for (int c = 0; c < 4; ++c) std::printf(" %+.3f", double(ev.m.e[r][c][0]));
            std::printf("\n");
        }
    }
    // sampling: average of f cos / pdf over random samples should be the directional albedo (around 0.8 for aluminum)
    {
        double sum = 0; int n = 0, fail = 0; uint32_t s = 12345u;
        auto rnd = [&]() { s = s * 1664525u + 1013904223u; return Real(s >> 8) / Real(1u << 24); };
        for (int i = 0; i < 20000; ++i) {
            auto bs = mat->sample(rec, wo, Point2(rnd(), rnd()), lambda);
            if (!bs) { ++fail; continue; }
            sum += double(bs->f[0]) * double(std::max(Real(0), bs->wi.z)) / double(bs->pdf); ++n;
        }
        std::printf("  sampling: E[f cos / pdf] = %.3f over %d samples (%d rejected)\n", sum / 20000, n, fail);
    }
    // speed
    {
        t0 = std::chrono::steady_clock::now(); double acc = 0; const int N = 200000;
        for (int i = 0; i < N; ++i) {
            const Real b = 0.1 + 1.4 * ((i / 7) % 199) / 199.0, p = 6.283 * (i % 977) / 977.0;
            const Vector3 wi(std::sin(b) * std::cos(p), std::sin(b) * std::sin(p), std::cos(b));
            acc += double(mat->eval(rec, wo, wi, lambda)[0]);
        }
        const double ns = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() / N * 1e9;
        std::printf("  eval: %.0f ns per call (checksum %.1f)\n", ns, acc);
    }
    return 0;
}
