// Cross-check driver: reads direction samples, writes decoded values for comparison with the Python reference (tests/).
//   test_decoder TABLE.hpbpca samples.txt out.txt
// samples.txt lines:  lx ly lz vx vy vz channel
// out.txt lines    :  phi_d theta_d theta_h m00  M00..M33 (16, table frame)  Md00..Md33 (16, default frames)  x_in(3) x_out(3)
#include "hpbrdf_pca.hpp"
#include <chrono>
#include <fstream>
#include <iostream>

int main(int argc, char **argv) {
    if (argc < 4) { std::cerr << "usage: test_decoder TABLE.hpbpca samples.txt out.txt\n"; return 1; }
    hpbrdf::PCATable T; std::string err;
    auto t0 = std::chrono::steady_clock::now();
    if (!T.load(argv[1], &err)) { std::cerr << "load failed: " << err << "\n"; return 1; }
    double load_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::cerr << "loaded " << argv[1] << "  grid " << T.n_phi << "x" << T.n_thd << "x" << T.n_thh << "  channels " << T.n_ch
              << (T.rgb ? " (rgb)" : " (nm)") << "  k " << T.k1 << "+" << T.k2 << (T.polarized ? "  polarized" : "  unpolarized")
              << "  coef " << T.coef.size() * 4 / 1048576 << " MB  in " << load_s << " s\n";

    std::ifstream in(argv[2]); std::ofstream out(argv[3]); out.precision(9);
    float lx, ly, lz, vx, vy, vz, ch; size_t n = 0;
    std::vector<float> scratch(T.k());
    while (in >> lx >> ly >> lz >> vx >> vy >> vz >> ch) {
        hpbrdf::Vec3 L(lx, ly, lz), V(vx, vy, vz);
        auto r = T.eval(L, V, ch, true, scratch.data());
        auto Md = T.polarized ? hpbrdf::rotate_mueller_basis(r.M, -L, r.g.x_in, hpbrdf::default_basis(-L), V, r.g.x_out, hpbrdf::default_basis(V)) : r.M;
        out << r.g.phi_d << ' ' << r.g.theta_d << ' ' << r.g.theta_h << ' ' << r.m00;
        for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) out << ' ' << r.M.m[i][j];
        for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) out << ' ' << Md.m[i][j];
        out << ' ' << r.g.x_in.x << ' ' << r.g.x_in.y << ' ' << r.g.x_in.z << ' ' << r.g.x_out.x << ' ' << r.g.x_out.y << ' ' << r.g.x_out.z << '\n';
        ++n;
    }
    const int N = 200000; float acc = 0; const float chx = T.rgb ? 1.f : 550.f;
    t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < N; ++i) {
        float a = 0.3f + 1.2f * (i % 1000) / 1000.f, b = 0.1f + 1.4f * ((i / 1000) % 200) / 200.f, p = 6.283f * (i % 977) / 977.f;
        hpbrdf::Vec3 L(std::sin(a), 0, std::cos(a)), V(std::sin(b) * std::cos(p), std::sin(b) * std::sin(p), std::cos(b));
        acc += T.eval(L, V, chx, false, scratch.data()).m00;
    }
    double t_unpol = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() / N;
    t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < N; ++i) {
        float a = 0.3f + 1.2f * (i % 1000) / 1000.f, b = 0.1f + 1.4f * ((i / 1000) % 200) / 200.f, p = 6.283f * (i % 977) / 977.f;
        hpbrdf::Vec3 L(std::sin(a), 0, std::cos(a)), V(std::sin(b) * std::cos(p), std::sin(b) * std::sin(p), std::cos(b));
        acc += T.eval(L, V, chx, true, scratch.data()).M.m[0][0];
    }
    double t_pol = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() / N;
    std::cerr << "wrote " << n << " samples;  eval time: unpolarized " << t_unpol * 1e9 << " ns, polarized " << t_pol * 1e9
              << " ns  (checksum " << acc << ")\n";
    return 0;
}
