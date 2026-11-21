"""
brdfkit command line.

  brdfkit info      TABLE
  brdfkit compress  TABLE OUT.npz [--k1 8] [--k2 32] [--dtype float16] [--step 1] [--floor FLOOR.npz] [--floor-keep 0.02] [--floor-max-gain 4] [--exchange] [--cosine light|camera]
  brdfkit eval      TABLE OUT.npz [--slices 0,30,60] [--sweep] [--floor FLOOR.npz] [--floor-keep 0.02] [--floor-max-gain 4]
  brdfkit floor     TABLE --out FLOOR.npz [--phi-step 12] [--theta-h-min 40] [--theta-d 25,75] [--stat p50]
  brdfkit edit      MODEL.npz --model MERL_PCA.npz --out EDITED.npz [--pc 1:-0.5 ...] [--toward nickel:0.5 ...] [--k 16] [--keep-albedo 45]
                    [--sym median|mean|none] [--calib CAL.npz [--calib-keep-level DEG]]   (no shift: just the cleaned model)
  brdfkit swap      SHAPE.npz (--source MODEL.npz | --nk NK.yml [--film-nm 20] [--n-film 1.65]) --out OUT.npz [--parts S,N,level] [--cutoff 30]
                    (shape of one material, colour / polarization / reflectance of another, docs/SWAP.md)
  brdfkit gloss     MODEL.npz [--out OUT.npz (--hwhm DEG | --alpha A | --kernel-alpha A)] [--kernel ggx|beckmann] [--ref-theta-d 30]
                    [--keep-albedo DEG] [--sym median|mean|none]   (roughen the lobe; without --out: report its width, docs/GLOSS.md)
  brdfkit calibrate GLASS_MODEL.npz --out CAL.npz [--cutoff 10] [--smooth 3] [--check MODEL.npz --refidx NK.yml [--film-nm 20] [--n-film 1.65]]
                    (radiometric error of the table from black glass, docs/CALIBRATION.md)
  brdfkit factor    TABLE --out MODEL.npz [--smooth 0,5] [--measured 20,80] [--floor FLOOR.npz] [--exchange] [--cosine light|camera] [--sym none|median|mean]
                    (A x S x N model; TABLE may then be MODEL.npz)
  brdfkit export    OUT.npz OUT.hpbpca [--dtype float16]
  brdfkit lu-chipman OUT.npz --refidx Al_Rakic.yml [Al_McPeak.yml ...] --out DIR [--raw TABLE] [--n-film 1.65]
  brdfkit merl-pca  DIR_OR_FILES... --out DIR [--stride 2] [--k 16] [--mapping log|nielsen] [--eps E] [--no-plot]
  brdfkit merl-project TABLE --model MERL_PCA.npz [MERL_PCA.npz ...] --out DIR [--k 16] [--sym full|none] [--reference MERL.binary] [--recompute] [--no-plot]
                    [--floor FLOOR.npz] [--floor-keep 0.02]
  brdfkit merl-subsets MERL_GRID.npz MERL_DIR --out DIR [--all-model MERL_PCA.npz] [--random 30] [--seed 0] [--ablate]
  brdfkit merl-compare RUN_DIR [RUN_DIR ...]      (summary of several merl-project outputs)
  brdfkit quirks    TABLE|REDUCED.npz --out REDUCED.npz [--dielectric | --nk NK.yml [--film-nm 20]] [--name NAME]
                    (inpainting, floor, radiometric, spectral, phi_d and polarization quirks of a measured table, docs/QUIRKS.md)
  brdfkit to-mitsuba TABLE OUT.pbsdf [--floor FLOOR.npz] [--exchange] [--cosine light|camera] [--calib CAL.npz] [--phi-step 1] [--bands LO,HI] [--band-step 1]
                    (a table or factor model, corrected, as a tensor_file for Mitsuba 3 `measured_polarized`, docs/MITSUBA.md)
  brdfkit umov      REDUCED.npz --reference GLASS_REDUCED.npz [--dark-max 540] [--region mid|far] [--n 1.505] [--out REPORT.txt]
                    (surface and body of a coloured dielectric split by polarization, on quirks reductions, docs/UMOV.md)
  brdfkit slots     TABLE|REDUCED.npz --out REDUCED.npz [--bands 550,700] [--phi-step 6] [--theta-h-min 15] [--name NAME]
                    (which of the table's two directions is the light, from the polarization of a diffuse body, docs/SLOTS.md)

TABLE: .hpbrdf / .pbsdf (Mitsuba tensor_file), .mat (KAIST pBRDF 2020), .binary (MERL), .npz (factor model of `brdfkit factor`).
"""
import argparse
import os


def _floor_args(p):
    p.add_argument('--floor', default=None, help='floor file of `brdfkit floor` (e.g. measured on black glass) to subtract from m00')
    p.add_argument('--floor-keep', type=float, default=0.02, help='never take m00 below this fraction of its measured value')
    p.add_argument('--floor-max-gain', type=float, default=4.0, help='cap on the rise of m_ij/m00 caused by the subtraction (0 = no cap)')


def _slot_args(p):
    p.add_argument('--exchange', action='store_true', help='exchange the two directions (hpBRDF tables are filed at phi_d + 180 deg, docs/SLOTS.md)')
    p.add_argument('--cosine', choices=('light', 'camera'), default=None,
                   help='divide by the cosine of this direction, as labelled after --exchange (hpBRDF: light)')
    p.add_argument('--min-cos', type=float, default=0.1, help='clamp of that cosine')


def _with_floor(table, a):
    from .analysis import clean as CL
    return CL.clean(table, a.floor or None, getattr(a, 'exchange', False), getattr(a, 'cosine', None), None,
                    a.floor_keep, a.floor_max_gain or None, getattr(a, 'min_cos', 0.1))


def main(argv=None):
    ap = argparse.ArgumentParser(prog='brdfkit', description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('info'); p.add_argument('table')
    p = sub.add_parser('compress'); p.add_argument('table'); p.add_argument('output')
    p.add_argument('--k1', type=int, default=8); p.add_argument('--k2', type=int, default=32)
    p.add_argument('--dtype', default='float16'); p.add_argument('--step', type=int, default=1, help='phi_d stride of the covariance pass')
    _floor_args(p); _slot_args(p)
    p = sub.add_parser('eval'); p.add_argument('table'); p.add_argument('npz')
    p.add_argument('--slices', default=''); p.add_argument('--sweep', action='store_true')
    _floor_args(p); _slot_args(p)
    p = sub.add_parser('floor', help='estimate the instrument floor (unpolarized offset) far from the specular lobe'); p.add_argument('table')
    p.add_argument('--out', required=True); p.add_argument('--phi-step', type=int, default=12); p.add_argument('--theta-h-min', type=float, default=40.0)
    p.add_argument('--theta-d', default='25,75', help='measured theta_d range in degrees'); p.add_argument('--stat', choices=('p10', 'p25', 'p50', 'p90'), default='p50')
    p = sub.add_parser('export'); p.add_argument('npz'); p.add_argument('output'); p.add_argument('--dtype', default='float16')
    p = sub.add_parser('lu-chipman'); p.add_argument('npz'); p.add_argument('--refidx', nargs='+', required=True)
    p.add_argument('--out', required=True); p.add_argument('--raw', default=None); p.add_argument('--n-film', type=float, default=1.65)
    p = sub.add_parser('factor', help='fit m00 = A(angles) S(theta_d, lambda), m_ij/m00 = N(theta_d, lambda)'); p.add_argument('table')
    p.add_argument('--out', required=True); p.add_argument('--smooth', default='0', help='smoothing widths of S over theta_d [deg], comma separated')
    p.add_argument('--measured', default='20,80', help='measured theta_d range [deg]; S and N are extended from its edges'); _floor_args(p); _slot_args(p)
    p.add_argument('--sym', choices=('none', 'median', 'mean'), default='none',
                   help='symmetrize A over the four equivalent phi_d after the fit (the reported errors are those of the fit)')
    p = sub.add_parser('edit', help='move the grey lobe A of a factor model along a MERL PCA (log mapping)'); p.add_argument('factor')
    p.add_argument('--model', required=True, help='merl_pca.npz (log mapping)'); p.add_argument('--out', required=True)
    p.add_argument('--pc', action='append', default=[], help='COMPONENT:T, shift by T times the spread of the PCA materials (repeatable)')
    p.add_argument('--toward', action='append', default=[], help='MATERIAL:T, move the fraction T of the way to that material (repeatable)')
    p.add_argument('--k', type=int, default=16); p.add_argument('--keep-albedo', type=float, default=None, help='rescale to keep the albedo at this angle of incidence [deg]')
    p.add_argument('--sym', choices=('median', 'mean', 'none'), default='median',
                   help='symmetrize A over the four equivalent phi_d before the edit (removes the one-sided arms of hpBRDF lobes)')
    p.add_argument('--calib', default=None, help='calibration file of `brdfkit calibrate`: divide A by its K(theta_d) first')
    p.add_argument('--calib-keep-level', type=float, default=None, help='with --calib: keep the level at this theta_d [deg] (fix the trend only)')
    p = sub.add_parser('swap', help='angular shape of one factor model with the colour, polarization and reflectance of another')
    p.add_argument('shape', help='factor model whose lobe shape is kept'); p.add_argument('--out', required=True)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--source', help='factor model giving the swapped parts'); g.add_argument('--nk', help='n, k table (.yml): Fresnel as the source')
    p.add_argument('--film-nm', type=float, default=0.0); p.add_argument('--n-film', type=float, default=1.65)
    p.add_argument('--parts', default='S,N,level', help='any of S (colour), N (polarization), level (reflectance), comma separated')
    p.add_argument('--cutoff', type=float, default=30.0, help='theta_h [deg] of the lobe energy that defines the level')
    p.add_argument('--smooth', type=float, default=3.0, help='smoothing of the level along theta_d [deg]')
    p = sub.add_parser('gloss', help='roughen the lobe of a factor model (convolve its microfacet slopes) to a target width')
    p.add_argument('factor'); p.add_argument('--out', default=None, help='output model; without it the width of the lobe is reported')
    g = p.add_mutually_exclusive_group()
    g.add_argument('--hwhm', type=float, default=None, help='target half width at half maximum of D(theta_h) [deg] at --ref-theta-d')
    g.add_argument('--alpha', type=float, default=None, help='target alpha_eq: the GGX roughness with that half width')
    g.add_argument('--kernel-alpha', type=float, default=None, help='convolve with a kernel of this roughness instead of a target')
    p.add_argument('--kernel', choices=('ggx', 'beckmann'), default='ggx', help='slope density of the added roughness')
    p.add_argument('--ref-theta-d', type=float, default=30.0, help='theta_d [deg] at which the width is measured')
    p.add_argument('--keep-albedo', type=float, default=None, help='rescale to keep the albedo at this angle of incidence [deg]')
    p.add_argument('--sym', choices=('median', 'mean', 'none'), default='median', help='symmetrize A over phi_d first (as edit)')
    p = sub.add_parser('calibrate', help='radiometric and polarimetric check against black glass (Fresnel of a dielectric)')
    p.add_argument('glass', help='factor model (.npz) of black glass, floor removed'); p.add_argument('--out', required=True)
    p.add_argument('--cutoff', type=float, default=10.0, help='theta_h [deg] up to which the lobe energy is integrated')
    p.add_argument('--smooth', type=float, default=3.0, help='smoothing of K along theta_d [deg]')
    p.add_argument('--check', default=None, help='factor model of a metal to check the correction on (albedo / Fresnel)')
    p.add_argument('--refidx', default=None, help='n, k table (.yml) of that metal'); p.add_argument('--film-nm', type=float, default=0.0)
    p.add_argument('--n-film', type=float, default=1.65); p.add_argument('--no-plot', action='store_true')
    p = sub.add_parser('merl-pca', help='inter-material PCA of a set of MERL tables'); p.add_argument('paths', nargs='+')
    p.add_argument('--out', required=True); p.add_argument('--stride', default='2', help='grid stride, one int or phi_d,theta_d,theta_h')
    p.add_argument('--k', type=int, default=16, help='components kept in the basis'); p.add_argument('--mapping', choices=('log', 'nielsen'), default='log')
    p.add_argument('--eps', type=float, default=None); p.add_argument('--no-plot', action='store_true')
    p = sub.add_parser('merl-project', help='project a table (e.g. hpBRDF aluminium) onto saved MERL PCA models'); p.add_argument('table')
    p.add_argument('--model', nargs='+', required=True, help='merl_pca.npz written by merl-pca'); p.add_argument('--out', required=True)
    p.add_argument('--k', type=int, default=16); p.add_argument('--sym', choices=('full', 'none'), default='full', help='average the phi_d-equivalent slices (full) or not')
    p.add_argument('--reference', default=None, help='MERL .binary of a comparison material'); p.add_argument('--recompute', action='store_true')
    p.add_argument('--no-plot', action='store_true'); _floor_args(p)
    p = sub.add_parser('merl-subsets', help='metal-only versus all-material PCA basis for a table on the MERL grid')
    p.add_argument('grid', help='<out>/<name>_merl_grid_full.npz written by merl-project'); p.add_argument('merl_dir'); p.add_argument('--out', required=True)
    p.add_argument('--all-model', default=None, help='saved 100-material merl_pca.npz (log), for leave-one-out numbers')
    p.add_argument('--random', type=int, default=30, help='random control subsets per size'); p.add_argument('--seed', type=int, default=0)
    p.add_argument('--ablate', action='store_true', help='refit without each material in turn (100 fits)')
    p = sub.add_parser('merl-compare', help='summarize several merl-project output folders'); p.add_argument('run_dirs', nargs='+')
    p = sub.add_parser('quirks', help='check a measured table for inpainting, floor, radiometric, spectral, phi_d and polarization quirks')
    p.add_argument('table', help='measured table, or the reduction (.npz) a previous run wrote'); p.add_argument('--out', required=True, help='reduction .npz (and .txt report)')
    p.add_argument('--name', default=None, help='name in the report (default: the file name)')
    g = p.add_mutually_exclusive_group()
    g.add_argument('--dielectric', action='store_true', help='fit the index from the specular polarization and compare the lobe energy with its Fresnel')
    g.add_argument('--nk', default=None, help='n, k table (.yml) of a metal for the Fresnel comparison')
    p.add_argument('--film-nm', type=float, default=0.0, help='dielectric film on the metal (nm)')
    p = sub.add_parser('slots', help='which of the two directions is the light: polarization of the diffuse body against the two angles')
    p.add_argument('table', help='measured polarized table, or the reduction (.npz) a previous run wrote'); p.add_argument('--out', required=True)
    p.add_argument('--bands', default='550,700', help='wavelength range [nm] to average (where the body is bright)')
    p.add_argument('--phi-step', type=int, default=6); p.add_argument('--theta-h-min', type=float, default=15.0, help='stay this far from the specular lobe [deg]')
    p.add_argument('--theta-d', default='22,78', help='measured theta_d range [deg]'); p.add_argument('--name', default=None)
    _floor_args(p); _slot_args(p)                                 # check a corrected view: after --exchange the verdict is 'as labelled'
    p = sub.add_parser('to-mitsuba', help='write a table or factor model, corrected, as a tensor_file for Mitsuba 3 measured_polarized')
    p.add_argument('table'); p.add_argument('output', help='.pbsdf / .hpbrdf'); _floor_args(p); _slot_args(p)
    p.add_argument('--calib', default=None, help='calibration of `brdfkit calibrate` (same --exchange / --cosine): divide by its K(theta_d)')
    p.add_argument('--phi-step', type=int, default=1, help='keep every n-th phi_d bin (must divide the number of intervals)')
    p.add_argument('--bands', default=None, help='wavelength range LO,HI [nm] to keep'); p.add_argument('--band-step', type=int, default=1)
    p = sub.add_parser('umov', help='surface and body of a coloured dielectric split by polarization, from quirks reductions')
    p.add_argument('reduced', help='reduction (.npz) of `brdfkit quirks` of the sample'); p.add_argument('--reference', required=True,
                   help='reduction of a dark sample (black glass): the spectrum of the floor')
    p.add_argument('--dark-max', type=float, default=540.0, help='the body is taken as black below this wavelength [nm]')
    p.add_argument('--region', choices=('mid', 'far'), default='mid'); p.add_argument('--n', type=float, default=None, help='index, for the Fresnel columns')
    p.add_argument('--out', default=None, help='report (.txt)'); p.add_argument('--name', default=None)
    a = ap.parse_args(argv)

    if a.cmd == 'info':
        from .io import load
        t = load(a.table); print(t.describe()); print(t.meta)
    elif a.cmd == 'compress':
        from .io import load
        from .compress import compress
        compress(_with_floor(load(a.table), a), a.output, a.k1, a.k2, a.dtype, a.step)
    elif a.cmd == 'eval':
        from .io import load
        from .compress import evaluate
        slices = [int(s) for s in a.slices.split(',')] if a.slices else None
        evaluate(_with_floor(load(a.table), a), a.npz, slices, a.sweep)
    elif a.cmd == 'floor':
        from .io import load
        from .analysis import floor as FL
        td = [float(s) for s in a.theta_d.split(',')]
        fl = FL.estimate(load(a.table), a.phi_step, a.theta_h_min, (td[0], td[1]), a.stat)
        FL.save(a.out, fl); FL.report(fl); print('wrote', a.out)
    elif a.cmd == 'export':
        from .export import write
        write(a.npz, a.output, a.dtype)
    elif a.cmd == 'factor':
        from .io import load
        from .analysis import factor as F
        sm = [float(x) for x in a.smooth.split(',')]; me = [float(x) for x in a.measured.split(',')]
        models = F.fit(_with_floor(load(a.table), a), (me[0], me[1]), sm)
        base = a.out[:-4] if a.out.endswith('.npz') else a.out
        for m in models:
            path = a.out if len(models) == 1 else f'{base}_s{m["smooth_deg"]:g}.npz'
            F.save(path, F.symmetrize(m, a.sym)); F.report(m); print('wrote', path, f'(A symmetrized over phi_d: {a.sym})')
    elif a.cmd == 'edit':
        from .analysis import factor as F, latent_edit as LE
        from .analysis.merl_pca import MaterialPCA
        pcs = [(int(x.split(':')[0]), float(x.split(':')[1])) for x in a.pc]
        tw = [(x.rsplit(':', 1)[0], float(x.rsplit(':', 1)[1])) for x in a.toward]
        base = F.load(a.factor)
        if a.calib:
            from .analysis import calibrate as CA
            base = CA.apply(base, CA.load(a.calib), a.calib_keep_level)
            print(f'radiometric calibration {a.calib} applied' + (f' (level kept at theta_d {a.calib_keep_level:g} deg)' if a.calib_keep_level is not None else ''))
        m, info = LE.edit(base, MaterialPCA.load(a.model), pcs, tw, a.k, a.keep_albedo, a.sym)
        F.save(a.out, m); LE.report(info); print('wrote', a.out)
    elif a.cmd == 'swap':
        from .analysis import factor as F, swap as SW
        shape = F.load(a.shape)
        src = F.load(a.source) if a.source else SW.fresnel_source(shape, a.nk, a.film_nm, a.n_film)
        m, info = SW.swap(shape, src, [p.strip() for p in a.parts.split(',')], a.cutoff, a.smooth)
        F.save(a.out, m); SW.report(info); print('wrote', a.out)
    elif a.cmd == 'gloss':
        from .analysis import factor as F, gloss as GL
        m = F.load(a.factor)
        if a.hwhm is None and a.alpha is None and a.kernel_alpha is None:
            if a.out:
                raise SystemExit('--out needs a target: --hwhm, --alpha or --kernel-alpha')
            m = F.symmetrize(m, a.sym); D = GL.microfacet_d(m)[0]
            w = {t: GL.width(m, t, D) for t in GL.REPORT_THETA_D}
            print('lobe hwhm [deg] at theta_d ' + ' / '.join(f'{t:g}' for t in w) + ': ' + ' / '.join(f'{v:.2f}' for v in w.values()))
            wr = GL.width(m, a.ref_theta_d, D)
            print(f'at theta_d {a.ref_theta_d:g} deg: hwhm {wr:.2f} deg, alpha_eq {GL.alpha_of_hwhm(wr):.3f}')
        else:
            if not a.out:
                raise SystemExit('a target needs --out')
            m, info = GL.set_gloss(m, a.hwhm, a.alpha, a.kernel_alpha, a.kernel, True, a.ref_theta_d, a.sym, a.keep_albedo)
            F.save(a.out, m); GL.report(info); print('wrote', a.out)
    elif a.cmd == 'calibrate':
        from .analysis import factor as F, calibrate as CA
        cal = CA.estimate(F.load(a.glass), a.cutoff, a.smooth)
        CA.save(a.out, cal)
        rows = None
        if a.check:
            if not a.refidx:
                raise SystemExit('--check needs --refidx (n, k of that metal)')
            rows = CA.check(F.load(a.check), cal, a.refidx, a.film_nm, a.n_film)
        lines = []
        CA.report(cal, rows, os.path.basename(a.check) if a.check else '', say=lambda s: (print(s), lines.append(s)))
        base = a.out[:-4] if a.out.endswith('.npz') else a.out
        open(base + '.txt', 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
        if not a.no_plot:
            CA.plot(cal, base + '.png', rows)
        print('wrote', a.out, base + '.txt', '' if a.no_plot else base + '.png')
    elif a.cmd == 'lu-chipman':
        from .analysis.lu_chipman import run_report
        run_report(a.npz, [(p, os.path.splitext(os.path.basename(p))[0]) for p in a.refidx], a.out, raw=a.raw, n_film=a.n_film)
    elif a.cmd == 'merl-pca':
        from .analysis.merl_pca import run_report
        stride = [int(s) for s in a.stride.split(',')]
        run_report(a.paths, a.out, stride[0] if len(stride) == 1 else stride, a.k, a.mapping, a.eps, plot=not a.no_plot)
    elif a.cmd == 'merl-project':
        from .analysis.merl_project import run_report
        run_report(a.table, a.model, a.out, a.k, a.sym, a.reference, a.recompute, plot=not a.no_plot, floor=a.floor, floor_keep=a.floor_keep)
    elif a.cmd == 'merl-subsets':
        from .analysis.merl_subsets import run_study
        run_study(a.grid, a.merl_dir, a.out, a.all_model, a.random, a.seed, a.ablate)

    elif a.cmd == 'merl-compare':
        from .analysis.merl_project import compare_runs
        compare_runs(a.run_dirs)
    elif a.cmd == 'quirks':
        from .analysis import quirks as QK
        if a.table.endswith('.npz'):
            red = QK.load(a.table)
        else:
            from .io import load
            red = QK.reduce(load(a.table)); QK.save(a.out, red); red = QK.load(a.out)
        optics = ('dielectric',) if a.dielectric else (('metal', a.nk, a.film_nm) if a.nk else None)
        lines = []
        QK.report(a.name or os.path.splitext(os.path.basename(a.table))[0], QK.measure(red, optics), say=lambda s: (print(s), lines.append(s)))
        base = a.out[:-4] if a.out.endswith('.npz') else a.out
        open(base + '.txt', 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
        print('wrote', base + '.txt')
    elif a.cmd == 'to-mitsuba':
        import numpy as np
        from .io import load, hpbrdf
        from .analysis import clean as CL
        t = CL.clean(load(a.table), a.floor or None, a.exchange, a.cosine, a.calib, a.floor_keep, a.floor_max_gain or None, a.min_cos)
        if (t.n_phi - 1) % a.phi_step:
            raise SystemExit(f'--phi-step must divide {t.n_phi - 1} (the first and last phi_d bins have to stay)')
        ch = np.arange(t.n_ch)
        if a.bands:
            lo, hi = (float(x) for x in a.bands.split(','))
            ch = ch[(t.channels >= lo) & (t.channels <= hi)]
        ch = ch[::a.band_step]
        if len(ch) == 0:
            raise SystemExit('no band left')
        phi = np.arange(0, t.n_phi, a.phi_step)
        print(f'{t.describe()}; corrections: {CL.describe(t.meta)}')
        hpbrdf.write_table(a.output, t, phi, ch, log=print)
        print(f'wrote {a.output}: {len(phi)}x{t.n_thd}x{t.n_thh} bins, {len(ch)} wavelengths, {os.path.getsize(a.output) / 2 ** 30:.2f} GB')
    elif a.cmd == 'umov':
        from .analysis import quirks as QK, umov as UM
        lam, rows = UM.split(QK.load(a.reduced), QK.load(a.reference), a.region, a.dark_max)
        lines = []
        UM.report(a.name or os.path.splitext(os.path.basename(a.reduced))[0], lam, rows, a.region, a.n, say=lambda s: (print(s), lines.append(s)))
        if a.out:
            open(a.out, 'w', encoding='utf-8').write('\n'.join(lines) + '\n'); print('wrote', a.out)
    elif a.cmd == 'slots':
        from .analysis import slots as SL
        if a.table.endswith('.npz'):
            red = SL.load(a.table)
        else:
            from .io import load
            b = [float(x) for x in a.bands.split(',')]; td = [float(x) for x in a.theta_d.split(',')]
            red = SL.reduce(_with_floor(load(a.table), a), (b[0], b[1]), a.phi_step, a.theta_h_min, (td[0], td[1])); SL.save(a.out, red); red = SL.load(a.out)
        lines = []
        SL.report(a.name or os.path.splitext(os.path.basename(a.table))[0], red, SL.measure(red), say=lambda s: (print(s), lines.append(s)))
        base = a.out[:-4] if a.out.endswith('.npz') else a.out
        open(base + '.txt', 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
        print('wrote', base + '.txt')


if __name__ == '__main__':
    main()
