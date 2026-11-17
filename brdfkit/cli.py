"""
brdfkit command line.

  brdfkit info      TABLE
  brdfkit compress  TABLE OUT.npz [--k1 8] [--k2 32] [--dtype float16] [--step 1] [--floor FLOOR.npz] [--floor-keep 0.02] [--floor-max-gain 4]
  brdfkit eval      TABLE OUT.npz [--slices 0,30,60] [--sweep] [--floor FLOOR.npz] [--floor-keep 0.02] [--floor-max-gain 4]
  brdfkit floor     TABLE --out FLOOR.npz [--phi-step 12] [--theta-h-min 40] [--theta-d 25,75] [--stat p50]
  brdfkit edit      MODEL.npz --model MERL_PCA.npz --out EDITED.npz [--pc 1:-0.5 ...] [--toward nickel:0.5 ...] [--k 16] [--keep-albedo 45]
                    [--sym median|mean|none] [--calib CAL.npz [--calib-keep-level DEG]]   (no shift: just the cleaned model)
  brdfkit calibrate GLASS_MODEL.npz --out CAL.npz [--cutoff 10] [--smooth 3] [--check MODEL.npz --refidx NK.yml [--film-nm 20] [--n-film 1.65]]
                    (radiometric error of the table from black glass, docs/CALIBRATION.md)
  brdfkit factor    TABLE --out MODEL.npz [--smooth 0,5] [--measured 20,80] [--floor FLOOR.npz] [--sym none|median|mean]
                    (A x S x N model; TABLE may then be MODEL.npz)
  brdfkit export    OUT.npz OUT.hpbpca [--dtype float16]
  brdfkit lu-chipman OUT.npz --refidx Al_Rakic.yml [Al_McPeak.yml ...] --out DIR [--raw TABLE] [--n-film 1.65]
  brdfkit merl-pca  DIR_OR_FILES... --out DIR [--stride 2] [--k 16] [--mapping log|nielsen] [--eps E] [--no-plot]
  brdfkit merl-project TABLE --model MERL_PCA.npz [MERL_PCA.npz ...] --out DIR [--k 16] [--sym full|none] [--reference MERL.binary] [--recompute] [--no-plot]
                    [--floor FLOOR.npz] [--floor-keep 0.02]
  brdfkit merl-subsets MERL_GRID.npz MERL_DIR --out DIR [--all-model MERL_PCA.npz] [--random 30] [--seed 0] [--ablate]
  brdfkit merl-compare RUN_DIR [RUN_DIR ...]      (summary of several merl-project outputs)

TABLE: .hpbrdf / .pbsdf (Mitsuba tensor_file), .mat (KAIST pBRDF 2020), .binary (MERL), .npz (factor model of `brdfkit factor`).
"""
import argparse
import os


def _floor_args(p):
    p.add_argument('--floor', default=None, help='floor file of `brdfkit floor` (e.g. measured on black glass) to subtract from m00')
    p.add_argument('--floor-keep', type=float, default=0.02, help='never take m00 below this fraction of its measured value')
    p.add_argument('--floor-max-gain', type=float, default=4.0, help='cap on the rise of m_ij/m00 caused by the subtraction (0 = no cap)')


def _with_floor(table, a):
    if not a.floor:
        return table
    from .analysis import floor as FL
    return FL.subtract(table, a.floor, a.floor_keep, a.floor_max_gain or None)


def main(argv=None):
    ap = argparse.ArgumentParser(prog='brdfkit', description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('info'); p.add_argument('table')
    p = sub.add_parser('compress'); p.add_argument('table'); p.add_argument('output')
    p.add_argument('--k1', type=int, default=8); p.add_argument('--k2', type=int, default=32)
    p.add_argument('--dtype', default='float16'); p.add_argument('--step', type=int, default=1, help='phi_d stride of the covariance pass')
    _floor_args(p)
    p = sub.add_parser('eval'); p.add_argument('table'); p.add_argument('npz')
    p.add_argument('--slices', default=''); p.add_argument('--sweep', action='store_true')
    _floor_args(p)
    p = sub.add_parser('floor', help='estimate the instrument floor (unpolarized offset) far from the specular lobe'); p.add_argument('table')
    p.add_argument('--out', required=True); p.add_argument('--phi-step', type=int, default=12); p.add_argument('--theta-h-min', type=float, default=40.0)
    p.add_argument('--theta-d', default='25,75', help='measured theta_d range in degrees'); p.add_argument('--stat', choices=('p10', 'p25', 'p50', 'p90'), default='p50')
    p = sub.add_parser('export'); p.add_argument('npz'); p.add_argument('output'); p.add_argument('--dtype', default='float16')
    p = sub.add_parser('lu-chipman'); p.add_argument('npz'); p.add_argument('--refidx', nargs='+', required=True)
    p.add_argument('--out', required=True); p.add_argument('--raw', default=None); p.add_argument('--n-film', type=float, default=1.65)
    p = sub.add_parser('factor', help='fit m00 = A(angles) S(theta_d, lambda), m_ij/m00 = N(theta_d, lambda)'); p.add_argument('table')
    p.add_argument('--out', required=True); p.add_argument('--smooth', default='0', help='smoothing widths of S over theta_d [deg], comma separated')
    p.add_argument('--measured', default='20,80', help='measured theta_d range [deg]; S and N are extended from its edges'); _floor_args(p)
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


if __name__ == '__main__':
    main()
