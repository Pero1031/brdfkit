"""
brdfkit command line.

  brdfkit info      TABLE
  brdfkit compress  TABLE OUT.npz [--k1 8] [--k2 32] [--dtype float16] [--step 1]
  brdfkit eval      TABLE OUT.npz [--slices 0,30,60] [--sweep]
  brdfkit export    OUT.npz OUT.hpbpca [--dtype float16]
  brdfkit lu-chipman OUT.npz --refidx Al_Rakic.yml [Al_McPeak.yml ...] --out DIR [--raw TABLE] [--n-film 1.65]
  brdfkit merl-pca  DIR_OR_FILES... --out DIR [--stride 2] [--k 16] [--mapping log|nielsen] [--eps E] [--no-plot]
  brdfkit merl-project TABLE --model MERL_PCA.npz [MERL_PCA.npz ...] --out DIR [--k 16] [--sym full|none] [--reference MERL.binary] [--recompute] [--no-plot]
  brdfkit merl-subsets MERL_GRID.npz MERL_DIR --out DIR [--all-model MERL_PCA.npz] [--random 30] [--seed 0] [--ablate]

TABLE: .hpbrdf / .pbsdf (Mitsuba tensor_file), .mat (KAIST pBRDF 2020), .binary (MERL).
"""
import argparse
import os


def main(argv=None):
    ap = argparse.ArgumentParser(prog='brdfkit', description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('info'); p.add_argument('table')
    p = sub.add_parser('compress'); p.add_argument('table'); p.add_argument('output')
    p.add_argument('--k1', type=int, default=8); p.add_argument('--k2', type=int, default=32)
    p.add_argument('--dtype', default='float16'); p.add_argument('--step', type=int, default=1, help='phi_d stride of the covariance pass')
    p = sub.add_parser('eval'); p.add_argument('table'); p.add_argument('npz')
    p.add_argument('--slices', default=''); p.add_argument('--sweep', action='store_true')
    p = sub.add_parser('export'); p.add_argument('npz'); p.add_argument('output'); p.add_argument('--dtype', default='float16')
    p = sub.add_parser('lu-chipman'); p.add_argument('npz'); p.add_argument('--refidx', nargs='+', required=True)
    p.add_argument('--out', required=True); p.add_argument('--raw', default=None); p.add_argument('--n-film', type=float, default=1.65)
    p = sub.add_parser('merl-pca', help='inter-material PCA of a set of MERL tables'); p.add_argument('paths', nargs='+')
    p.add_argument('--out', required=True); p.add_argument('--stride', default='2', help='grid stride, one int or phi_d,theta_d,theta_h')
    p.add_argument('--k', type=int, default=16, help='components kept in the basis'); p.add_argument('--mapping', choices=('log', 'nielsen'), default='log')
    p.add_argument('--eps', type=float, default=None); p.add_argument('--no-plot', action='store_true')
    p = sub.add_parser('merl-project', help='project a table (e.g. hpBRDF aluminium) onto saved MERL PCA models'); p.add_argument('table')
    p.add_argument('--model', nargs='+', required=True, help='merl_pca.npz written by merl-pca'); p.add_argument('--out', required=True)
    p.add_argument('--k', type=int, default=16); p.add_argument('--sym', choices=('full', 'none'), default='full', help='average the phi_d-equivalent slices (full) or not')
    p.add_argument('--reference', default=None, help='MERL .binary of a comparison material'); p.add_argument('--recompute', action='store_true')
    p.add_argument('--no-plot', action='store_true')
    p = sub.add_parser('merl-subsets', help='metal-only versus all-material PCA basis for a table on the MERL grid')
    p.add_argument('grid', help='<out>/<name>_merl_grid_full.npz written by merl-project'); p.add_argument('merl_dir'); p.add_argument('--out', required=True)
    p.add_argument('--all-model', default=None, help='saved 100-material merl_pca.npz (log), for leave-one-out numbers')
    p.add_argument('--random', type=int, default=30, help='random control subsets per size'); p.add_argument('--seed', type=int, default=0)
    p.add_argument('--ablate', action='store_true', help='refit without each material in turn (100 fits)')
    a = ap.parse_args(argv)

    if a.cmd == 'info':
        from .io import load
        t = load(a.table); print(t.describe()); print(t.meta)
    elif a.cmd == 'compress':
        from .io import load
        from .compress import compress
        compress(load(a.table), a.output, a.k1, a.k2, a.dtype, a.step)
    elif a.cmd == 'eval':
        from .io import load
        from .compress import evaluate
        slices = [int(s) for s in a.slices.split(',')] if a.slices else None
        evaluate(load(a.table), a.npz, slices, a.sweep)
    elif a.cmd == 'export':
        from .export import write
        write(a.npz, a.output, a.dtype)
    elif a.cmd == 'lu-chipman':
        from .analysis.lu_chipman import run_report
        run_report(a.npz, [(p, os.path.splitext(os.path.basename(p))[0]) for p in a.refidx], a.out, raw=a.raw, n_film=a.n_film)
    elif a.cmd == 'merl-pca':
        from .analysis.merl_pca import run_report
        stride = [int(s) for s in a.stride.split(',')]
        run_report(a.paths, a.out, stride[0] if len(stride) == 1 else stride, a.k, a.mapping, a.eps, plot=not a.no_plot)
    elif a.cmd == 'merl-project':
        from .analysis.merl_project import run_report
        run_report(a.table, a.model, a.out, a.k, a.sym, a.reference, a.recompute, plot=not a.no_plot)
    elif a.cmd == 'merl-subsets':
        from .analysis.merl_subsets import run_study
        run_study(a.grid, a.merl_dir, a.out, a.all_model, a.random, a.seed, a.ablate)


if __name__ == '__main__':
    main()
