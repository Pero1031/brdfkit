"""
brdfkit command line.

  brdfkit info      TABLE
  brdfkit compress  TABLE OUT.npz [--k1 8] [--k2 32] [--dtype float16] [--step 1]
  brdfkit eval      TABLE OUT.npz [--slices 0,30,60] [--sweep]
  brdfkit export    OUT.npz OUT.hpbpca [--dtype float16]
  brdfkit lu-chipman OUT.npz --refidx Al_Rakic.yml [Al_McPeak.yml ...] --out DIR [--raw TABLE] [--n-film 1.65]
  brdfkit merl-pca  DIR_OR_FILES... --out DIR [--stride 2] [--k 16] [--mapping log|nielsen] [--eps E] [--no-plot]

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


if __name__ == '__main__':
    main()
