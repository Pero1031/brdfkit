"""
End-to-end demo: from the downloaded datasets to edited hpBRDF materials ready for a renderer.

    python examples/latent_demo.py --hpbrdf-dir data/kaist/hpbrdf --merl-dir data/merl --out out/demo
    python examples/latent_demo.py ... --renderer D:/GoLD_rayt/out/brdfkit-spec/GoLD_rayt.exe --env path/to/env.hdr   (also render)

Needs (see docs/WALKTHROUGH.md for where to get them):
  <hpbrdf-dir>/Aluminum.hpbrdf      the material to edit (any hpBRDF material works, --material)
  <hpbrdf-dir>/Black_glass.hpbrdf   the reference for the instrument floor (optional: without it the floor is not removed)
  <merl-dir>/*.binary               the 100 MERL materials, for the latent space

Steps (each is skipped when its output already exists, unless --force):
  1. MERL inter-material PCA, log mapping           -> <out>/merl_pca/merl_pca.npz
  2. instrument floor measured on black glass       -> <out>/floor.npz
  3. factor model of the material, floor removed    -> <out>/<material>_factor.npz      (A x S x N, docs/FACTOR.md)
     and its lobe symmetrized over phi_d             -> <out>/<material>_factor_sym-median.npz   (--sym; removes the streak)
     (--calib CAL.npz: first corrected with the black-glass calibration -> ..._sym-median_cal.npz, docs/CALIBRATION.md)
  4. edits of its grey lobe along MERL components   -> <out>/edits/<name>.npz          (docs/LATENT_EDIT.md)
  5. each model baked into a .hpbpca                -> <out>/hpbpca/<name>.hpbpca
  6. GoLD_rayt scene files, one row of spheres per component -> <out>/scenes/row_pc<N>_sym-median.json
  7. (with --renderer) the renders and a contact sheet -> <out>/renders/
Times on an 8-thread desktop: step 2 about 20 s, step 3 7-10 min (two passes over 13 GB), step 5 2-3 min per model,
step 7 3-4 min per row (1280 px, 64 spp). Step 1 reads the 100 MERL files once.
"""
import argparse
import glob
import json
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from brdfkit.io import load                                    # noqa: E402
from brdfkit.compress import compress                          # noqa: E402
from brdfkit.export import write as export_hpbpca              # noqa: E402
from brdfkit.analysis import floor as FL, factor as FA, latent_edit as LE, merl_pca as MP, calibrate as CA   # noqa: E402

DEFAULT_EDITS = {1: (-1.0, -0.5, 0.5, 1.0), 3: (-1.5, -0.75, 0.75, 1.5)}


def step(path, force, what):
    if os.path.exists(path) and not force:
        print(f'[skip] {what}: {path} exists')
        return False
    print(f'[run ] {what} -> {path}', flush=True)
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--hpbrdf-dir', default=os.path.join(ROOT, 'data', 'kaist', 'hpbrdf'))
    ap.add_argument('--merl-dir', default=os.path.join(ROOT, 'data', 'merl'))
    ap.add_argument('--out', default=os.path.join(ROOT, 'out', 'demo'))
    ap.add_argument('--material', default='Aluminum')
    ap.add_argument('--floor-ref', default='Black_glass', help='hpBRDF material whose far field gives the instrument floor')
    ap.add_argument('--edits', default='', help='e.g. "1:-1,-0.5,0.5,1;3:-1.5,1.5" (component: shifts in units of the MERL spread)')
    ap.add_argument('--keep-albedo', type=float, default=45.0, help='angle of incidence [deg] whose albedo every edit keeps; <0 = none')
    ap.add_argument('--sym', choices=('median', 'mean', 'none'), default='median',
                    help='phi_d symmetrization of the lobe, applied to the unedited material and to every edit')
    ap.add_argument('--calib', default=None, help='calibration file of `brdfkit calibrate` (black glass): corrects the theta_d '
                    'dependence and scale of the lobe energy (docs/CALIBRATION.md); output names get _cal')
    ap.add_argument('--renderer', default=None, help='GoLD_rayt.exe (spectral build) to render the rows')
    ap.add_argument('--renderer-cwd', default=None, help='working directory of the renderer (default: the folder of the exe)')
    ap.add_argument('--env', default=None, help='environment map (.hdr) for the scenes; default: the renderer\'s own')
    ap.add_argument('--spp', type=int, default=64)
    ap.add_argument('--force', action='store_true')
    a = ap.parse_args()

    out = os.path.abspath(a.out)
    for d in ('merl_pca', 'edits', 'hpbpca', 'scenes', 'renders'):
        os.makedirs(os.path.join(out, d), exist_ok=True)
    src = os.path.join(a.hpbrdf_dir, f'{a.material}.hpbrdf')
    if not os.path.exists(src):
        sys.exit(f'missing {src}: download the hpBRDF material first (docs/WALKTHROUGH.md)')
    edits = DEFAULT_EDITS
    if a.edits:
        edits = {int(part.split(':')[0]): tuple(float(v) for v in part.split(':')[1].split(',')) for part in a.edits.split(';')}
    keep = None if a.keep_albedo < 0 else a.keep_albedo

    # 1. MERL PCA
    pca_path = os.path.join(out, 'merl_pca', 'merl_pca.npz')
    if step(pca_path, a.force, 'MERL PCA (log)'):
        files = sorted(glob.glob(os.path.join(a.merl_dir, '*.binary')))
        if len(files) < 20:
            sys.exit(f'found {len(files)} MERL .binary files in {a.merl_dir}; the latent space needs the database (100 materials)')
        MP.run_report(files, os.path.dirname(pca_path), stride=2, k=16, mapping='log', plot=False)
    pca = MP.MaterialPCA.load(pca_path)

    # 2. floor
    floor_path = os.path.join(out, 'floor.npz'); ref = os.path.join(a.hpbrdf_dir, f'{a.floor_ref}.hpbrdf')
    fl = None
    if os.path.exists(floor_path) and not a.force:
        print(f'[skip] floor: {floor_path} exists'); fl = FL.load(floor_path)
    elif os.path.exists(ref):
        print(f'[run ] floor from {ref}', flush=True)
        fl = FL.estimate(load(ref)); FL.save(floor_path, fl); FL.report(fl)
    else:
        print(f'[warn] {ref} not found: the instrument floor is NOT removed (the material will look bluish; docs/FLOOR.md)')

    # 3. factor model
    fac_path = os.path.join(out, f'{a.material}_factor.npz')
    if step(fac_path, a.force, 'factor model'):
        t = load(src)
        if fl is not None:
            t = FL.subtract(t, fl)
        m = FA.fit(t, smooth_deg=(0.0,))[0]
        FA.save(fac_path, m); FA.report(m)
    base = FA.load(fac_path)
    tag = ('' if a.sym == 'none' else f'_sym-{a.sym}') + ('_cal' if a.calib else '')   # output names say how the lobe was treated
    if a.calib:
        base = CA.apply(base, CA.load(a.calib))
    if tag:
        clean_path = os.path.join(out, f'{a.material}_factor{tag}.npz')
        base = FA.symmetrize(base, a.sym)
        what = ', '.join(([f'lobe symmetrized over phi_d ({a.sym})'] if a.sym != 'none' else []) + (['radiometric calibration'] if a.calib else []))
        if step(clean_path, a.force, what):
            FA.save(clean_path, base)
        fac_path = clean_path

    # 4. edits (the base is already symmetrized, so the edit itself does not symmetrize again)
    names = {}
    for pc, ts in edits.items():
        names[pc] = []
        for t in ts:
            name = f'{a.material}_pc{pc}_{t:+g}{tag}'; path = os.path.join(out, 'edits', name + '.npz')
            if step(path, a.force, f'edit PC{pc} {t:+g}'):
                m, info = LE.edit(base, pca, pcs=[(pc, t)], keep_albedo_deg=keep, sym='none')
                info['phi_symmetry'] = a.sym; m['meta']['edit']['phi_symmetry'] = a.sym
                FA.save(path, m); LE.report(info)
            names[pc].append((t, name, path))

    # 5. bake
    def bake(npz_model, name):
        hp = os.path.join(out, 'hpbpca', name + '.hpbpca'); pc_npz = os.path.join(out, 'hpbpca', name + '_pca.npz')
        if step(hp, a.force, f'bake {name}'):
            compress(load(npz_model), pc_npz, 8, 32, 'float16', step=4)
            export_hpbpca(pc_npz, hp)
            for f in (pc_npz, pc_npz[:-4] + '_basis_full.npz'):
                if os.path.exists(f):
                    os.remove(f)
        return hp
    hp_base = bake(fac_path, f'{a.material}_factor{tag}')
    rows = {pc: [(t, bake(p, n)) for t, n, p in lst] for pc, lst in names.items()}

    # 6. scenes: one row per component, the unedited material in the middle
    scenes = {}
    for pc, lst in rows.items():
        items = sorted(lst + [(0.0, hp_base)], key=lambda x: x[0])
        n = len(items); x0 = -0.1 * (n - 1) / 2
        mats, objs = {}, []
        for i, (t, hp) in enumerate(items):
            mats[f'm{i}'] = {'type': 'measured', 'file': hp.replace('\\', '/'), 'alphaSample': 0.1}
            objs.append({'shape': 'sphere', 'center': [x0 + 0.1 * i, 0.04, 0.0], 'radius': 0.04, 'material': f'm{i}'})
        mats['ground'] = {'type': 'diffuse', 'albedo': 0.3}
        objs.append({'shape': 'quad', 'corner': [-1, 0, 1], 'edge1': [2, 0, 0], 'edge2': [0, 0, -2], 'material': 'ground'})
        env = {'type': 'hdr', 'file': os.path.abspath(a.env).replace('\\', '/'), 'scale': 1.0} if a.env else 'hdr'
        sc = {'options': {'width': 1280, 'spp': a.spp}, 'camera': {'from': [0.0, 0.17, 0.12 * n + 0.02], 'at': [0.0, 0.035, 0.0], 'fov': 30},
              'environment': env, 'materials': mats, 'objects': objs}
        path = os.path.join(out, 'scenes', f'row_pc{pc}{tag}.json'); json.dump(sc, open(path, 'w'), indent=1)
        scenes[pc] = (path, [t for t, _ in items])
        print(f'[done] scene {path}  (shifts {", ".join(f"{t:+g}" for t, _ in items)})')

    # 7. render
    if not a.renderer:
        print('\nRender the scenes with GoLD_rayt (spectral build):  GoLD_rayt.exe <scene.json> --denoise --stokes --out row.png')
        return
    cwd = a.renderer_cwd or os.path.dirname(os.path.abspath(a.renderer))
    if not os.path.isdir(cwd):
        sys.exit(f'renderer working directory {cwd} does not exist (--renderer-cwd)')
    pngs = {}
    for pc, (path, ts) in scenes.items():
        png = os.path.join(out, 'renders', f'row_pc{pc}{tag}.png')
        if step(png, a.force, f'render PC{pc}'):
            subprocess.run([a.renderer, path, '--seed', '1', '--denoise', '--stokes', '--out', png], cwd=cwd, check=True,
                           stdout=subprocess.DEVNULL)
        den = png[:-4] + '_denoised.png'
        pngs[pc] = (den if os.path.exists(den) else png, ts)
    try:
        import matplotlib; matplotlib.use('Agg')
        import matplotlib.pyplot as plt, matplotlib.image as mpimg
    except ImportError:
        return
    fig, axs = plt.subplots(len(pngs), 1, figsize=(12.8, 3.2 * len(pngs)), squeeze=False)
    for ax, (pc, (png, ts)) in zip(axs[:, 0], pngs.items()):
        ax.imshow(mpimg.imread(png)[200:480]); ax.axis('off')
        ax.set_title(f'{a.material}: grey lobe moved along MERL PC{pc} by ' + ', '.join(f'{t:+g}' for t in ts) + ' (spread units)'
                     + (f', lobe symmetrized over phi_d ({a.sym})' if a.sym != 'none' else '')
                     + (', calibrated on black glass' if a.calib else ''), loc='left')
    plt.tight_layout(); sheet = os.path.join(out, 'renders', f'sheet{tag}.png'); plt.savefig(sheet, dpi=90)
    print(f'[done] {sheet}')


if __name__ == '__main__':
    main()
