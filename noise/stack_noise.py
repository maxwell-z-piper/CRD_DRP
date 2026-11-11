#!/usr/bin/env python3
"""Make science and Gaussian-noise stacks with one fixed KcwiKit recipe.

Run this after detector reduction and sky/telluric calibration. There are two
stages:

1. ``--science-only`` stacks the chosen reduced exposures and writes a recipe
   pointing back to that new science stack.
2. The normal run first replays the science stack to verify the recipe, then
   replaces each exposure's flux with Gaussian draws from its UNCERT extension
   and stacks those draws using the same alignment, masks, and settings.

The input exposures are read only. These mocks measure what fixed stacking does
to reduced-exposure noise; they do not repeat the detector or ZAP reductions.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import contextmanager
import hashlib
import json
import multiprocessing
import os
from pathlib import Path
import shutil

import numpy as np
from astropy.io import fits


def sha256(path):
    """Fingerprint a file so downstream calibration can detect changed inputs."""
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path, value):
    """Replace a manifest only after its complete JSON has reached disk."""
    temp = Path(path).with_suffix('.tmp')
    temp.write_text(json.dumps(value, indent=2) + '\n')
    temp.replace(path)


@contextmanager
def working_directory(path):
    """KcwiKit finds its .list, .par, and shift files in the working directory."""
    previous = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


def make_noise_exposure(source, destination, rng):
    """Write one mock exposure, changing only the primary (flux) image.

    Every valid pixel is drawn from N(0, UNCERT²). We copy the original FITS
    headers, uncertainty, mask, flags, and other extensions unchanged.
    """
    with fits.open(source, memmap=False) as hdus:
        sigma = np.asarray(hdus['UNCERT'].data, dtype=float)
        original = hdus[0].data
        if sigma.shape != original.shape:
            raise ValueError(f'UNCERT shape mismatch: {source}')
        # Some real UNCERT entries are negative. KcwiKit squares them, so their
        # Gaussian sigma is the absolute value. Zeros and nonfinite entries
        # must keep the same valid/invalid pattern as the science exposure.
        valid = np.isfinite(original) & np.isfinite(sigma) & (sigma != 0)
        draw = rng.normal(size=sigma.shape) * np.where(valid, np.abs(sigma), 0.)
        # KcwiKit treats exact input zeros as invalid before projecting the
        # variance. Preserve that support as well as non-finite locations.
        draw[original == 0] = 0.
        draw[~np.isfinite(original)] = original[~np.isfinite(original)]
        hdus[0].data = draw
        hdus[0].header['NOISEMC'] = (True, 'Gaussian draw from reduced UNCERT')
        hdus.writeto(destination, overwrite=False)


def stack_products(directory):
    """Return KcwiKit's i=flux, v=variance, m=mask, e=exposure files."""
    return {letter: Path(directory) / f'stack_{letter}cubes.fits' for letter in 'ivme'}


def exposure_key(path):
    """Match an exposure across KSkyWizard ZAP and KcwiKit shift filenames."""
    # KcwiKit may record either kr..._icubes or kr..._zap_icubes in a shift
    # list. Both identify the same numbered exposure for this order check.
    name = Path(path).name
    if name.endswith('.fits'):
        name = name[:-5]
    if name.endswith('_icubes'):
        name = name[:-7]
    if name.endswith('_zap'):
        name = name[:-4]
    return name


def compare_products(actual, reference, letters='ivme', rtol=2e-5):
    """Require matching grids and pixel values, one wavelength plane at a time.

    Science replay checks i/v/m/e. Gaussian replay checks v/m/e: random flux
    should change, while variance, valid coverage, and exposure time should not.
    """
    for letter in letters:
        with fits.open(actual[letter], memmap=True) as a, fits.open(reference[letter], memmap=True) as b:
            if a[0].shape != b[0].shape:
                raise ValueError(f'{letter}: replay/reference shapes differ: {a[0].shape}, {b[0].shape}')
            for key in ('CRVAL1','CRVAL2','CRVAL3','CRPIX1','CRPIX2','CRPIX3',
                        'CD1_1','CD1_2','CD2_1','CD2_2','CD3_3'):
                if not np.isclose(a[0].header[key], b[0].header[key], rtol=0, atol=1e-10):
                    raise ValueError(f'{letter}: replay/reference WCS mismatch at {key}')
            for plane in range(a[0].shape[0]):
                x, y = a[0].data[plane], b[0].data[plane]
                atol = 0. if letter == 'm' else max(1e-30, np.nanmax(np.abs(y)) * 1e-8)
                if not np.allclose(x, y, rtol=0 if letter == 'm' else rtol, atol=atol, equal_nan=True):
                    raise ValueError(f'{letter}: replay/reference values differ at plane {plane}; '
                                     'check exposure choice, shifts, masks, flux factors and KcwiKit version')


def validate_recipe(recipe, require_reference=True):
    """Check inputs before a costly stack, and return saved shifts in order.

    A first ``--science-only`` run has no reference yet. Noise generation must
    use the recipe emitted by that run, which names the reference science cube.
    """
    if recipe.get('schema') != 'crd-kcwikit-replay-v1':
        raise ValueError('Expected schema crd-kcwikit-replay-v1')
    required = [recipe['par_file'], recipe['shift_file']]
    if require_reference:
        if not recipe.get('reference_icube'):
            raise ValueError('Build this science stack with --science-only first; then use recipe_for_noise.json')
        required.append(recipe['reference_icube'])
    inputs = recipe['exposures']
    if not inputs or len({x['path'] for x in inputs}) != len(inputs):
        raise ValueError('Exposure list must be nonempty and unique')
    required += [x['path'] for x in inputs]
    required += list(recipe.get('sidecars', {}).values())
    for path in required:
        if not Path(path).is_file():
            raise FileNotFoundError(path)
    for row in inputs:
        if len(row['trim']) != 2 or any(int(n) < 0 for n in row['trim']):
            raise ValueError('Each trim must be two nonnegative integers')
        with fits.open(row['path']) as hdus:
            for name in ('UNCERT', 'MASK', 'FLAGS'):
                if name not in hdus:
                    raise ValueError(f"Missing {name}: {row['path']}")
    shifts = np.loadtxt(recipe['shift_file'], dtype=str, ndmin=2)
    if shifts.shape != (len(inputs), 3):
        raise ValueError('Shift file must have one path/x/y row per exposure, in the same order')
    for row, shift in zip(inputs, shifts):
        if exposure_key(row['path']) != exposure_key(shift[0]):
            raise ValueError('Shift row/exposure order mismatch')
    kwargs = recipe.get('stack_kwargs', {})
    if kwargs.get('crr', False) or kwargs.get('medcube', False) or kwargs.get('use_astrom', False):
        raise ValueError('This replay supports the documented fixed-weight stack, without CR refitting, medcube or astrom refitting')
    if kwargs.get('method', 'drizzle') != 'drizzle' or kwargs.get('cubed', False):
        raise ValueError('This recipe expects flux-calibrated icubes and drizzle stacking')
    # KcwiKit can silently discover sidecar files beside its .par. An omitted
    # sidecar could change the science stack without changing this recipe.
    prefix = Path(recipe['par_file']).with_suffix('')
    for suffix in ('.flx.list', '.pre.list', '.preshift.list', '.astrom.list'):
        p = Path(str(prefix) + suffix)
        if p.exists() and suffix not in recipe.get('sidecars', {}):
            raise ValueError(f'Include existing sidecar explicitly in recipe: {p}')
    return shifts


def run_stack(recipe, shifts, directory, kcwi, rng=None):
    """Build a private KcwiKit workspace and stack science or one noise draw.

    ``rng=None`` means science: symlink each original exposure, which KcwiKit
    only reads. An RNG means noise: write a fresh mock exposure for every row.
    The saved shift numbers are reused in both cases; no new alignment is fit.
    """
    directory.mkdir(parents=True, exist_ok=False)
    rows = []
    for i, exposure in enumerate(recipe['exposures']):
        source = Path(exposure['path'])
        dest = directory / f'exposure_{i:03d}_icubes.fits'
        if rng is None:
            # KcwiKit only reads its inputs. Symlinks save disk for science replay.
            dest.symlink_to(source.resolve())
        else:
            make_noise_exposure(source, dest, rng)
        region = source.with_suffix('.thum.reg')
        if region.exists():
            shutil.copy2(region, dest.with_suffix('.thum.reg'))
        rows.append(f"{dest.with_name(dest.name.replace('_icubes.fits',''))} {exposure['trim'][0]} {exposure['trim'][1]}\n")
    # KcwiKit locates these three files by the shared basename ``stack``.
    # The new local exposure names keep the mock and science calls identical.
    (directory/'stack.list').write_text(''.join(rows))
    shutil.copy2(recipe['par_file'], directory/'stack.par')
    (directory/'stack.shift.list').write_text(''.join(
        f'exposure_{i:03d}_icubes {row[1]} {row[2]}\n' for i,row in enumerate(shifts)))
    for suffix, source in recipe.get('sidecars', {}).items():
        if suffix != '.flx.list':
            raise ValueError('Only flux-normalization sidecars are currently supported; use an explicit fixed-grid recipe')
        shutil.copy2(source, directory/('stack'+suffix))
    kwargs = dict(recipe.get('stack_kwargs', {}))
    kwargs.update(overwrite=True)
    with working_directory(directory):
        kcwi.kcwi_stack('stack.list', **kwargs)
    products = stack_products(directory)
    for path in products.values():
        if not path.is_file():
            raise RuntimeError(f'KcwiKit did not produce {path}')
    # Only remove our generated exposure copies after a successful stack.
    for path in directory.glob('exposure_*_icubes.fits'):
        path.unlink()
    # KcwiKit keeps large projected/trimmed intermediates in this private cache.
    # The four final cubes and the exact recipe remain outside it.
    cache = directory / 'kcwi_stack'
    if cache.is_dir():
        shutil.rmtree(cache)
    return products


def run_noise_realization(recipe, shifts, root, index, seed, replay):
    """Run one draw in its own process and private KcwiKit working directory."""
    from kcwikit.kcwi import kcwi

    products = run_stack(recipe, shifts, root/f'noise_{index:04d}', kcwi,
                         np.random.default_rng(seed))
    compare_products(products, replay, letters='vme')
    digest = sha256(products['i'])
    # Script 2 needs only mock flux. Keep the checked v/m/e products out of
    # the final ensemble so parallel jobs do not multiply disk use needlessly.
    for letter in 'vme':
        products[letter].unlink()
    return index, str(products['i']), digest


def main(argv=None):
    """Run a quick input audit, build a reference, or generate an ensemble."""
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument('--recipe', required=True, type=Path,
                   help='JSON with exposure paths, saved shifts, and KcwiKit settings')
    p.add_argument('--output-dir', required=True, type=Path,
                   help='New directory for outputs; existing directories are not overwritten')
    p.add_argument('--realizations', type=int, default=16,
                   help='Independent noise stacks to make (minimum 4; default 16)')
    p.add_argument('--seed', type=int, default=20261021,
                   help='Seed for reproducible, independent Gaussian draws')
    p.add_argument('--workers', type=int, default=1,
                   help='Independent KcwiKit stacks to run at once (default 1; each needs substantial RAM and temporary disk)')
    p.add_argument('--audit-only', action='store_true',
                   help='Check recipe inputs without running KcwiKit or creating outputs')
    p.add_argument('--science-only', action='store_true',
                   help='Build a NEW science stack and emit its noise-replay recipe')
    args = p.parse_args(argv)
    recipe = json.loads(args.recipe.read_text())
    shifts = validate_recipe(recipe, require_reference=not args.science_only)
    print(f"Validated {len(recipe['exposures'])} input exposures for {recipe['arm']}", flush=True)
    if args.audit_only:
        print('Input audit only; stack equivalence still requires the automatic science replay.')
        return 0
    if args.realizations < 4:
        raise ValueError('Use at least four independent realizations')
    if args.workers < 1:
        raise ValueError('--workers must be at least 1')
    from kcwikit.kcwi import kcwi
    root = args.output_dir.expanduser().resolve()
    root.mkdir(parents=True, exist_ok=False)
    if args.science_only:
        # This reference stack is the anchor for every later mock. Keeping its
        # path in a new recipe prevents accidentally pairing mocks with an old
        # stack made from different exposures or settings.
        products = run_stack(recipe, shifts, root/'science_stack', kcwi)
        recipe['reference_icube'] = str(products['i'])
        write_json(root/'recipe_for_noise.json', recipe)
        print(f"Science stack: {products['i']}\nNoise recipe: {root/'recipe_for_noise.json'}")
        return 0
    ref_i = Path(recipe['reference_icube']).resolve()
    reference = {x: ref_i.with_name(ref_i.name.replace('_icubes.fits', f'_{x}cubes.fits')) for x in 'ivme'}
    # Record the exact recipe, code, science products, and exposure bytes used.
    # Script 2 later checks these hashes before calculating noise multipliers.
    recipe_sha = sha256(args.recipe)
    reference_sha = {x: sha256(path) for x,path in reference.items()}
    source_sha = {row['path']: sha256(row['path']) for row in recipe['exposures']}
    replay = run_stack(recipe, shifts, root/'science_replay', kcwi)
    # Stop before drawing noise if the supposedly identical science stack
    # cannot be reproduced. A mismatch invalidates the calibration reference.
    compare_products(replay, reference)
    manifest = {
        'schema': 'crd-noise-stack-ensemble-v1',
        'arm': recipe['arm'],
        'seed': args.seed,
        'recipe': recipe,
        'recipe_sha256': recipe_sha,
        'source_sha256': source_sha,
        'kcwikit_source_sha256': sha256(kcwi.__file__),
        'reference_sha256': reference_sha,
        'reference_products': {key: str(path) for key, path in reference.items()},
        'science_replay_passed': True,
        'requested_realizations': args.realizations,
        'workers': args.workers,
        'noise_icubes': [],
        'noise_sha256': {},
        'complete': False,
        'scope': (
            'Independent reduced-exposure noise through fixed KcwiKit stacking; '
            'upstream DRP and sky-model uncertainty not simulated'
        ),
    }
    manifest_path = root/'noise_manifest.json'
    write_json(manifest_path, manifest)
    # SeedSequence fixes the draw for each index regardless of completion order.
    seeds = np.random.SeedSequence(args.seed).spawn(args.realizations)
    def record(result):
        index, path, digest = result
        print(f'Noise stack {index+1}/{args.realizations} validated', flush=True)
        manifest['noise_sha256'][path] = digest
        manifest['noise_icubes'] = sorted(manifest['noise_sha256'])
        write_json(manifest_path, manifest)

    if args.workers == 1:
        for i, seed in enumerate(seeds):
            record(run_noise_realization(recipe, shifts, root, i, seed, replay))
    else:
        # KcwiKit changes its process working directory. Separate spawned
        # processes, not threads, keep the per-draw .list and cache isolated.
        context = multiprocessing.get_context('spawn')
        with ProcessPoolExecutor(max_workers=min(args.workers, args.realizations),
                                 mp_context=context) as pool:
            futures = [pool.submit(run_noise_realization, recipe, shifts, root,
                                   i, seed, replay) for i, seed in enumerate(seeds)]
            try:
                for future in as_completed(futures):
                    record(future.result())
            except BaseException:
                for future in futures:
                    future.cancel()
                raise
    # ``complete`` stays false until every requested realization has passed.
    manifest['complete'] = True
    write_json(manifest_path, manifest)
    print(manifest_path)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
