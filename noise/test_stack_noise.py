"""Checks the fixed science/noise replay recipe without a costly real stack."""
from __future__ import annotations

from pathlib import Path
import json
import os
import subprocess
import sys
import numpy as np
from astropy.io import fits

from stack_noise import make_noise_exposure, run_stack, stack_products, validate_recipe
from make_recipe import recipe_from_list


def test_noise_exposure_preserves_uncertainty_and_flags(tmp_path):
    shape = (32, 3, 2)
    primary = fits.PrimaryHDU(np.full(shape, 7., dtype=np.float32))
    primary.data[0,0,0] = 0.
    primary.header['CRVAL3'] = 5000.
    sigma = np.full(shape, 2., dtype=np.float32)
    sigma[0,0,1] = -2.
    fits.HDUList([primary,
        fits.ImageHDU(np.zeros(shape, dtype=np.uint8), name='MASK'),
        fits.ImageHDU(sigma, name='UNCERT'),
        fits.ImageHDU(np.zeros(shape, dtype=np.uint8), name='FLAGS')
    ]).writeto(tmp_path/'source.fits')
    dest = tmp_path/'noise.fits'
    make_noise_exposure(tmp_path/'source.fits', dest, np.random.default_rng(7))
    with fits.open(dest) as hdus:
        assert hdus[0].header['NOISEMC']
        assert hdus[0].header['CRVAL3'] == 5000.
        assert abs(np.mean(hdus[0].data)) < .5
        assert 1.4 < np.std(hdus[0].data) < 2.6
        assert hdus[0].data[0,0,0] == 0.
        assert hdus['UNCERT'].data[0,0,1] == -2.
        assert np.all(np.abs(hdus['UNCERT'].data) == 2.)
        assert not np.any(hdus['FLAGS'].data)


def test_science_and_mock_use_same_saved_shift_and_trim(tmp_path):
    class FakeKCWI:
        def __init__(self): self.calls = []
        def kcwi_stack(self, filename, **kwargs):
            self.calls.append((Path('stack.list').read_text(),
                               Path('stack.shift.list').read_text(), kwargs))
            for path in stack_products(Path.cwd()).values():
                fits.PrimaryHDU(np.zeros((2,2,2))).writeto(path)
    source = tmp_path/'kr241226_00172_zap_icubes.fits'
    fits.HDUList([fits.PrimaryHDU(np.ones((2,2,2))),
        fits.ImageHDU(np.zeros((2,2,2)),name='MASK'),
        fits.ImageHDU(np.ones((2,2,2)),name='UNCERT'),
        fits.ImageHDU(np.zeros((2,2,2)),name='FLAGS')]).writeto(source)
    par = tmp_path/'science.par'
    par.write_text('stack_orientation -45\n')
    recipe = {'exposures':[{'path':str(source),'trim':[3,3]}],
              'par_file':str(par), 'stack_kwargs':{'drizzle':.9}}
    shifts = np.asarray([['kr241226_00172_icubes','0.125','-0.25']])
    fake = FakeKCWI()
    run_stack(recipe, shifts, tmp_path/'science', fake)
    run_stack(recipe, shifts, tmp_path/'mock', fake, np.random.default_rng(5))
    assert len(fake.calls) == 2
    assert ' 3 3\n' in fake.calls[0][0]
    assert fake.calls[0][1] == fake.calls[1][1]
    assert '0.125 -0.25' in fake.calls[0][1]
    assert fake.calls[0][2]['drizzle'] == .9


def test_recipe_accepts_shift_file_from_zap_alignment(tmp_path):
    source = tmp_path/'kr241226_00172_zap_icubes.fits'
    fits.HDUList([fits.PrimaryHDU(np.ones((2,2,2))),
        fits.ImageHDU(np.zeros((2,2,2)),name='MASK'),
        fits.ImageHDU(np.ones((2,2,2)),name='UNCERT'),
        fits.ImageHDU(np.zeros((2,2,2)),name='FLAGS')]).writeto(source)
    par = tmp_path/'red.par'
    par.write_text('dimension 67 112\n')
    shifts = tmp_path/'red.shift.list'
    recipe = {'schema':'crd-kcwikit-replay-v1', 'arm':'RH3',
              'par_file':str(par), 'shift_file':str(shifts),
              'reference_icube':None,
              'exposures':[{'path':str(source),'trim':[3,3]}],
              'stack_kwargs':{'drizzle':.9,'method':'drizzle'}, 'sidecars':{}}
    for root in ('kr241226_00172_zap_icubes','kr241226_00172_icubes'):
        shifts.write_text(f'{tmp_path/root} 0.0 0.0\n')
        assert validate_recipe(recipe, require_reference=False).shape == (1,3)


def test_recipe_uses_every_listed_final_cube_and_trim(tmp_path):
    roots = [tmp_path/'kb300101_00001_zap', tmp_path/'kb300101_00002_zap']
    for root in roots:
        fits.HDUList([fits.PrimaryHDU(np.ones((2,2,2))),
            fits.ImageHDU(np.zeros((2,2,2)), name='MASK'),
            fits.ImageHDU(np.ones((2,2,2)), name='UNCERT'),
            fits.ImageHDU(np.zeros((2,2,2)), name='FLAGS')
        ]).writeto(str(root) + '_icubes.fits')
    stack_list = tmp_path/'blue_kskywizard.list'
    stack_list.write_text(f'{roots[0]} 3 4\n{roots[1]} 2 3\n')
    (tmp_path/'blue_kskywizard.par').write_text('dimension 67 112\n')
    (tmp_path/'blue_kskywizard.shift.list').write_text(
        f'{roots[0]}_icubes 0 0\n{roots[1]}_icubes 1 -1\n')
    recipe = recipe_from_list(stack_list, 'BL')
    assert recipe['arm'] == 'BL'
    assert len(recipe['exposures']) == 2
    assert [row['trim'] for row in recipe['exposures']] == [[3, 4], [2, 3]]
    assert [row['path'] for row in recipe['exposures']] == [str(root) + '_icubes.fits' for root in roots]


def test_parallel_draws_finish_with_sorted_complete_manifest(tmp_path):
    """Exercise spawned workers with a tiny stand-in for KcwiKit's file API."""
    package = tmp_path/'kcwikit'/'kcwi'
    package.mkdir(parents=True)
    (package.parent/'__init__.py').write_text('')
    (package/'__init__.py').write_text('')
    (package/'kcwi.py').write_text('''from pathlib import Path
import numpy as np
from astropy.io import fits
def kcwi_stack(filename, **kwargs):
    root = Path(Path(filename).read_text().split()[0])
    flux = fits.getdata(str(root) + '_icubes.fits')
    header = fits.Header()
    for key in ('CRVAL1','CRVAL2','CRVAL3','CRPIX1','CRPIX2','CRPIX3',
                'CD1_1','CD1_2','CD2_1','CD2_2','CD3_3'):
        header[key] = 1.0
    for letter in 'ivme':
        data = flux if letter == 'i' else np.ones_like(flux)
        fits.PrimaryHDU(data, header=header).writeto('stack_' + letter + 'cubes.fits')
''')
    source = tmp_path/'kb300101_00001_icubes.fits'
    fits.HDUList([fits.PrimaryHDU(np.ones((3,2,2))),
        fits.ImageHDU(np.zeros((3,2,2)), name='MASK'),
        fits.ImageHDU(np.ones((3,2,2)), name='UNCERT'),
        fits.ImageHDU(np.zeros((3,2,2)), name='FLAGS')]).writeto(source)
    par = tmp_path/'blue.par'
    par.write_text('dimension 2 2\n')
    shift = tmp_path/'blue.shift.list'
    shift.write_text(f'{source} 0 0\n')
    recipe_path = tmp_path/'recipe.json'
    recipe_path.write_text(json.dumps({
        'schema':'crd-kcwikit-replay-v1', 'arm':'BL',
        'par_file':str(par), 'shift_file':str(shift), 'reference_icube':None,
        'exposures':[{'path':str(source), 'trim':[0,0]}],
        'stack_kwargs':{'drizzle':0.9, 'method':'drizzle'}, 'sidecars':{},
    }))
    env = dict(os.environ, PYTHONPATH=str(tmp_path))
    script = str(Path(__file__).with_name('stack_noise.py'))
    science = tmp_path/'science'
    subprocess.run([sys.executable, script, '--recipe', str(recipe_path),
                    '--output-dir', str(science), '--science-only'],
                   env=env, check=True, capture_output=True, text=True)
    noise = tmp_path/'noise'
    subprocess.run([sys.executable, script, '--recipe', str(science/'recipe_for_noise.json'),
                    '--output-dir', str(noise), '--realizations', '4', '--workers', '2'],
                   env=env, check=True, capture_output=True, text=True)
    manifest = json.loads((noise/'noise_manifest.json').read_text())
    assert manifest['complete'] and manifest['workers'] == 2
    assert len(manifest['noise_icubes']) == 4
    assert manifest['noise_icubes'] == sorted(manifest['noise_icubes'])
    serial = tmp_path/'serial_noise'
    subprocess.run([sys.executable, script, '--recipe', str(science/'recipe_for_noise.json'),
                    '--output-dir', str(serial), '--realizations', '4', '--workers', '1'],
                   env=env, check=True, capture_output=True, text=True)
    for parallel_path in manifest['noise_icubes']:
        serial_path = serial/Path(parallel_path).parent.name/'stack_icubes.fits'
        np.testing.assert_array_equal(fits.getdata(parallel_path), fits.getdata(serial_path))
