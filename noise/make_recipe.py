#!/usr/bin/env python3
"""Create a fixed science/noise recipe from the exact KcwiKit stack list.

Each non-comment row in the list is ``cube_root lower_trim upper_trim``.
``cube_root`` must name the final exposure, including ``_zap`` when that is
the product to stack. KcwiKit appends ``_icubes.fits`` to this root.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from stack_noise import validate_recipe


def recipe_from_list(list_path: Path, arm: str, *, drizzle: float = 0.9,
                     nsigma_clip: float = 10, npix_trim: int = 3) -> dict:
    list_path = list_path.expanduser().resolve()
    if not list_path.is_file() or not list_path.name.endswith('.list'):
        raise ValueError(f'Expected an existing KcwiKit .list file: {list_path}')
    stem = str(list_path)[:-5]
    exposures = []
    for line_number, line in enumerate(list_path.read_text().splitlines(), start=1):
        line = line.split('#', 1)[0].strip()
        if not line:
            continue
        fields = line.split()
        if len(fields) != 3:
            raise ValueError(f'{list_path}:{line_number}: expected cube root and two trim integers')
        root, lower, upper = fields
        if root.endswith('_icubes.fits'):
            raise ValueError(f'{list_path}:{line_number}: give the cube root without _icubes.fits')
        try:
            trim = [int(lower), int(upper)]
        except ValueError as exc:
            raise ValueError(f'{list_path}:{line_number}: trims must be integers') from exc
        root_path = Path(root).expanduser()
        if not root_path.is_absolute():
            root_path = list_path.parent / root_path
        exposures.append({'path': str(Path(str(root_path.resolve()) + '_icubes.fits')),
                          'trim': trim})
    flux_sidecar = Path(stem + '.flx.list')
    recipe = {
        'schema': 'crd-kcwikit-replay-v1',
        'arm': arm,
        'par_file': stem + '.par',
        'shift_file': stem + '.shift.list',
        'reference_icube': None,
        'exposures': exposures,
        'stack_kwargs': {'drizzle': drizzle, 'method': 'drizzle',
                         'nsigma_clip': nsigma_clip, 'npix_trim': npix_trim},
        'sidecars': {'.flx.list': str(flux_sidecar)} if flux_sidecar.is_file() else {},
    }
    # Also checks that every final exposure exists, has uncertainty/flags,
    # and matches a saved shift row in the same order.
    validate_recipe(recipe, require_reference=False)
    return recipe


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--list', required=True, type=Path, dest='list_path')
    parser.add_argument('--arm', required=True, choices=('BL', 'RH3'))
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--drizzle', type=float, default=0.9)
    parser.add_argument('--nsigma-clip', type=float, default=10)
    parser.add_argument('--npix-trim', type=int, default=3)
    args = parser.parse_args(argv)
    if not 0 < args.drizzle <= 1 or args.nsigma_clip <= 0 or args.npix_trim < 0:
        parser.error('Require 0 < drizzle <= 1, nsigma-clip > 0, and npix-trim >= 0')
    recipe = recipe_from_list(args.list_path, args.arm, drizzle=args.drizzle,
                              nsigma_clip=args.nsigma_clip, npix_trim=args.npix_trim)
    output = args.output.expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('x') as stream:
        json.dump(recipe, stream, indent=2)
        stream.write('\n')
    print(f"Wrote {output} from {len(recipe['exposures'])} listed {args.arm} exposures")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
