
# Final KcwiKit stack validation

Run this immediately after KcwiKit stacking and before CRD_DAP Script 01.

```bash
python validate_kcwikit_stacks.py \
    --blue-icube /path/to/blue_icubes.fits \
    --red-icube  /path/to/red_icubes.fits \
    --output-dir stack_validation
```

With standard KcwiKit filenames, matching `_vcubes.fits`, `_mcubes.fits`, and
`_ecubes.fits` files are inferred automatically.

The script is specifically intended to catch the finite, orders-of-magnitude
negative artifacts previously caused by non-finite trimmed samples entering the
KcwiKit int32/reprojection path.

A negative value is **not** automatically a failure. By default, a catastrophic
candidate must satisfy both:

- local robust z <= -50
- negative depth >= 25 times the local absolute flux scale

These thresholds are deliberately conservative relative to ordinary reduction
residuals.

Outputs:

- `stack_validation_summary.json`
- `stack_validation_candidates.ecsv`
- `blue_minimum_flux_map.png`
- `red_minimum_flux_map.png`
- `blue_catastrophic_count_map.png`
- `red_catastrophic_count_map.png`
- `*_worst_negative_candidates.png`

Exit code `1` means the reduction should not proceed to CRD_DAP science
analysis until the stack is investigated.
