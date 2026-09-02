#!/usr/bin/env python3
"""Validate final CRD_DRP package and write CRD_DRP_reduction_manifest.json."""

from __future__ import annotations
import argparse,json,hashlib
from datetime import datetime,timezone
from pathlib import Path
from astropy.table import Table
from crd_drp_utils import resolve_stack_paths,load_stack,load_1d_mask_fits,scan_negative_candidates

VERSION="1.0.0"

def parser():
    p=argparse.ArgumentParser()
    for arm in ("blue","red"):
        p.add_argument(f"--{arm}-icube",required=True)
        p.add_argument(f"--{arm}-vcube"); p.add_argument(f"--{arm}-mcube"); p.add_argument(f"--{arm}-ecube")
        p.add_argument(f"--{arm}-atmospheric-mask",required=True)
        p.add_argument(f"--{arm}-atmospheric-intervals",required=True)
    p.add_argument("--stack-integrity-summary")
    p.add_argument("--atmospheric-mask-summary")
    p.add_argument("--max-mask-fraction",type=float,default=0.45)
    p.add_argument("--output-dir",default="reduction_package_validation")
    return p

def sha256(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""):h.update(chunk)
    return h.hexdigest()

def validate_arm(arm,args,out):
    paths=resolve_stack_paths(getattr(args,f"{arm}_icube"),getattr(args,f"{arm}_vcube"),getattr(args,f"{arm}_mcube"),getattr(args,f"{arm}_ecube"))
    cube=load_stack(paths)
    mp=Path(getattr(args,f"{arm}_atmospheric_mask")).expanduser().resolve()
    ip=Path(getattr(args,f"{arm}_atmospheric_intervals")).expanduser().resolve()
    if not ip.is_file():raise FileNotFoundError(ip)
    mask,_,_=load_1d_mask_fits(mp,cube.wavelength)
    frac=float(mask.mean())
    if frac>args.max_mask_fraction:raise ValueError(f"{arm.upper()} mask fraction {frac:.3f} exceeds {args.max_mask_fraction:.3f}")
    rows,struct=scan_negative_candidates(cube)
    annotated=[]; unexplained=[]; explained=[]
    for row in rows:
        r=dict(row); inside=bool(mask[r["GROUP_START_INDEX"]:r["GROUP_END_INDEX"]+1].any()); r["INSIDE_ATMOSPHERIC_MASK"]=inside
        if r["RAW_SEVERITY"]=="FAIL":(explained if inside else unexplained).append(r)
        r["ARM"]=arm.upper(); annotated.append(r)
    if annotated:Table(rows=annotated).write(out/f"{arm}_postmask_negative_candidates.ecsv",format="ascii.ecsv",overwrite=True)
    structural_fail=struct["n_positive_exposure_nonfinite_science_voxels"]>0 or struct["n_positive_exposure_negative_variance_voxels"]>0
    status="FAIL" if structural_fail or unexplained else "PASS"
    files={"icube":paths.icube,"vcube":paths.vcube,"mcube":paths.mcube,"ecube":paths.ecube,"atmospheric_mask_fits":mp,"atmospheric_intervals_ecsv":ip}
    rec={}
    for k,p in files.items():
        if p is None:rec[k]=None
        else:
            pp=Path(p).resolve(); rec[k]={"path":str(pp),"sha256":sha256(pp),"size_bytes":pp.stat().st_size}
    return {"arm":arm.upper(),"status":status,"shape_wave_y_x":list(map(int,cube.flux.shape)),"wavelength_min_angstrom":float(cube.wavelength[0]),"wavelength_max_angstrom":float(cube.wavelength[-1]),"atmospheric_mask_fraction":frac,"atmospheric_masked_native_pixels":int(mask.sum()),"n_extreme_fail_groups_inside_atmospheric_mask":len(explained),"n_unexplained_extreme_fail_groups_outside_atmospheric_mask":len(unexplained),"files":rec}

def main():
    args=parser().parse_args(); out=Path(args.output_dir).expanduser().resolve(); out.mkdir(parents=True,exist_ok=True)
    upstream={}
    upstream_fail=False
    if args.stack_integrity_summary:
        p=Path(args.stack_integrity_summary).expanduser().resolve(); d=json.loads(p.read_text()); upstream["stack_integrity_summary"]={"path":str(p),"overall_status":d.get("overall_status")}; upstream_fail=d.get("overall_status") in {"FAIL","ERROR"}
    if args.atmospheric_mask_summary:
        p=Path(args.atmospheric_mask_summary).expanduser().resolve(); d=json.loads(p.read_text()); upstream["atmospheric_mask_summary"]={"path":str(p)}
    arms={}
    for arm in ("blue","red"):
        try:arms[arm]=validate_arm(arm,args,out)
        except Exception as exc:arms[arm]={"arm":arm.upper(),"status":"ERROR","error":{"type":type(exc).__name__,"message":str(exc)}}
    states=[arms[a]["status"] for a in ("blue","red")]
    overall="FAIL" if upstream_fail or "FAIL" in states or "ERROR" in states else "PASS"
    summary={"pipeline":"CRD_DRP","script":Path(__file__).name,"version":VERSION,"created_utc":datetime.now(timezone.utc).isoformat(),"overall_status":overall,"upstream_validation":upstream,"arms":arms}
    sp=out/"reduction_package_validation.json"; sp.write_text(json.dumps(summary,indent=2)+"\n")
    manifest={"schema":"CRD_DRP_reduction_manifest_v1","pipeline":"CRD_DRP","created_utc":datetime.now(timezone.utc).isoformat(),"validation_status":overall,"validation_summary":str(sp.resolve()),"wavelength_frame_contract":"Atmospheric masks are 1-D boolean masks sampled on each final stack's native observed wavelength grid.","mask_semantics":{"true":"exclude from stellar analysis","false":"not excluded by atmospheric mask"},"blue":arms["blue"],"red":arms["red"]}
    mp=out/"CRD_DRP_reduction_manifest.json"; mp.write_text(json.dumps(manifest,indent=2)+"\n")
    print(f"\nCRD_DRP REDUCTION PACKAGE: {overall}\nBLUE: {arms['blue']['status']}\nRED:  {arms['red']['status']}\nManifest: {mp}")
    return 1 if overall=="FAIL" else 0

if __name__=="__main__":
    raise SystemExit(main())
