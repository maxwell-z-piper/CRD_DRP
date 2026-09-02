#!/usr/bin/env python3
"""Build fixed observed-frame BLUE/RED atmospheric masks for CRD_DRP."""

from __future__ import annotations
import argparse, json
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.table import Table
from crd_drp_utils import (
    resolve_stack_paths, load_stack, extract_spatial_median_spectrum,
    robust_highpass_z, infer_wavelength_medium, convert_medium,
    bridge_small_false_gaps, dilate_1d, wavelength_edges,
    save_1d_mask_fits, resolve_path
)

VERSION="1.0.0"

def parser():
    p=argparse.ArgumentParser()
    for arm in ("blue","red"):
        p.add_argument(f"--{arm}-icube",required=True)
        p.add_argument(f"--{arm}-vcube"); p.add_argument(f"--{arm}-mcube"); p.add_argument(f"--{arm}-ecube")
        p.add_argument(f"--{arm}-reference-cube",action="append",default=[])
        p.add_argument(f"--{arm}-reference-object-mask",action="append",default=[])
        p.add_argument(f"--{arm}-science-medium",choices=("auto","air","vacuum"),default="auto")
        p.add_argument(f"--{arm}-reference-medium",choices=("auto","air","vacuum"),default="auto")
    p.add_argument("--known-line-list",default=None)
    p.add_argument("--sky-line-sigma",type=float,default=8.0)
    p.add_argument("--sky-reference-min-fraction",type=float,default=0.50)
    p.add_argument("--continuum-window-angstrom",type=float,default=75.0)
    p.add_argument("--sky-padding-pixels",type=int,default=2)
    p.add_argument("--bridge-pixels",type=int,default=1)
    p.add_argument("--known-line-half-width-angstrom",type=float,default=3.5)
    p.add_argument("--output-dir",default="atmospheric_masks")
    return p

def built_in_lines():
    return [(5577.338,"[O I] 5577","air"),(6300.304,"[O I] 6300","air"),(6363.776,"[O I] 6364","air")]

def load_extra_lines(path):
    if path is None: return []
    t=Table.read(path); low={n.lower():n for n in t.colnames}
    wcol=next((v for k,v in low.items() if "wave" in k or "lambda" in k),None)
    lcol=next((v for k,v in low.items() if "label" in k or "name" in k),None)
    mcol=next((v for k,v in low.items() if "medium" in k),None)
    if wcol is None: raise ValueError("Known-line table needs a wavelength column.")
    return [(float(r[wcol]),str(r[lcol]) if lcol else "user sky line",str(r[mcol]).lower() if mcol else "air") for r in t]

def resolve_masks(refs,masks):
    if not masks:return [None]*len(refs)
    if len(masks)==1:return [resolve_path(masks[0])]*len(refs)
    if len(masks)!=len(refs):raise ValueError("Supply zero, one shared object mask, or one mask per reference cube.")
    return [resolve_path(m) for m in masks]

def intervals(wave,mask,empirical,known_labels):
    edges=wavelength_edges(wave); idx=np.flatnonzero(mask); rows=[]
    if idx.size:
        starts=np.r_[idx[0],idx[1:][np.diff(idx)>1]]; ends=np.r_[idx[:-1][np.diff(idx)>1],idx[-1]]
        for iid,(i0,i1) in enumerate(zip(starts,ends),1):
            labels=sorted({lab for labs in known_labels[i0:i1+1] for lab in labs})
            rows.append({"INTERVAL_ID":iid,"INDEX_LO":int(i0),"INDEX_HI":int(i1),"N_PIXELS":int(i1-i0+1),"WAVE_LO_ANGSTROM":float(edges[i0]),"WAVE_HI_ANGSTROM":float(edges[i1+1]),"REASON":"empirical_sky+known_sky_line" if labels and np.any(empirical[i0:i1+1]) else ("known_sky_line" if labels else "empirical_sky"),"KNOWN_LINES":"; ".join(labels)})
    return Table(rows=rows)

def build_arm(arm,args,lines,out):
    paths=resolve_stack_paths(getattr(args,f"{arm}_icube"),getattr(args,f"{arm}_vcube"),getattr(args,f"{arm}_mcube"),getattr(args,f"{arm}_ecube"))
    cube=load_stack(paths); wave=cube.wavelength
    sm=getattr(args,f"{arm}_science_medium"); sm=infer_wavelength_medium(cube.header,"vacuum") if sm=="auto" else sm
    refs=[resolve_path(p) for p in getattr(args,f"{arm}_reference_cube")]
    omasks=resolve_masks(refs,getattr(args,f"{arm}_reference_object_mask"))
    zrefs=[]; meta=[]
    for p,om in zip(refs,omasks):
        rw,rs,hdr,nsky=extract_spatial_median_spectrum(p,om)
        rm=getattr(args,f"{arm}_reference_medium"); rm=infer_wavelength_medium(hdr,sm) if rm=="auto" else rm
        rw=convert_medium(rw,rm,sm)
        _,z,sig=robust_highpass_z(rw,rs,args.continuum_window_angstrom)
        zrefs.append(np.interp(wave,rw,z,left=np.nan,right=np.nan))
        meta.append({"path":str(p),"object_mask":str(om) if om else None,"n_sky_spaxels":nsky,"medium":rm,"robust_sigma":sig})
    if zrefs:
        zarr=np.vstack(zrefs); detections=zarr>=args.sky_line_sigma
        frac=np.nanmean(detections.astype(float),axis=0); medz=np.nanmedian(zarr,axis=0)
        empirical=frac>=args.sky_reference_min_fraction
    else:
        frac=np.zeros_like(wave); medz=np.full_like(wave,np.nan); empirical=np.zeros_like(wave,bool)
    empirical=dilate_1d(bridge_small_false_gaps(empirical,args.bridge_pixels),args.sky_padding_pixels)
    known=np.zeros_like(wave,bool); labels=np.empty(wave.size,dtype=object)
    for i in range(wave.size):labels[i]=[]
    for lw,label,medium in lines:
        center=float(convert_medium([lw],medium,sm)[0]); sel=np.abs(wave-center)<=args.known_line_half_width_angstrom
        known|=sel
        for i in np.flatnonzero(sel):labels[i].append(label)
    final=empirical|known
    save_1d_mask_fits(out/f"{arm}_atmospheric_mask.fits",final,wave,arm,sm,paths.icube)
    tab=intervals(wave,final,empirical,labels); tab.meta["PIPELINE"]="CRD_DRP"; tab.meta["ARM"]=arm.upper(); tab.write(out/f"{arm}_atmospheric_mask.ecsv",format="ascii.ecsv",overwrite=True)
    ref=Table(); ref["WAVELENGTH_ANGSTROM"]=wave; ref["EMPIRICAL_SKY_Z"]=medz; ref["REFERENCE_DETECTION_FRACTION"]=frac; ref["EMPIRICAL_MASK"]=empirical; ref["KNOWN_LINE_MASK"]=known; ref["FINAL_MASK"]=final; ref.write(out/f"{arm}_atmospheric_reference.ecsv",format="ascii.ecsv",overwrite=True)
    fig,ax=plt.subplots(figsize=(13,5)); ax.plot(wave,medz,lw=.8); ax.fill_between(wave,0,1,where=final,alpha=.2,transform=ax.get_xaxis_transform()); ax.set_xlabel("Observed wavelength [A]"); ax.set_ylabel("Empirical sky z"); ax.set_title(f"{arm.upper()} atmospheric mask: {100*np.mean(final):.2f}% native pixels masked"); fig.tight_layout(); fig.savefig(out/f"{arm}_atmospheric_mask_qc.png",dpi=160); plt.close(fig)
    return {"arm":arm.upper(),"science_medium":sm,"reference_cubes":meta,"mask_fits":str((out/f"{arm}_atmospheric_mask.fits").resolve()),"mask_intervals_ecsv":str((out/f"{arm}_atmospheric_mask.ecsv").resolve()),"n_native_pixels":int(final.size),"n_masked_pixels":int(np.count_nonzero(final)),"masked_fraction":float(np.mean(final)),"n_intervals":len(tab)}

def main():
    args=parser().parse_args(); out=Path(args.output_dir).expanduser().resolve(); out.mkdir(parents=True,exist_ok=True)
    lines=built_in_lines()+load_extra_lines(args.known_line_list)
    arms={arm:build_arm(arm,args,lines,out) for arm in ("blue","red")}
    payload={"pipeline":"CRD_DRP","script":Path(__file__).name,"version":VERSION,"created_utc":datetime.now(timezone.utc).isoformat(),"selection_principle":"Fixed observed-frame masks constructed independently of CRD_DAP/pPXF residuals.","arms":arms,"parameters":vars(args)}
    (out/"atmospheric_mask_summary.json").write_text(json.dumps(payload,indent=2)+"\n")
    for arm in ("blue","red"):
        s=arms[arm]; print(f"{arm.upper()}: {s['n_masked_pixels']}/{s['n_native_pixels']} ({100*s['masked_fraction']:.2f}%) masked")
    return 0

if __name__=="__main__":
    raise SystemExit(main())
