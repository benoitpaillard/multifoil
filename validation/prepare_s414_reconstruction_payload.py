#!/usr/bin/env python3
"""Package all inputs and verify preparation in an isolated checkout."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile

import run_s414_reconstruction as run

LAUNCH=run.ROOT/'validation/results/s414_reconstruction_launch_20261002'
SHELL='''#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
exec 9>launch.lock
flock -n 9
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
.venv/bin/python -m unittest -v test_s414_reconstruction test_fixed_lift_policy test_optimization test_paneling > preflight.log 2>&1
exec .venv/bin/python validation/run_s414_reconstruction.py --workers 4
'''


def geometry_preflight():
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from s414_geometry import elements, SCENARIOS
    from viscous import new_model
    from multifoil import validate_geometry
    variants={s:dict(elements(s)) for s in SCENARIOS}
    rows=[]
    for scenario,geometry in variants.items():
        models=[new_model(name,xy,1e6,.1,9.,(.25,0.),1e-5) for name,xy in geometry.items()]
        validate_geometry([(m.geom.name,m.foil.x.T[::-1]) for m in models])
        for name,xy in geometry.items():
            model=next(m for m in models if m.geom.name==name)
            if model.regularization!=0. or model.geom.chord!=1.:
                raise ValueError('Unexpected edge opening or changed reference chord')
            rows.append(dict(scenario=scenario,element=name,nodes=len(xy),
                finite_TE_thickness_c=float(np.linalg.norm(xy[0]-xy[-1])),
                solver_TE_thickness_c=float(np.linalg.norm(model.foil.x[:,0]-model.foil.x[:,-1])),
                solver_TE_regularization=model.regularization,solver_reference_chord=model.geom.chord,
                maximum_corresponding_sample_displacement_c=float(np.linalg.norm(xy-variants['nominal'][name],axis=1).max()),
                nonintersecting=True,solver_geometry_nonintersecting=True))
    run.write(LAUNCH/'geometry_preflight.json',dict(scenarios=rows,
        note='Source interpolation and placement sensitivities; not statistical uncertainty bounds. No flow solve.'))
    fig,axes=plt.subplots(1,3,figsize=(13,4),layout='constrained')
    for name,xy in variants['nominal'].items():axes[0].plot(xy[:,0],xy[:,1],label=name)
    axes[0].set(title='Nominal modified test outline',xlabel='x/c',ylabel='z/c');axes[0].axis('equal');axes[0].legend()
    for name,xy in variants['source_cubic'].items():
        delta=np.linalg.norm(xy-variants['nominal'][name],axis=1)*457.2
        axes[1].plot(xy[:,0],delta,'.',ms=2,label=name)
    axes[1].set(title='Cubic versus PCHIP samples',xlabel='x/c',ylabel='Corresponding-point displacement (mm)');axes[1].legend()
    for scenario in ('nominal','aft_z_minus','aft_z_plus'):
        xy=variants[scenario]['aft'];axes[2].plot(xy[:,0],xy[:,1],label=scenario)
    aft=variants['nominal']['aft'];nose=aft[np.argmin(aft[:,0])]
    axes[2].set(title='Aft nose: ±0.13 mm shifts',xlabel='x/c',ylabel='z/c',
        xlim=(nose[0]-.0008,nose[0]+.013),ylim=(nose[1]-.006,nose[1]+.006));axes[2].legend(fontsize=8)
    for ax in axes:ax.grid(alpha=.2)
    fig.suptitle('S414 declared input sensitivities · finite trailing edges retained')
    fig.savefig(LAUNCH/'geometry_preflight.png',dpi=160);plt.close(fig)


def main():
    LAUNCH.mkdir(exist_ok=True)
    if (LAUNCH/'launch_receipt.json').exists():
        raise ValueError('Preserve launched payload; receipt already exists')
    geometry_preflight()
    spec=run.specification();run.prepare(run.ROOT/run.STUDY,spec)
    names=list(spec['sha256'])+['vendor/__init__.py','LICENSE','THIRD_PARTY.md',
        'validation/tn808_geometry.py','validation/data/nlr7301/NLR7301MainECARP.dat',
        'validation/data/nlr7301/NLR7301FlapECARP.dat','validation/data/nlr7301/LICENSE',
        'validation/monitor_s414_reconstruction.py',
        'validation/prepare_s414_reconstruction_payload.py',
        run.QUAL+'/FINDINGS.md',run.STUDY+'/study.json']
    files={name:run.ROOT/name for name in names}
    shell=LAUNCH/'run_s414_reconstruction.sh';shell.write_text(SHELL)
    files['run_s414_reconstruction.sh']=shell
    manifest={name:run.sha(path) for name,path in sorted(files.items())}
    run.write(LAUNCH/'payload_manifest.json',manifest)
    with tempfile.TemporaryDirectory(prefix='s414_reconstruction_payload_') as tmp:
        root=Path(tmp)
        for name,path in files.items():
            dest=root/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,dest)
        check=subprocess.run([sys.executable,'validation/run_s414_reconstruction.py','--prepare-only'],
            cwd=root,text=True,capture_output=True)
        (LAUNCH/'isolated_payload_check.log').write_text(check.stdout+check.stderr)
        if check.returncode:raise RuntimeError('Isolated payload preparation failed')
        expected=run.sha(run.ROOT/run.STUDY/'study.json')
        if run.sha(root/run.STUDY/'study.json')!=expected:raise ValueError('Isolated manifest differs')
        tests=subprocess.run([sys.executable,'-m','unittest','-v','test_s414_reconstruction',
            'test_fixed_lift_policy','test_optimization','test_paneling'],
            cwd=root,text=True,capture_output=True)
        (LAUNCH/'isolated_payload_tests.log').write_text(tests.stdout+tests.stderr)
        if tests.returncode:raise RuntimeError('Isolated payload tests failed')
    archive=LAUNCH/'s414_reconstruction_payload.tar.gz'
    with tarfile.open(archive,'w:gz') as tar:
        for name,path in sorted(files.items()):tar.add(path,arcname=name)
        tar.add(LAUNCH/'payload_manifest.json',arcname='payload_manifest.json')
    receipt=dict(files=len(files),archive_bytes=archive.stat().st_size,archive_sha256=run.sha(archive),
        study_sha256=expected,isolated_payload_preparation_pass=True,isolated_payload_tests_pass=True)
    run.write(LAUNCH/'preparation.json',receipt);print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
