from pathlib import Path
import os,subprocess,json,datetime,hashlib,sys
old=Path("/home/jony/multifoil-runs/s702-repeatability-20260916")
rt=Path("/home/jony/multifoil-runs/s702-policy-20260921/handover/policy_runtime")
env=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",MKL_NUM_THREADS="1",NUMEXPR_NUM_THREADS="1")
command=[str(old/".venv/bin/python"),"-m","unittest","-v","test_viscous.ViscousPhysics.test_multiple_elements_order_and_far_limit"]
r=subprocess.run(command,cwd=old,env=env,capture_output=True,text=True)
(rt/"original_far_field.log").write_text(r.stdout+r.stderr)
receipt={"completed_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"command":command,"cwd":str(old),"exit_code":r.returncode,"log":"original_far_field.log","sha256":{n:hashlib.sha256((old/n).read_bytes()).hexdigest() for n in ("multifoil.py","viscous.py","optimization.py","paneling.py","vendor/mfoil.py","test_viscous.py")}}
(rt/"original_far_field.json").write_text(json.dumps(receipt,indent=2)+"\n")
