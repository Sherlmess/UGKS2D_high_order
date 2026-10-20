#!/usr/bin/env python3
"""Release-only two-step cold-start checks; outputs go to an ignored fresh directory."""
import os, pathlib, re, shutil, subprocess, json, time
root=pathlib.Path(__file__).resolve().parents[1]
base=root/'runs'/time.strftime('release_smoke_%Y%m%dT%H%M%S')
base.mkdir(parents=True,exist_ok=False)
results={}
for method in ['high_implicit','low_implicit','low_explicit']:
 d=base/method;shutil.copytree(root/'cases'/method,d)
 p=d/'Config.ini';text=p.read_text()
 for key in ['MAX_ITER','Implicit_MaxIter']:
  text,n=re.subn(r'(?m)^(\s*'+key+r'\s*=).*$',r'\1 2',text);assert n==1
 p.write_text(text)
 env=os.environ.copy();env.update(MPI_OVERSUBSCRIBE='1',R5_WALL_LIMIT='60')
 start=time.monotonic()
 try:
  cp=subprocess.run(['bash',str(root/'run.sh'),str(d)],cwd=root,env=env,timeout=90,capture_output=True,text=True)
  result={'exit_code':cp.returncode,'wall_seconds':time.monotonic()-start,'stderr':cp.stderr}
 except subprocess.TimeoutExpired:
  result={'exit_code':'TIMEOUT','wall_seconds':time.monotonic()-start}
 results[method]=result
 (base/'summary.json').write_text(json.dumps(results,indent=2)+'\n')
 print(method,result,flush=True)
 if result['exit_code']!=0:raise SystemExit(1)
print(base)
