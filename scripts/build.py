#!/usr/bin/env python3
"""Build only this snapshot; standard-library-only driver."""
from pathlib import Path
import hashlib
import json
import os
import platform
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]

def build():
    env = os.environ.copy()
    if platform.system() == 'Darwin':
        sdk = subprocess.check_output(['xcrun', '--show-sdk-path'], text=True).strip()
        env['SDKROOT'] = sdk
        env['LIBRARY_PATH'] = sdk + '/usr/lib' + (':' + env['LIBRARY_PATH'] if env.get('LIBRARY_PATH') else '')
    out = ROOT / 'verification'
    out.mkdir(exist_ok=True)
    start = time.monotonic()
    with (out / 'build.log').open('w') as log:
        proc = subprocess.run(['make', 'release', 'FC=mpifort'], cwd=ROOT, env=env,
                              stdout=log, stderr=subprocess.STDOUT, timeout=180)
    report = {'exit_code': proc.returncode, 'wall_s': time.monotonic()-start,
              'compiler': subprocess.check_output(['mpifort', '--version'], text=True),
              'platform': platform.platform(), 'SDKROOT': env.get('SDKROOT')}
    if proc.returncode == 0:
        report['binary_sha256'] = hashlib.sha256((ROOT/'bin/ugks_solver').read_bytes()).hexdigest()
        test = ROOT / 'bin/test_s2_termination'
        subprocess.run(['mpifort', '-I'+str(ROOT/'obj'), str(ROOT/'tests/test_s2_termination.f90'),
                        str(ROOT/'obj/mod_S2SteadyTermination.o'), '-o', str(test)],
                       cwd=ROOT, env=env, check=True, timeout=30)
        result = subprocess.check_output([str(test)], cwd=ROOT, env=env, text=True, timeout=10)
        report['termination_test'] = result.strip()
    (out/'build.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))
    if proc.returncode:
        raise SystemExit(proc.returncode)

if __name__ == '__main__':
    build()
