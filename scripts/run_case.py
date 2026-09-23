#!/usr/bin/env python3
"""Run one case in a NEW directory. No hidden restart or launch environment."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import shutil
import signal
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
CASES = json.loads((ROOT/'examples/cases.json').read_text())

def steady_finished_in_nine(text):
    rows = [line.split() for line in text.splitlines()
            if line.split() and line.split()[0] == 'S2Steady']
    return len(rows) == 1 and rows[0][:3] == ['S2Steady', 'T', '9']

def run_case(case, output):
    spec = CASES[case]
    source = ROOT/'examples'/case
    output = Path(output).resolve()
    binary = ROOT/'bin/ugks_solver'
    if not binary.is_file():
        raise RuntimeError('Run python3 scripts/build.py first')
    output.mkdir(parents=True, exist_ok=False)
    for name in ('Config.ini', 'Mesh2d.dat'):
        shutil.copy2(source/name, output/name)
    env = os.environ.copy()
    for key in list(env):
        if key.startswith(('S2I1_', 'S2DO1_', 'S2D3_')):
            del env[key]
    env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    if spec['restart']:
        env['S2DO1_RESTART_PATH'] = str(source/'restart.bin')
    command = ['mpirun', '-np', '1', str(binary)]
    info = {'case':case, 'argv':command, 'cwd':str(output),
            'restart':env.get('S2DO1_RESTART_PATH'), 'timeout_s':spec['max_seconds'],
            'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest()}
    (output/'launch.json').write_text(json.dumps(info,indent=2)+'\n')
    start = time.monotonic()
    with (output/'stdout.log').open('w') as log:
        proc = subprocess.Popen(command, cwd=output, env=env, stdout=log,
                                stderr=subprocess.STDOUT, start_new_session=True)
        try:
            rc = proc.wait(timeout=spec['max_seconds'])
            info['timed_out'] = False
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid, signal.SIGKILL)
            rc = proc.wait()
            info['timed_out'] = True
    info.update(exit_code=rc,wall_s=time.monotonic()-start)
    passed = rc == 0 and not info['timed_out']
    if passed:
        final = output/'native_state_final.bin'
        passed = final.is_file()
        if passed:
            info['final_sha256'] = hashlib.sha256(final.read_bytes()).hexdigest()
        if spec['expected_final'] and passed:
            info['final_bytewise_equal'] = final.read_bytes() == (source/'expected_final.bin').read_bytes()
            passed = info['final_bytewise_equal']
        if case == 'steady_finish':
            status = output/'ConvergenceStatus.txt'
            info['convergence_status'] = status.read_text().strip() if status.exists() else None
            info['nine_corrections_and_converged'] = steady_finished_in_nine(info['convergence_status'] or '')
            passed = passed and info['nine_corrections_and_converged']
    info['PASS'] = passed
    (output/'result.json').write_text(json.dumps(info,indent=2)+'\n')
    print(json.dumps(info,indent=2))
    return info

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('case', choices=sorted(CASES))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if not run_case(args.case, args.output)['PASS']:
        raise SystemExit(1)
