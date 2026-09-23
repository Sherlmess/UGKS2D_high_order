#!/usr/bin/env python3
from pathlib import Path

boundary = Path('src/mod_Boundary.f90').read_text()
start = boundary.index('subroutine bc_freestream_DIS()')
end = boundary.index('end subroutine bc_freestream_DIS', start)
body = boundary[start:end].lower()
if 'elseif(ip_y==np_y-1)' in body.replace(' ', ''):
    raise SystemExit('free-wave left and upper boundaries must use independent IF blocks')

solve = Path('src/mod_Solve.f90').read_text().replace(' ', '').lower()
if 'time_dis=dble(iter)*dt_global' in solve:
    raise SystemExit('physical disturbance time must be accumulated, not reconstructed from the current dt')

print('PASS: unsteady baseline source guards')
