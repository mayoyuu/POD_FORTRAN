#!/usr/bin/env python3
"""FPM runner: execute a compiled test in an explicitly prepared private fixture tree."""
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def main():
    executable=Path(sys.argv[1]).resolve(strict=True)
    sandbox=Path(os.environ['POD_TEST_SANDBOX']).resolve(strict=True)
    if sandbox.name.startswith('p4-native-') is False or not (sandbox/'fpm.toml').is_file():
        raise ValueError('expected explicitly prepared p4-native fixture tree')
    if (sandbox/'SRP').is_symlink() or (sandbox/'OPM').is_symlink():
        raise ValueError('test outputs must not target historical SRP/OPM symlinks')
    logs=sandbox/'regression_logs'; logs.mkdir(exist_ok=True)
    start=time.monotonic()
    log_path=logs/(executable.name+'.log')
    with log_path.open('x') as log:
        try:
            result=subprocess.run(['timeout','--kill-after=10',os.environ.get('POD_TEST_TIMEOUT','300'),
                str(executable),*sys.argv[2:]],cwd=sandbox,stdout=log,stderr=subprocess.STDOUT)
            status=result.returncode
        except OSError as exc:
            log.write(str(exc)); status=127
    record={'test':executable.name,'returncode':status,'elapsed_seconds':time.monotonic()-start,
            'cwd':str(sandbox),'log':str(log_path),'timeout_seconds':os.environ.get('POD_TEST_TIMEOUT','300')}
    with (logs/'results.jsonl').open('a') as stream: stream.write(json.dumps(record)+'\n')
    print(json.dumps(record),flush=True)
    if status: print(log_path.read_text()[-2500:],flush=True)
    return status


if __name__=='__main__': sys.exit(main())
