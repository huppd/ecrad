#!/usr/bin/env python3
#-*- coding: utf-8 -*-

__author__ = "Mikhail Zhigun"

import config
from typing import Optional, NamedTuple
import subprocess, os, tempfile, argparse, sys
from contextlib import contextmanager
assert sys.version_info[0] >= 3 and sys.version_info[1] >= 6, 'Python >= 3.6 is required'

TMPFS_DIR = '/dev/shm'


@contextmanager
def prepare_dir(dir_path=None, parent_dir=None):
    if dir_path is not None:
        os.makedirs(dir_path, exist_ok=True)
        yield dir_path
    else:
        d = tempfile.TemporaryDirectory(dir=parent_dir)
        yield d.name
        d.cleanup()


class TestcaseCFG(NamedTuple):
    openmp_enabled: bool
    block_size: int

DEFAULT_TEST_CASE_CFG = TestcaseCFG(True, 32)


def run_testcase(name: str, solver_type: str, cfg: TestcaseCFG, data_multiplier: int = 1, output_file_path: Optional[str] = None,
                 output_report_path: Optional[str] = None, reference_file_path: Optional[str] = None,
                 output_cfg_path: Optional[str] = None, working_dir: Optional[str] = None) -> bool:
    in_working_dir = working_dir
    in_file_path = os.path.join(config.TEST_DATA_INPUT_DIR, name + '.nc')
    with prepare_dir(in_working_dir, TMPFS_DIR) as working_dir:
        basename = f'{name}_{solver_type}_{data_multiplier}'
        if output_file_path is None:
            output_file_path = os.path.join(working_dir, f'{basename}.nc')
        if output_report_path is None:
            output_report_path = os.path.join(working_dir, f'{basename}_report.xml')
        if output_cfg_path is None:
            output_cfg_path = os.path.join(working_dir, f'{basename}_config.nam')
        assert os.path.exists(config.INPUT_DRIVER_CFG_TEMPLATE_FILEPATH), \
            f'Input cfg template "{config.INPUT_DRIVER_CFG_TEMPLATE_FILEPATH}" not found'
        cfg_file_txt = None
        with open(config.INPUT_DRIVER_CFG_TEMPLATE_FILEPATH, 'r') as f:
            cfg_file_txt = f.read()
        cfg_file_txt = cfg_file_txt.format(openmp_enabled=cfg.openmp_enabled,
                                           block_size=cfg.block_size,
                                           output_report_path=output_report_path,
                                           output_cfg_path=output_cfg_path,
                                           cfg_dir_path=config.RADIATION_CFG_DIR,
                                           sw_solver_type=solver_type,
                                           lw_solver_type=solver_type)
        with open(output_cfg_path, 'w') as f:
            f.write(cfg_file_txt)
        args = (config.DRIVER_BIN, output_cfg_path, in_file_path, output_file_path)
        p = subprocess.run(args=args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=working_dir)

        def on_fail(reason: str, p=p, args=args):
            print(f'ECRAD driver run failed. Reason: {reason}')
            print(f'args: {args}')
            stderr = p.stderr.decode('utf-8').lower()
            print(f'stderr: {stderr}')

        if p.returncode != 0:
            on_fail('return code not zero')
            return False
        if not os.path.exists(output_file_path):
            on_fail('output file missing')
            return False
        if not os.path.exists(output_file_path):
            on_fail('report missing')
            return False
        return True


def test_case_ref_file_path(test_case_name: str, solver_type: str):
    return os.path.join(config.TEST_DATA_REFERENCE_DIR, f'{test_case_name}_{solver_type}_result.nc')


def reset_ci_tests():
    for test_case_name in config.TEST_CASES:
        for solver_type in config.SOLVER_TYPES:
            out_file_path = test_case_ref_file_path(test_case_name, solver_type)
            print(f'Generating reference file for test case: {test_case_name} solver: {solver_type}')
            run_testcase(test_case_name, solver_type, DEFAULT_TEST_CASE_CFG, output_file_path=out_file_path)


def main():
    parser = argparse.ArgumentParser(description="ECRAD automated tests runner")
    parser.add_argument("--reset-ci-tests", help="Generate reference results for ci testcases",
                        action='store_true')
    args = parser.parse_args()
    if args.reset_ci_tests:
        reset_ci_tests()


if __name__ == '__main__':
    main()