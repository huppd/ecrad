#!/usr/bin/env python3
#-*- coding: utf-8 -*-

__author__ = "Mikhail Zhigun"

import sys
assert sys.version_info[0] >= 3 and sys.version_info[1] >= 6, 'Python >= 3.6 is required'
import config
from typing import Optional, NamedTuple
import subprocess, os, tempfile, argparse
from contextlib import contextmanager
import xarray as xr
import numpy as np

TMPFS_DIR = '/dev/shm'

def file_exists(path: str):
    return os.path.exists(path) and os.path.isfile(path)

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


def run_testcase(name: str, solver_type: str, cfg: TestcaseCFG, data_multiplier: int = 1,
                 output_file_path: Optional[str] = None, output_report_path: Optional[str] = None,
                 reference_file_path: Optional[str] = None, cmp_abs_tolerance: float=0.0001,
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
        assert file_exists(config.INPUT_DRIVER_CFG_TEMPLATE_FILEPATH), \
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
            print(f'ECRAD driver run failed! Reason: {reason}')
            print(f'args: {args}')
            stderr = p.stderr.decode('utf-8').lower()
            print(f'stderr: {stderr}')

        if p.returncode != 0:
            on_fail('return code not zero')
            return False
        if not file_exists(output_file_path):
            on_fail('output file missing')
            return False
        if not file_exists(output_file_path):
            on_fail('report missing')
            return False
        if reference_file_path is not None:
            if not file_exists(reference_file_path):
                print(f'Reference file {reference_file_path} not found')
            if not nc_equal(reference_file_path, output_file_path, cmp_abs_tolerance):
                print(f'Output file {output_file_path} is not equal to {reference_file_path} within absolute tolerance'+
                      f' {cmp_abs_tolerance}')
                return False
        return True


def tol_test_case_ref_file_path(test_case_name: str, solver_type: str):
    return os.path.join(config.TOLERANCE_TEST_DATA_REFERENCE_DIR, f'{test_case_name}_{solver_type}_result.nc')


def reset_ci_tests():
    for test_case_name in config.TEST_CASES:
        for solver_type in config.SOLVER_TYPES:
            out_file_path = tol_test_case_ref_file_path(test_case_name, solver_type)
            print(f'Generating reference file for test case: {test_case_name} solver: {solver_type}')
            run_testcase(test_case_name, solver_type, DEFAULT_TEST_CASE_CFG, output_file_path=out_file_path)


def run_ci_tests():
    run_ci_tol_tests()


def run_ci_tol_tests():
    print('Tolerance tests:')
    for test_case_name in config.TEST_CASES:
        for solver_type in config.SOLVER_TYPES:
            ref_file_path = tol_test_case_ref_file_path(test_case_name, solver_type)
            print(f'\t{test_case_name} solver: {solver_type}')
            run_testcase(test_case_name, solver_type, DEFAULT_TEST_CASE_CFG, reference_file_path=ref_file_path)


def nc_equal(ref_filename: str, in_filename: str, abs_threshold: float=0.0001) -> bool:
    ref = xr.open_dataset(ref_filename)
    input = xr.open_dataset(in_filename)
    # Compare dimensions
    ref_dim_names = sorted(ref.dims.keys())
    in_dim_names = sorted(input.dims.keys())
    if ref_dim_names != in_dim_names:
        return False
    for dim_name in ref_dim_names:
        if ref.dims[dim_name] != input.dims[dim_name]:
            return False
    # Compare data
    ref_var_names = sorted(ref.data_vars.keys())
    in_var_names = sorted(input.data_vars.keys())
    if ref_var_names != in_var_names:
        return False
    for var_name in ref_var_names:
        ref_var_data = ref.data_vars[var_name]
        in_var_data = input.data_vars[var_name]
        if ref_var_data.name != in_var_data.name:
            return False
        if ref_var_data.dims != in_var_data.dims:
            return False
        if ref_var_data.sizes != in_var_data.sizes:
            return False
        if not np.allclose(ref_var_data, in_var_data, atol=abs_threshold, equal_nan=True):
            return False
    return True


def multiply_dataset(in_filename: str, out_filename: str, n_times: int):
    """Enlarges dataset by copying it n_times"""
    assert False, 'Implementation not finished'
    assert n_times > 0
    data = xr.open_dataset(in_filename)
    out = xr.Dataset()

    new_dims = {key: val for key, val in data.dims.items() if key != "column"}

    out = out.expand_dims(new_dims)
    out = out.expand_dims({"column": n})

    lon = data["longitude"][:]
    lat = data["latitude"][:]

    out["longitude"] = (("column",), np.linspace(np.min(lon), np.max(lon), n))
    out["latitude"] = (("column",), np.linspace(np.min(lat), np.max(lat), n))

    n_repeat = n // data.dims["column"] + 1

    for vn in data.data_vars:
        var = data[vn]
        if vn not in ["longitude", "latitude"] and "column" in var.dims:
            colax = var.dims.index("column")
            out[vn] = (var.dims, np.repeat(var[:], n_repeat, axis=colax)[:n])

    out.to_netcdf(out_filename)


def main():
    parser = argparse.ArgumentParser(description="ECRAD automated tests runner")
    parser.add_argument("--reset-ci-tests", help="Generate reference results for ci testcases", action='store_true')
    parser.add_argument("--run-ci-tests", help="Run continuous integration tests", action='store_true')
    parser.add_argument("--run-ci-tol-tests", help="Run continuous integration tolerance tests", action='store_true')
    args = parser.parse_args()
    if args.reset_ci_tests:
        reset_ci_tests()
    if args.run_ci_tests:
        run_ci_tests()
    if args.run_ci_tol_tests:
        run_ci_tol_tests()


if __name__ == '__main__':
    main()