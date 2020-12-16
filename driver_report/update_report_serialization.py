import argparse
import os
import sys
from os.path import join as join_path

assert sys.version_info[0] >= 3 and sys.version_info[1] >= 7, 'Python >= 3.7 is required'
import xsdtools

THIS_DIR_PATH = os.path.dirname(os.path.realpath(__file__))
DRIVER_REPORT_RELATIVE_PATH = 'driver_report.xsd'
DRIVER_REPORT_FORTRAN_RELATIVE_PATH = './report'
DRIVER_REPORT_FORTRAN_SRC_DIR = join_path(THIS_DIR_PATH, DRIVER_REPORT_FORTRAN_RELATIVE_PATH)
DRIVER_REPORT_SCHEMA_FILE = join_path(THIS_DIR_PATH, DRIVER_REPORT_RELATIVE_PATH)


def file_exists(path: str):
    return os.path.exists(path) and os.path.isfile(path)


def dir_exists(path: str):
    return os.path.exists(path) and os.path.isdir(path)


if __name__ == "__main__":
    assert file_exists(DRIVER_REPORT_SCHEMA_FILE), f'Report schema file "{DRIVER_REPORT_SCHEMA_FILE}" not found'
    assert dir_exists(DRIVER_REPORT_FORTRAN_SRC_DIR), f'Report Fortran source dir "{DRIVER_REPORT_FORTRAN_SRC_DIR}" not found'
    codegen = xsdtools.FortranGenerator(DRIVER_REPORT_SCHEMA_FILE)
    codegen.render_to_files('*', output_dir=DRIVER_REPORT_FORTRAN_SRC_DIR)
