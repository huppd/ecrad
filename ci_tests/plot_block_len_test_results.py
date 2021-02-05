#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Mikhail Zhigun"

import argparse
import pickle
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from tests_runner import BlockLengthTestRes, BlockLengthTestsRes
from collections import OrderedDict as odict


def main():
    parser = argparse.ArgumentParser(description="ECRAD block tests results plotter")
    parser.add_argument("--block-len-test-file-path", type=str, default='block-len-test-results.bin', required=True)
    args = parser.parse_args()
    f_path = args.block_len_test_file_path
    with open(f_path, 'rb') as f:
        all_res = pickle.load(f)
    sizes = [sz for sz in range(all_res.min_block_size, all_res.max_block_size + 1)]
    n_sizes = len(sizes)

    def size_idx(sz: int, min_size: int = all_res.min_block_size):
        return sz - min_size

    timings_by_type = odict()
    for solver_type in all_res.solver_types:
        timings_by_type[solver_type] = [None] * n_sizes
    for res in all_res.results:
        if not res.omp_enabled:
            ave_time = sum(res.timings) / len(res.timings)
            timings_by_type[res.solver_name][size_idx(res.block_size)] = float(ave_time) / 1000.

    fig, ax = plt.subplots()  # Create a figure and an axes.
    for solver_type, timings in timings_by_type.items():
        ax.plot(sizes, timings, label=solver_type)
    ax.set_xlabel('# columns per block')
    ax.set_ylabel('runtime, ms')
    ax.set_title(f'{all_res.cpu_model_name}, single thread')
    ax.legend()
    ax.set_ylim(ymin=0)
    fig.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.grid()

    plt.show()
    return True


if __name__ == '__main__':
    sys.exit(0 if main() else 1)
