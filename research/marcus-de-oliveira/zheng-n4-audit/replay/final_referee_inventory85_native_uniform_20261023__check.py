"""Own replay for final_referee_inventory85_native_uniform_20261023.

This document carries the same proof as root_inventory85_native_uniform_corrected_20261023
(same displays, same numbering); the only mathematical difference is the order of the
second and third entries of C' in (19), which check L5 of the shared replay tests.
This script runs the shared replay (own code, in this directory) unchanged.
Run:  python -I -B <this file> [--no-numeric]
"""
import os
import runpy
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
runpy.run_path(os.path.join(HERE, "root_inventory85_native_uniform_corrected_20261023__check.py"), run_name="__main__")
