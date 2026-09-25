"""Include the A1 coordinate contract in lightweight repository CI."""
import importlib.util
from pathlib import Path

path = (Path(__file__).resolve().parents[2] /
        'RQ_Specified/A1_transitional_epithelial_state_distinction/scripts/tests/test_direct_intervals.py')
spec = importlib.util.spec_from_file_location('a1_interval_tests', path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
IntervalTests = module.IntervalTests
