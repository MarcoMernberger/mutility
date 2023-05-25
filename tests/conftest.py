# -*- coding: utf-8 -*-
"""
    Dummy conftest.py for mutility.

    If you don't know what this is for, just leave it empty.
    Read more about conftest.py under:
    https://pytest.org/latest/plugins.html
"""
import sys
import pathlib

root = pathlib.Path(__file__).parent.parent
sys.path.append(str(root / "src"))
