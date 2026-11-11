import re
from pathlib import Path

import numpy as np
import pandas as pd

FILENAME_PATTERN = re.compile(
    r"grid_Ra(?P<Ra>[\d.eE+-]+)_Ha(?P<Ha>[\d.eE+-]+)_N(?P<N>[\d.eE+-]+)_Pr(?P<Pr>[\d._eE+-]+)\.csv$"
)


def parse_filename(name):
    match = FILENAME_PATTERN.search(name)
    if match is None:
        return None
    values = match.groupdict()
    values["Pr"] = values["Pr"].replace("_", ".")
    return values


def list_grid_files(directory):
    files = []
    for path in sorted(Path(directory).glob("grid_*.csv")):
        params = parse_filename(path.name)
        if params is not None:
            files.append((path, params))
    return files


def load_fields(path):
    frame = pd.read_csv(path)
    n = frame["i"].nunique()
    shape = (frame["j"].nunique(), n)
    fields = {name: frame[name].to_numpy().reshape(shape) for name in ("i", "j", "T", "C", "u", "v")}
    fields["X"] = fields.pop("i")
    fields["Y"] = fields.pop("j")
    return fields


def stream_function(u):
    return np.cumsum(u, axis=0)
