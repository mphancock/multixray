from pathlib import Path
import pickle
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import seaborn as sns
sns.set_theme()
import math
import pandas as pd
import random
import numpy as np
import sys

import IMP
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/src")))
from refine import refine_hs_to_max_ff


if __name__ == "__main__":
    exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/183_test")