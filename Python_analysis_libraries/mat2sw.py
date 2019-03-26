# import required modules
import numpy as np
import scipy
import pandas as pd
import seaborn as sns
import matplotlib.pylab as plt
import matplotlib.ticker as ticker
import Tkinter
from tkFileDialog import askdirectory
import glob
import mdlab as mdl

cnd = "blocked"
dur = ""
fn = "_".join([cnd, dur]) + ".npz"
# pdb.set_trace()
mdl.mat2spikewaveforms(fn)

print()
print()
print()

# load the data from matlab files using gui
root = Tkinter.Tk()
path_name = askdirectory(initialdir="c:/", title='select where matrix files are located')
root.withdraw()  # close the main root window
kw = "*.npz"

files = glob.glob(path_name + "/" + kw)
path_length = len(path_name) + 1

condtion = []
id = []
spikewaveform = []
for fn in files:
    # add recording names to the header for backup use
    recording_fn = fn[path_length:]
    print(recording_fn)
    data = np.load(fn)
    tmp_id = [data["birdid"][i] + mdl.separator + str(data["unit"][i]) for i in range(len(data["unit"]))]
    tmp_spikewaveform = data["spikewaveform"]
    tmp_condition = [recording_fn] * len(tmp_id)

    condtion.extend(tmp_condition)
    id.extend(tmp_id)
    spikewaveform.extend(tmp_spikewaveform)
