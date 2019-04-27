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

# cnd = "blocked"
# dur = ""
# fn = "_".join([cnd, dur]) + ".npz"
# # pdb.set_trace()
# mdl.mat2spikewaveforms(fn)

print("middle...")
print("middle...")
print("middle...")

# load the data from matlab files using gui
root = Tkinter.Tk()
path_name = askdirectory(initialdir="c:/", title='select where matrix files are located')
root.withdraw()  # close the main root window
kw = "*.npz"

files = glob.glob(path_name + "/" + kw)
path_length = len(path_name) + 1

condition = []
id = []
spikewaveform = []
nan_label = []
for fn in files:
    # add recording names to the header for backup use
    recording_fn = fn[path_length:]
    print(recording_fn)
    data = np.load(fn)
    tmp_id = [data["birdid"][i] + mdl.separator + str(int(data["unit"][i][0])) for i in range(len(data["unit"]))]
    tmp_spikewaveform = data["spikewaveform"]
    tmp_condition = [recording_fn] * len(tmp_id)

    tmp_nan_label = np.any(np.isnan(tmp_spikewaveform))
    condition.extend(tmp_condition)
    id.extend(tmp_id)
    spikewaveform.extend(tmp_spikewaveform)
    nan_label.append(tmp_nan_label)

df = pd.DataFrame({"condition": condition, "id": id, "spikewaveform": np.asarray(spikewaveform)})

print("end....")
print("end....")
print("end....")

import spike_waveform_clustering as swc

spike_waveform_list = spikewaveform / np.max(spikewaveform, axis=1)[:, np.newaxis]

new_representations = swc.pca_spike_waveform(spikewaveform, num_components=4)