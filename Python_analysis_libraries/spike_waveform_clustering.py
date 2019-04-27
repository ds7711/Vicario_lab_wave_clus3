import numpy as np
import pandas as pd
import scipy.io
import matplotlib.pylab as plt
from sklearn.decomposition import PCA
from os import walk, path
import sys
sys.path.append('C:\Users\md\Dropbox\MD_scripts\Python_lab_scripts\core_scripts')
sys.path.append('/home/md/Dropbox/MD_scripts/Python_lab_scripts/core_scripts')
sys.path.append('C:\Users\Mingwen\Dropbox\MD_scripts\Python_lab_scripts\core_scripts')
import mdlab as mdl


# constants

def list_files():
    """
    find files in the selected directories and return their absolute path
    :return:
    """
    # from Tkinter import *
    import Tkinter
    from tkFileDialog import askopenfilename, askdirectory

    # load the data from matlab files using gui
    root = Tkinter.Tk()
    path_name = askdirectory(parent=root, initialdir="./")
    root.withdraw()  # close the main root window
    # path = "C:\\Users\\md\\Dropbox\\Lab_Data\\2015_NCM_syllable_surprisal\\Raw_data\\YW570\\"
    fn_list = []
    for (dirpath, dirnames, filenames) in walk(path_name):
        # print dirpath, dirnames, filenames
        filenames = [dirpath + "/" + item for item in filenames]
        fn_list.extend(filenames)
    return (fn_list)


def read_mat_sua_spikes(fn_path):
    """
    read .mat file and extract the unit name and its average waveform
    :param fn_path:
    :return:
        header_list: unit name
        spike_waveform_list: average waveform
    """
    trig_kw = "_trig"
    IDstim_kw = "_IDstim"
    data_kw = "_ass_n"
    spike_waveform_idx = 3
    data_dict = scipy.io.loadmat(fn_path)
    keys = data_dict.keys()

    header_list = []
    spike_waveform_list = []
    for key in keys:
        if trig_kw in key or IDstim_kw in key:
            continue
        elif data_kw in key:
            data_array = data_dict[key][0, 0]
            spike_waveforms = data_array[spike_waveform_idx]
            avg_waveform = spike_waveforms.mean(axis=0)
            header_list.append(key)
            spike_waveform_list.append(avg_waveform)
        else:
            continue
    return (header_list, spike_waveform_list)


def pca_spike_waveform(spike_waveforms, num_components=None):
    """
    use PCA to reduce the dimension of the spikewaveform
    :param spike_waveforms:
    :param num_components:
    :return:
    """
    pca = PCA(n_components=num_components)
    pca.fit(spike_waveforms)
    representations = pca.transform(spike_waveforms)
    print(pca.explained_variance_ratio_)
    return (representations)


def af_clustering(spike_waveforms, preference=None):
    af = AffinityPropagation(preference=preference)
    labels = af.fit_predict(spike_waveforms)
    return(labels)



def extract_unit_identity(unit_name):
    """
    extract bird and unit name from the filename and return a string that contains them
    :param unit_name:
    :return:
    """
    fn_type_kw = "MD"
    unit_kw = "_ass_nw_"
    unit_kw_length = len(unit_kw)
    birdname_ln = 5
    if fn_type_kw in unit_name:
        unit_name = unit_name[2:]
    birdname = unit_name[:birdname_ln]
    st_idx = unit_name.find(unit_kw)
    unit = unit_name[st_idx+unit_kw_length:]
    unit = "".join(unit.split("_"))
    extracted_name = birdname + mdl.separator + unit
    return(extracted_name)

def class_scatters(X, labels):
    """
    scatter plot with color coding
    :param X:
    :param labels:
    :return:
    """
    unique_labels = np.unique(labels)
    n_clusters_ = len(unique_labels)
    colors = list('bgrcmykbgrcmykbgrcmykbgrcmyk')
    for k, col in zip(range(n_clusters_), colors):
        class_members = labels == k
        plt.plot(X[class_members, 0], X[class_members, 1], col + '.')


def code2label(labels, code_label_dict):
    new_labels = [code_label_dict[item] for item in labels]
    return(np.asarray(new_labels))