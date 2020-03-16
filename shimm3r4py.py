import numpy as np
import matplotlib.pyplot as plt
import csv
import math

from scipy.signal import butter, lfilter  #for signal filtering
import scipy

"""
    returns a nx2 data array foreach ADC channel
"""
def read_data(file_path):    
    with open(file_path) as tsvfile:
        reader = csv.reader(tsvfile, dialect='excel-tab')
        row_idx = 0
        data_list = [[], []]
        for row in reader:
            if(row_idx <= 3):
                row_idx += 1
                continue
            ch_1, ch_2 = float(row[0]), float(row[1])
            data_list[0].append(ch_1)
            data_list[1].append(ch_2)

    array = np.asarray(data_list, dtype=np.float32)
    return np.reshape(array, (len(data_list[0]), -1))

def get_gsr(data, range_setting):
    rf = [40200, 287000, 1000000, 3300000]
    
    rf = lambda b: range_setting/(((b*(3/4095))/0.5) -1)
    
    rows, cols = data.shape
    data = np.apply_along_axis(rf, 1, data)

    return data


def butterworth_lowp(data, cutt_off, srate, order=2):
    nyq = 0.5 * srate 
    normal_cutoff = cut_off / nyq 
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = lfilter(b, a, data)
    return y

def butterworth_highp(data, cutt_off, srate, order=2):
    nyq = 0.5 * srate
    normal_cutoff = cut_off / nyq 
    b, a = butter(order, normal_cutoff, btype='high', analog=False)
    y = lfilter(b, a, data)
    return y

def find_peaks_onset_and_offset(data, onset = 0.01, offset = 0):
    peaks = [] # [onset, max_idx, offset_idx] x #peaks

    is_onset = false
    last_peak = 0

    rows, cols = data.shape
    for point in range(rows):
        x = data[point]
        if is_onset:
            if x <= offset:
                peak_onset = np.argmax(data[last_peak:point])
                if peak_onset >= onset:
                    peaks.append([last_peak, peak_onset+last_peak, point])
                is_onset = false
            else:
                if x > offset:
                    last_peak = point
                    is_onset = True
    return peaks


def phasic_filter(data, srate, seconds=4):
    phasic_signal = []
    rows, cols = data.shape
    for sample in range(rows):
        smin = sample - seconds * srate
        smax = sample + seconds * srate

        if smin < 0:
            smin = sample
        if smax > rows:
            smax = sample

        newsample = data[sample] - np.mean(data[smin:smax])
        phasic_signal.append(newsample)
    return  phasic_signal

def get_gsr_features(filtered_gsr, srate, peaks):
    results = {}
    results["peak"]  = {"start":peak[0], "max": peak[1], "end": peak[2]}
    results["rise_time"] = (peak[1] - peak[0])/srate
    results["latency"] = (peak[0])/srate
    results["amplitude"] = filtered_gsr[peak[1]] - filtered_gsr[peak[0]]
    results["eda_peak"] = filtered_gsr[peak[1]]

    half_amplitude_idx = filtered_gsr.index(
        min(
            filtered_gsr[peak[1]:peak[2]], 
            key = lambda a: abs(a - results["amplitude"]/2.0)
            )
        )
    results["decay_time"] = (half_amplitude_idx - peak[1])/srate
    results["scr_width"] = (half_amplitude_idx - peak[0])/srate
    
    return results

def analyze_gsr(raw_gsr, srate, low_pass = 1.0, high_pass = 0.05, phasic_seconds = 10):

    results = {}

    filtered_gsr = butterworth_lowp(raw_gsr, low_pass, srate, order=2)
    filtered_gsr = butterworth_highp(raw_gsr, low_pass, srate, order=2)

    downsampling_fac = int(srate/10)
    nsamples = int(len(filtered_gsr) / downsampling_fac)
    filtered_gsr = scipy.signal.resample(filtered_gsr, nsamples)
    filtered_gsr = phasic_gsr_Filter(filtered_gsr, srate, seconds=phasic_seconds)

    peaks = find_peaks_onset_and_offset(filtered_gsr)
    
    for peak in peaks:
        results[peaks.index(peak)] = get_gsr_features(filtered_gsr, 10, peak)
    return(results)




def main():
    data = read_data("data/data_device4B78.csv")
    gsr = get_gsr(data, 1)


if __name__ == "__main__":
    main()