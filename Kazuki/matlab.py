#########################################################################
#                                                                       #
# This routine will take alpha-omega data that has been                 #
# converted to matlab format and convert it to a python                 #
# dictionary for analysis.                                              #
#                                                                       #
# The data format is: data[matlab variable name][number                 #
#      of epochs][data values]                                          #
#                                                                       #
# First a list of all the mat files in a directory are                  #
# brought into the program using the (get_mat_files)                    #
# function                                                              #
#                                                                       #
# Then since the AO system only stores the data in 2                    #
# minute segments all of these files need to be combined.               #
# This is done in the (combine_variables) function.                     #
#                                                                       #
# The next function is the get_variable_names function                  #
# This functions pulls the specific variables to be analyzed            #
# from the Full combined matlab files via a variable list               #
# and is stored in the function vars.                                   #
#                                                                       #
# The function modify_variable_names calls the function                 #
# the function combine_variables. The purpose of this                   #
# function is to create two dictionaries. The first                     #
# dictionary holds the values of the variable names with                #
# the word combine added to the end as the key and the                  #
# values of the variable as the data of that key.                       #
# The second dictionary holds the variable name as the key              #
# and the sample rate as the value for that key.                        #
# Once in that format there are routines that will                      #
# determine the start of the movement of either the                     #
# accelerometer data or the EMG data.                                   #
#                                                                       #
# The start_end function determines the star and end of the movement by #
# determining the start of the movement from either the EMG signal or   #
# the accelerometer signal. The actual variable name is passed to this  #
# function through the vars input variable.                             #
#                                                                       #
#                                                                       #
#                                                                       #
#                                                                       #
#                                                                       #
# Variables:                                                            #
#                                                                       #
#     Directory_path:   The path to the directory where the             #
#           variables are held                                          #
#     file_extension:   The file extension to use - for                 #
#           purposes of this program these will all be                  #
#           matlab files so the extension will be .mat                  #
#     variable_list:    The list of variables to search for.            #
#     vars:             The vector that holds variable names to         #
#           be analyzed from get_variable_names.                        #
#     vars from the start_end function:   This is the actual variable   #
#           used to determine the start and end of the movement         #
#           and is passed to the start_end function                     #
#     combined:         The data dictionary with the key being the name #
#           of the variable and the value of the key being the data     #
#           for that variable.                                          #
#     samp_freq:        The sample frequency dictionary with the key    #
#           being the name of the variable and the value being the      #
#           the sample frequency for that variable.                     #
#     abs_val:          This is the rectified accelerometer of EMG      #
#           in the start_end function used to determine the start of the#
#           movement.                                                   #
#     ts_times:         This is a vector of the start movement times    #
#           from the start_end function.                                #
#     te_times:         This is a vector of the end of the movement     #
#           times from the start end function.                          #



import math
from operator import truediv

import time

import scipy.io as sio
from scipy import signal
import matplotlib.pyplot as plt
import numpy as np
import os
import neurokit2 as nk
from scipy import signal
from scipy.fft import fftshift

from matplotlib.cbook import flatten
from tensorpac import Pac, EventRelatedPac, PreferredPhase
from tensorpac.utils import PeakLockedTF, PSD, ITC, BinAmplitude


# Get a list of all the mat files in the directory path

# directory_path = input(r"Input a directory path:")
directory_path = r'C:\Users\biost\Desktop\Desktop 05212024\Kazuki Data\Gaubatz_Michael\Gaubatz_LBrain_RBody'

# file_extension = input(r"Input a file extension (.xxx): ")
file_extension = '.mat'

#variable_list = input(r'Enter the variable list file (*.txt)')
variable_list = 'v.txt'

# Get all the mat file names so they can be combined
def get_mat_files(directory_path, file_extension):
      print('Get Matlab File Names')
      # maf_files is the list of all the mat files in the directory
      mat_files = []
      for file in os.listdir(directory_path):
            if file.endswith(file_extension):
                  mat_files.append(file)
      return mat_files

# Combined all the variables files together and get sample frequency for each variable
def combine_variables(directory_path, mat_files, variable):
      print('Combining the files to one long data stream')
      print(variable)
      accl = [0]
      sfcv = []
      # print(variable)
      for dp in mat_files:
            dir_path = fr"{directory_path}\{dp}"
            main_mat = sio.loadmat(dir_path)
            hold = main_mat[variable]
            accl.extend(hold.flatten())
      sfcv = (main_mat[variable + '_KHz'] * 1000.0).flatten()
      return accl, sfcv

# Get the specific variable names from a file that will be analyzed
def get_variable_names(variable_list):
      print('Get Variable Names')
      vars=[]

      # Ask the user to enter the file with the variable list in it
      # variable list is a test case
      # Get the list of variables that are going to be analyzed


      # variable_list = input('Enter the file that contains the variable names: ')

      file = open(variable_list, 'r')
      for line in file:
            vars.append(line.strip())
      file.close()
      # vars is a list of the variable names
      return vars


def modify_var_names(vars):
      # Create a dictionary to hold all the combined data and then
      # another dictionary to hold the sample rate data.
      print('Modify Variable Names - creating a dictionary holding each variable name and then the data')
      combined = {}
      samp_freq = {}

      for name in vars:
            #print(name)
            hold, sf = combine_variables(directory_path, mat_files, name)
            combined[name + '_combine'] = hold
            samp_freq[name] = sf
            print(f"samp_freq:{samp_freq[name]}")
      #print(f"samp_freq:{samp_freq}")
      return combined, samp_freq


def start_end(vars, combined, samp_freq, sensitivity):
# Calculate the start and end times for each trial based on the trigger being
# either the accl or EMG start point. Go 2 sec before and 2 sec after so the full
# data will be 2 sec before, 2 second of movement, and 2 sec after.

# Sensitivity is the value that the function looks for to calculate the start point of the trial

# Get the average value of the accl or EMG variable

      # accl, sf = combine_variables(directory_path, mat_files, variable)
#      for i in range(len(combined)):
      l = [0]
      for i in range(len(l)):
            accl = combined[vars + '_combine']

            # Rectify and square of values greater than 250
            print(f"Rectify and Normalizing {vars} data.")

            abs_val = np.absolute(accl)
            abs_val_hold = abs_val
            abs_val_mean = np.mean(abs_val)
            abs_val = abs_val - abs_val_mean

            print(f'Sensitivity {sensitivity}')

            #k = 0
            j = 0
            while j < len(abs_val):
                  if abs_val[j] > sensitivity:
                        abs_val[j] = 100000.0

                        #k = k + 1
                        #print(f'k: {k}')

                        j = j + int((4 * samp_freq[0]))
                        if j >= len(abs_val):
                              break
                        abs_val[j] = 0.0
                  else:
                        abs_val[j] = 0.0
                        j = j + 1

            while j < len(abs_val):
                  if abs_val[j] < 100000.0:
                        abs_val[j] = 0.0
                  j = j + 1




            ###for j in range(len(abs_val)):
            ####      #if abs_val[j] > 1000.0:
            ###      if abs_val[j] > sensitivity:
            ###            abs_val[j] = 10000.0
            ###      else:
            ###            abs_val[j] = 0.0

            # Get the start and end times of each motor movement
            print(f"Get the start and end times of each motor movement from the index variable "
                  f"{vars}")
            ts_times = []
            te_times = []
            flag = 1
            j = 0
            num_trials = 0
            while j < len(abs_val):
                  j = j + 1
                  # Need this break in case the index j gets to big, and we get an index error
                  if j >= len(abs_val):
                        break
                  if abs_val[j] == 100000.0 and flag == 1:
                        ts_times.append(j / samp_freq)
                        # Add 4 times the sample rate to the start point to get the full EMG
                        # signal - note there
                        # are 5 seconds between movements so 4 sample rates is 4 seconds.
                        j = j + int((4 * samp_freq[0]))
                        num_trials += 1
                        te_times.append(j / samp_freq)

            #print(f'num_trials: {num_trials}')
            ts_samples = ts_times * samp_freq
            te_samples = te_times * samp_freq

            # Create a x-axis for the time plot
            times = np.arange(len(abs_val)) * 1 / samp_freq

            print(f'Plot the {vars} data (both absolute and raw)')

            x1 = times

            plt.subplot(311)
            plt.title(f"{vars} data squared off and raw")
            plt.plot(x1, abs_val)
            plt.subplot(312)
            plt.plot(x1, combined[vars + '_combine'])
            plt.subplot(313)
            plt.plot(x1, abs_val_hold)
            plt.show()

            # Plot accl (which is either the accl or EMG data depending on what is being
            # used as the trigger in the time and frequency domains
            #print('Plot')

            j = 0

            hold0e = (te_times[j] * samp_freq) + (samp_freq)
            hold0s = (ts_times[j] * samp_freq)
            hold = len(accl[int(hold0s[0]):int(hold0e[0])])

            yavg = np.zeros(hold)
            while j < num_trials - 1:
                  hold0e = (te_times[j] * samp_freq) + (samp_freq)
                  hold0s = (ts_times[j] * samp_freq)

                  hold0s = int(hold0s[0])
                  hold0e = int(hold0e[0])

                  # Rectify the values
                  yy1 = np.absolute(accl[hold0s:hold0e])

                  # Check to make sure the array size are the
                  # same and if not then adjust the last values of that array
                  if len(yy1) < len(yavg):
                        sub = (abs(len(yy1) - len(yavg)))
                        yhold = np.zeros(sub)
                        yy1 = np.append(yy1, yhold)
                  elif len(yy1) > len(yavg):
                        sub = (abs(len(yy1) - len(yavg))) * -1
                        yy1 = yy1[:sub]
                  yavg = yavg + yy1

                  x1 = times[(hold0s - int(2 * samp_freq[0])):hold0e]
                  #y1 = np.absolute(accl[(hold0s - int(2 * samp_freq[0])):hold0e])

                  ##print(f'Plot the segmented data for trial {j}.')

                  ##plt.subplot(211)
                  ##plt.title(f"accl data for trial: {j}")
                  ##plt.plot(x1, y1)
                  ##plt.show()

                  j = j + 1
      yavg = yavg / num_trials

      # Check to make sure the array size are the
      # same and if not then adjust the last values of that array
      if len(yavg) < len(x1):
            sub = (abs(len(yavg) - len(x1)))
            yhold = np.zeros(sub)
            yavg = np.append(yavg, yhold)
      elif len(yavg) > len(x1):
            sub = (abs(len(yavg) - len(x1))) * -1
            yavg = yavg[:sub]
      print("Average Plot")
      plt.subplot(211)
      plt.title(f"Power Spectral Density (PSD): {vars}")
      plt.plot(x1, yavg)
      plt.subplot(212)
      plt.psd(yavg, NFFT=1024, Fs=samp_freq)
      plt.show()

      print(f'len of ts_samples: {len(ts_samples)}')
      print(f'len of ts_times: {len(ts_times)}')


      # Using neurokit2

      #print(f'Using neurokit2 on {vars}')

      #accl = combined[vars + '_combine']

      #plt.plot(accl)

      #eog_signals, info = nk.eog_process(accl, sampling_rate=samp_freq[0])
      #nk.eog_plot(eog_signals, info)

      #signals, info = nk.emg_process(accl, sampling_rate=samp_freq[0])
      #nk.emg_plot(signals, info)
      # print(f'info {info}')


      #emg_cleaned = nk.emg_clean(accl, sampling_rate=samp_freq[0])
      #emg_amplitude = nk.emg_amplitude(emg_cleaned)
      #activity, info = nk.emg_activation(emg_amplitude=emg_amplitude,
#                                                     method="threshold")

      #nk.events_plot([info["EMG_Offsets"], info["EMG_Onsets"]], emg_cleaned)

      return yavg, hold0s, hold0e, ts_times, te_times, ts_samples, te_samples



def segment_data_variable(variable, accl_start_samples, sf_accl, sf):
      # Segment the specific variable according to the start times of the
      #                  accl or EMG movement
      # Segment the data to the following:
      #     2 sec before the movement
      #     2 sec during the movement
      #     2 sec after the movement

      # First convert the start times from the Accl variable to the
      #               appropriate sample points since sample rates
      #                  are different between variables.
      # vst = variable start times
      #                        (Accl_Start_time)*(Variable_Sample_rate)
      #  Variable_Start_time = -------------------------------------------
      #                             Accl_Sample_rate

      t_samp_var = np.round((accl_start_samples * sf) / sf_accl)
      vst = t_samp_var.flatten()
      before = {}
      after = {}
      during = {}
      all = {}

      for i in range(len(accl_start_samples)):
          during[i] = variable[int(vst[i]):(int(vst[i]) + int(2 * sf[0]) - 1)]
          before[i] = variable[(int(vst[i]) - int(sf[0] * 2)):(int(vst[i]) - 1)]
          after[i] = variable[(int(vst[i]) + int(sf[0] * 2)):(int(vst[i]) +
                                                       (2 * int(sf[0] * 2)) - 1)]
          all[i] = variable[(int(vst[i]) - int(sf[0] * 2)):(int(vst[i]) +
                                                      (2 * int(sf[0] * 2)) - 1)]
      return(all, before, during, after)


def add_motor_condition(y_text, fontsize=14, color='k', ax=None):
    # The values of x_time are the center coordinates of the text on the plot
    x_times = [-1.0, 1.0, 3.0]
    #x_times = [-2.0, 0.0, 2.0]
    x_conditions = ['Pre', 'Movement', 'Post']
    if ax is None: ax = plt.gca()  # noqa
    plt.sca(ax)
    plt.axvline(0., lw=2, color=color)
    plt.axvline(2.0, lw=2, color=color)
    for x_t, t_t in zip(x_times, x_conditions):
        # print(f'x_t {x_t}')
        plt.text(x_t, y_text, t_t, color=color, fontsize=fontsize, ha='center',
                 va='center', fontweight='bold')
    return

def get_avg_seg_variables(data, noe, sf, vars):
      # Calculate the average value for all the variables.

      result = np.zeros(len(data[1]))
      single = np.zeros(len(data[1]))

      #total = data - np.mean(data)
      #result = []
      for i in range(noe - 1):
            # Sum up all the data in the various epochs
            # result = []
            for j in range(len(data[i])):
                  #result.append(total[j] + data[i][j])
                  result[j] = result[j] + data[i][j]
                  single[j] = single[j] + data[5][j]
            #total = result
      # Divide each value in the total by the number of epochs to get the average
      result = np.divide(result, (noe - 1))
      single = np.divide(single, (noe - 1))

      # Find the maximum amplitude of the data for placement of the text
      # pre/movement/post in the graph
      max_amplitude = max(result)
      max_amplitudes = max(single)

      x = []
      for j in range(len(result)):
            x.append((j / sf) - 2.0) # The -2.0 is so time starts at 2 sec before the movement.

      print(f'Segmented data average for: {vars}')

      plt.plot(x, result)
      plt.xlabel("Time(sec)")
      plt.ylabel("uV")
      plt.title(f"Average over segments for: {vars}")
      add_motor_condition(1.0 * max_amplitude)
      plt.show()

      plt.plot(x, single)
      plt.xlabel("Time(sec)")
      plt.ylabel("uV")
      plt.title(f"Single epoch (5) from: {vars}")
      add_motor_condition(1.0 * max_amplitudes)
      plt.show()

      return

def calcualte_PSD(data, sf, noe, nfft, hold2, k, extra_text=''):
      # noe = Number of Epochs
      print(f'Calculating the PSD for variable {hold2}' + extra_text)
      #window_size = int(sf[0])
      window_size = nfft
      # Calculate PSD using Welch's method
      psd = np.zeros(nfft // 2 + 1)
      # Calculate the PSD and sum them
      for i in range(noe - 1): # Sum over the epochs
            n = 0
            print(f'PSD Epoch #: {noe - (i + 1)}')
            #print(f'int(len(data[i])/sf[0]): {int(len(data[i]) / sf[0])}')
            #print(f'sf[0]: {sf[0]}')
            #for j in range(int(len(data[i]) / sf[0])): # Sum over the data in the epochs
            for j in range(int(len(data[i]) / nfft)):  # Sum over the data in the epochs
                  print(f'PSD subepoch #: {j + 1}')
                  # (i.e. if 6 sec then would do six FFTs)
                  windowed_data = data[i][n:(n + window_size)] * np.hanning(window_size)
                  n = n + window_size
                  fft_result = np.fft.fft(windowed_data, n=nfft)
                  psd += np.abs(fft_result[:nfft // 2 + 1]) ** 2  # Take the positive frequencies

      psd /= len(data[k] * noe)  # Normalize by the total length of the signal
      frequency = np.fft.fftfreq(nfft, 1 / sf)[:nfft // 2 + 1]

      # Plot the PSD
      plt.figure(figsize=(10, 6))
      plt.plot(frequency, psd)
      plt.xlabel("Frequency (Hz)")
      plt.ylabel("Power Spectral Density")
      plt.title(f"Power Spectral Density (PSD): {hold2} - " + extra_text)
      plt.yscale('log')
      plt.xlim(0, 100)
      plt.grid(True)
      plt.show()

      #f, Pxx_den = signal.periodogram(data[k], sf[0], window='hamming', nfft=256,
      #                                scaling='spectrum')
      #plt.semilogy(f, Pxx_den)
      #plt.ylim([1e-7, 1e2])
      #plt.xlabel('frequency [Hz]')
      #plt.ylabel('PSD [V**2/Hz]')
      #plt.xlim(0, 100)
      #plt.yscale('log')
      #plt.title(f"Periodogram: {hold2} - " + extra_text)
      #plt.show()

      return

# Calculate the spectrogram
def calcualte_spectrogram(data, sf, hold2, k, extra_text=''):
      print(f'Calculating the spectrogram for variable {hold2}' + extra_text)
      f, t, Sxx = signal.spectrogram(data[k], sf[0])
      plt.pcolormesh(t, f, Sxx, shading='gouraud')
      plt.ylabel('Frequency [Hz]')
      plt.xlabel('Time [sec]')
      plt.title(f"Spectrogram for: {hold2} - " + extra_text)
      plt.ylim(0, 100)
      plt.show()

      return


def line_filter(variable, sf, line):

      # Remove 60 Hz and the DC offset from the signals.

      # First subtract out the DC component

      variable = variable - np.mean(variable)

      f0 = line
      Q = 20.0
      b, a = signal.iirnotch(f0, Q, sf)
      filtered_data = signal.filtfilt(b, a, variable)
      return filtered_data


def citaf(data, noe, hold2):
      # Combine the data into the format that tensorpac needs
      # Format is an array that is (noe, data)
      # Data is the data in dictionary format for the analysis
            # data[number of epochs][values]
      # tp_array is the array in tensorpac format

      #print(f'keys of data: {data.keys()}')

      print(f'Converting data to a format tensorpac needs for variable: {hold2}')

      tp_array = np.zeros((noe, len(data[0])))

      #print(f'tp_array shape: {tp_array.shape}')
      #print(f'data shape: {data[29].shape}')
      #print(f'noe: {noe}')

      for i in range(noe):
            #print(f'i:{i}')
            for j in range(len(data[i])):
                  tp_array[i, j] = data[i][j]

      return tp_array

def calcualte_bispectra(x, y, z, sf, fft_size, noe, hold2, k):
      # noe = Number of Epochs
      # x, y, and z are the three data variables - they can all be the same or
      #   different depending on if you want to do a cross spectra

      tbstart = time.time()

      global x_fft_result, y_fft_result, z_fft_result
      print(f'Calculating the bispectrum/bi-phase/bicoherance for variable {hold2}')
      #window_size = int(sf[0])
      window_size = fft_size
      nfft = window_size
      sample_freq = sf[0]

      # Test Data
      # First signal consisting of a one second 10 <-> 100hz coupling
      # sf = 1000.
      # x1, tvec = pac_signals_wavelet(f_pha=10, f_amp=100, n_epochs=n_epochs, noise=2,
      #                             n_times=n_times, sf=sf)
      # Second signal : one second of random noise
      # x2 = np.random.rand(n_epochs, 1000)
      # now, concatenate the two signals across the time axis
      # x = np.concatenate((x1, x2), axis=1)
      # time = np.arange(x.shape[1]) / sf
      #################################################
      #           Test                                #
      #################################################
      ##sample_freq = 1000
      ##window_size = 1000
      ##nfft = 1000
      ##noe = 10
      ##duration = 10
      ##frequency1 = 11
      ##frequency2 = 7
      ##t = np.linspace(0, duration, int(sample_freq * 1), endpoint=False)
      ##sine_wave1 = np.sin(2 * np.pi * frequency1 * t)
      ##sine_wave2 = np.sin(2 * np.pi * frequency2 * t)
      ##x1 = (sine_wave1 * sine_wave2)
      ##plt.figure(figsize=(10, 6))
      ##plt.plot(t, x1)
      ##plt.xlabel("Time (s)")
      ##plt.ylabel("Amplitude")
      ##plt.title("Multiplied Sine Waves")
      ##plt.xlim(0, 1)
      ##plt.grid(True)
      ##plt.show()


      bisp = np.zeros((nfft, nfft), dtype=complex)
      bisp_mag = np.zeros((nfft, nfft))
      bisp_phase = np.zeros((nfft, nfft))
      sum1 = np.zeros((nfft, nfft))
      sum2 = np.zeros((nfft, nfft))
      sum3 = np.zeros((nfft, nfft))



      # Calculate the FFTs
      for i in range(noe - 1): # Sum over the epochs
            print(f'Epoch: {noe - 1 - i}')
            n = 0
            # print(f'int(len(x[i]) / samp_freq {int(len(x[i]) / sample_freq)}')
            #for j in range(int(len(x[i]) / sample_freq)): # Sum over the data in the epochs (sub-epochs)
            for j in range(int(len(x[i]) / nfft)):  # Sum over the data in the epochs (sub-epochs)
                  # (i.e. if 6 sec then would do six FFTs)
                  # print(f'j: {j}')
                  print(f'SubEpoch: {(int(len(x[i]) / nfft)) - j}')

                  # Window the data using a Hanning window
                  x_windowed_data = x[i][n:(n + window_size)] * np.hanning(window_size)
                  y_windowed_data = y[i][n:(n + window_size)] * np.hanning(window_size)
                  z_windowed_data = z[i][n:(n + window_size)] * np.hanning(window_size)
                  n = n + window_size

                  # Calculate the FFT for the triple product
                  x_fft_result = np.fft.fft(x_windowed_data, n=nfft)
                  y_fft_result = np.fft.fft(y_windowed_data, n=nfft)
                  z_fft_result = np.fft.fft(z_windowed_data, n=nfft)

                  #for m in range(int((sample_freq / 2) - 1)):
                  #      for n in range(int((sample_freq / 2) - 1)):

                  for m in range(int((nfft / 2) - 1)):
                        for n in range(int((nfft / 2) - 1)):
                              # Calculate the triple product at all frequencies to 1/2 the sample freq.
                              bisp[m, n] = bisp[m, n] + ((x_fft_result[m] *
                                          y_fft_result[n]) *
                                          np.conj(z_fft_result[m + n]))
                              # Calculate the denominator for the bicoherance
                              sum1[m, n] = sum1[m, n] + ((np.absolute(x_fft_result[m]*
                                                                      y_fft_result[n]))**2)
                              sum2[m, n] = sum2[m, n] + ((np.absolute(z_fft_result[m + n]))**2)


      # Calculate the magnitude and phase of the bispectrum

      bisp_mag = np.absolute(bisp) / (noe + (int(len(x[1]) / sample_freq)))
      bisp_phase = np.angle(bisp) / (noe + (int(len(x[1]) / sample_freq)))

      #print(f'bisp_mag[100, 100]: {bisp_mag[100, 100]}')

      # Calculate the bicoherence
      sum3 = ((sum1 / (noe + (int(len(x[1]) / sample_freq)))) *
              (sum2 / (noe + (int(len(x[1]) / sample_freq)))))
      bic = bisp_mag / np.sqrt(sum3)

      #print(f'sum3[100, 100]: {sum3[100, 100]}')

      # significance value is from table 3.5 in JLS thesis
      sig = 2.131 / np.sqrt((noe + (int(len(x[1]) / sample_freq))))

      for m in range(int((sample_freq / 2) - 1)):
            for n in range(int((sample_freq / 2) - 1)):
                  if bic[m, n] < sig:
                        bic[m, n] = 0.0

      # Get a frequency list for the x-axis of the plots
      frequency = np.fft.fftfreq(nfft, d=(1 / sample_freq))
      plt.figure(figsize=(10, 6))
      plt.xlabel("Frequency 1 (Hz)")
      plt.ylabel("Frequency 2 (Hz)")
      plt.title(f"Bispectrum of: {hold2}")
      max_value = bisp_mag.max()
      #plt.contour(frequency, frequency, bisp_mag, levels=np.linspace(0, 1e9, 50),
      #            cmap='viridis', extend='both')
      plt.contour(frequency, frequency, bisp_mag, levels=np.linspace(0, max_value, 50),
                  cmap='viridis', extend='both')
      plt.colorbar()
      plt.ylim(0, 80)
      plt.xlim(0, 80)
      # plt.grid(True)
      plt.show()

      # Plot the biphase
      plt.xlabel("Frequency 1 (Hz)")
      plt.ylabel("Frequency 2 (Hz)")
      plt.title(f"Bisphase of: {hold2}")

      plt.contourf(frequency, frequency, bisp_phase, cmap='viridis')
      plt.colorbar()
      plt.ylim(0, 80)
      plt.xlim(0, 80)
      # plt.grid(True)
      plt.show()

      # Plot the bicoherance
      plt.figure(figsize=(10, 6))
      # plt.plot(frequency, psd)
      plt.xlabel("Frequency 1 (Hz)")
      plt.ylabel("Frequency 2 (Hz)")
      plt.title(f"Bcoherance of: {hold2}")


      plt.contourf(frequency, frequency, bic, levels=np.linspace(0, 1, 25), cmap='viridis')
      plt.colorbar()
      plt.ylim(0, 80)
      plt.xlim(0, 80)
      #plt.grid(True)
      plt.show()

      print(f'Bispectral running time in seconds = {time.time() - tbstart}')

      return bisp_mag, bisp_phase, bic

def erpac_sub (data, sf, time_x, v_text):

      #################################################
      #                                               #
      # The code below are modified from examples     #
      # taken from:                                   #
      # https://etiennecmb.github.io/tensorpac        #
      # /auto_examples/tuto/plot_real_data.html       #
      #                                               #
      #################################################

      print(f'Calculating ERPAC for variable: {v_text}')

      # define an ERPAC object
      p = EventRelatedPac(f_pha=[8, 30], f_amp=(50, 200, 1, 1))
      # p = EventRelatedPac(f_pha=[15, 30], f_amp='hres')

      # extract phases and amplitudes
      pha = p.filter(sf, data, ftype='phase', n_jobs=1)
      amp = p.filter(sf, data, ftype='amplitude', n_jobs=1)

      ###############################################################################
      # Compute the ERPAC using the two implemented methods and plot it
      ###############################################################################

      print(f'Plotting ERPAC for variable: {v_text}')

      # implemented ERPAC methods
      methods = ['circular', 'gc']

      plt.figure(figsize=(14, 8))
      for n_m, m in enumerate(methods):
            # compute the erpac
            erpac = p.fit(pha, amp, method=m, smooth=100, n_jobs=-1).squeeze()

            # plot
            plt.subplot(len(methods), 1, n_m + 1)
            p.pacplot(erpac, time_x, p.yvec, xlabel='Time (second)' * n_m,
                      cmap='Spectral_r', ylabel='Amplitude frequency', title=p.method + v_text,
                      cblabel='ERPAC', vmin=0., rmaxis=True, polar=False)
            plt.axvline(0.0, linestyle='--', color='r', linewidth=2)
            plt.axvline(2.0, linestyle='--', color='r', linewidth=2)

      plt.tight_layout()
      p.show()

      print(f'Plotting the Inter-trial Coherence for: {v_text}')

      itc = ITC(data, sf, f_pha=(8, 30, 2, 0.2))
      itc.plot(times=time_x, cmap='plasma', fz_labels=15, fz_title=10)
      add_motor_condition(25, color='white')
      plt.suptitle(f'{v_text}')
      plt.show()

      print(f'Plotting ERPPAC for alpha-phase for variable: {v_text}')

      rp_obj = EventRelatedPac(f_pha=[8, 12], f_amp=(30, 160, 30, 2))
      erpac = rp_obj.filterfit(sf, data, method='gc', smooth=100)
      rp_obj.pacplot(erpac.squeeze(), time_x, rp_obj.yvec, xlabel='Time (seconds)',
                     ylabel='Amplitude frequency (Hz)',
                     title='Event-Related PAC occurring for alpha phase for' + v_text,
                     fz_labels=15, fz_title=8)
      add_motor_condition(135, color='white')
      plt.show()

      print(f'Plotting ERPPAC for beta-phase for variable: {v_text}')

      rp_obj = EventRelatedPac(f_pha=[13, 30], f_amp=(30, 160, 30, 2))
      erpac = rp_obj.filterfit(sf, data, method='gc', smooth=100)
      rp_obj.pacplot(erpac.squeeze(), time_x, rp_obj.yvec, xlabel='Time (seconds)',
                     ylabel='Amplitude frequency (Hz)',
                     title='Event-Related PAC occurring for beta phase for' + v_text,
                     fz_labels=15, fz_title=8)
      add_motor_condition(135, color='white')
      plt.show()


      print(f'Plotting the PAC for the three times in movement for: {v_text}')

      p_obj = Pac(idpac=(6, 0, 0), f_pha=(8, 30, 4, .2), f_amp=(60, 200, 20, 2))
      # extract all the phases and amplitudes
      pha_p = p_obj.filter(sf, data, ftype='phase')
      amp_p = p_obj.filter(sf, data, ftype='amplitude')
      # define time indices where rest, planning and execution are defined
      time_rest = slice(0, int(sf * 2))
      time_prep = slice(int(sf * 2), int(sf * 4))
      time_exec = slice(int(sf * 4), int(sf * 6))
      # define phase / amplitude during rest / planning / execution
      pha_rest, amp_rest = pha_p[..., time_rest], amp_p[..., time_rest]
      pha_prep, amp_prep = pha_p[..., time_prep], amp_p[..., time_prep]
      pha_exec, amp_exec = pha_p[..., time_exec], amp_p[..., time_exec]
      # compute PAC inside rest, planning, and execution
      pac_rest = p_obj.fit(pha_rest, amp_rest).mean(-1)
      pac_prep = p_obj.fit(pha_prep, amp_prep).mean(-1)
      pac_exec = p_obj.fit(pha_exec, amp_exec).mean(-1)

      vmax = np.max([pac_rest.max(), pac_prep.max(), pac_exec.max()])
      kw = dict(vmax=vmax, vmin=.04, cmap='viridis')
      plt.figure(figsize=(14, 4))
      plt.suptitle(f'PAC for the three movement times for: {v_text}')
      plt.subplot(131)
      p_obj.comodulogram(pac_rest, title="PAC Pre-Movement [-2, 0]s", **kw)
      plt.subplot(132)
      p_obj.comodulogram(pac_prep, title="PAC Movement [0, 2]s", **kw)
      plt.ylabel('')
      plt.subplot(133)
      p_obj.comodulogram(pac_exec, title="PAC Post Movement [2, 4]s", **kw)
      plt.ylabel('')
      plt.tight_layout()
      plt.suptitle(v_text)
      plt.show()

      #########################################################################
      #                                                                       #
      #                                                                       #
      # To align time-frequency maps in Python, techniques often involve      #
      # adjusting the time or frequency axes to match a reference or to       #
      # compensate for variations. Here's a breakdown of common methods:      #
      # Peak-locking                                                          #
      # If the time-frequency maps are related to events with specific        #
      # peaks, aligning them to the peak time can be effective.               #
      #                                                                       #
      #                                                                       #
      #########################################################################


      print(f'Align time-frequency map for (alpha phase): {v_text}')

      peak = PeakLockedTF(data, sf, 0., times=time_x, f_pha=[8, 12],
                          f_amp=(5, 160, 30, 2))
      plt.figure(figsize=(8, 8))
      plt.suptitle(f'Align time-frequency map for (alpha phase): {v_text}')
      ax_1, ax_2 = peak.plot(zscore=True, baseline=(250, 750), cmap='Spectral_r',
                             vmin=-1, vmax=2)
      add_motor_condition(135, color='black', ax=ax_1)
      plt.tight_layout()
      plt.show()

      print(f'Align time-frequency map for (beta phase): {v_text}')

      peak = PeakLockedTF(data, sf, 0., times=time_x, f_pha=[13, 30],
                          f_amp=(5, 160, 30, 2))
      plt.figure(figsize=(8, 8))
      plt.suptitle(f'Align time-frequency map for (beta phase): {v_text}')
      ax_1, ax_2 = peak.plot(zscore=True, baseline=(250, 750), cmap='Spectral_r',
                             vmin=-1, vmax=2)
      add_motor_condition(135, color='black', ax=ax_1)
      plt.tight_layout()
      plt.show()


      print(f'Binning gamma amplitude according to alpha for: {v_text}')

      # define phase and amplitude filtering properties
      kw_filt = dict(f_pha=[8, 12], f_amp=[75, 105], n_bins=20)
      # bin the rest, planning and execution periods. Note that ideally, the entire
      # trial should be filtered and then binning should be performed
      bin_rest = BinAmplitude(data[:, time_rest], sf, **kw_filt)
      bin_prep = BinAmplitude(data[:, time_prep], sf, **kw_filt)
      bin_exec = BinAmplitude(data[:, time_exec], sf, **kw_filt)

      plt.figure(figsize=(16, 5))
      plt.suptitle(f'Binning gamma amplitude according to alpha for: {v_text}')
      # bin rest period
      plt.subplot(1, 3, 1)
      bin_rest.plot(normalize=True, color='gray', unit='deg')
      plt.ylim(0.5, 1.2), plt.title("Pre", fontsize=18)
      # bin planning period
      plt.subplot(1, 3, 2)
      bin_prep.plot(normalize=True, unit='deg')
      plt.ylim(0.5, 1.2), plt.ylabel(''), plt.title("Movement", fontsize=18)
      # bin execution period
      plt.subplot(1, 3, 3)
      bin_exec.plot(normalize=True, color='red', unit='deg')
      plt.ylim(0.5, 1.2), plt.ylabel(''), plt.title("Post", fontsize=18)
      plt.tight_layout()
      plt.show()

      print(f'Binning gamma amplitude according to beta for: {v_text}')

      # define phase and amplitude filtering properties
      kw_filt = dict(f_pha=[13, 30], f_amp=[75, 105], n_bins=20)
      # bin the rest, planning and execution periods. Note that ideally, the entire
      # trial should be filtered and then binning should be performed
      bin_rest = BinAmplitude(data[:, time_rest], sf, **kw_filt)
      bin_prep = BinAmplitude(data[:, time_prep], sf, **kw_filt)
      bin_exec = BinAmplitude(data[:, time_exec], sf, **kw_filt)

      plt.figure(figsize=(16, 5))
      plt.suptitle(f'Binning gamma amplitude according to beta for: {v_text}')
      # bin rest period
      plt.subplot(1, 3, 1)
      bin_rest.plot(normalize=True, color='gray', unit='deg')
      plt.ylim(0.5, 1.2), plt.title("Pre", fontsize=18)
      # bin planning period
      plt.subplot(1, 3, 2)
      bin_prep.plot(normalize=True, unit='deg')
      plt.ylim(0.5, 1.2), plt.ylabel(''), plt.title("Movement", fontsize=18)
      # bin execution period
      plt.subplot(1, 3, 3)
      bin_exec.plot(normalize=True, color='red', unit='deg')
      plt.ylim(0.5, 1.2), plt.ylabel(''), plt.title("Post", fontsize=18)
      plt.tight_layout()
      plt.show()


      print(f'Identifying the preferred alpha phase for: {v_text}')
      # define the preferred phase object
      pp_obj = PreferredPhase(f_pha=[8, 12])
      # only extract the alpha phase
      pp_pha = pp_obj.filter(sf, data, ftype='phase')
      pp_pha_rest = pp_pha[..., time_rest]
      pp_pha_prep = pp_pha[..., time_prep]
      pp_pha_exec = pp_pha[..., time_exec]
      # compute the preferred phase (reuse the amplitude computed above)
      ampbin_rest, _, vecbin = pp_obj.fit(pp_pha_rest, amp_rest, n_bins=72)
      ampbin_prep, _, vecbin = pp_obj.fit(pp_pha_prep, amp_prep, n_bins=72)
      ampbin_exec, _, vecbin = pp_obj.fit(pp_pha_exec, amp_exec, n_bins=72)
      # mean binned amplitude across trials
      ampbin_rest = np.squeeze(ampbin_rest).mean(-1).T
      ampbin_prep = np.squeeze(ampbin_prep).mean(-1).T
      ampbin_exec = np.squeeze(ampbin_exec).mean(-1).T

      plt.figure(figsize=(18, 5.2))
      plt.suptitle(f'Preferred Alpha phase for: {v_text}')
      kw_plt = dict(cmap='Spectral_r', interp=.1, cblabel='Amplitude bins',
                    vmin=0.012, vmax=0.016, colorbar=True, y=1.05, fz_title=18)
      pp_obj.polar(ampbin_rest, vecbin, p_obj.yvec, subplot=131, title='Pre',
                   **kw_plt)
      pp_obj.polar(ampbin_prep, vecbin, p_obj.yvec, subplot=132,
                   title='Movement', **kw_plt)
      pp_obj.polar(ampbin_exec, vecbin, p_obj.yvec, subplot=133,
                   title='Post', **kw_plt)
      plt.tight_layout()
      plt.show()

      print(f'Identifying the preferred beta phase for: {v_text}')
      # define the preferred phase object
      pp_obj = PreferredPhase(f_pha=[13, 30])
      # only extract the alpha phase
      pp_pha = pp_obj.filter(sf, data, ftype='phase')
      pp_pha_rest = pp_pha[..., time_rest]
      pp_pha_prep = pp_pha[..., time_prep]
      pp_pha_exec = pp_pha[..., time_exec]
      # compute the preferred phase (reuse the amplitude computed above)
      ampbin_rest, _, vecbin = pp_obj.fit(pp_pha_rest, amp_rest, n_bins=72)
      ampbin_prep, _, vecbin = pp_obj.fit(pp_pha_prep, amp_prep, n_bins=72)
      ampbin_exec, _, vecbin = pp_obj.fit(pp_pha_exec, amp_exec, n_bins=72)
      # mean binned amplitude across trials
      ampbin_rest = np.squeeze(ampbin_rest).mean(-1).T
      ampbin_prep = np.squeeze(ampbin_prep).mean(-1).T
      ampbin_exec = np.squeeze(ampbin_exec).mean(-1).T

      plt.figure(figsize=(18, 5.2))
      plt.suptitle(f'Preferred Beta phase for: {v_text}')
      kw_plt = dict(cmap='Spectral_r', interp=.1, cblabel='Amplitude bins',
                    vmin=0.012, vmax=0.016, colorbar=True, y=1.05, fz_title=18)
      pp_obj.polar(ampbin_rest, vecbin, p_obj.yvec, subplot=131, title='Pre',
                   **kw_plt)
      pp_obj.polar(ampbin_prep, vecbin, p_obj.yvec, subplot=132,
                   title='Movement', **kw_plt)
      pp_obj.polar(ampbin_exec, vecbin, p_obj.yvec, subplot=133,
                   title='Post', **kw_plt)
      plt.tight_layout()
      plt.show()

if __name__ == '__main__':

      # If flag_for_bispectra_plots is 0 then do not calculate the bispectral plots
      # If flag_for_bispectral_plots is 1 then calculate and print the bispectral plots
      flag_for_bispectral_plots = 0

      # If flag_for_spectral_plots is 0 then do not calculate and print the PSD plots
      # If flag_for_spectral_plots is 1 then calculate and print the PSD plots
      flag_for_spectral_plots = 0

      # If flag_for_spectrogram_plots is 0 then do not calculate and print the spectrogram plots
      # If flag_for_spectrogram_plots is 1 then calculate and print the spectrogram plots
      flag_for_spectrogram_plots = 0

      # If flag_for_ERPAC is 0 then do not print the ERPAC
      # If flag_for_ERPAC is 1 then print the ERPAC
      flag_for_ERPAC_plots = 0

      # The sensitivity variable is used to select the start points for the
      # accl or EMG triggered data amplitudes.
      sensitivity = 250.0

      # This tells the code that we are running data from Kazuki's motor test
      # and it should only have 30 epochs - no more.
      # flag = 1 Kazuki
      # flag = 0 Other.
      flag2 = '1'

      # This is the fft epoch length for the PSD and the bispectra
      nfft = 1024

      # directory_path and file extension are entered by the user
      mat_files = get_mat_files(directory_path, file_extension)
      # variable_list is a list of the variables in the mat file that are to be analyzed
      vars = get_variable_names(variable_list)
      combined, samp_freq = modify_var_names(vars)

      # Add a 60 Hz filter to the data

      for k in range(len(vars)):
            print(f'Removing line noise from: {vars[k]}')
            sf = samp_freq[vars[k]]
            variable = combined[vars[k] + '_combine']
            filtered_data = line_filter(variable, sf[0], 60)

            combined[vars[k] + '_combine'] = filtered_data

      # Get the starting and ending times for the movements based on the accl
      # data
      for j in range(len(vars)):
            #if vars[j] == 'CANALOG_IN_1___Accl':
            #if vars[j] == 'CEMG_2___02___Flx_Ext':
            if vars[j] == 'CEMG_2___01___Bi_Tri':
                  # Store the variable sample rate that the segments are built on
                  sfreq = samp_freq[vars[j]]

                  # sensitivity is the variable the determines the cutoff for the detection point
                  # of the above variable to calculate the start of each trial
                  # Sensitivity variable is used here
                  print(f'Calculating the start and end times for each trial based on the {vars[j]} variable.')
                  yavg, hold0s, hold0e, ts_times, te_times, ts_samples, te_samples = \
                        (start_end(vars[j], combined, sfreq, sensitivity))

      # Build a dictionary of the data for the segments before, during
      # after and for all the data
      before_seg = {}
      during_seg = {}
      after_seg = {}
      all_seg = {}
      for k in range(len(vars)):
            #print(f'k: {k}')
            hold = vars[k] + '_combine'
            hold2 = vars[k] + '_segmented'
            # Just for testing purposes - for plot to make data is correct
            # hold3 is the variable that is being segmented data label
            # k=2 is the accelerometer variable
            # k=0 is the bi_tri variable
            # k=1 is the flex_ext variable
            if k == 0:
                  hold3 = vars[k] + '_segmented'

            all, before, during, after = segment_data_variable(combined[hold],
                                    ts_samples, sfreq, samp_freq[vars[k]])

            # Fill the dictionaries with the segmented data
            before_seg[hold2] = before
            during_seg[hold2] = during
            after_seg[hold2] = after
            all_seg[hold2] = all

      print("Plot the complete segmented data for one segment")

      times = np.arange(len(all_seg[hold3][1]))

      x1 = all_seg[hold3][1]
      plt.plot(times, x1)
      plt.title(f"Segmented data for first trial of: {hold3}")
      plt.show()

      for k in range(len(vars)):

            # combined variables are the long dictionaries that hold all the data from one
            # variable across all mat files
            # segments variables are the variables that have been chopped at the start
            # of the triggered movement variable and then two seconds in front and 2 seconds
            # behind

            if (vars[k][0:7] == 'CANALOG') or (vars[k][0:4] == 'CEMG'):
                  print(f'Skipping {vars[k]}')
            else:
                  hold = vars[k] + '_combine'
                  hold2 = vars[k] + '_segmented'
                  x = after_seg[hold2]
                  y = after_seg[hold2]
                  z = after_seg[hold2]

                  if flag_for_spectral_plots == 1:
                        t_spect_start = time.time()
                        print('Before FFT Plots')
                        calcualte_PSD(before_seg[hold2], samp_freq[vars[k]],
                                      len(ts_samples), nfft, hold2, k, 'Before')
                        print('During FFT Plots')
                        calcualte_PSD(during_seg[hold2], samp_freq[vars[k]],
                                      len(ts_samples), nfft, hold2, k, 'During')
                        print('After FFT Plots')
                        calcualte_PSD(after_seg[hold2], samp_freq[vars[k]],
                                      len(ts_samples), nfft, hold2, k, 'After')
                        print(f'Spectral running time in seconds = {time.time() - t_spect_start}')
                  else:
                        print('Skipping PSD plots')

                  if flag_for_bispectral_plots == 1:
                        t_bis_start = time.time()
                        mag, phase, bic = calcualte_bispectra(x, y, z,samp_freq[vars[k]],
                                                              nfft, len(ts_samples), hold2, k)

                        print(f'Bispectral running time in seconds = {time.time() - t_bis_start}')

                  else:
                        print('Skipping bispectral plots')

      if flag_for_spectrogram_plots == 1:

            for k in range(len(vars)):
                  if (vars[k][0:7] == 'CANALOG') or (vars[k][0:4] == 'CEMG'):
                        print(f'Skipping {vars[k]}')
                  else:
                        hold = vars[k] + '_combine'
                        hold2 = vars[k] + '_segmented'
                  calcualte_spectrogram(all_seg[hold2], samp_freq[vars[k]], hold2, k, extra_text='All')
      else:
            print(f'Skipping spectrogram plots')


      for k in range(len(vars)):
            hold2 = vars[k] + '_segmented'
            # Number of epochs
            noe = len(ts_samples)
            get_avg_seg_variables(all_seg[hold2], noe,
                            samp_freq[vars[k]], hold2)

      # Place the data in the format the tensorpac needs:
      #     Array of signals of shape (n_epochs, n_times).

      # Calculate the ERPAC

      if flag_for_ERPAC_plots == 1:
            for k in range(len(vars)):

                  # Skip the trigger variables (accelerometer and emg)
                  if (vars[k][0:7] == 'CANALOG') or (vars[k][0:4] == 'CEMG'):
                        print(f'Skipping {vars[k]}')
                  else:
                        t_epac_start = time.time()
                        hold2 = vars[k] + '_segmented'
                        noe = len(ts_samples)

                        # Convert to tensorpac format
                        if flag2 == 1:
                              if noe > 30:
                                    noe = 30
                        tp_data_array = citaf(all_seg[hold2], (noe), hold2)

                        # Get the sample frequency and the time values
                        #sfq = sf[0] * 1.0
                        sfq = samp_freq[vars[k]][0] * 1.0

                        # Get the time vector - note the -2 is so we get 2 sec before, 2 sec
                        # during and 2 sec after, and they line up at zero as the start trigger.
                        time_1 = (np.arange(tp_data_array.shape[1]) / sfq) - 2.0

                        # Calculate the Event Related PAC for variable vars[k]
                        erpac_sub(tp_data_array, sfq, time_1, vars[k])

                        print(f'PAC running time in seconds = {time.time() - t_epac_start}')
      else:
            print('Skipping the ERPAC Plots')
