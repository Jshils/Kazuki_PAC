from operator import truediv

import scipy.io as sio
from scipy import signal
import matplotlib.pyplot as plt
import numpy as np
import os

from matplotlib.cbook import flatten
from tensorpac import Pac, EventRelatedPac, PreferredPhase
from tensorpac.utils import PeakLockedTF, PSD, ITC, BinAmplitude

# Get a list of all the mat files in the directory path

# directory_path = input(r"Input a directory path:")
directory_path = r'C:\Users\biost\Desktop\Desktop 05212024\Kazuki Data\Gaubatz_Michael\Gaubatz_LBrain_RBody'

# Get the list of variables that are going to be analyzed
variable_list = 'v.txt'

# file_extension = input(r"Input a file extension (.xxx): ")
file_extension = '.mat'

# Get all the mat file names so they can be combined
def get_mat_files(directory_path, file_extension):
      print('Get Matlab File Names')
      # maf_files is the list of all the mat files in the directory
      mat_files = []
      for file in os.listdir(directory_path):
            if file.endswith(file_extension):
                  mat_files.append(file)
      return mat_files

# Combined all the variables files together and get sample frequency
def combine_variables(directory_path, mat_files, variable):
      print('Combining the files to one long data stream')
      print(variable)
      accl = [0]
      sf = []
      # print(variable)
      for dp in mat_files:
            dir_path = fr"{directory_path}\{dp}"
            main_mat = sio.loadmat(dir_path)
            hold = main_mat[variable]
            accl.extend(hold.flatten())
      sf = (main_mat[variable + '_KHz'] * 1000.0).flatten()
            #sf.extend((main_mat[variable + '_KHz'] * 1000.0).flatten())
      #print(f"sf:{sf}")
      return accl, sf

# Get the specific variable names from a file that will be analyzed
def get_variable_names():
      print('Get Variable Names')
      vars=[]

      # Ask the user to enter the file with the variable list in it
      # variable list is a test case
      name_file = variable_list

      #name_file = input('Enter the file that contains the variable names: ')

      file = open(name_file, 'r')
      for line in file:
            vars.append(line.strip())
      file.close()
      #print(vars)
      # vars is a list of the variable names
      return vars


def modify_var_names(vars):
      print('Modify Variable Names - create a dictionary')
      #print('key=combine names:vale=string of values')
      #print('Key=sf: values=sampfreq')
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

      # print(combined.keys())

def print_greater_than_x(data, x):
    # Prints numbers greater than x and their indices in a list.

    for index, number in enumerate(data):
        if number > x:
            print(f"Number: {number}, Index: {index}")


def start_end(vars, combined, samp_freq):

# Get the average value of the accl variable

      # accl, sf = combine_variables(directory_path, mat_files, variable)
#      for i in range(len(combined)):
      l = [0]
      for i in range(len(l)):
            #print(f"vars: {vars}")
            print('Set time axis')
            print(f'combine: {vars}')
            accl = combined[vars + '_combine']

            # Rectify and square of values greater than 250

            print("Rectify and Normalizing Data")

            abs_val = np.absolute(accl)
            for j in range(len(abs_val)):
                  if abs_val[j] > 500.0:
                        abs_val[j] = 10000.0
                  else:
                        abs_val[j] = 0.0

            # Get the start and end times of each motor movement
            print("Get the start and end times of each motor movement from the accl variable")
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
                  if abs_val[j] == 10000.0 and flag == 1:
                        ts_times.append(j / samp_freq[vars])
                        flag = 0
                        # Add 10000 to the start point to get the full EMG signal
                        j = j + 10000
                        num_trials += 1
                  elif abs_val[j] == 0.0 and flag == 0:
                        te_times.append(j / samp_freq[vars])
                        flag = 1
                        j = j + 1

            print(f'num_trials: {num_trials}')
            #print(f'Length of ts_times: {len(ts_times)}')
            #print(f'Length of te_times: {len(te_times)}')
            #print(f'ts_times: {ts_times}')
            #ts_times = ts
            #te_times = te
            ts_samples = ts_times * samp_freq[vars]
            te_samples = te_times * samp_freq[vars]
            #print(f'ts_samples: {ts_samples}')


            # Create a x-axis for the time plot
            times = np.arange(len(abs_val)) * 1 / samp_freq[vars]

            print('Plot the accl data (both absolute and raw)')

            x1 = times

            ##plt.subplot(211)
            ##plt.title(f"Accelerometer data squared off and raw")
            ##plt.plot(x1, abs_val)
            ##plt.subplot(212)
            ##plt.plot(x1, combined[vars + '_combine'])
            ##plt.show()


            # Plot accl in the time and frequency domains
            print('Plot')

            j = 0

            hold0e = (te_times[j] * samp_freq[vars]) + (samp_freq[vars])
            hold0s = (ts_times[j] * samp_freq[vars])
            hold = len(accl[int(hold0s[0]):int(hold0e[0])])

            yavg = np.zeros(hold)
            while j < num_trials - 1:
                  #hold0s = (ts[j] * samp_freq[vars[i]]) - (samp_freq[vars[i]])
                  hold0e = (te_times[j] * samp_freq[vars]) + (samp_freq[vars])

                  hold0s = (ts_times[j] * samp_freq[vars])


                  hold0s = int(hold0s[0])
                  hold0e = int(hold0e[0])


                  #print(j)

                  # Rectify the values

                  #y1 = np.absolute(accl[hold0s:hold0e])
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

                  #x1 = times[hold0s:hold0e]
                  x1 = times[(hold0s - int(2 * samp_freq[vars][0])):hold0e]
                  y1 = np.absolute(accl[(hold0s - int(2 * samp_freq[vars][0])):hold0e])

                  ##print(f'Plot the segmented data for trial {j}.')

                  ##plt.subplot(211)
                  ##plt.title(f"Accl data for trial: {j}")
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
      plt.psd(yavg, NFFT=1024, Fs=samp_freq[vars])
      plt.show()

      return yavg, hold0s, hold0e, ts_times, te_times, ts_samples, te_samples



def segment_data_variable(variable, accl_start_samples, sf_accl, sf):
      # Segment the specific variable according to the start times of the accl movmemnt
      #print(f'len(variable): {len(variable)}')
      #print(f'accl_start_samples: {accl_start_samples}')
      # Segment the data to the following:
      #     2 sec before the movement
      #     2 sec during the movement
      #     2 sec after the movement

      # First convert the start times from the Accl variable to the appropriate sample points
      # vst = variable start times
      #                        (Accl_Start_time)*(Variable_Sample_rate)
      #  Variable_Start_time = -------------------------------------------
      #                             Accl_Sample_rate

      t_samp_var = np.round((accl_start_samples * sf) / sf_accl)
      vst = t_samp_var.flatten()
      #print(f'vst: {vst}')
      #print(f'vst[1]: {vst[1]}')
      #print(f'vst[2]: {vst[2]}')
      #print(f'vst[3]: {vst[3]}')
      #print(f'sf: {sf}')
      #print(f'Variable_Sample_rate (sf[0]): {sf[0]}')
      #print(f'Accl_Sample_rate (sf_accl): {sf_accl}')
      #print(f'Accl_Start_sample(accl_start_samples): {accl_start_samples}')
      #print(f'Variable Start Sample (vst): {vst}')
      before = {}
      after = {}
      during = {}
      all = {}

      #print(f'Len accl_start_samples: {len(accl_start_samples)}')
      #print(f'len(accl_start_samples): {len(accl_start_samples)}')
      for i in range(len(accl_start_samples)):
          during[i] = variable[int(vst[i]):(int(vst[i]) + int(2 * sf[0]) - 1)]
          before[i] = variable[(int(vst[i]) - int(sf[0] * 2)):(int(vst[i]) - 1)]
          after[i] = variable[(int(vst[i]) + int(sf[0] * 2)):(int(vst[i]) +
                                                       (2 * int(sf[0] * 2)) - 1)]
          all[i] = variable[(int(vst[i]) - int(sf[0] * 2)):(int(vst[i]) +
                                                      (2 * int(sf[0] * 2)) - 1)]
          #print(f'i:{i}')
          #print(f'Len of all[i]: {len(all[i])}')
          #plt.title(f"all[i]")
          #plt.plot(all[i])
          #plt.show()

      #print(f"Keys before: {before.keys()}")
      #print(f"Keys during: {during.keys()}")
      #print(f"Keys after: {after.keys()}")
      #print(f'len(all[1]): {len(all[1])}')
      #print(f'len(before[1]): {len(before[1])}')
      #print(f'len(during[1]): {len(during[1])}')
      #print(f'len(after[1]): {len(after[1])}')
      return(all, before, during, after)


def add_motor_condition(y_text, fontsize=14, color='k', ax=None):
    x_times = [-1.0, 1.0, 3.0]
    x_conditions = ['Pre', 'Movement', 'Post']
    if ax is None: ax = plt.gca()  # noqa
    plt.sca(ax)
    plt.axvline(0., lw=2, color=color)
    plt.axvline(2.0, lw=2, color=color)
    for x_t, t_t in zip(x_times, x_conditions):
        plt.text(x_t, y_text, t_t, color=color, fontsize=fontsize, ha='center',
                 va='center', fontweight='bold')

def get_avg_seg_variables(data, noe, sf, vars):
      # Calculate the average value for all the variables.

      #print(f'len(data[1]): {len(data[1])}')
      #print(f'len data: {len(data)}')

      total = np.zeros(len(data[1]))
      #print(f'len(total) after np.zeros: {len(total)}')
      for i in range(noe - 1):
            # Sum up all the data in the various epochs
            result = []
            for j in range(len(data[i])):
                  result.append(total[j] + data[i][j])
            #print(f'len(result): {len(result)}')
            total = result
      # Divide each value in the total by the number of epochs to get the average
      total = np.divide(total, noe)

      #print(f'len(total): {len(total)}')

      # Find the maximum amplitude of the data for placement of the text
      # pre/movement/post in the graph
      max_amplitude = max(total)

      #print(f'len(total): {len(total)}')

      x = []
      #print(f'sf: {sf}')
      for j in range(len(total)):
            x.append((j / sf) - 2.0)
            #if j == 1376:
                  #print(f'x[j]: {x[j]}')

      print(f'Segmented data average for: {vars}')

      plt.plot(x, total)
      plt.xlabel("Time(sec)")
      plt.ylabel("uV")
      plt.title(f"Average over segments for: {vars}")
      add_motor_condition(0.9 * max_amplitude)
      plt.show()


def calcualte_PSD(data, sf, noe, hold2, k):
      # noe = Number of Epochs
      print(f'Calculating the PSD for variable {hold2}')
      window_size = int(sf[0])
      nfft = window_size
      # Calculate PSD using Welch's method
      psd = np.zeros(nfft // 2 + 1)
      # Calculate the PSD and sum them
      for i in range(noe - 1): # Sum over the epochs
            #print(f'i:{i}')
            #print(f'Data[j][1]: {data[j][1]}')
            #windowed_data = data[k][i:(i + window_size)] * np.hanning(window_size)  # Apply Hanning window
            n = 0
            for j in range(int(len(data[i]) / sf[0])): # Sum over the data in the epochs (i.e. if 6 sec then would do six FFTs)
                  #windowed_data = data[i][j:(j + window_size)] * np.hanning(window_size)
                  #print(f'n: {n}')
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
      plt.title(f"Power Spectral Density (PSD): {hold2}")
      plt.yscale('log')
      # plt.ylim(0, 500)
      plt.xlim(0, 100)
      plt.grid(True)
      plt.show()
      return

def line_filter(variable, sf, line):
      f0 = line
      Q = 10.0

      b, a = signal.iirnotch(f0, Q, sf)
      filtered_data = signal.filtfilt(b, a, variable)
      return filtered_data


def citaf(data, noe, hold2, k ):
      # Combine the data into the format that tensorpac needs to read the data in
      # Format is an array that is (noe, data)

      print(f'Converting data to a format tensorpac needs for variable: {hold2}')

      #temp_array = np.full((noe, int(len(data[1]))), data)

      #print(f'data[1]: {data[1]}')
      #p = data[1]
      #print(f'p.shape: {p.shape}')
      #q = p[0]
      #r = p[0]
      #print(f'r: {r}')
      #print(f'q.shape: {q.shape}')
      #print(f'len(data): {len(data)}')
      #print(f'len(data[0]): {len(data[0])}')

      hold = len(data[0])
      hold_noe = noe - 1

      x = np.zeros((hold_noe, hold))



      for i in range(noe - 1):
            for j in range(len(data[0])):
                  x[i, j] = data[i][j]
                  hold_key = data.keys()
                  hold_value = data.values()



      #print(f'temp_array.shape: {temp_array.shape}')
      #tp_data_array = temp_array
      return x



if __name__ == '__main__':
      # directory_path and file extension are entered by the user
      mat_files = get_mat_files(directory_path, file_extension)
      # variable_list is a list of the variables in the mat file that are to be analyzed
      vars = get_variable_names()
      combined, samp_freq = modify_var_names(vars)

      # Add a 60 Hz filter to the data

      for k in range(len(vars)):
            print(f'Removing line noise from: {vars[k]}')
            sf = samp_freq[vars[k]]
            variable = combined[vars[k] + '_combine']
            filtered_data = line_filter(variable, sf[0], 60)
            print(f'Len of filtered data: {len(filtered_data)}')
      #f0 = 60.0
      #Q = 9.0
      #fs = samp_freq['CMacro_LFP_03']
      #b, a = signal.iirnotch(f0, Q, fs[0])
      #filtered_data = signal.filtfilt(b, a, combined['CMacro_LFP_03_combine'])

            combined[vars[k] + '_combine'] = filtered_data

      # Get the starting and ending times for the movements based on the accl
      # data
      for j in range(len(vars)):
            print(f'j: {j}')
            if vars[j] == 'CANALOG_IN_1___Accl':
                  # Store the variable sample rate that the segments are built on
                  sf_accl = samp_freq[vars[j]]
                  yavg, hold0s, hold0e, ts_times, te_times, ts_samples, te_samples = \
                        (start_end(vars[j], combined, samp_freq))

      # Build a dictionary of the data for the segments before, during
      # after and for all the data  - note we are really only using all the data
      # dictionary
      before_seg = {}
      during_seg = {}
      after_seg = {}
      all_seg = {}
      for k in range(len(vars)):
            print(f'k: {k}')
            hold = vars[k] + '_combine'
            hold2 = vars[k] + '_segmented'
            # Just for testing purposes - for plot to make data is correct
            # hold3 is the accelerometer segmented data label
            if k == 2:
                  hold3 = vars[k] + '_segmented'

            all, before, during, after = segment_data_variable(combined[hold],
                                    ts_samples, sf_accl, samp_freq[vars[k]])

            # Fill the dictionaries with the segmented data
            before_seg[hold2] = before
            during_seg[hold2] = during
            after_seg[hold2] = after
            all_seg[hold2] = all
            print(f'len(all): {len(all)}')
            print(f'len(all[1]): {len(all[1])}')
            #print(f'len(all_seg[hold2][k]):{len(all_seg[hold2][k])}')
            #print(f'len all[k]: {len(all[k])}')
      #print(f'Keys all_seg: {all_seg.keys()}')
      #print(f'Keys all_seg[hold2]: {all_seg[hold2].keys()}')
      print("Plot the complete segmented data for one segment")
      #print(f'hold0e: {hold0e}')
      #print(f'all_seg[hold2][16]: {all_seg[hold2][16]}')
      #print(f'len(all[hold2][2]): {len(all_seg[hold2][2])}')
      #print(f'samp_freq: {samp_freq}')
      #times = np.arange(len(all_seg[hold2])) * 1 / samp_freq[hold2]
      #times = np.arange(len(all_seg[hold2][10]))
      #print(f'hold3: {hold3}')
      times = np.arange(len(all_seg[hold3][1]))
      #print(f'len(times): {len(times)}')
      #x1 = all_seg[hold2][10]
      x1 = all_seg[hold3][1]
      plt.plot(times, x1)
      plt.title(f"Segmented data for first trial of: {hold3}")
      plt.show()

      for k in range(len(vars)):
            #print(f'samp_freq[vars[0]]: {samp_freq[vars[k]]}')

            # combined variables are the long dictionaries that hold all the data from one
            # variable across all mat files
            # segments variables are the variables that have been chopped at the start
            # of the accl movement and then two seconds in front and 2 seconds behind

            hold = vars[k] + '_combine'
            hold2 = vars[k] + '_segmented'
            calcualte_PSD(all_seg[hold2], samp_freq[vars[k]],
                          len(ts_samples), hold2, k)

      for k in range(len(vars)):
            hold2 = vars[k] + '_segmented'
            #print(f'all_seg[hold2][k].shape: {all_seg[hold2][k].shape}')
            # Number of epochs
            noe = len(ts_samples)
            #get_avg_seg_variables(all_seg[hold2][k], noe,
            #                      samp_freq[vars[k]], hold2)
            get_avg_seg_variables(all_seg[hold2], noe,
                            samp_freq[vars[k]], hold2)


      # Place the data in the format the tensorpac needs:
      #     Array of signals of shape (n_epochs, n_times).
      #  Will do for one value at first

      for k in range(len(vars)):
            #print(f'vars[k][0:7]: {vars[k][0:7]}')
            #print(f'vars[k][0:4]: {vars[k][0:4]}')
            if (vars[k][0:7] == 'CANALOG') or (vars[k][0:4] == 'CEMG'):
                  print(f'Skipping {vars[k]}')
            else:
                  hold2 = vars[k] + '_segmented'
                  noe = len(ts_samples)
                  tp_data_array = citaf(all_seg[hold2], noe, hold2, k)

                  # Get the time values
                  x = []
                  sf = samp_freq[vars[k]]
                  for j in range(len(all_seg[hold2][k])):
                        x.append((j / sf) - 2.0)

                  sf = sf[0] * 1.0

                  time = (np.arange(tp_data_array.shape[1]) / sf) - 2.0

                  # define an ERPAC object
                  p = EventRelatedPac(f_pha=[10, 30], f_amp=(50, 200, 1, 1))

                  #print(f'tp_data_array.shape: {tp_data_array.shape}')
                  #print(f'type(tp_data_array): {type(tp_data_array)}')
                  #print(tp_data_array[1, 1])

                  # extract phases and amplitudes
                  pha = p.filter(sf, tp_data_array, ftype='phase', n_jobs=1)
                  amp = p.filter(sf, tp_data_array, ftype='amplitude', n_jobs=1)

                  ###############################################################################
                  # Compute the ERPAC using the two implemented methods and plot it
                  ###############################################################################

                  # implemented ERPAC methods
                  methods = ['circular', 'gc']

                  plt.figure(figsize=(14, 8))
                  for n_m, m in enumerate(methods):
                        # compute the erpac
                        erpac = p.fit(pha, amp, method=m, smooth=100, n_jobs=-1).squeeze()

                        # plot
                        plt.subplot(len(methods), 1, n_m + 1)
                        p.pacplot(erpac, time, p.yvec, xlabel='Time (second)' * n_m,
                                  cmap='Spectral_r', ylabel='Amplitude frequency', title=p.method + vars[k],
                                  cblabel='ERPAC', vmin=0., rmaxis=True)
                        plt.axvline(0.0, linestyle='--', color='r', linewidth=2)
                        plt.axvline(2.0, linestyle='--', color='r', linewidth=2)

                  plt.tight_layout()
                  p.show()