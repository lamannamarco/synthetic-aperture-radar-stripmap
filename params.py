""" Class containing the FMCW radar parameters """

import numpy as np
from tabulate import tabulate
from scipy.constants import speed_of_light as c

class Params:


    def __init__(self, fc=10e9, fs=10e6, R_max=3e3, bandwidth=150e6):
        """ Initializing the class

        The chirp length is not an input: it is derived so that the maximum beat
        frequency (target at R_max) falls exactly at the edge of the ADC band.
        """

        # Primary parameters
        self.fc = fc  # [Hz] Center frequency (X band)
        self.fs = fs  # [Hz] Sampling frequency of the beat signal
        self.R_max = R_max  # [m] Maximum range
        self.bandwidth = bandwidth  # [Hz] Chirp bandwidth

        # Derived parameters
        self.lc = c / self.fc  # [m] Operating wavelength
        self.alpha = self.fs * (2 * self.R_max / c)**-1  # [Hz/s] Chirp rate such that f_beat(R_max) = fs (complex I/Q samples)
        self.Tc = self.bandwidth / self.alpha  # [s] Chirp length
        self.N_ADC = int(np.ceil(self.Tc * self.fs))  # ADC samples per chirp
        self.Tc_full = self.N_ADC / self.fs  # [s] Chirp length rounded to an integer number of ADC samples
        self.t_full = np.arange(self.N_ADC) / self.fs  # [s] Fast time (sampled at fs)
        self.R_delta = c / (2 * self.bandwidth)  # [m] Range resolution

        Params.print_params_to_console(self)


    def print_params_to_console(self):
        """ Printing to console all the radar parameters """

        data = [
                ['fc', f'{self.fc*1e-9:.3f}', 'GHz', ''],
                ['Wavelength', f'{self.lc*1e3:.3f}', 'mm', ''],
                ['B', f'{self.bandwidth*1e-6:.3f}', 'MHz', ''],
                ['fs', f'{self.fs*1e-6:.3f}', 'MHz', ''],
                ['Nadc', f'{self.N_ADC}', '/', 'Samples per chirp'],
                ['Tc', f'{self.Tc_full*1e6:.3f}', 'us', 'Nadc/fs'],
                ['B/Tc', f'{self.alpha*1e-6/1e6:.3f}', 'MHz/us', ''],
                ['Maximum range', f'{self.R_max:.3f}', 'm', ''],
                ['Range resolution', f'{self.R_delta:.3f}', 'm', ''],
        ]
        headers = ['Item', 'Value', 'Unit', 'Comments']
        table = tabulate(data, headers=headers, disable_numparse=True)
        print('\nRADAR PARAMETERS TABLE')
        print(table)
