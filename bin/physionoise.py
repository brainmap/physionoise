import string, sys
import logging

import numpy, scipy
from numpy import *
from scipy import *
from scipy import linalg, linspace, polyval, polyfit, sqrt, stats, randn, signal, interpolate
from scipy.signal import cspline1d, cspline1d_eval, iirdesign, lfilter, lp2lp
from scipy.signal.signaltools import hilbert
from scipy.ndimage import spline_filter1d
from scikits.timeseries.lib.moving_funcs import cmov_average
from pylab import *
from matplotlib import *
from matplotlib import pyplot as plt
from matplotlib.mlab import psd
from matplotlib.mlab import slopes



usage = """Usage: physionoise -c cardiac_signal -o cardiac_trigger -r respiratory_signal [options]"""

about = """Author: Dan Kelley
    Waisman Laboratory for Brain Imaging and Behavior
    University of Wisconsin-Madison
    Last Update: November 1, 2007

This program performs the following:

1  Removes noise with a low pass Butterworth filter with zero phase shift
   
2  Calculates the residual = Raw signal - filtered signal
   
3  Idenitifies dirty peaks with big residuals that need cleaning
   
4  Dirty peaks +/- a window are removed to produce filtered, cleaned data
   
5  Calculates the best fitting spline through the cleaned data
   
6  Spline interpolates over the removed dirty peaks to produce a clean
   respiratory waveform (RW) and cardiac waveform (CW)
   This is fast on Mac OS. Other platforms use the --speed option.
   
7  Finds the respiratory peaks (RPpd) and cardiac peaks (CPpd).
   Peak finding was slightly modified from the open source peakdet.m algorithm 
   to include magnitude and time thresholds and ported into python.
   
8  The top and bottom envelopes are generated using the RPpd. The 
   respiratory volumes over time (RVT) is their difference as in:
      Birn RM, Diamond JB, Smith MA, Bandettini PA. Separating 
      respiratory-variation-related fluctuations from 
      neuronal-activity-related fluctuations in fMRI. Neuroimage. 2006 Jul 
      15;31(4):1536-48. Epub 2006 Apr 24. PMID: 16632379 
   Same is done for the cardiac waveforms using CPpd to produce CVT curves.
   
9  Cardiac rate time (CRT) courses based on the TTL (CRTttl) are generated as the 
   inverse change in time between CP centered points on the initial peak as in:
      Shmueli K, van Gelderen P, de Zwart JA, Horovitz SG, Fukunaga M, Jansma
      JM, Duyn JH. Low-frequency fluctuations in the cardiac rate as a source
      of variance in the resting-state fMRI BOLD signal. Neuroimage. 2007
      Nov 1;38(2):306-20. Epub 2007 Aug 9. PMID: 17869543
   An RRT is also generated in the same manner using the RPpd.
   
10 Cardiac rate time (CRT) courses based on the third derivative (CRTd3)
   and the third derivative's R-wave estimate (CRTd3R) are generated as in:
      Chan GS, Middleton PM, Celler BG, Wang L, Lovell NH. Automatic detection
      of left ventricular ejection time from a finger photoplethysmographic
      pulse oximetry waveform: comparison with Doppler aortic measurement.
      Physiol Meas. 2007 Apr;28(4):439-52. Epub 2007 Mar 20. PMID: 17395998
   
11 PhysioNoise outputs the RW (downsampled to 40Hz by default), CPd3, CPd3R,
   and CPttl time courses as options for 3dretroicor and the RVT, RRT, CVT,
   CRTd3, CRTd3R, and CRTttl curves sampled on the TR and half TR as
   covariate options for neuroanalysis. 
"""


# this is the logger object that will be available gloabally, careful!
logger = None


def parseargs():
	from optparse import OptionParser, OptionGroup

	parser = OptionParser(usage=usage)
	parser.add_option("-r","--resp",action="store",type="string",dest="sigfile",
		help="read respiratory signal from FILE",
		metavar="FILE")
	parser.add_option("-c","--card",action="store",type="string",dest="csigfile",
		help="read cardiac signal from FILE",
		metavar="FILE")
	parser.add_option("-o","--ox",action="store",type="string",dest="ctrigfile",
		help="read cardiac pulse ox trigger from FILE",
		metavar="FILE")
	parser.add_option("--TR",action="store",type="float",dest="TR",
		help="TR is TR seconds [2.0]",
		metavar="TR", default=2.0)
	parser.add_option("--numTR",action="store",type="int",dest="NumTR",
		help="Number of TRs [480]",
		metavar="TR", default=480)
	parser.add_option("--inputHz",action="store",type="int",dest="OrigSamp",
		help="input frequency in Hz[40]",
		metavar="Hz", default=40)
	parser.add_option("--outputHz",action="store",type="int",dest="NewSamp",
		help="desired output frequency in Hz[40]",
		metavar="Hz", default=40)
	parser.add_option("--fpass",action="store",type="float",dest="fpass",
		help="Butterworth filter:Stop frequency of passband in Hz[2.0]",
		metavar="Hz", default=2.0)
	parser.add_option("--fstop",action="store",type="float",dest="fstop",
		help="Butterworth filter:Start frequency of the stopband in Hz[10.0]",
		metavar="Hz", default=10.0)
	parser.add_option("--trimwindow",action="store",type="int",dest="TrimWindow",
		help="Number of points PTS to clean on either side of dirty residual peaks[10]",
		metavar="PTS", default=10)
	parser.add_option("--plr",action="store",type="float",dest="FwdRm",
		help="Percentage of points on the left of dirty peaks to clean on the right [1.2]",
		metavar="Percentage", default=1.2)
	parser.add_option("--statthresh",action="store",type="int",dest="FoldSTD",
		help="Mark points beyond FOLD standard deviation(s) as dirty\
		residual peaks for cleaning  [6]",
		metavar="FOLD", default=6)
	parser.add_option("--rmagthresh",action="store",type="int",dest="RMagThresh",
		help="Respiratory Peak must be RMAG units away from nearest extrema  [100]",
		metavar="RMAG", default=100)
	parser.add_option("--rtimethresh",action="store",type="float",dest="RTimeThresh",
		help="Respratory Peak must be RTIME seconds away from the nearest extrema  [0]",
		metavar="RTIME", default=0)
	parser.add_option("--cmagthresh",action="store",type="int",dest="CMagThresh",
		help="Cardiac Peak must be CMAG units away from nearest extrema  [20]",
		metavar="CMAG", default=20)
	parser.add_option("--ctimethresh",action="store",type="float",dest="CTimeThresh",
		help="Cardiac Peak must be CTIME seconds away from the nearest extrema  [0]",
		metavar="CTIME", default=0)
	parser.add_option("-l","--plots",action="store_true",dest="PlotAll",
		help="this turns on interactive plotting",
		default=False)
	parser.add_option("-s","--speed",action="store_true",dest="Speed",
		help="This skips artifact detection and spline interpolation.\
		Instead a cubic spline filter is used and does a decent job removing\
		artifact. [artifact detect]",
		default=False)
	parser.add_option("-p","--prefix",action="store",type="string",dest="Prefix",
		help="Output prefix [PreIcor]",
		default='PreIcor')
	parser.add_option("-v","--verbose",action="count",dest="verbosity",
		help="Indicate wordiness by adding more v's",
		default=0)
	parser.add_option("-a","--about",action="store_true",dest="about",
		help="Displays information about the program's authors and how it works.")
	parser.add_option("--saveFiles", action="store_true",dest="saveFiles",
		help="Write a gamut regressors out to separate text files.",
		default=False)
	parser.add_option("-t","--truncate", action="store", type="string", dest="truncate",
		help="Truncate the signals to expected length based on TRs. Truncation can \
		be either at the front, [end], or none.",
		metavar="SIDE", default="end")

	options, args = parser.parse_args()

	if options.about:
		print(about)
		sys.exit(0)

	if not (options.sigfile and options.csigfile and options.ctrigfile):
		print "Error: Respiratory signal, cardiac signal, and cardiac\
		trigger files are all required."
		parser.print_help()
		sys.exit(0)
	
	return (options, args)


def initialize_logger():
	"""docstring for initialize_logger"""
	global logger
	logger = logging.getLogger("physionoise")
	console_handler = logging.StreamHandler()
	console_formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
	console_handler.setFormatter(console_formatter)
	logger.addHandler(console_handler)
	logger.setLevel(logging.WARNING)


def chart_summary_of_signals(signals, labels, figure):
	plt.figure(figure)
	position_base = 321
	for i in range(len(signals)):
		plt.subplot(position_base + i)
		plt.plot(signals[i])
		plt.ylabel(labels[i])


def chart_envelopes(signal_spline, maxima, topenv, minima, botenv, svt, title, legend, figure, row=1):
	position = 210 + row
	plt.figure(figure)
	plt.subplot(position)
	plt.hold(True)
	plt.plot(signal_spline, 'b')
	plt.plot(topenv, 'g')
	plt.plot(botenv, 'g')
	plt.plot(numpy.where(maxima)[0], signal_spline[maxima], 'ro')
	plt.plot(numpy.where(minima)[0], signal_spline[minima], 'ro')
	plt.plot(svt, 'c')
	plt.plot(numpy.zeros(len(signal_spline)), 'k')
	plt.hold(False)
	plt.legend(legend, shadow = False, loc = 1)
	ltext = plt.gca().get_legend().get_texts()
	plt.setp(ltext[0], fontsize = 8, color = 'b')
	plt.setp(ltext[1], fontsize = 8, color = 'g')
	plt.setp(ltext[2], fontsize = 8, color = 'g')
	plt.setp(ltext[3], fontsize = 8, color = 'r')
	plt.setp(ltext[4], fontsize = 8, color = 'r')
	plt.setp(ltext[5], fontsize = 8, color = 'c')
	plt.title(title)


def extrema_indices(timeseries, ds, dt):
	"""Finds indices of all extrema in the given timeseries.

	Modified from peakdet.m by Eli Billauer, 3.4.05 (Explicitly not copyrighted).
	This function is released to the public domain; Any use is allowed.
	Two paramters control the algorithm.
	ds: signal change necessary to identify an extremum.
	dt: the minimum distance in timepoints between any two extrema
	"""
	# lists to hold the found extrema indices and magnitudes
	e_indices, e_magnitudes = [], []
	# Definitions:
	# cem -- current_extremum_magnitude
	# cei -- current_extremum_index
	cem = float("-infinity")
	cei = -2 * dt
	ascending = True

	for i in xrange(len(timeseries)):

		signal = timeseries[i]

		if ascending and signal > cem:
			cem = signal
			cei = i
		elif not ascending and signal < cem:
			cem = signal
			cei = i

		beyond_last_extremum = len(e_indices) == 0 or (cei - e_indices[-1]) > dt
		beginning_to_descend = signal < (cem - ds)
		beginning_to_ascend  = signal > (cem + ds)

		if beyond_last_extremum:
			if ascending and beginning_to_descend:
				e_indices.append(cei)
				e_magnitudes.append(cem)
				ascending = False
			elif not ascending and beginning_to_ascend:
				e_indices.append(cei)
				e_magnitudes.append(cem)
				ascending = True

	return e_magnitudes, e_indices


def extrema(timeseries, ds, dt):
	"""docstring for extrema"""
	e_magnitudes, e_indices = extrema_indices(timeseries, ds, dt)
	return numpy.array([ i in e_indices for i in range(len(timeseries)) ])


def minima(timeseries, ds, dt):
	e_magnitudes, e_indices = extrema_indices(timeseries, ds, dt)
	offset = 0 if e_magnitudes[0] < e_magnitudes[1] else 1
	e_magnitudes, e_indices = e_magnitudes[offset::2], e_indices[offset::2]
	return numpy.array([ i in e_indices for i in range(len(timeseries)) ])


def maxima(timeseries, ds, dt):
	e_magnitudes, e_indices = extrema_indices(timeseries, ds, dt)
	offset = 0 if e_magnitudes[0] > e_magnitudes[1] else 1
	e_magnitudes, e_indices = e_magnitudes[offset::2], e_indices[offset::2]
	return numpy.array([ i in e_indices for i in range(len(timeseries)) ])


def zeromin(signal):
	"""simply adjusts a signal so the min is zero"""
	return signal - numpy.nanmin(signal)

def dump(signal, filebase, prefix=None, samplerate=None):
	"""dumps a signal to a text file"""
	filename = filebase
	if samplerate: filename = "%s_%s" % (filename, samplerate)
	if prefix:     filename = "%s_%s" % (prefix, filename)
	filename += ".txt"
	numpy.savetxt(filename, signal, fmt='%8.6f')


def stat_summary(an_array):
	"""returns a tuple (N, mean, standard_deviation)"""
	return (len(an_array), mean(an_array), std(an_array))


def report_value(message, value):
	"""Displays a formatted message / value pair"""
	print("%-40s\t%s" % (message, value))


def report_values(msg_val_pairs):
	"""Displays a list of formatted message / value pairs
	
	feed it a list of message / value tuples"""
	for kv in msg_val_pairs:
		report_value(*kv)



def lfilter_zi(b,a):
	"""compute the zi state from the filter parameters. see [Gust96].
	
	Based on:
	[Gust96] Fredrik Gustafsson, Determining the initial states in forward-backward 
	filtering, IEEE Transactions on Signal Processing, pp. 988--992, April 1996, 
	Volume 44, Issue 4
	"""
	n = max(len(a),len(b))
	zin = ( eye(n-1) - hstack( (-a[1:n,newaxis], vstack((eye(n-2),zeros(n-2)))) ) )
	zid = b[1:n] - a[1:n]*b[0]
	zi_matrix = linalg.inv(zin)*(matrix(zid).transpose())
	zi_return = []
	
	#convert the result into a regular array (not a matrix)
	for i in range(len(zi_matrix)):
		zi_return.append(float(zi_matrix[i][0]))
		
	return array(zi_return)



def filtfilt(b, a, x):
	"""performs a filter on a signal"""
	# b: numerator, a: denominator, x: signal
	# For now only accepting 1d arrays
	ntaps = max(len(a), len(b))
	edge = ntaps*3
	if x.ndim != 1:
		raise ValueError, "filtfilt is only accepting 1 dimension arrays."
	
	#x must be bigger than edge
	if x.size < edge:
		raise ValueError, "Input vector needs to be bigger than 3 * max(len(a),len(b)."
		
	if len(a) < ntaps: a = r_[a, zeros(len(b)-len(a))]
	if len(b) < ntaps: b = r_[b, zeros(len(a)-len(b))]
		
	zi = lfilter_zi(b, a)
	
	# Grow the signal to have edges for stabilizing 
	# the filter with inverted replicas of the signal
	s = r_[ 2*x[0] - x[edge:1:-1], x, 2*x[-1] - x[-1:-edge:-1] ]
	
	# in the case of one go we only need one of the extrems
	# both are needed for filtfilt
	(y, zf) = lfilter(b, a, s, -1, zi*s[0])
	(y, zf) = lfilter(b, a, flipud(y), -1, zi*y[-1])
	
	return flipud(y[edge-1:-edge+1])


def downsample(signal, ratio):
	"""downsamples a signal by a given sampling ratio"""
	resample_size = int(len(signal) / ratio)
	return scipy.signal.resample(signal, resample_size)

def resample_to_tr(signal, input_sample_rate, num_trs, tr_time, tr_multiplier=1.0):
	"""resamples a signal in tr space.
	
	Checks whether the given tr information
	is consistent with the length of the signal.  Optionally you can specify to
	resample in some multiple of tr's as well."""
	run_duration = num_trs * tr_time
	expected_signal_length = int(run_duration * input_sample_rate)
	resample_size = int(num_trs / tr_multiplier)
	
	if len(signal) != expected_signal_length:
		truncated_signal = signal[0:expected_signal_length]
		return scipy.signal.resample(truncated_signal, resample_size)
	else:
		return scipy.signal.resample(signal, resample_size)


def connect_dots(x, y, xmin=None, xmax=None):
	"""given a sequence of x,y pairs connects them with line segments and returns
	an entire filled in signal series"""
	if (xmin and x[0] != xmin):
		x = numpy.append([0],x)
	if (xmax and x[-1] != xmax):
		x = numpy.append(x,[0])
	spline = scipy.interpolate.interp1d(x, y)
	return spline( range(x[0], x[-1]+1) )

def envelope(signal, extrema):
	"""finds the given locations (usually extrema) on a signal and
	connects the dots"""
	extrema[0] = True; extrema[-1] = True
	return connect_dots(indices_of(extrema), signal[extrema])

def event_rate(events, sample_rate):
	"""Given some events in a series (booleans), returns the rate of the events
	over time"""
	event_times = indices_of(events)
	dt = diff(event_times)
	dt = float(sample_rate) / dt
	event_times[0] = 0
	event_times = numpy.append(event_times, len(events)-1)
	dt = numpy.append(dt[0], dt)
	dt = numpy.append(dt, dt[-1])
	return connect_dots(event_times, dt)


def spectral_analysis(filtered_signal, sample_rate, signal_name="unknown"):
	"""Runs a quick spectral analysis of a filtered signal, returns a bunch of
	information about it."""
	x           = int( (log(sample_rate) - log(0.01)) / log(2) )
	NFFT        = int(2**x)
	nooverlap   = NFFT / 2
	pxx, freqs  = psd( filtered_signal,
		NFFT=NFFT, Fs=sample_rate, noverlap=nooverlap,
		detrend=detrend_linear, window=window_hanning
	)
	i = pxx.argmax() # index of the maximum power frequency
	max_frequency = freqs[i]
	max_power     = pxx[i]
	period        = 1.0 / max_frequency
	
	if max_frequency == 0:
		logger.warning("%s power peak at zero. Threshold at %f Hz" % (signal_name, freqs[4]))
		i = pxx[5:].argmax()
		max_frequency = freqs[i+5]
		max_power     = pxx[i+5]
		period        = 1.0 / max_frequency
	
	print "Butterworth Filter: %s signal" % signal_name
	report_value("Spectral peak (Hz)", max_frequency)
	report_value("Peak power density", max_power)
	report_value("Spectral period (s)", period)
	
	return(pxx, freqs, max_frequency, max_power, period)


def build_butterworth_filter(passband_frequency, stopband_frequency):
	"""Builds a butterworth filter and returns the paramters needed to use the
	filter on a signal."""
	# filter_b and filter_a are numerator and denominator of the low-pass filter
	# original: gpass=5 gstop=60
	filter_b, filter_a = iirdesign(
		passband_frequency, stopband_frequency,
		gpass=5, gstop=20, analog=0, ftype="butter"
	)
	
	w, h = scipy.signal.freqz(filter_b, filter_a)
	response_magnitude = abs(h)
	sample_frequencies = w / pi
	
	return (filter_b, filter_a, response_magnitude, sample_frequencies)


def censor_signal(signal, filtered_signal, fold, left_trim, right_trim_percent):
	"""Finds points where residual is greater than the given fold, and trims
	given number of point from the left, and some multiple of that from the right.
	"""
	ltrim = int(left_trim)
	rtrim = int(left_trim * right_trim_percent)

	indexes = numpy.arange(len(signal))
	residual = signal - filtered_signal
	threshold = float(fold) * numpy.std(residual)

	censored = numpy.absolute(residual) > threshold
	accepted = numpy.zeros(len(residual)) == 0
	
	logger.info( "Censoring n dirty peaks: %d" % (len(numpy.where(censored)[0])) )
	
	for i in numpy.where(censored)[0]:
		for j in range(i - ltrim, i + rtrim):
			if (j>=0 and j<len(accepted)): accepted[j] = False
			
	return (signal[accepted], filtered_signal[accepted], indexes[accepted])


def load_signal(filename):
	"""loads float valued signal from a file, one point per line"""
	return numpy.fromfile(filename, 'f', sep='\n')


def mean_centered(signal):
	"""centers a signal on its mean"""
	signal_mean = numpy.mean(signal[0:len(signal)/2])
	return signal - signal_mean


def indices_of(bool_array):
	"""returns an array of indexes of each event in a boolean array"""
	return numpy.where(bool_array)[0]


def spikewave(signal, locations):
	"""creates a spikewave given a signal and an array of events (booleans)"""
	spike = numpy.zeros(len(signal))
	spike[locations] = signal[locations]
	return spike



def main(options, args):
	"""The main procedure of physionoise"""
	initialize_logger()
	if (options.verbosity == 1): logger.setLevel(logging.INFO)
	if (options.verbosity == 2): logger.setLevel(logging.DEBUG)
	if (options.verbosity > 2):  logger.setLevel(0)
	
	respiratory_signal     =  mean_centered( load_signal(options.sigfile) )
	respiratory_triggers   =  load_signal(options.ctrigfile).astype('int')
	cardiac_signal         =  mean_centered( load_signal(options.csigfile) )
	cardiac_triggers       =  load_signal(options.ctrigfile).astype('int')
	
	respiratory_xaxis      = range(len(respiratory_signal))
	cardiac_xaxis          = range(len(cardiac_signal))
	
	input_sample_rate      = options.OrigSamp
	output_sample_rate     = options.NewSamp
	sample_ratio           = float(input_sample_rate) / output_sample_rate
	nyquist                = input_sample_rate / 2.0
	passband_frequency     = options.fpass / nyquist
	stopband_frequency     = options.fstop / nyquist
	
	tr_time                = options.TR
	num_tr                 = options.NumTR
	input_signal_length    = int(num_tr * tr_time * input_sample_rate)
	output_signal_length   = int(num_tr * tr_time * output_sample_rate)
	
	# these options are for tuning the peakdet algorithm
	respiratory_dm         = options.RMagThresh
	respiratory_dt         = options.RTimeThresh  # in seconds!
	cardiac_dm             = options.CMagThresh
	cardiac_dt             = options.CTimeThresh  # in seconds!
	
	plots_requested        = options.PlotAll
	use_fast_spline        = options.Speed
	prefix                 = options.Prefix
	truncate               = options.truncate
	
	# converting absolute event times from the trigger file to a boolean array
	CPttl = (numpy.zeros(len(cardiac_signal)) != 0) # array of Falses same size as signal
	CPttl[cardiac_triggers] = True
	
	
	### NOTE on encoding events ###
	# in all cases where we want to denote a series of events in a signal, we will
	# use an array of booleans.  For example, if we want to know where signal
	# peaks are we create an array of the same size as the signal with True entries
	# where each of the peaks occured and False everywhere else.
	# These boolean arrays can be used conveniently to index into the signals
	# e.g. signal[peaks] => produces the magnitude of the signal at the peaks.
	#
	
	
	####  Butterworth Filter Data  ###
	filter_b, filter_a, butter_response_magnitude, butter_sample_frequences = build_butterworth_filter(
		passband_frequency, stopband_frequency
	)
	
	filtered_respiratory = filtfilt(filter_b, filter_a, respiratory_signal)
	filtered_cardiac     = filtfilt(filter_b, filter_a, cardiac_signal)
	
	# the spectral analyses are just for plotting, they are currently not saved
	# to file at the end
	pxx, freqs, maxfreq, PeakPower, TPeriod = spectral_analysis(
		filtered_respiratory, input_sample_rate, "Respiratory"
	)
	
	Cpxx, Cfreqs, Cmaxfreq, CPeakPower, CTPeriod = spectral_analysis(
		filtered_cardiac, input_sample_rate, "Cardiac"
	)
	
	
	
	####  Spline interpolation  ####
	respiratory_residual = respiratory_signal - filtered_respiratory
	cardiac_residual     = cardiac_signal - filtered_cardiac

	if use_fast_spline:
		logger.info( "Speedy Spline Filtering" )
		respiratory_spline  = scipy.ndimage.spline_filter1d(filtered_respiratory, order=3)
		cardiac_spline      = scipy.ndimage.spline_filter1d(filtered_cardiac, order=3)
	else:
		logger.info( "Slower (better?) Spline Filtering" )
		logger.info( "Censoring regions of high variance" )
		
		censored_resp, resp_included, resp_included_indexes = censor_signal(
			respiratory_signal, filtered_respiratory,
			options.FoldSTD, options.TrimWindow, options.FwdRm
		)
		censored_cardiac, card_included, card_included_indexes = censor_signal(
			cardiac_signal, filtered_cardiac,
			options.FoldSTD, options.TrimWindow, options.FwdRm
		)
		
		logger.debug("resp len(%d) censored resp len(%d)" % (len(respiratory_signal), len(censored_resp)))
		logger.debug("card len(%d) censored card len(%d)" % (len(cardiac_signal), len(censored_cardiac)))
		
		logger.info( "Spline interpolation" )
		# Generate cubic spline....Too slow when not on Mac OS X but looks nice
		resp_spline_estimate  = scipy.interpolate.splrep(resp_included_indexes, resp_included, k=3)
		respiratory_spline    = scipy.interpolate.splev(respiratory_xaxis, resp_spline_estimate)
		card_spline_estimate  = scipy.interpolate.splrep(card_included_indexes, card_included, k=3)
		cardiac_spline        = scipy.interpolate.splev(cardiac_xaxis, card_spline_estimate)
	
	
	
	#----------------------------------------------------------  Finding peaks
	respiratory_extrema = extrema(respiratory_spline, respiratory_dm, respiratory_dt)
	respiratory_maxima  = maxima(respiratory_spline, respiratory_dm, respiratory_dt)
	respiratory_minima  = minima(respiratory_spline, respiratory_dm, respiratory_dt)
	
	cardiac_extrema     = extrema(cardiac_spline, cardiac_dm, cardiac_dt)
	cardiac_maxima      = maxima(cardiac_spline, cardiac_dm, cardiac_dt)
	cardiac_minima      = minima(cardiac_spline, cardiac_dm, cardiac_dt)
	
	
	
	#------------------------------------------------------- Find Volume over Time
	logger.info( "Calculating Envelopes" )
	logger.info( "Respiratory" )
	
	resp_top_envelope    = envelope(respiratory_spline, respiratory_maxima)
	resp_bottom_envelope = envelope(respiratory_spline, respiratory_minima)
	resp_abs_volume      = resp_top_envelope - resp_bottom_envelope
	resp_rate            = event_rate(respiratory_maxima, input_sample_rate) 
	rvt                  = resp_abs_volume * resp_rate
	
	logger.info( "Cardiac" )
	
	card_top_envelope    = envelope(cardiac_spline, cardiac_maxima)
	card_bottom_envelope = envelope(cardiac_spline, cardiac_minima)
	card_abs_volume      = card_top_envelope - card_bottom_envelope
	card_rate            = event_rate(cardiac_maxima, input_sample_rate)
	cvt                  = card_abs_volume * card_rate
	
	
	
	#------ Find third derivatives of Cardiac signal to create new trigger file
	logger.info( "Finding Cardiac Derivatives" )
	
	cardiac_d1 = slopes(cardiac_xaxis, cardiac_spline) #* input_sample_rate
	cardiac_d2 = slopes(cardiac_xaxis, cardiac_d1) #* input_sample_rate
	cardiac_d3 = slopes(cardiac_xaxis, cardiac_d2) #* input_sample_rate
	smoothed_caradiac_d3 = cmov_average(cardiac_d3, input_sample_rate / 2)
	
	# caution: using 1.0 for the amplitude threshold for peak finding, this might
	# not always work well.
	card_d3_maxima  = maxima(smoothed_caradiac_d3, 1.0, cardiac_dt)
	card_d3_minima  = minima(smoothed_caradiac_d3, 1.0, cardiac_dt)
	
	logger.info( "Thresholding CPd3" )
	CPd3 = card_d3_minima
	card_d3_bottom_envelope = envelope(smoothed_caradiac_d3, card_d3_minima)
	card_d3_top_envelope    = envelope(smoothed_caradiac_d3, card_d3_maxima)
	# There have been several approaches to filtering acceptable d3 minima 
	# to include in CPd3:
	# 1. the min must be in a region where d1 is negative
	# 2. the min must be less than the top envelope at that point
	# 3. the min must be less than 0.5 * the mean of the top envelope
	#
	# The current approach, calculate the bottom envelope and compute a moving
	# average of it at 1.5 * the sampling rate, accept any minima that
	# are <= the moving average.  The motivation for this is that in most cases
	# there are two classes of minima that are relatively seperable locally but 
	# not globally.  The moving average provides a consistent margin between the 
	# classes at the local level.
	
	CPd3 &= (card_d3_bottom_envelope <= cmov_average(card_d3_bottom_envelope, input_sample_rate*1.5))
	
	# Older minima filtering approaches
	#
	# CPd3 &= (cardiac_d1 >= -0.5)
	# CPd3 &= (card_d3_bottom_envelope <= mean(card_d3_bottom_envelope))
	# CPd3 &= (abs(card_d3_bottom_envelope) >= abs(card_d3_top_envelope))
	# CPd3 &= (abs(card_d3_bottom_envelope) >= 0.5*mean(abs(card_d3_top_envelope)))
	# CPd3 = -1.0 * spikewave(smoothed_caradiac_d3, CPd3)
	
	# R waves are d3 peaks that immediately precede CPd3 blips
	CPd3Rwave = (numpy.zeros(len(CPd3)) != 0)  # boolean array of all falses.
	for i in indices_of(CPd3):
		previous_d3_maxima = card_d3_maxima[0:i]
		if (len(indices_of(previous_d3_maxima)) > 0):
			index_of_prior_d3_maximum = indices_of(previous_d3_maxima)[-1]
			CPd3Rwave[index_of_prior_d3_maximum] = True
	
	
	
	#----------- Calculate respiratory (RRT) and cardiac rate per time (CRT)
	logger.info( "Calculating Rates over Time" )
	
	RRT    = event_rate(respiratory_maxima, input_sample_rate)
	CRTttl = event_rate(CPttl, input_sample_rate)
	CRTd3  = event_rate(CPd3, input_sample_rate)
	CRTd3R = event_rate(CPd3Rwave, input_sample_rate)
	
	
	
	
	#---------------------------------------------------------------- DownSampling
	logger.info( "Downsampling" )
	
	downsampled_respiratory = downsample(respiratory_spline, sample_ratio)
	downsampled_cardiac     = downsample(cardiac_spline,     sample_ratio)
	
	
	
	#-------------------------------------------------------------------- Plotting
	if plots_requested:
		# Butterworth filter analysis summary
		plt.figure(1)
		plt.title("Butterworth filter spectral analysis")
		plt.subplot(311)
		plt.semilogy(butter_sample_frequences * nyquist, butter_response_magnitude)
		plt.ylabel('Magnitude')
		plt.title('Bode Plot')
		plt.ylim((10e-4, 10e0))
		plt.xlim((0.0, options.fstop))
		plt.figure(1)
		plt.subplot(312) 
		plt.plot(freqs, pxx)
		plt.ylabel('Respiratory Power')
		plt.xlim((0.0, options.fstop))
		plt.figure(1)
		plt.subplot(313) 
		plt.plot(Cfreqs, Cpxx)
		plt.xlabel('frequency(Hz)')
		plt.ylabel('Cardiac Power')
		plt.xlim((0.0, options.fstop))
		
		# overview of the filter and spline interpolated signals
		plt.figure(2)
		plt.subplot(211)
		plt.plot(respiratory_xaxis, filtered_respiratory, 'b-',
		         respiratory_xaxis, respiratory_spline,   'go',
		         respiratory_xaxis, respiratory_signal,   'k',
		         respiratory_xaxis, respiratory_residual, 'r'
		)
		plt.legend(('Filtered', 'Spline', 'Raw', 'Residual'), shadow=False, loc=1)
		ltext = plt.gca().get_legend().get_texts()
		plt.setp(ltext[0], fontsize = 8, color = 'b')
		plt.setp(ltext[1], fontsize = 8, color = 'g')
		plt.setp(ltext[2], fontsize = 8, color = 'k')
		plt.setp(ltext[3], fontsize = 8, color = 'r')
		plt.title('Respiratory')
		plt.subplot(212)
		plt.plot(cardiac_xaxis, filtered_cardiac, 'b-',
		         cardiac_xaxis, cardiac_spline,   'go',
		         cardiac_xaxis, cardiac_signal,   'k',
		         cardiac_xaxis, cardiac_residual, 'r'
		)
		plt.legend(('Filtered','Spline','Raw','Residual'), shadow=False, loc=1)
		ltext = plt.gca().get_legend().get_texts()
		plt.setp(ltext[0], fontsize = 8, color = 'b')
		plt.setp(ltext[1], fontsize = 8, color = 'g')
		plt.setp(ltext[2], fontsize = 8, color = 'k')
		plt.setp(ltext[3], fontsize = 8, color = 'r')
		plt.title('Cardiac')
		
		# top and bottom peak finding and envelopes
		chart_envelopes(respiratory_spline,
			respiratory_maxima, resp_top_envelope,
			respiratory_minima, resp_bottom_envelope,
			rvt, 'Respiratory Data',
			('Respiratory Spline','Top Envelope','Bottom Envelope','Maxima','Minima','RVT'),
			figure=3, row=1
		)
		chart_envelopes(cardiac_spline,
			cardiac_maxima, card_top_envelope,
			cardiac_minima, card_bottom_envelope,
			cvt, 'Cardiac Data',
			('Cardiac Spline','Top Envelope','Bottom Envelope','Maxima','Minima','CVT'),
			figure=3, row=2
		)
		
		# summary of the volume over time waves and  cardiac wave event finding results
		signals_to_plot = (CRTd3R, CRTttl, RRT, CRTd3, rvt, cvt)
		labels = ('CRTd3R', 'CRTttl', 'RRT', 'CRTd3', 'RVT', 'CVT')
		chart_summary_of_signals(signals_to_plot, labels, figure=5)
		
		# cardiac wave events viewed on the spline 
		# make cardiac_d3 in roughly the same range as the spline
		adjustment_ratio = numpy.nanmax(cardiac_spline) / numpy.nanmax(smoothed_caradiac_d3)
		adjusted_d3      = smoothed_caradiac_d3 * 1.5 * adjustment_ratio
		plt.figure(6)
		plt.title('Event Locations on the Cardiac Waveform')
		plt.hold()
		plt.plot(cardiac_signal, '#aaaaaa')
		plt.plot(cardiac_spline, 'k')
		plt.plot(adjusted_d3, '#ff5555')
		plt.plot(indices_of(CPttl), cardiac_spline[CPttl], 'yo')
		plt.plot(indices_of(CPd3), cardiac_spline[CPd3], 'bo')
		plt.plot(indices_of(CPd3Rwave), cardiac_spline[CPd3Rwave], 'go')
		plt.legend(('Raw Cardiac','Cardiac Spline','Cardiac d3','CPttl','CPd3','CPd3Rwave'), shadow=False, loc=1)
		ltext = plt.gca().get_legend().get_texts()
		plt.setp(ltext[0], fontsize = 8, color = '#aaaaaa')
		plt.setp(ltext[1], fontsize = 8, color = 'k')
		plt.setp(ltext[2], fontsize = 8, color = '#ff5555')
		plt.setp(ltext[3], fontsize = 8, color = 'y')
		plt.setp(ltext[4], fontsize = 8, color = 'b')
		plt.setp(ltext[5], fontsize = 8, color = 'g')
		
		
	
	
	#--------------------------------------------------- Trim length for retroicor
	logger.info( "Trim length" )
	
	# expected length at input sample rate
	ei = tr_time * num_tr * input_sample_rate
	
	# expected length at downsample rate
	ed =   tr_time * num_tr * output_sample_rate
	
	if (truncate == 'end'):
		logger.info( "Truncating at end of signal" )
		respiratory_spline      = respiratory_spline[0:ei]
		cardiac_spline          = cardiac_spline[0:ei]
		CPd3                    = CPd3[0:ei]
		CPd3Rwave               = CPd3Rwave[0:ei]
		CPttl                   = CPttl[0:ei]
		RRT                     = RRT[0:ei]
		CRTd3                   = CRTd3[0:ei]
		CRTttl                  = CRTttl[0:ei]
		CRTd3R                  = CRTd3R[0:ei]
		resp_top_envelope       = resp_top_envelope[0:ei]
		card_top_envelope       = card_top_envelope[0:ei]
		downsampled_cardiac     = downsampled_cardiac[0:ed]
		downsampled_respiratory = downsampled_respiratory[0:ed]
	elif (truncate == 'front'):
		logger.info( "Truncating at front of signal" )
		ei_front = len(respiratory_spline) - ei
		ed_front = len(downsampled_respiratory) - ed
		
		respiratory_spline      = respiratory_spline[ei_front:]
		cardiac_spline          = cardiac_spline[ei_front:]
		CPd3                    = CPd3[ei_front:]
		CPd3Rwave               = CPd3Rwave[ei_front:]
		CPttl                   = CPttl[ei_front:]
		RRT                     = RRT[ei_front:]
		CRTd3                   = CRTd3[ei_front:]
		CRTttl                  = CRTttl[ei_front:]
		CRTd3R                  = CRTd3R[ei_front:]
		resp_top_envelope       = resp_top_envelope[ei_front:]
		card_top_envelope       = card_top_envelope[ei_front:]
		downsampled_cardiac     = downsampled_cardiac[ed_front:]
		downsampled_respiratory = downsampled_respiratory[ei_front:]
	else:
		logger.info( "Not truncating signal to expected TR size." )
	
	
	#---------------------------------------------------------------- Reporting
	peaktypes = {
		"CPttl": CPttl, "CPd3": CPd3, "CPd3Rwave": CPd3Rwave
	}
	
	rpt = report_value # just aliasing the report_value function
	tblhdr_fmt = "%-14s%-14s%-14s"
	tblrow_fmt = "%-14d%-14.4f%-14.4f"
	lbl_fmt = "%-14s%-14s"
	divider_width = 90
	divider = "=" * divider_width
	
	print
	rpt( lbl_fmt % ("Peak type","Aspect"), tblhdr_fmt % ('N', 'mean', 'std') )
	print divider
	
	for peak_label, peakseries in peaktypes.iteritems():
		mags    = cardiac_spline[peakseries]
		times   = indices_of(peakseries)
		dt      = diff(times)
		gt3std  = mags[mags > mean(mags) + 3*std(mags)]
		gt4std  = mags[mags > mean(mags) + 4*std(mags)]
		gtmean  = mags[mags > mean(mags) + 1.0]
		gtbelow = mags[mags > mean(mags) - 0.4]
		
		reports = {
			"all timepoints": peakseries, "peak magnitudes": mags, "peak times": times,
			"interpeak times": dt,        "peaks > 3std": gt3std,  "peaks > 4std": gt4std,
			"peaks > mean + 1": gtmean,   "peaks > mean - 0.4": gtbelow }
		for report_label, report_series in reports.iteritems():
			n, m, s = stat_summary(report_series)
			rpt(lbl_fmt % (peak_label, report_label), tblrow_fmt % (n,m,s))
			
	
	# Dumping results to files
	if (options.saveFiles):
		logger.info( "Saving files" )
		
		signals  = { 'resp_spline': respiratory_spline,    'card_spline': cardiac_spline,
		             'RVT'      : rvt,
		             'RRT'      : RRT,                     'CRTd3' : CRTd3,
		             'CRTd3R'   : CRTd3R,                  'CRTttl': CRTttl }
		for fbase, signal in signals.iteritems():
			dump(zeromin(signal), fbase, prefix, input_sample_rate)
			for tr_label, tr_multiple in { 'TR_': 1.0, 'HalfTR_': 0.5 }.iteritems():
				resampled_signal = resample_to_tr(signal, input_sample_rate, num_tr, tr_time, tr_multiple)
				dump(zeromin(resampled_signal), tr_label + fbase, prefix, input_sample_rate)
				
		peaksigs = { 'CPd3': CPd3, 'CPd3R': CPd3Rwave, 'CPttl': CPttl }
		for fbase, peak in peaksigs.iteritems():
			dump(peak, fbase, prefix, input_sample_rate)
		
		# probably not needed, left over from original physionoise.py
		dump(zeromin(downsampled_respiratory), 'resp_spline_downsampled', prefix, output_sample_rate)
	
	
	if plots_requested:
		plt.show()
	
	
	sys.exit(0)


if __name__ == '__main__':
	options, args = parseargs()
	main(options, args)