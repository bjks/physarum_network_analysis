from __future__ import division
import numpy as np
from scipy import integrate, interpolate
from scipy.signal import butter, lfilter, spectrogram

slack_l, slack = 0.1, 1
cutoff = 50
L = 25

from scipy.io import wavfile
x = 
sr, x = wavfile.read('capriccio.wav')
x = x[:(L + slack) * sr, 0]
x = x

# sr = 44100
# x = np.random.normal(size=((L + slack) * sr,))

b, a = butter(2, 2 * cutoff / sr, btype='low')  # Butterworth

# cutoff function
def f(t):
    return (10000 - 1000 * np.clip(t, 0, 9) - 1000 * np.clip(t-19, 0, 0.8)) \
            / cutoff

# and its reciprocal
def fr(_, t):
    return cutoff / (10000 - 1000 * t.clip(0, 9) - 1000 * (t-19).clip(0, 0.8))

# modulate time
# calculate upper end of td first
tdmax = integrate.quad(f, 0, L + slack_l, points=[9, 19, 19.8])[0]
span = (0, tdmax)
t = np.arange(x.size) / sr
tdinfo = integrate.solve_ivp(fr, span, np.zeros((1,)),
                             t_eval=np.arange(0, span[-1], 1 / sr),
                             vectorized=True)
td = tdinfo.y.ravel()
# modulate signal
xd = interpolate.interp1d(t, x)(td)
# and linearly filter
yd = lfilter(b, a, xd)
# modulate signal back to linear time
