def predictAmp(r, a0=1, freq=3, wavetype='body', Q=None, beta=None):
    """
    Calculates amplitude of body or surface wave given a radius from
    the source and several optional attenuation parameters.

    :param: radius: Distance between source and station (km)
    :param: a0: initial amplitude (default 1)
    :param: freq: center frequency (hz, default 3)
    :param: Q: quality factor (default 66 for body wave, 50 for surface wave)
    :param: beta: wave speed km/s (default 2.3 for body wave, 1 for surface wave)
    :return: amplitude in m/s

    .. note::
        Formula from Battaglia, 2003.
    """

    if wavetype == 'body':
        if not Q:
            Q = 66
        if not beta:
            beta = 2.3
        # Calculate B
        B = (np.pi * freq) / (Q * beta)
        # Calculate Amplitude
        amp = a0 * (np.exp(-B * r) / r)

    elif wavetype == 'surface':
        if not Q:
            Q = 50
        if not beta:
            beta = 1
        # Calculate B
        B = (np.pi * freq) / (Q * beta)
        # Calculate Amplitude
        amp = a0 * (np.exp(-B * r) / np.sqrt(r))

    return amp