def plot_fft(t, signal, reconstruct=True, **kwargs):
    """This function will compute the Fast Fourier Transform of a given signal
    and plot the results as a function of period.  If desired, it will also
    re-compute the signal based on the FFT results (a summation of a series of 
    sines and cosines with various frequencies) to verify the FFT results.
    
    Inputs:
    t           = Signal's x-axis inputs in seconds
    signal      = Signal, same length as t, any unit
    reconstruct = Option to re-compute the signal based on FFT results 
                  (adds time, default is True)
    title       = Title for plots
    period_xlim = X-limits to use for period subplot
    period_ylim = Y-limits to use for period subplot
    """
    
    # Compute number of data points
    n = len(t)
    
    # Compute FFT
    FFT = fft.fft(signal)  
    FFT = FFT[:n/2]  #extract the positive-frequency terms
    
    # Compute frequencies
    freqs = fft.fftfreq(signal.size, t[1]-t[0])
    freqs = freqs[:n/2]  #extract the positive-frequency terms
    
    # Compute periods
    T = 1 / freqs / 3600 #hours
        
    # Compute amplitudes 
    amp = abs(FFT) * 2 / n

    # Compute amplitudes 
    phase = angle(FFT) * 2 / n
    
    # Compute Fourier Coefficients
    a = real(FFT) * 2 / n
    b = -imag(FFT) * 2 / n
    
    # Re-compute original signal based on FFT results (for verification)
    def reconstruct_signal(t, freqs, a, b):
        comp_sig = zeros(len(t))
        for j in arange(len(t)):
            comp_sig[j] = a[0]/2 + sum([(
                          a[i] * cos(2 * pi * freqs[i] * (t[j] - t[0])) + 
                          b[i] * sin(2 * pi * freqs[i] * (t[j] - t[0]))) 
                          for i in arange(1, len(a))])     
        return comp_sig
    if reconstruct == True:
        signal_fft = reconstruct_signal(t, freqs, a, b)
    
    # Plot results
    # Original vs Computed signal
    figure(figsize=[8,11])
    subplot(211)
    if kwargs.has_key('cxctime'):
        pl = plot_cxctime
    else:  
        pl = plot
        xlabel('Time [s]')
    pl(t, signal, 'b', label='orig')
    if reconstruct == True:
        pl(t, signal_fft, 'r', label='fft')
        legend()
    ylabel('Signal')
    if kwargs.has_key('title'):
        title(kwargs['title'])
    # Amplitude vs Period
    subplot(212)
    plot(T, amp,'x-')
    xlabel('Period [hrs]')
    ylabel('Signal Amplitude')
    if kwargs.has_key('period_xlim'):
        xlim(kwargs['period_xlim'])
    if kwargs.has_key('period_ylim'):
        xlim(kwargs['period_ylim'])
    tight_layout()


close('all')    

x = fetch.Msidset(['pitch', 'pftank1t', 'pftank2t', 'pline01t', 'pline02t', 
                   'pline03t', 'pline04t', 'pline05t', 'pline06t', 
                   'pline07t'], '2011:001', '2013:001', stat='5min')

# IPS Tank
plot_fft(x['pftank1t'].times, x['pftank1t'].vals, reconstruct=False, title='IPS Tank:  PFTANK1T', cxctime=True, period_xlim=[0,50], period_ylim=[0,1])
savefig('fft_pftank1t.png')
plot_fft(x['pftank2t'].times, x['pftank2t'].vals, reconstruct=False, title='IPS Tank:  PFTANK2T', cxctime=True, period_xlim=[0,50], period_ylim=[0,1])
savefig('fft_pftank2t.png')

# Circuit 5105 - Not cycling
plot_fft(x['pline01t'].times, x['pline01t'].vals, reconstruct=False, title='Circuit 5105:  PLINE01T', cxctime=True, period_xlim=[0,50], period_ylim=[0,4])
savefig('fft_pline01t.png')
plot_fft(x['pline02t'].times, x['pline02t'].vals, reconstruct=False, title='Circuit 5105:  PLINE02T', cxctime=True, period_xlim=[0,50], period_ylim=[0,4])
savefig('fft_pline02t.png')
plot_fft(x['pline03t'].times, x['pline03t'].vals, reconstruct=False, title='Circuit 5105:  PLINE03T', cxctime=True, period_xlim=[0,50], period_ylim=[0,4])
savefig('fft_pline03t.png')
plot_fft(x['pline04t'].times, x['pline04t'].vals, reconstruct=False, title='Circuit 5105:  PLINE04T', cxctime=True, period_xlim=[0,50], period_ylim=[0,4])
savefig('fft_pline04t.png')

# Circuit 5107
plot_fft(x['pline05t'].times, x['pline05t'].vals, reconstruct=False, title='Circuit 5107:  PLINE05T', cxctime=True, period_xlim=[0,50], period_ylim=[0,2])
savefig('fft_pline05t.png')
plot_fft(x['pline06t'].times, x['pline06t'].vals, reconstruct=False, title='Circuit 5107:  PLINE06T', cxctime=True, period_xlim=[0,50], period_ylim=[0,2])
savefig('fft_pline06t.png')
plot_fft(x['pline07t'].times, x['pline07t'].vals, reconstruct=False, title='Circuit 5107:  PLINE07T', cxctime=True, period_xlim=[0,50], period_ylim=[0,2])    
savefig('fft_pline07t.png')    