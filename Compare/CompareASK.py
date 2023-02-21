import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc

# Define the signal-to-noise ratio (SNR) in dB
snr_dB = np.arange(0, 15, 1)
snr = 10 ** (snr_dB / 10)

# Define the number of bits
# N = 10**6

# Calculate the BER for ASK modulation
ber1 = 0.5 * (erfc(np.sqrt(0.1*snr))) # BASK 50% depth
ber2 = 0.5 * (erfc(np.sqrt(0.5*snr))) # BASK 100% depth

# Plot the BER diagram
plt.semilogy(snr_dB, ber1, 'o-', color='r',label = 'BASK 50%')
plt.semilogy(snr_dB, ber2, 'o-', color='g',label = 'BASK 100%')
plt.xlabel('Signal-to-Noise Ratio (dB)')
plt.ylabel('Bit Error Rate (BER)')
plt.title('BER Diagram for ASK Modulation')
plt.grid(True)
plt.yticks(fontsize=5)
plt.xticks(snr_dB,fontsize=12)
plt.legend()
plt.show()