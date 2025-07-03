import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import sys
sys.path.append('/workspaces/affective-polarization')
from functions import *

# Configuration
m = 3
n = 100
num_networks = 500

# Generate networks and measure gamma values
gamma_values = []
for i in range(num_networks):
    G = nx.barabasi_albert_graph(n, m)
    
    # Estimate gamma using powerlaw package
    degrees = [deg for _, deg in G.degree()]
    fit = powerlaw.Fit(degrees, discrete=True, xmin=m, verbose=False)
    gamma = fit.alpha
    
    gamma_values.append(gamma)

gamma_array = np.array(gamma_values)

# Plot histogram
plt.hist(gamma_array, bins=50, color='steelblue', edgecolor='black')
plt.axvline(2.3, color='olive', linestyle='--', label='2.3')
plt.axvline(3.0, color='indianred', linestyle='--', label='3.0')
plt.title(f"Histogram of Fitted γ Values (n={n}, m={m})")
plt.xlabel("γ (scaling exponent)")
plt.ylabel("Frequency")
plt.legend()
plt.tight_layout()
plt.savefig('network_generation/plots/gamma_histogram.png', dpi=300, bbox_inches='tight')
print("Plot saved as gamma_histogram.png")

# Print summary
mask = (gamma_array >= 2.3) & (gamma_array <= 3.0)
pct = 100 * np.sum(mask) / len(gamma_array)
print(f"Percentage of gamma values between 2.3 and 3.0: {pct:.2f}% ({np.sum(mask)}/{len(gamma_array)})")