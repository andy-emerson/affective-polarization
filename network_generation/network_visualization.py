import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from functions import *

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Generate networks
n = 40
m = None
ER, BA = generate_networks(n, m)

# Calculate what the number of edges should be
if m == None:
    m = 1
    while True:
        size = m * (n - m)
        if m + 1 <= n <= np.exp(2 * m):
            break
        m += 1
size = m * (n - m)

# Estimate gamma
degrees = [deg for _, deg in BA.degree()]
fit = powerlaw.Fit(degrees, discrete=True, xmin=m, verbose=False)
gamma = fit.alpha

# Print graph properties
print(f"Returned size:     {size}")
print(f"Actual size of ER: {ER.number_of_edges()}")
print(f"Actual size of BA: {BA.number_of_edges()}")
print(f"γ estimate of BA: {gamma:.2f}")

# Create figure and two subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Plot first graph
plt.sca(ax1)
nx.draw(ER, node_color='steelblue', edge_color='gray', node_size=100)
ax1.set_title(f'Erdős–Rényi (n={n}, e={ER.number_of_edges()})')

# Plot second graph
plt.sca(ax2)
nx.draw(BA, node_color='indianred', edge_color='gray', node_size=100)
ax2.set_title(f'Barabási–Albert (n={n}, e={BA.number_of_edges()}, γ={gamma:.2f})')

plt.savefig('network_generation/plots/network_visualization.png', dpi=300, bbox_inches='tight')
print("Plot saved as network_visualization.png")