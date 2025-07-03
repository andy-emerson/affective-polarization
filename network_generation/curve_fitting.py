import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from functions import *

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Parameters
n = 300      # Number of nodes
k = 30       # Number of graphs to generate
m = 3        # Set manually

# --- Collect data ---
degree_counts_list = []  # list of pandas Series for each graph
all_degrees = []         # flattened degrees across all graphs
gamma_list = []          # list of gamma estimates

for _ in range(k):
    _, G = generate_networks(n, m)
    
    degrees = np.array([deg for _, deg in G.degree()])
    degree_counts = pd.Series(degrees).value_counts().sort_index()
    degree_counts_list.append(degree_counts)

    all_degrees.extend(degrees.tolist())

    # Estimate gamma directly from the graph
    # Estimate gamma directly from the graph
    degrees_for_gamma = [deg for _, deg in G.degree()]
    fit = powerlaw.Fit(degrees_for_gamma, discrete=True, xmin=m, verbose=False)
    gamma_estimate = fit.alpha
    gamma_list.append(gamma_estimate)

# --- Aggregate median counts ---
degree_counts_df = pd.concat(degree_counts_list, axis=1).fillna(0)
median_counts = degree_counts_df.median(axis=1)

# --- Log-log linear regression on median counts ---
degrees = median_counts.index.values
counts = median_counts.values

mask = counts > 0  # Only use degrees with positive count
log_degrees = np.log10(degrees[mask])
log_counts = np.log10(counts[mask])

slope, intercept, r_value, p_value, std_err = linregress(log_degrees, log_counts)

# --- Final overall gamma estimate (average of all graph gammas) ---
overall_gamma = np.mean(gamma_list)

# --- Prepare data for plotting ---
all_degrees = np.array(all_degrees)

# --- Plotting ---

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# First plot: Degree histogram (fixed bins of width 1)
bins = np.arange(min(all_degrees), max(all_degrees) + 1, 1)
hist_vals, bin_edges, _ = axes[0].hist(all_degrees, bins=bins, edgecolor="black", color='steelblue', alpha = 0.8, label='Degree Histogram')

# Lock y-axis limit based on histogram
ymin, ymax = axes[0].get_ylim()

# Compute bin centers
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

# Overlay shifted power-law fit line, with gamma label
x_fit = bin_centers
valid_mask = (x_fit - m) > 0
x_fit_valid = x_fit[valid_mask]
y_fit = ((x_fit_valid - m) ** (-overall_gamma)) * np.max(hist_vals)

axes[0].plot(x_fit_valid, y_fit, color='indianred', linestyle='-', label=f'Power-law Fit (γ ≈ {overall_gamma:.2f})')

# Restore original y-limits
axes[0].set_ylim(ymin, ymax)

axes[0].set_xlabel('Degree')
axes[0].set_ylabel('Frequency')
axes[0].set_title('Degree Histogram with Power-Law Fit')
axes[0].legend()

# Second plot: Linear regression fit and scatter
axes[1].scatter(degrees, counts, edgecolor="black",color='steelblue', label='Median Degree Counts', zorder=2)
axes[1].plot(degrees, 10**(intercept + slope * np.log10(degrees)), color='indianred', linestyle='-', label=f'Linear Fit (slope = {slope:.2f})', zorder=1)
axes[1].set_xscale('log')
axes[1].set_yscale('log')
axes[1].set_xlabel('Degree')
axes[1].set_ylabel('Median Count')
axes[1].set_title('Log-Log Linear Regression on Median Counts')
axes[1].legend()

plt.tight_layout()
plt.savefig('network_generation/plots/curve_fitting.png', dpi=300, bbox_inches='tight')
print("Plot saved as curve_fitting.png")