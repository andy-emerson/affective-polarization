import matplotlib.pyplot as plt

n_values = list(range(3, 25))  # n from 3 to 24 inclusive
m = 3

# Compute E(n), average degree, and density using isolated seed math
E_values = [m*(n - m) for n in n_values]
avg_degree = [2*m - (2*m**2)/n for n in n_values]
density = [2*E/(n*(n-1)) for E, n in zip(E_values, n_values)]

# Set up two side-by-side plots
fig, axes = plt.subplots(1, 2, figsize=(14,5))  # 1 row, 2 columns

# --- Left Plot: Average Degree ---
axes[0].plot(n_values, avg_degree, marker='o', color='steelblue', linestyle='-', label='Average Degree')
limit = 2 * m
axes[0].axhline(
    y=limit,
    color='indianred',
    linestyle='--',
    label=fr'$\lim_{{n\to\infty}} \langle k \rangle = {limit}$'
)

axes[0].set_xlim(3, 24)
axes[0].set_ylim(0, 8)
axes[0].set_xticks(range(3, 25))
axes[0].set_yticks(range(0, 9))
axes[0].set_xlabel('Number of nodes (n)')
axes[0].set_ylabel('Average Degree')
axes[0].set_title('Average Degree vs Number of Nodes (m=3)')
axes[0].legend()

# --- Right Plot: Density ---
axes[1].plot(n_values, density, marker='o', color='steelblue', linestyle='-', label='Density')

axes[1].set_xlim(3, 24)
axes[1].set_ylim(0, 1)
axes[1].set_xticks(range(3, 25))
axes[1].set_yticks([i/10 for i in range(0, 11)])
axes[1].set_xlabel('Number of nodes (n)')
axes[1].set_ylabel('Density')
axes[1].set_title('Density vs Number of Nodes (m=3)')
axes[1].legend()

# Layout and save
plt.tight_layout()
plt.savefig('Comparing ER and BA/plots/degree_and_density.png', dpi=300, bbox_inches='tight')
print("Plot saved as degree_and_density.png")