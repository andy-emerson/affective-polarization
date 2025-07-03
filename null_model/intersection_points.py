import numpy as np
import matplotlib.pyplot as plt

# --- First Plot: Threshold vs BA Edges ---
def plot_threshold_vs_edges(ax, m=3):
    n_values = np.arange(m, 401)
    y1 = (n_values * np.log(n_values)) / 2
    y2 = m * (n_values - m)

    diff = y1 - y2
    sign_changes = np.where(np.diff(np.sign(diff)))[0]

    line1, = ax.plot(n_values, y1, color='steelblue', label='Connectivity Threshold for Erdös-Renyí')
    line2, = ax.plot(n_values, y2, color='indianred', linestyle='--', label='Edges in Barabási-Albert')

    dots = []
    for idx in sign_changes:
        n_cross = n_values[idx]
        y_cross = y1[idx]
        ax.axvline(x=n_cross, color='gray', linestyle='dashed')
        dot, = ax.plot(n_cross, y_cross, 'o', color='black', markersize=6, label='Intersection')
        dots.append(dot)

    ax.set_xlabel('Number of Nodes')
    ax.set_ylabel('Number of Edges')
    ax.set_title(f'Comparison for $m={m}$')
    ax.set_xlim(m, 400)
    ax.set_ylim(0, 1200)

    handles = [line1, line2, dots[0]]
    labels = [h.get_label() for h in handles]
    ax.legend(handles, labels)


# --- Second Plot: Exponential Fit to Intersections ---
def plot_exponential_fit(ax):
    m_points = np.array([1, 2, 3, 4, 5])
    n2_points = np.array([4.92155, 45.8579, 385.001, 2948.78, 21976.4])
    m_curve = np.linspace(1, 5, 200)
    n_curve = np.exp(2 * m_curve)

    ax.plot(m_curve, n_curve, '-', label=r'$n = e^{2m}$', color='steelblue')
    ax.plot(m_points, n2_points, 'o', label='Intersection Points', color='black')
    ax.fill_between(m_curve, n_curve, color='steelblue', alpha=0.3)

    ax.set_xlabel('$m$')
    ax.set_ylabel('$n$')
    ax.set_title('Exponential Fit to Intersection Points')
    ax.legend()


# --- Combine both plots ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
plot_threshold_vs_edges(ax1)
plot_exponential_fit(ax2)
plt.tight_layout()
plt.savefig('null_model/plots/intersection_points.png', dpi=300, bbox_inches='tight')
print("Plot saved as intersection_points.png")