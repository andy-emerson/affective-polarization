# Confirmed Required
import networkx as nx
import numpy as np
import powerlaw
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sqlite3
import pandas as pd
from scipy.stats import pearsonr, linregress
from scipy.stats import ConstantInputWarning
from scipy.sparse.linalg import eigsh
import warnings

######################
### NETWORK MODELS ###
######################

def generate_networks(n, m=None, gamma_bounds=[2.3, 3.0]):
    """
    Generate a connected Erdős–Rényi graph and a Barabási–Albert graph, 
    ensuring the BA graph has a gamma within given bounds.

    Parameters
    ----------
    n : int
        Number of nodes.
    m : int or None, optional
        Edges to attach from a new node in BA model. 
        If None, the function automatically selects the smallest valid m.
    gamma_bounds : list of float, optional
        Lower and upper bounds for acceptable gamma [lower, upper].

    Returns
    -------
    er_graph : networkx.Graph
        A connected Erdős–Rényi graph.
    ba_graph : networkx.Graph
        A Barabási–Albert graph with acceptable gamma.
    gamma : float
        Power-law exponent of the BA graph.
    size : int
        Number of edges in the ER graph (and approximately the BA graph).
    """
    assert isinstance(gamma_bounds, list) and len(gamma_bounds) == 2, "gamma_bounds must be a list of two numbers"
    lower, upper = gamma_bounds

    if m is None:
        assert n >= 2, "n must be at least 2"
        # Find smallest m satisfying conditions
        m = 1
        while True:
            size = m * (n - m)
            if m + 1 <= n <= np.exp(2 * m):
                break
            m += 1
    else:
        # Manual m provided: input validation
        assert m >= 1, "m must be at least 1"
        assert m + 1 <= n <= np.exp(2 * m), f"n must satisfy {m+1} ≤ n ≤ e^(2m)"
        size = m * (n - m)

    # --- MAIN LOOP ---
    
    # Generate a connected ER graph
    while True:
        er_graph = nx.gnm_random_graph(n, size)
        if nx.is_connected(er_graph):
            break

    # Generate BA graph and validate gamma
    while True:
        ba_graph = nx.barabasi_albert_graph(n, m)
        degrees = [deg for _, deg in ba_graph.degree()]
        fit = powerlaw.Fit(degrees, discrete=True, xmin=m, verbose=False)
        gamma = fit.alpha
        if lower <= gamma <= upper:
            break

    return er_graph, ba_graph

###################
### AFFILIATION ###
###################

def random_subgraphs(graph, ratios):
    """
    Partition a graph's nodes into disjoint subgraph views (not copies),
    according to given relative ratios using the Hamilton (largest remainder) method.

    Parameters
    ----------
    graph : networkx.Graph
        Input graph (unchanged).
    ratios : list of numbers
        Relative proportions of nodes assigned to each subgraph.

    Returns
    -------
    subgraphs : list of networkx.Graph
        List of new independent subgraph objects.
    """

    ratios = np.array(ratios, dtype=float)
    proportions = ratios / ratios.sum()

    n_nodes = graph.number_of_nodes()

    # --- Hamilton Method for group sizes ---
    raw_sizes = proportions * n_nodes
    base_sizes = np.floor(raw_sizes).astype(int)
    remainders = raw_sizes - base_sizes

    missing = n_nodes - base_sizes.sum()

    if missing > 0:
        extra_indices = np.argsort(-remainders)[:missing]
        for idx in extra_indices:
            base_sizes[idx] += 1

    # --- Shuffle and slice ---
    nodes = list(graph.nodes())
    np.random.shuffle(nodes)

    subgraphs = []
    idx = 0
    for size in base_sizes:
        subgraph = graph.subgraph(nodes[idx:idx+size])
        subgraphs.append(subgraph)
        idx += size

    return subgraphs

#----------------------------------------------

def spectral_subgraphs(graph):
    """
    Partition a connected graph into two communities using the sign of the Fiedler vector 
    (second-smallest eigenvector of the normalized Laplacian matrix). This corresponds to 
    basic spectral clustering via a binary cut.

    Parameters
    ----------
    graph : networkx.Graph
        An undirected, connected graph.

    Returns
    -------
    subgraphs : list of networkx.Graph
        Two disjoint subgraphs corresponding to the sign-based partition.
    """

    assert nx.is_connected(graph), "Graph must be connected."

    # Compute normalized Laplacian
    L = nx.normalized_laplacian_matrix(graph).astype(float)

    # Compute the two smallest eigenvalues and corresponding eigenvectors
    _, vecs = eigsh(L, k=2, which="SM")

    # The second eigenvector is the Fiedler vector
    fiedler_vector = vecs[:, 1].flatten()

    # Partition nodes by sign of Fiedler vector
    nodes = np.array(graph.nodes())
    group1 = nodes[fiedler_vector >= 0]
    group2 = nodes[fiedler_vector < 0]

    return [graph.subgraph(group1), graph.subgraph(group2)]


#----------------------------------------------

def modular_subgraphs(graph):
    """
    Partition a graph into two communities using the sign of the leading eigenvector
    of the modularity matrix (Newman's spectral modularity method).

    Parameters
    ----------
    graph : networkx.Graph
        An undirected graph.

    Returns
    -------
    subgraphs : list of networkx.Graph
        Two disjoint subgraphs corresponding to the sign-based partition.
    """

    # Adjacency matrix
    A = nx.to_numpy_array(graph)

    # Degree vector and total edge count
    k = A.sum(axis=1)
    m = k.sum() / 2

    # Construct modularity matrix: B = A - (k_i * k_j) / (2m)
    B = A - np.outer(k, k) / (2 * m)

    # Compute leading eigenvector of modularity matrix
    _, vecs = eigsh(B, k=1, which="LA")
    leading_vector = vecs[:, 0]

    # Partition nodes by sign of the leading eigenvector
    nodes = np.array(graph.nodes())
    group1 = nodes[leading_vector >= 0]
    group2 = nodes[leading_vector < 0]

    return [graph.subgraph(group1), graph.subgraph(group2)]

###################
### ATTRIBUTION ###
###################

def random_opinions(G, subgraphs, opinion_ratios):
    """
    Assigns binary opinions randomly within each subgraph.

    Parameters
    ----------
    G : networkx.Graph
        The full graph.
    subgraphs : list of networkx.Graph
        List of disjoint subgraph views from G.
    opinion_ratios : list of float
        Fraction of nodes in each subgraph to randomly assign the elite opinion (1).

    Returns
    -------
    networkx.Graph
        The same graph G with node attribute "opinion" assigned.
    """

    assert len(subgraphs) == len(opinion_ratios), \
        "Each subgraph must have a corresponding opinion_ratio."

    opinions = {}

    for sg, ratio in zip(subgraphs, opinion_ratios):
        nodes = list(sg.nodes())
        np.random.shuffle(nodes)

        n_select = int(np.ceil(ratio * len(nodes)))
        selected = set(nodes[:n_select])

        for node in nodes:
            opinions[node] = 1 if node in selected else 0

    nx.set_node_attributes(G, opinions, "opinion")
    return G

# ---------------------------------------------------------------------

def centrality_opinions(G, subgraphs, elite_opinions, opinion_ratios):
    """
    Assigns binary opinions to nodes based on subgraph membership and eigenvector centrality.

    Parameters
    ----------
    G : networkx.Graph
        The full graph (not the subgraphs).
    subgraphs : list of networkx.Graph
        List of disjoint subgraph views from G (not copies).
    elite_opinions : list of int (0 or 1)
        The elite opinion assigned to each subgraph.
    opinion_ratios : list of float (0.0 to 1.0)
        Fraction of top-centrality nodes in each subgraph to receive the elite opinion.

    Returns
    -------
    networkx.Graph
        The same graph G with node attribute "opinion" assigned.
    """

    assert len(subgraphs) == len(elite_opinions) == len(opinion_ratios), \
        "Each subgraph must have a corresponding elite_opinion and opinion_ratio."

    # Compute eigenvector centrality for the full graph
    centrality = nx.eigenvector_centrality_numpy(G)

    opinions = {}

    for sg, elite_opinion, ratio in zip(subgraphs, elite_opinions, opinion_ratios):
        nodes = list(sg.nodes())
        scores = [(node, centrality[node]) for node in nodes]
        scores.sort(key=lambda x: x[1], reverse=True)  # sort by centrality descending

        n_select = int(np.ceil(ratio * len(nodes)))
        selected = {node for node, _ in scores[:n_select]}

        for node in nodes:
            opinions[node] = elite_opinion if node in selected else 1 - elite_opinion

    # Assign opinions to the full graph
    nx.set_node_attributes(G, opinions, "opinion")
    return G


####################
### POLARIZATION ###
####################

def neighborhood_update(G, node_list, alpha=1, beta=1, delta=0):
    """
    Selects a random node from a given list and updates its opinion based on the polarization formula:
        α(b - a) - β(d - c)

    Parameters:
    G (networkx.Graph): Input graph with "partition" and "opinion" attributes.
    node_list (list): List of node IDs to choose from.
    alpha (float): Weight for same-partition influence.
    beta (float): Weight for different-partition influence.
    delta (float): Decision threshold.

    Returns:
    None: Modifies G in-place.
    """
    assert all(node in G for node in node_list), "All nodes in node_list must be in G."

    # Select a random node from the given list
    node = np.random.choice(node_list)

    # Get node attributes
    node_partition = G.nodes[node]["affiliation"]
    node_opinion = G.nodes[node]["opinion"]

    # Initialize counts
    a = b = c = d = 0

    # Count neighbors based on partition and opinion
    for neighbor in G.neighbors(node):
        neighbor_partition = G.nodes[neighbor]["affiliation"]
        neighbor_opinion = G.nodes[neighbor]["opinion"]

        # Same partition
        if neighbor_partition == node_partition:
            if neighbor_opinion == 0:
                a += 1
            else:
                b += 1
        # Different partition
        else:  
            if neighbor_opinion == 0:
                c += 1
            else:
                d += 1

    # Compute polarization score
    polarization_score = alpha * (b - a) - beta * (d - c)

    # Update opinion based on threshold (keep as boolean)
    if polarization_score < -delta:
        G.nodes[node]["opinion"] = 0
    elif polarization_score > delta:
        G.nodes[node]["opinion"] = 1

def polarizing_neighborhoods(G, alpha=1, beta=1, delta=0):
    """
    Loops through all nodes and identifies which would change their opinions
    based on the polarization formula:
        α(b - a) - β(d - c)

    Parameters:
    G (networkx.Graph): Input graph with "affiliation" and "opinion" attributes.
    alpha (float): Weight for same-partition influence.
    beta (float): Weight for different-partition influence.
    delta (float): Decision threshold.

    Returns:
    List[int]: A list of node IDs that would change their opinions.
    """
    polarizing_nodes = []

    for node in G.nodes():
        node_partition = G.nodes[node]["affiliation"]
        node_opinion = G.nodes[node]["opinion"]

        a = b = c = d = 0

        # Count neighbors based on partition and opinion
        for neighbor in G.neighbors(node):
            neighbor_partition = G.nodes[neighbor]["affiliation"]
            neighbor_opinion = G.nodes[neighbor]["opinion"]

            # Same partition
            if neighbor_partition == node_partition:
                if neighbor_opinion == 0:
                    a += 1
                else:
                    b += 1
            # Different partition
            else:
                if neighbor_opinion == 0:
                    c += 1
                else:
                    d += 1

        polarization_score = alpha * (b - a) - beta * (d - c)

        # Determine if the node would change its opinion
        if (polarization_score < -delta and node_opinion == 1) or (polarization_score > delta and node_opinion == 0):
            polarizing_nodes.append(node)

    return polarizing_nodes

##################
### SIMULATION ###
##################

def time_sim(G: nx.Graph, t: int, alpha=1, beta=1, delta=0) -> None:
    """
    Runs the polarization simulation for a fixed number of steps.
    At each step, considers all nodes in G for update.

    Parameters:
    - G: A NetworkX graph where each node has 'opinion' and 'affiliation' attributes.
    - t: Number of simulation steps to run.
    - alpha, beta, delta: Parameters for polarization dynamics.
    """
    all_nodes = list(G.nodes)
    for _ in range(t):
        neighborhood_update(G, all_nodes, alpha, beta, delta)

def event_sim(G: nx.Graph, alpha=1, beta=1, delta=0) -> None:
    """
    Runs the affective polarization simulation on the graph G.
    Modifies G in-place by updating node 'opinion' attributes.

    Parameters:
    - G: A NetworkX graph where each node has 'opinion' and 'affiliation' attributes.
    - alpha, beta, delta: Parameters for polarization dynamics.
    """
    while True:
        polarizing_nodes = polarizing_neighborhoods(G, alpha, beta, delta)
        if not polarizing_nodes:
            break
        neighborhood_update(G, polarizing_nodes, alpha, beta, delta)

def record_events(G, alpha=1, beta=1, delta=0):
    """
    Runs the affective polarization (AP) simulation on the given graph G.

    Returns:
    pd.DataFrame: "history" where:
                  - Rows (index) represent time steps (including initial conditions).
                  - Columns include:
                      - nodes that could change per step
                      - count of opinion=1 nodes in each partition
    """
    # Get partitions and initialize storage
    partitions = nx.get_node_attributes(G, "affiliation")
    unique_partitions = sorted(set(partitions.values()))

    # Initialize empty DataFrame
    history = pd.DataFrame(columns=["possible_changes"] + [f"partition_{p}" for p in unique_partitions])

    step = 0
    while True:  # Loop until no more changes are possible
        # Retrieve current state
        opinions = nx.get_node_attributes(G, "opinion")
        partition_counts = {p: 0 for p in unique_partitions}
        for node, part in partitions.items():
            if opinions[node] == 1:
                partition_counts[part] += 1

        # Track number of possible changes at this step
        polarizing_nodes = polarizing_neighborhoods(G, alpha, beta, delta)
        possible_changes = len(polarizing_nodes)

        # Create a new row as a dictionary
        row = {"possible_changes": possible_changes}
        row.update({f"partition_{p}": partition_counts[p] for p in unique_partitions})

        # Append the new row to the DataFrame using pd.concat() (more efficient than df.append())
        history = pd.concat([history, pd.DataFrame([row])], ignore_index=True)

        # If no possible changes, stop the simulation
        if not polarizing_nodes:
            break

        # Apply opinion updates to the graph
        neighborhood_update(G, polarizing_nodes, alpha, beta, delta)

        step += 1

    # Set time steps as the index
    history.index.name = "time_step"

    return history  # Returning a DataFrame that was constructed incrementally

# ------------------------------------------------------

def record_stabilization_trials(n, m, k, partition_ratios, opinion_ratios, elite_opinions, alpha, beta, delta):
    """
    Run opinion dynamics simulations on ER and BA graphs with both random and centrality-based opinion assignments.
    Store summary metrics in a SQLite database named 'Stabilization Trials.db'.

    Parameters
    ----------
    n : int
        Number of nodes in the graph.
    m : int
        Parameter for BA graph (used only in gamma estimation).
    k : int
        Number of trials per graph/opinion variant.
    partition_ratios : list of float
        Relative sizes of the two partitions (e.g. [1, 1]).
    opinion_ratios : list of float
        Initial proportion of opinion-1 nodes per partition (e.g. [0.6, 0.4]).
    elite_opinions : list of int
        Elite opinion value per partition for centrality-based assignment (e.g. [1, 1]).
    alpha, beta, delta : float
        Parameters controlling opinion change dynamics.
    """

    assert len(partition_ratios) == 2
    assert len(opinion_ratios) == 2
    assert len(elite_opinions) == 2

    # --- Setup database ---
    db_path = "Stabilization Trials.db"
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute("""
    CREATE TABLE IF NOT EXISTS trial_results (
        trial_id INTEGER PRIMARY KEY AUTOINCREMENT,
        graph_type TEXT,
        opinion_type TEXT,
        "order" INTEGER,
        size INTEGER,
        imbalance_ratio REAL,
        opinion_ratio1 REAL,
        opinion_ratio2 REAL,
        alpha REAL,
        beta REAL,
        delta REAL,
        timesteps INTEGER,
        C0 INTEGER,
        rho REAL,
        gamma REAL
    )
    """)
    conn.commit()

    graph_types = ["ER", "BA"]
    opinion_types = ["random", "centrality"]
    p_ratios = np.array(partition_ratios) / sum(partition_ratios)
    imbalance_ratio = max(p_ratios) / min(p_ratios)

    for graph_type in graph_types:
        if graph_type == "ER":
            G, _ = generate_networks(n)
        else:
            _, G = generate_networks(n, m)

        subgraphs = random_subgraphs(G, p_ratios.tolist())

        for i, sg in enumerate(subgraphs):
            for node in sg.nodes:
                G.nodes[node]["affiliation"] = i

        opinion_variants = {}
        G_random = G.copy()
        G_random = random_opinions(G_random, subgraphs, opinion_ratios)
        opinion_variants["random"] = G_random

        G_central = G.copy()
        G_central = centrality_opinions(G_central, subgraphs, elite_opinions, opinion_ratios)
        opinion_variants["centrality"] = G_central

        for opinion_type, graph_variant in opinion_variants.items():
            for _ in range(k):
                G_sim = graph_variant.copy()

                gamma = estimate_gamma(G_sim, m=m)
                C0 = len(polarizing_neighborhoods(G_sim, alpha, beta, delta))

                timesteps = 0
                while True:
                    changeable = polarizing_neighborhoods(G_sim, alpha, beta, delta)
                    if not changeable:
                        break
                    neighborhood_update(G_sim, changeable, alpha, beta, delta)
                    timesteps += 1

                rho = (C0 / timesteps) if timesteps > 0 else float("nan")

                cursor.execute("""
                INSERT INTO trial_results (
                    graph_type, opinion_type, "order", size,
                    imbalance_ratio,
                    opinion_ratio1, opinion_ratio2,
                    alpha, beta, delta,
                    timesteps, C0, rho, gamma
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """, (
                    graph_type, opinion_type,
                    G_sim.number_of_nodes(),
                    G_sim.number_of_edges(),
                    imbalance_ratio,
                    opinion_ratios[0], opinion_ratios[1],
                    alpha, beta, delta,
                    timesteps, C0, rho, gamma
                ))

    conn.commit()
    conn.close()

# --------------------------------------------

def analyze_trial_correlations(graph_type=None, opinion_type=None, **filters):
    """
    Compute Pearson correlation p-values between continuous input variables and outcomes.
    Filters on graph_type and opinion_type before analysis.
    
    Parameters:
        graph_type: 'ER' or 'BA' (optional)
        opinion_type: 'random' or 'centrality' (optional)
        **filters: any additional filters to apply to trial_results
        
    Returns:
        DataFrame of p-values: rows = input variables, columns = outcome variables
    """

    # --- Load data ---
    conn = sqlite3.connect("Stabilization Trials.db")
    df = pd.read_sql_query("SELECT * FROM trial_results", conn)
    conn.close()

    # --- Apply filters ---
    if graph_type is not None:
        df = df[df["graph_type"] == graph_type]
    if opinion_type is not None:
        df = df[df["opinion_type"] == opinion_type]

    for key, value in filters.items():
        if key not in df.columns:
            raise ValueError(f"Filter key '{key}' not found in trial_results columns.")
        df = df[df[key] == value]

    # --- Add derived column: S0 = n - C0 ---
    df["S0"] = df["order"] - df["C0"]

    # --- Define variables ---
    outcome_vars = ["rho", "S0", "timesteps"]
    input_vars = [
        "order", "size", "imbalance_ratio", "opinion_ratio1", "opinion_ratio2",
        "alpha", "beta", "delta", "gamma"
    ]

    results = pd.DataFrame(index=input_vars, columns=outcome_vars)

    for x in input_vars:
        for y in outcome_vars:
            if x == y:
                results.loc[x, y] = None
                continue

            valid = df[[x, y]].dropna()
            if len(valid) < 5:
                results.loc[x, y] = None
                continue

            try:
                with warnings.catch_warnings():
                    warnings.filterwarnings("error", category=ConstantInputWarning)
                    _, pval = pearsonr(valid[x], valid[y])
                    results.loc[x, y] = pval
            except ConstantInputWarning:
                results.loc[x, y] = None
            except Exception:
                results.loc[x, y] = None

    return results.astype(float)

# ------------------------------------------------

def record_polarization_trials(n, m, k, partition_ratios, opinion_ratios, elite_opinions, alpha, beta, delta, t):
    """
    Run opinion dynamics simulations on ER and BA graphs with both random and centrality-based opinion assignments.
    Store polarization metrics in a SQLite database named 'Polarization Trials.db'.
    """

    assert len(partition_ratios) == 2
    assert len(opinion_ratios) == 2
    assert len(elite_opinions) == 2

    db_path = "Polarization Trials.db"
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute("""
    CREATE TABLE IF NOT EXISTS trial_results (
        trial_id INTEGER PRIMARY KEY AUTOINCREMENT,
        graph_type TEXT,
        opinion_type TEXT,
        "order" INTEGER,
        size INTEGER,
        imbalance_ratio REAL,
        opinion_ratio1 REAL,
        opinion_ratio2 REAL,
        alpha REAL,
        beta REAL,
        delta REAL,
        gamma REAL,
        initial_intergroup_diff REAL,
        final_opinion_ratio1 REAL,
        final_opinion_ratio2 REAL,
        final_intergroup_diff REAL,
        change_from_initial_ratio1 REAL,
        change_from_initial_ratio2 REAL,
        change_in_intergroup_diff REAL
    )
    """)
    conn.commit()

    graph_types = ["ER", "BA"]
    p_ratios = np.array(partition_ratios) / sum(partition_ratios)
    imbalance_ratio = max(p_ratios) / min(p_ratios)
    initial_intergroup_diff = abs(opinion_ratios[0] - opinion_ratios[1])

    for graph_type in graph_types:
        if graph_type == "ER":
            G, _ = generate_networks(n)
        else:
            _, G = generate_networks(n, m)

        subgraphs = random_subgraphs(G, p_ratios.tolist())

        for i, sg in enumerate(subgraphs):
            for node in sg.nodes:
                G.nodes[node]["affiliation"] = i

        opinion_variants = {
            "random": random_opinions(G.copy(), subgraphs, opinion_ratios),
            "centrality": centrality_opinions(G.copy(), subgraphs, elite_opinions, opinion_ratios)
        }

        for opinion_type, graph_variant in opinion_variants.items():
            for _ in range(k):  # Removed tqdm
                G_sim = graph_variant.copy()
                gamma = estimate_gamma(G_sim, m=m)
                time_sim(G_sim, t, alpha, beta, delta)

                affiliations = nx.get_node_attributes(G_sim, "affiliation")
                opinions = nx.get_node_attributes(G_sim, "opinion")

                group_0 = [n for n in G_sim if affiliations[n] == 0]
                group_1 = [n for n in G_sim if affiliations[n] == 1]

                final_ratio1 = sum(opinions[n] == 1 for n in group_0) / len(group_0)
                final_ratio2 = sum(opinions[n] == 1 for n in group_1) / len(group_1)
                final_intergroup_diff = abs(final_ratio1 - final_ratio2)

                change_from_initial_ratio1 = abs(opinion_ratios[0] - final_ratio1)
                change_from_initial_ratio2 = abs(opinion_ratios[1] - final_ratio2)
                change_in_intergroup_diff = abs(initial_intergroup_diff - final_intergroup_diff)

                cursor.execute("""
                INSERT INTO trial_results (
                    graph_type, opinion_type, "order", size,
                    imbalance_ratio,
                    opinion_ratio1, opinion_ratio2,
                    alpha, beta, delta, gamma,
                    initial_intergroup_diff,
                    final_opinion_ratio1, final_opinion_ratio2,
                    final_intergroup_diff,
                    change_from_initial_ratio1,
                    change_from_initial_ratio2,
                    change_in_intergroup_diff
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """, (
                    graph_type, opinion_type,
                    G_sim.number_of_nodes(),
                    G_sim.number_of_edges(),
                    imbalance_ratio,
                    opinion_ratios[0], opinion_ratios[1],
                    alpha, beta, delta, gamma,
                    initial_intergroup_diff,
                    final_ratio1, final_ratio2,
                    final_intergroup_diff,
                    change_from_initial_ratio1,
                    change_from_initial_ratio2,
                    change_in_intergroup_diff
                ))

    conn.commit()
    conn.close()

# ----------------------------------------------------

def analyze_polarization_correlations(**filters):
    """
    Compute p-values and coefficients for correlations between input variables and polarization outcomes
    using linear regression. Accepts filter kwargs to subset trials.

    Returns
    -------
    pvals : pandas.DataFrame
        Rows = input variables, Columns = outcome variables; entries = p-values.
    coefs : pandas.DataFrame
        Same shape as pvals; entries = estimated coefficients (sign + magnitude).
    """

    # --- Load data ---
    conn = sqlite3.connect("Polarization Trials.db")
    df = pd.read_sql_query("SELECT * FROM trial_results", conn)
    conn.close()

    # --- Apply filters ---
    for key, value in filters.items():
        if key not in df.columns:
            raise ValueError(f"Filter key '{key}' not found in trial_results columns.")
        df = df[df[key] == value]

    # --- Define variables specific to polarization trials ---
    outcome_vars = [
        "final_opinion_ratio1",
        "final_opinion_ratio2",
        "final_intergroup_diff",
        "change_from_initial_ratio1",
        "change_from_initial_ratio2",
        "change_in_intergroup_diff"
    ]

    all_input_vars = [
        "graph_type", "opinion_type", "order", "size", "imbalance_ratio",
        "opinion_ratio1", "opinion_ratio2", "alpha", "beta", "delta", "gamma",
        "initial_intergroup_diff"
    ]

    categorical_vars = {"graph_type", "opinion_type"}

    # --- Prepare result tables ---
    pvals = pd.DataFrame(index=all_input_vars, columns=outcome_vars, dtype=float)
    coefs = pd.DataFrame(index=all_input_vars, columns=outcome_vars, dtype=float)

    # --- Run regression for each outcome/predictor pair ---
    for y in outcome_vars:
        for x in all_input_vars:
            sub = df[[x, y]].dropna()
            if sub.shape[0] < 5:
                pvals.loc[x, y] = None
                coefs.loc[x, y] = None
                continue

            # Build formula
            if x in categorical_vars:
                formula = f"{y} ~ C({x})"
            else:
                formula = f"{y} ~ {x}"

            # Fit OLS model
            model = smf.ols(formula=formula, data=sub).fit()

            # Extract predictor’s p-value and coefficient
            pvals.loc[x, y] = model.pvalues.iloc[1]
            coefs.loc[x, y] = model.params.iloc[1]

    return pvals.astype(float), coefs.astype(float)

##################
### VALIDATION ###
##################

def facebook_graph(npz_file):
    """
    Generates an unweighted, undirected NetworkX graph from a sparse adjacency matrix (.npz file).

    Parameters:
        npz_file (str): Path to the .npz file containing the sparse upper-triangular adjacency matrix.

    Returns:
        G (networkx.Graph): A NetworkX graph object.
    """
    # Load the stored sparse upper triangular adjacency matrix
    A_triu = sp.load_npz(npz_file)

    # Symmetrize the adjacency matrix efficiently to restore the full undirected graph
    A_full = A_triu + A_triu.T

    # Convert the sparse adjacency matrix into an unweighted NetworkX graph
    return nx.from_scipy_sparse_array(A_full, create_using=nx.Graph())

#####################
### VISUALIZATION ###
#####################

def generate_colors(num_groups):
    """
    Generate a list of colors based on the number of groups.

    Parameters
    ----------
    num_groups : int
        Number of groups (partitions, subgraphs, etc.)

    Returns
    -------
    colors : list of str
        List of hex color codes, one for each group.
    """
    # Predefined good colors
    base_colors = [
        "#cd5c5c",  # muted red
        "#4682b4",  # steel blue
        "#6ba384",  # soft green
        "#cfcf55"   # muted yellow
    ]

    colors = base_colors.copy()

    # Add muted random colors if needed
    if num_groups > len(colors):
        np.random.seed(42)  # Optional for consistent results
        for _ in range(num_groups - len(colors)):
            random_color = "#{:06x}".format(np.random.randint(0x777777, 0xCCCCCC))
            colors.append(random_color)

    return colors[:num_groups]

# --------------------------------

def gamma_histogram(m=None, n=None, generate=False):
    """
    Analyze and optionally generate gamma values stored in the database.

    Parameters
    ----------
    m : int, tuple, or None
        If int: exact m. 
        If tuple: ('<', int) or ('<=', int).
        If None: ignore m in filtering.
    n : int, tuple, or None
        Same logic as m, but for n (number of nodes).
    generate : bool
        If True, generate new gamma values according to provided m and/or n.
    """
    
    # --- Helper to parse m/n filters ---
    def parse_filter(value):
        if value is None:
            return None
        elif isinstance(value, tuple) and value[0] in ('<', '<='):
            return value
        elif isinstance(value, int):
            return ('=', value)
        else:
            raise ValueError("Only '=', '<', '<=' allowed. Must be int or ('<', int), ('<=', int).")

    m_filter = parse_filter(m)
    n_filter = parse_filter(n)
    
    # --- Generate if appropriate ---
    if generate and m_filter:
    
        m_val = m_filter[1]
    
        # Open SQLite connection
        conn = sqlite3.connect("Scale Free Test.db")
        cursor = conn.cursor()
    
        # Ensure table exists
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS gamma_results (
            m INTEGER,
            n INTEGER,
            gamma REAL
        )
        """)
    
        # Find candidate n values for this m
        n_min = 2 * m_val
        n_max = int(np.floor(np.exp(2 * m_val)))
        n_candidates = range(n_min, n_max + 1)
    
        # Apply n filter if given
        if n_filter:
            if n_filter[0] == '<':
                n_candidates = [n_val for n_val in n_candidates if n_val < n_filter[1]]
            elif n_filter[0] == '<=':
                n_candidates = [n_val for n_val in n_candidates if n_val <= n_filter[1]]
            elif n_filter[0] == '=':
                n_candidates = [n_val for n_val in n_candidates if n_val == n_filter[1]]
    
        for n_val in tqdm(n_candidates, desc=f"Generating for m={m_val}"):
            G = nx.barabasi_albert_graph(n_val, m_val)
            gamma = estimate_gamma(G, m_val)
    
            cursor.execute(
                "INSERT INTO gamma_results (m, n, gamma) VALUES (?, ?, ?)",
                (m_val, n_val, gamma)
            )
    
        conn.commit()
        conn.close()

    # --- Connect to database to read ---
    conn = sqlite3.connect("Scale Free Test.db")
    cursor = conn.cursor()

    # --- Build WHERE clause dynamically ---
    where_clauses = []
    params = []

    if m_filter:
        where_clauses.append(f"m {m_filter[0]} ?")
        params.append(m_filter[1])
    if n_filter:
        where_clauses.append(f"n {n_filter[0]} ?")
        params.append(n_filter[1])

    # Finalize query
    if where_clauses:
        where_clause = " AND ".join(where_clauses)
        query = f"SELECT m, n, gamma FROM gamma_results WHERE {where_clause}"
    else:
        query = "SELECT m, n, gamma FROM gamma_results"

    cursor.execute(query, params)
    rows = cursor.fetchall()
    conn.close()

    if not rows:
        print("No gamma values found for specified m and n.")
        return

    # --- Organize data ---
    gamma_list = []
    for _, _, gamma in rows:
        if gamma is not None and not np.isnan(gamma):
            gamma_list.append(gamma)

    gamma_array = np.array(gamma_list)

    # --- Plot histogram ---
    plt.hist(gamma_array, bins=50, color='steelblue', edgecolor='black')
    plt.axvline(2.3, color='olive', linestyle='--', label='2.3')
    plt.axvline(3.0, color='indianred', linestyle='--', label='3.0')
    plt.title("Histogram of Fitted γ Values")
    plt.xlabel("γ (scaling exponent)")
    plt.ylabel("Frequency")
    plt.legend()
    plt.tight_layout()
    plt.show()

    # --- Print percentage summary ---
    mask = (gamma_array >= 2.3) & (gamma_array <= 3.0)
    pct = 100 * np.sum(mask) / len(gamma_array)
    print(f"Percentage of gamma values between 2.3 and 3.0: {pct:.2f}% ({np.sum(mask)}/{len(gamma_array)})")

# -------------------------------------------------------

def subgraph_plot(G, subgraphs):
    """
    Plot a full graph and its given subgraphs side-by-side.

    Parameters
    ----------
    G : networkx.Graph
        The original full graph.
    subgraphs : list of networkx.Graph
        List of already-partitioned subgraphs.
    """

    # --- Generate colors ---
    colors = generate_colors(len(subgraphs))

    # --- Build node-to-color mapping ---
    node_color_mapping = {}
    for i, sg in enumerate(subgraphs):
        for node in sg.nodes():
            node_color_mapping[node] = colors[i]

    # --- Create subplots ---
    fig, axes = plt.subplots(1, len(subgraphs) + 1, figsize=(6 * (len(subgraphs) + 1), 6))

    # --- Plot full graph first ---
    axes[0].set_title("Full Graph", fontsize=14)
    nx.draw(
        G,
        ax=axes[0],
        node_color=[node_color_mapping.get(node, "gray") for node in G.nodes()],
        edge_color="gray",
        node_size=100
    )

    # --- Plot each subgraph separately ---
    for i, sg in enumerate(subgraphs):
        axes[i+1].set_title(f"Subgraph {i+1}", fontsize=14)
        nx.draw(
            sg,
            ax=axes[i+1],
            node_color=colors[i],
            edge_color="gray",
            node_size=100
        )

    plt.tight_layout()
    plt.show()