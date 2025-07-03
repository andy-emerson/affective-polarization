# Introduction

Welcome to the Jupyter notebook for the 2024 **Mathematics Clinic** at Claremont Graduate University. This year, the research topic is **opinion dynamics**. In particular, we will be modeling **affective polarization**, broadening the application of previous research. Affective polarization refers to the tendency of individuals to view their own political group positively and opposing groups negatively, in ways that not only reflect partisan identity but also influence the formation and interpretation of political beliefs.

There have been many different attempts to model affective polarization. Some of the most common have used Bayesian updating. Others model it as a coordination game. Our research builds on the network algorithm model developed by Buddhika Nettasinghe, Ashwin Rao, Bohan Jiang, Allon G. Percus, and Kristina Lerman, in 2025 and expanded upon by Buddhika Nettasinghe, Allon G. Percus, and Kristina Lerman, also in 2025. The authors examine affective polarization, where individuals exhibit **in-group favoritism** (trusting and supporting like-minded individuals) and **out-group animosity** (distrusting and disliking those with opposing views). To analyze the evolution of polarization, they use an affective polarization algorithm on both **complete graphs** and **stochastic block models (SBM)** to simulate the evolution of beliefs in social networks.

SBM is a natural choice for a network model to simulate affective polarization, because affective polarization operates on group identies, and SBM generates networks with well-defined group structures. Complete graphs are a natural choice as a baseline for comparison. By utilizing SBM and complete graphs as the network models on which to apply their polarization algorithm, their paper provides a theoretically grounded yet practically relevant model to understand the persistence and consequences of affective polarization. However, it leaves open the question of how other network structures might influence the evolution of opinions on a network. 

Our study will apply the same polarization algorithm used by Nettasinghe *et al* to the **Erdős-Rényi (ER)** and **Barabási-Albert (BA)** network models to capture different structural characteristics of real-world social networks. The BA model provides us with a scale-free network, where well-connected individuals may have disproportionate influence over the evolution of the opinions both within and between groups. The ER model will serve as our baseline of comparison.

# Generation

The affective polarization algorithm is a network algorithm that modifies node attributes representing individual opinions. We begin, therefore, by generating a network. Specifically, we will be using the Erdős-Rényi model and the Barabási-Albert model. The latter provies us with a **scale-free random graph** to study the effects the scale-free property has on opinion dynamics, the former provides us a **uniform random graph** as a baseline for comparison.

If we wish to compare the influence that structure has on the behavior of our network, the two models should be otherwise similar. This includes the size (the number of edges) and order (the number of nodes), but in our case, it will also be important that both graphs will be connected. The polarization algorithm only applies to connected components. Therefore, disjoint subgraphs would essentially confuse the behavior of multiple small networks as a single large model. Measuring the dynamics of a a disconnected graph is as simple as aggregating the dynamics of its connected components. Therefore, our research will focus on connected graphs, and so it will be important during network generation to generate ER and BA networks that both connected.

The challenge, therefore, becomes, how do we generate both graphs such that they are same size, order and connectivity?

## The Erdős-Rényi Model

The Erdős-Rényi (ER) model was introduced by **Paul Erdös** and **Alfréd Rényi** in 1959. It provides a simple yet powerful framework for studying the properties of networks where edges form randomly between nodes. The ER model comes in two main variants: $G(n, m)$, where a graph is selected uniformly at random from the set of all graphs with $n$ nodes (i.e. order $n$) and exactly $m$ edges (i.e. size $m$), and $G(n, p)$, where each of the $\binom{n}{2}$ possible edges between $n$ nodes is included independently with probability $p$. The latter, $G(n, p)$, is perhaps the more commonly used formulation; however, because we need to maintain precise control of the size of our networks, we will be using the $G(n, m)$ model.

### Minimum Size

ER graphs are not always connected. The number of edges $m$ determines whether the graph is likely to be connected. Erdës and Rényi proved that, as the number of nodes $n$ increases towards infinity, the probability that the graph is connected undergoes a sharp transition at the **connectivity threshold**. Specifically,
$$
\lim_{n \to \infty} \mathbb{P}(G(n, m) \text{ is connected}) =
\begin{cases}
0 & \text{if } m < \displaystyle{\frac{n \ln n}{2}} \\
1 & \text{if } m > \displaystyle{\frac{n \ln n}{2}}
\end{cases}
$$
Therefore, our ER graphs will need to have *at least* $m=\displaystyle{\frac{n \ln n}{2}}$ number of edges in order to serve as a useful baseline for comparison with BA graphs. Remember that this is not a study of ER graphs generally. It is instead a baseline for comparison for the connected components of a BA graph of a given size and order. Therefore, a simple search for connected ER graphs of the same size and order will suffice. 

## The Barabási–Albert Model

The Barabási–Albert (BA) model was introduced by **Albert-László Barabási** and **Réka Albert** in 1999. Like the Erdős–Rényi model, BA networks are generated by specifying the number of nodes in the graph; however, unlike the ER model, which begins with $n$ isolated nodes and adds $m$ edges randomly, BA model begins with only $m$ nodes, and then grows the graph by adding $n-m$ additional nodes, each with $m$ new edges. Therefore, in the BA model, $m$ represents the number of new edges added per new node, rather than the total number of edges in the graph as in the ER model. The initial $m$ nodes are called the **seed** of teh BA graph. When a new node is added to the seed, the probability of attaching to an existing node is proportional to its degree. This process, called **preferential attachment**, produces a connected graph with a **power-law** (or "rich-get-richer") degree distribution:

$$p(k) \sim k^{-\gamma},\quad\mathbb{P}(2.5 \leq \gamma \leq 3) \approx 1$$

The scaling exponent $\gamma$ will be one of the most important measures for determining how network structure influences polarization.

The seed could be any graph on $m$ nodes. It is often a complete graph $K_m$. However, for this project, we use the NetworkX implementation of the BA model, which initializes the graph with $m$ *isolated* nodes. The first new node then adds $m$ edges, one for each of the isolated nodes, forming a connected graph on which preferential attachment can begin. If, for example, $m=3$ and $n>3$, then the initial graph will be 3 isolated nodes. Then the fourth node is added, the degree distribution will be $\{3, 1, 1, 1\}$, immediately giving us a degree distribution that appears to follow a power-law.

The total number of edges for a BA graph generated this way is $|E| = m(n-m)$, where $n$ is the number of nodes. The average degree is therefore

$$
\langle k \rangle = \frac{2|E|}{n} = 2m - \frac{2m^2}{n}.
$$

The average degree asymptotically approaches $2m$ as $n$ grows large:

$$
\lim_{n\to\infty} \left(2m - \frac{2m^2}{n}\right) = 2m.
$$

### Maximum Order

As $n$ increases relative to $m$, our BA graph becomes less dense, approaching zero. Therefore, it must eventually drop below the connectivity threshold for ER graphs. The question becomes: where do these two functions intersect for a given $m$?

$$m(n-m) \leq \frac{n \ln n}{2}$$

For this, let's begin by finding the points of intersection.

Note that the functions intersect at two points; however, our graphs will be significantly larger than the minimum threshold, so we will focus on the larger of the two numbers.

| $m$ | $n_1$ | $n_2$ |
|:---|:-----|:-----|
| 1 | 1.00000 | 4.92155 |
| 2 | 2.64141 | 45.8579 |
| 3 | 3.87469 | 385.001 |
| 4 | 5.00876 | 2948.78 |
| 5 | 6.10425 | 21976.4 |

Curve fitting the larger points of intersection strongly suggests that the maximum order of a BA graph that will have a sufficient number of edges for an ER graph of the same size and order to be connected follows closely (but not exactly) to a simple exponential function.
$$n\leq e^{2m}$$

## Scale-Free Networks

So far, we know that the order $n$ of our network should satisfy $m \leq n \leq e^{2m}$. This allows us to easily find connected Erdös-Renyí graphs of the same size and order of a given Barabási-Albert graph. However, the focus of this research is to study how structural properties (specifically, the scale-free property) of networks influence opinion dynamics.

A graph is said to have the scale-free property when its structure looks statistically similar at different scales of observation. In such graphs, a small number of nodes have very high degrees, acting as hubs, while the vast majority of nodes have relatively few connections. This uneven distribution of connections remains consistent even if one examines only a subset of the graph or zooms out to view the entire network.

The scale-free property is closely associated with degree distributions that follow a power law. Formally, a graph is said to exhibit a power-law degree distribution if the probability $P(k)$ that a randomly chosen node has degree $k$ satisfies $P(k) \sim k^{-\gamma}$ for some constant exponent $\gamma > 1$. This relationship implies that while most nodes have small degrees, a few nodes have very large degrees, and the probability of encountering a node of degree $k$ decays polynomially rather than exponentially. The scale-free property thus reflects the invariance of the degree distribution under rescaling, with the functional form $k^{-\gamma}$ preserved across different levels of observation.

The Barabási–Albert model is commonly used to study scale-free networks because it generates graphs whose degree distributions approximate a power law. As the network grows according to its attachment mechanism, the resulting structure exhibits the key feature of scale-free graphs: a small number of highly connected nodes alongside many nodes with few connections. The Barabási–Albert model thus provides a simple generative framework that reproduces the characteristic statistical self-similarity of scale-free degree distributions observed in many real-world networks.

While the Barabási–Albert model provides a convenient way to generate scale-free networks, our interest is not in the model itself, but in the structural properties of the networks it produces. We use the BA model as a tool to generate graphs that approximate real-world systems, particularly those exhibiting scale-free degree distributions. In particular, we focus on graphs where the degree distribution exponent $\gamma$ falls between 2.3 and 3.0, consistent with the range observed in many empirical studies of social networks. When $n$ is sufficiently large, the mean value of the power-law exponent $\gamma$ for Barabási–Albert graphs approaches 3. However, our $n$ is bounded above by $e^{2m}$ and may be significantly less than that. Therefore, it is important to determine what $\gamma$ values BA graphs will have within the range of $n$ that we will be working with.

So, now that we know the range of values for which it is possible for ER and BA graphs to be of the same size and order and also that the vast majority of ER graphs will be connected while the vast majority of BA graphs (which are always connected) will have $2.3\leq\gamma\leq 3.0$; we can easily generate such graphs.

Our `generate_networks` function takes the number of nodes $n$, calculates the smallest value of $m$ that will generate a number of edges in a BA above teh connectivity threshold for an ER graph, then generate those graphs, keeping only those that are connected (in the case of ER) with a $\gamma$ between 2.3 and 3.0 (in the case of BA). Alternatively, an $m$ can be specified and the boundaries for $\gamma$ can be adjusted.

A power-law distribution appears linear when plotted on logarithmic scales for both the x-axis and y-axis.

$$P(k) \sim k^{-\gamma} \longrightarrow \log P(k) \sim -\gamma \log k$$

As a result, the data points form a straight line on a log-log plot, with a slope equal to $-\gamma$, making it easy to visually identify power-law behavior. In the case of Barabási–Albert graphs, this implies that their degree distributions, when plotted on log-log axes, should display an approximately straight line, confirming the presence of scale-free structure.

# Affiliation

Now that we have our graphs, we need to provide our nodes with group identity. Affective polarization is dependent on in-group and out-group identities, so we need *at least* two subgraphs. Regardless of how many subgraphs are identified, their union should be the entire set of nodes (i.e. all nodes belong to a group) and their intersection should be the empty set (i.e. no nodes belong to more than one group). The question remains, how should we allocate nodes to subgraphs?

Real-world social networks have community structures. These structures have recognizable intergroup and intragroup patterns of connection. Specifically, intragroup connections tend to be dense, which we call clustering, and intergroup connections tend to be sparse, which we call modularity. While in the real-world we often begin with group identity and then measure its clustering and modularity coefficients, we can also discover communities within a population by partitioning the population in ways that either maximize intragroup connections or minimizing intergroup connections.

The network generator most directly related to community detection is the stochastic block model, making it an excellent choice for studying affective polarization. However, our research is primarily focused on *other* network structures, namely the scale-free property. Nevertheless, if we wish to distinguished between the influence of the scale-free property and the influence of community structures, there will need to be community structures within our scale-free network. We accomplished this in the way nodes are assigned to subgraphs.

We use both spectral partitioning, maximizing intergroup connections by sorting nodes according to their Fiedler vector and partitioning the set. The Fiedler vector is associated with the second-smallest eigenvalue on the Laplacian matrix, $L = D - A$, where $A$ is the adjacency matrix and $D$ is the diagonalmatrix of node degrees. We also use modularity maximization, which is also a spectral method using the leading eigenvector of the modularity matrix

$$B_{ij} = A_{ij} - \frac{k_i k_j}{2m}$$

where $k_i$ is the degree of node $i$ and $m$ is the number of edges in $G$. These two method are compared both to each other as well as to a random subgraph assignment as a baseline of comparison. In each case, the sizes of the partitions are maintained, and only random assignment is not deterministic. You can see the visual differences below for each type of graph (ER or BA) and each method of partitioning (random, spectral, modular).

Given that real-world networks tend to have a high level of clustering and modularity simultaneously yet our networks have been generated to model community structure only as a secondary property, it is natural to ask to what extent our two deterministic methods of partition produce similar partitions.

# Simulation

Now that we have constructed our network, we can use a polarization algorithm to simulate opinion dynamics on the network. we are primarily interested in studying two distinct ways in which network stucture influences behavior: first, the time it takes to reach a stable state (i.e. stabilization; second, what that stable state is (i.e. polarization). The best method of implementing the polarization algorithm depends on the question being asked.

In order to understand stabilization time, we must know when the stable state is reached. Therefore, we wll use an event-based simualtion, counting the number of nodes that are capable of changing each iteration of the simulation and focusing on those nodes, which are collectively referred to as teh polarizing neighborhood. The simulation ends when the number of nodes in the polarizing neighborhood reaches zero. In this way, we will know exactly when a given trial reaches a stable state.

However, the event-based simulation is significantly more complex then a time-based simulation, which will terminate after a predetermined number of time steps, regardless of whether or not a stable state has been reached. A time-based simulation is more appropriate for understanding how behavior evolves over time because we can run many more simulations due to the improved processing speed. We can use what we learn about stabilization time to estimate a worst-case scenario number of time steps, which will still be faster in the time-based simulation than the best-case scenario in the event-based simulation.

Nevertheless, because previous research focused on the later, our code will begin with the time-based simulation for validation purposes. Both simulations use the same affective polarization algorithm.

The **affective polarization algorithm** takes a given network along with an $\alpha , \beta , \delta$, which respectively represent in-group love, out-group hate, and threshold parameters. Then it will pick a vertex at random and count the following:

$a = $ the number of neighbors in the same partition with opinion 0.

$b = $ the number of neighbors in the same partition with opinion 1.

$c = $ the number of neighbors in a different partition with opinion 0.

$d = $ the number of neighbors in a different partition with opinion 1.

These counts are then used in the following calculation:
$$\alpha(b-a)-\beta(d-c)$$
If this value is less than $-\delta$, it sets the opinion attribute of that vertex to 0. If it is greater than $\delta$, sets the opinion attribute to 1, and does nothing otherwise. The function `neighborhood_update` simulates opinion dynamics in a social network by selecting a random node from a given list and updating its opinion based on the influence of its neighbors. The update follows a polarization formula that weighs the influence of like-minded and opposing neighbors within and across affiliations. If the computed polarization score exceeds a given threshold ($\delta$), the node's opinion is updated accordingly. The function modifies the graph in place.

## Stabilizations Trials

**Motivation**

We are primarily interested in the system’s **stable state**, which is the point at which no further opinion changes occur. However, the **time-driven simulation** (time-sim) selects nodes at random and may select already stable nodes, making it difficult to determine whether a lack of change indicates convergence or just bad sampling.

To address this, we introduce an **event-driven simulation** (event-sim), which only applies the polarization update rule to nodes that are known to be unstable. This gives us valuable data about the system’s dynamics that we can use to estimate when the time-sim is likely to have stabilized.

**Why Not Use Event-Sim Exclusively?**

While the event-sim is more efficient in terms of convergence (fewer steps to stabilization), it is computationally expensive. Each step requires identifying all unstable nodes, which may involve scanning the entire network.

By contrast, the time-sim is much cheaper per step—especially on large graphs—but offers no built-in signal that stabilization has been reached. The goal is to use data from event-sim to define a **stabilization threshold** $t*$, at which point we can be confident the time-sim has converged.

**Definitions and Assumptions**

We define the following quantities:

- $n$: Total number of nodes in the network  
- $S_t$: Number of stable nodes at time $t$  
- $T$: Number of events required for stabilization in event-sim  
- $t^*$: Expected number of time steps to stabilization in time-sim  
- $p_t$: Probability that the time-sim selects an unstable node at time $t$

In the time-sim, where nodes are selected uniformly at random, the probability of selecting an unstable node is:
$$p_t = 1 - \frac{S_t}{n}$$

**Deriving the Stabilization Threshold**

Although one event may stabilize multiple nodes, we can estimate the average rate of stabilization as:
$$\rho = \frac{n - S_0}{T}$$

Here, $S_0$ is the number of stable nodes at the start of the simulation, and $T$ is the total number of stabilization events in the event-sim. The value $\rho$ captures the average number of nodes stabilized per event.

We introduce an **instability tolerance** $\varepsilon \in (0,1)$, which defines how close we want to be to a fully stable system. For example, $\varepsilon = 0.01$ corresponds to 99% stabilization.

The **stabilization threshold** is then given by:
$$t^* = -\frac{n}{\rho} \ln(\varepsilon)$$

This expression tells us how many time steps the time-sim should run in order to reach the desired level of convergence, based on estimates from the event-sim.

**Derivation**

Assuming stabilization progresses smoothly, we approximate the number of unstable nodes at any point in time-sim as decreasing exponentially with the number of events $E$):
$$U(E) = n - \rho E$$

Then the probability of selecting an unstable node at time $t$, when $E(t)$ events have occurred, is:
$$p_t = \frac{U(E)}{n} = \frac{n - \rho E(t)}{n}$$

Each time step contributes to event progress with probability $p_t$. So the expected change in the number of events is:
$$\frac{dE}{dt} = p_t = \frac{n - \rho E}{n}$$

This is a linear ordinary differential equation. We solve:
$$\frac{dE}{dt} = \frac{n - \rho E}{n}$$

Separate variables:
$$\frac{dE}{n - \rho E} = \frac{dt}{n}$$

Integrate both sides:
$$-\frac{1}{\rho} \ln(n - \rho E) = \frac{t}{n} + C$$

Apply the initial condition \( E(0) = 0 \):
$$C = -\frac{1}{\rho} \ln(n)$$

Substitute and simplify:
$$-\frac{1}{\rho} \ln(n - \rho E) + \frac{1}{\rho} \ln(n) = \frac{t}{n}$$
$$\ln\left( \frac{n}{n - \rho E} \right) = \frac{\rho t}{n}$$

Solving for $E(t)$:
$$E(t) = \frac{n}{\rho} \left(1 - e^{-\rho t / n} \right)$$

**Step 4: Solve for $t^*$**

We define a stabilization threshold by asking: how long until a fraction $1 - \varepsilon$ of the total stabilization events $T$ has occurred?**

So set:
$$E(t^*) = (1 - \varepsilon) T$$

Substitute into the expression for $E(t)$:
$$(1 - \varepsilon) T = \frac{n}{\rho} \left(1 - e^{-\rho t^* / n} \right)$$

Solve for $t^*$:
$$1 - \frac{(1 - \varepsilon) T \rho}{n} = e^{-\rho t^* / n}$$
$$t^* = -\frac{n}{\rho} \ln\left(1 - \frac{(1 - \varepsilon) T \rho}{n} \right)$$

Now plug in $\rho = \frac{n - S_0}{T}$. After simplification, the threshold becomes:
$$t^* = -\frac{n}{\rho} \ln(\varepsilon)$$

This is the expected number of time steps in the time-driven simulation needed to reach a state where only an $\varepsilon$-fraction of events remain.
