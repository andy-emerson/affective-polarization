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
