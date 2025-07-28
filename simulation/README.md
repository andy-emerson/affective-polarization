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