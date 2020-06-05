# SpinGlasses_MscThesis

## Spin Glass simulation on small complete graphs with different setups for my MSc Thesis.
**Code comments and Thesis in Hungarian.**

## Summary from the thesis:
In the computer simulations of spin glasses the aim is to create as large systems as possible, which can be brought to thermodinamic equilibrium. These simulations usually need advanced (sometimes specific) hardware and software, and long runtimes. The results of these simulations carry important information about the system, but some details can be lost, such as the past of the system or the importance of individual samples.

In my work, I studied spin glasses on a complete graph (with J_ij = ±1 random coupling) with two main methods: counting the whole phase space and Metropolis – Monte Carlo method. Mainly, I studied smaller systems (N = 20), where the distribution of the correlations is usually not symmetric for individual samples. The symmetry in the correlations occurs when we average them over the samples. (It also occurs in larger systems: N = 128.)

I have run spin glass simulations containing homogenous magnetic field for individual samples, which I later averaged over the samples. The results show that even a smaller magnetic field can change the behavior of the system. Likewise, I studied how one spin’s perturbation affects the whole system. Due to the topology, one spin’s behavior directly affects all the other spins, like a magnetic field. I found that the system responds almost the same for fixing the orientation of a single spin as for the presence of a smaller magnetic field.

In the end, I studied the importance of the degree of frustration. On an individual sample and averaged over frustrated samples also I found that a system with a lower degree of frustration (“frustrated system”) shows weaker response to the presence of a smaller magnetic field. I also found that systems with a negative degree of frustration show strongly different behavior from systems with a positive degree of frustration. It seems that systems with zero degree of frustration behave “regularly”: the degrees of frustration of random sampled systems are mostly zero.
