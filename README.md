# Volume constrained barostats in i-Pi.

Convential barostats are known to induce large volume fluctuations in the system, which can cause unwanted phase transitions in flexible systems. For example MIL-53(Al) has two known phases: an large-pored (LP) phase and a close-pored (CP) phase. If the fluctuation on the cell volume caused by the barostat is sufficiently, a transition from the LP to the CP phase can occur spontaneously. In these situations mechanical properties, like the pressure-vs-volume curve, can not be reliably computed and one is forced to go to much larger unit cells, which is computationally often not feasible. To tackle this problem the (N, V, $\sigma_a = 0$, T)-ensemble is proposed in [1], which allows cell shape fluctuations while keeping the cell volume fixed. MD simulations in this ensemble ensures that the system relaxes under the sole constraint of a constant volume. This enables one to measure the ensemble average of the instantaneous hydrostatic pressure $\langle P_i \rangle$, essentially compmuting the pressure as a function of the volume. If done for different volumes a pressure-vs-volume curve is constructed, which can even be integrated to obtain the free energy as a function of volume.

In the MTTK barostat the cell tensor $\textbf{h}$ is updated using the cell momentum tensor $\textbf{p}_g$:

```math
\dot{\textbf{h}} = \frac{\textbf{p}_g \textbf{h}}{W},
```
with W the barostat mass.
The volume of the unit cell is given by $V = \text{det}\left(\textbf{h}\right)$, from this we can compute the time derivative of the volume as:

```math
\begin{aligned}
  \dot{V} &= \frac{\text{d}}{\text{dt}} \det\left(\mathbf{h}\right) \\
          &= \det\left(\mathbf{h}\right) \text{Tr}\left(\mathbf{h}^{-1} \dot{\mathbf{h}}\right) \\
          &= \det\left(\mathbf{h}\right) \text{Tr}\left(\mathbf{h}^{-1} \frac{\mathbf{p}_g \mathbf{h}}{W}\right) \\
          &= \frac{\det\left(\mathbf{h}\right)}{W} \text{Tr}\left(\mathbf{p}_g\right)
\end{aligned}
```
where we have used the Jacobi's formula for the derivative of a determinant and the cyclic property of traces.
We now find that a constant volume gives:

```math
\dot{V} = 0 \iff \text{Tr}\left(\textbf{p}_g\right) = 0.
```
Thus we find that the cell momentum tensor must be traceless in order for the cell volume to remain fixed.
In code format this simply amounts to inserting the line 
```
p_g -= np.eye(3) * np.trace(self.p) / 3.0
```
before `p_g` is used to update the unit cell tensor. 

Probably the simplest way to implement this is by adding a `vol_constraint` boolean argument to the `Barostat`-class. If `vol_constraint` is set to `True` the above line of code is executed at each update step, ensuring that the volume is kept constant.


# References
[1. Rogge, S. M. J. et al. A Comparison of Barostats for the Mechanical Characterization of Metal–Organic Frameworks. J. Chem. Theory Comput. 11, 5583–5597 (2015).](https://pubs.acs.org/doi/10.1021/acs.jctc.5b00748)
