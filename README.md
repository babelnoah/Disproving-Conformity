# Disproving-Conformity
Code identifies discrepancies in the mathematical model proposed by Smaldino and Epstein in "Social conformity despite individual preferences for distinctiveness" published in Royal Society (used as evidence).

## Relevant paper disproven:

The paper in question can be found through the following link:

[Social conformity despite individual preferences for distinctiveness](https://royalsocietypublishing.org/doi/10.1098/rsos.140437)

This peer-reviewed article was scrutinized through the application of various statistical and linear algebra methods of the original study's analysis. We identify a potential discrepancy within the model that arises under the assumption of stable equilibrium. Our analysis suggests that the positions of conformists and non-conformists tend towards infinity for the distance between them to approach a constant.Contrary to what might be expected, the average positions of conformists and non-conformists always diverge. This continual divergence can be explained by the inherent structure of the preference functions and the adjustment mechanism. We can analyze the stability of this system by calculating the spectral radius of the matrix. A spectral radius equal to 1 indicates a neutral equilibrium where the system does not move once it reaches the equilibrium state. Solving the characteristic equation of the matrix yields eigenvalues λ₁ = 0.6 and λ₂ = 1. Indeed, this aligns expected stabilization of the xC and xN  growth rate, causing the relative distance between them to remain constant, converging to a value of 0.25.

![conformity_rate_2](https://github.com/babelnoah/Disproving-Conformity/assets/114769700/b280fd59-ac12-41ce-bd13-46268d749a6c)

However, the dominant eigenvalue λ₂ = 1 with its corresponding symmetric eigenvector ~ [-0.707, -0.707] prescribes a linear escalation of the average positions of both subpopulations. This leads to a problem of infinite growth, a counter-intuitive and non-physical result. Such an assertion is counter-intuitive and non-physical, as it implies an unbounded linear increase in the positions of the agents in both subpopulations.

![what](https://github.com/babelnoah/Disproving-Conformity/assets/114769700/ff7da68e-0ff4-4682-9a57-7da617020b98)

To unravel the system's behavior computationally, we commenced with a mixed population of conformists and non-conformists. Iterative position updates were executed until 10^6 iterations were reached.

These results have been sent to The Royal Society, and the paper is being considered for retraction. 

## Acknowledgements

Some of the computing for this project was performed on the Sherlock cluster. We would like to thank Stanford University and the Stanford Research Computing Center for providing computational resources and support that contributed to these research results.
