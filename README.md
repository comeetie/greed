# greed : Bayesian greedy clustering 

Greed enable model based clustering of networks, counts data matrix and much more with different type of generative models. Model selection and clustering is performed in combination by optimizing the Integrated Classification Likelihood

(which is equivalent to minimizing the description length). 

Their are four models availables : 

- sbm : Stochastick Block Models (directed), 
- dcsbm : degree corected Block Models (directed), 
- mm: Mixture of Multinomials, 
- mreg : Mixture of regression. 

With the Integrated Classification Likelihood the parameters of the models are integrated out leaving $p(x,z;\alpha)$ whiwh can be optimized in $z$. The optimization is performed in this thanks to a combination of greedy local search and a genetic algorithm, several optimization algorithms are available.

Since the Integrated Classification Likelihood introduces a natural regularisation for complex models such strategie automaticaly find a good value for the number of cluster, the user needs only to provide an initial guess.

Eventually, the whole path of solutions from K* to 1 cluster is extracted this enable a partial ordering of the clusters, and the evaluation of simpler clustering. The package also provides some ploting functionality.

