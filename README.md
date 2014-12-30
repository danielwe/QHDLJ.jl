# QHDLJ

This package can be used for simulating large-scale quantum optical circuits according to the scheme outlined in [this paper] [1].

## Key facts about QHDLJ

1. It includes routines for defining component models as well as a basic set of $\chi_2$ and $\chi_3$ non-linear models.  
2. It interfaces with my RK4.jl ODE/SDE solver package for performing stochastic simulations that model the effect of quantum shot noise.
3. It provides methods for composing basic or composite components into circuits.
4. It features methods for finding fixpoints
5. It allows for computing linearized models and transfer functions.
6. It has methods for computing the power spectral density of a system linearized at some point. This allows to analyze whether a system produces squeezed output states.
7. It's very fast, comparable to C/C++. The non-linear part of the ODEs are defined symbolically for each component and then concatenated and compiled for a composite circuit.

## Installation

Install [Julia][2], then open a Julia shell and run:

	Pkg.clone("https://github.com/ntezak/RK4.jl.git")
	Pkg.clone("https://bitbucket.org/ntezak/QHDLJ.jl.git")


## How to get started

See the [example notebook](http://nbviewer.ipython.org/urls/bitbucket.org/ntezak/qhdlj.jl/raw/master/examples/QHDLJ%20Example.ipynb) for how to get started.



## Roadmap

1. [QHDL][3] parser: Allow for importing/exporting QHDL files
2. Interface with [QNET][4]
3. Improve/vectorize symbolic composition of non-linear ODE parts.


[1]: http://arxiv.org/abs/1402.5983 "Santori, C., Pelc, J. S., Beausoleil, R. G., Tezak, N., Hamerly, R., & Mabuchi, H. (2014). Quantum noise in large-scale coherent nonlinear photonic circuits."

[2]: http://julialang.org

[3]: http://rsta.royalsocietypublishing.org/content/370/1979/5270 "Tezak, N., Niederberger, A., Pavlichin, D. S., Sarma, G., & Mabuchi, H. (2012). Specification of photonic circuits using quantum hardware description language. Philosophical Transactions. Series A, 370(1979), 5270–90."

[4]: http://mabuchilab.github.io/QNET/
