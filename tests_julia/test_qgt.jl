using PauliPropagation
include("../julia_functions/qgt_compute.jl")
using Random

nl = 4
nq = 3
topology = bricklayertopology(nq; periodic=false)
circuit = hardwareefficientcircuit(nq, nl; topology=topology)

nparams = countparameters(circuit)
Random.seed!(42)
thetas = randn(nparams);

print(qgt_element(nq, circuit, thetas, 5, 3, max_weight=5))

