import pathlib

import matplotlib.pyplot as plt
import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit.library import TwoLocal
from qiskit.primitives import StatevectorEstimator
from qiskit.quantum_info import SparsePauliOp
from qiskit_algorithms.gradients import LinCombQGT
from rustworkx.generators import grid_graph
from windsurfer.hamiltonians.graph import entanglement_from_graph

directory = pathlib.Path(__file__).parent.resolve()


class PPEstimator:
    def run(
        qc: list[QuantumCircuit],
        obs: list[SparsePauliOp],
        parameter_values: list[list[float]],
        **options
    ):
        qc_to_pp(c)


graph = grid_graph(4, 4)
qc = TwoLocal(
    graph.num_nodes(), ["rz"], ["cx"], entanglement_from_graph(graph), reps=1
).decompose()

rng = np.random.default_rng(0)
parameter_values = rng.uniform(-1, 1, qc.num_parameters)


qgt = LinCombQGT(StatevectorEstimator()).run([qc], [parameter_values]).result().qgts[0]
print(qgt)

# fig, axs = plt.subplots(2, 1, figsize=(12, 5))
#
# axs[0].imshow(qgt.real)
#
# plt.show()
