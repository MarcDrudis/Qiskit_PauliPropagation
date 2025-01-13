import numpy as np
from qiskit.circuit import Parameter, QuantumCircuit
from qiskit.circuit.library import PauliGate, TwoLocal, XGate, YGate, ZGate
from qiskit.circuit.library.n_local.n_local import get_entangler_map
from qiskit.compiler import transpile
from qiskit.converters import circuit_to_dag, dag_to_circuit
from qiskit.dagcircuit import DAGCircuit
from qiskit.quantum_info import SparsePauliOp, Statevector
from rustworkx.generators import grid_graph
from tqdm import tqdm
from windsurfer.gradient_samplers.QGT_circuits import (QGT_Circuits,
                                                       get_qgt_circuit)
from windsurfer.hamiltonians.graph import entanglement_from_graph

from qiskit_paulipropagation.julia import (pp, pp_estimator, propagation,
                                           qc_to_pp, sparsepauliop_to_pp,
                                           supported_gates)

# from qiskit_paulipropagation.julia import pp_estimator

if __name__ == "__main__":
    graph = grid_graph(4, 4)
    qc = TwoLocal(
        graph.num_nodes(), ["ry"], ["cx"], entanglement_from_graph(graph), reps=1
    ).decompose()
    print(qc)
    # qc = transpile(qc, basis_gates=supported_gates)
    # print(qc)
    observable = SparsePauliOp.from_sparse_list([("Z", [-1], 1)], qc.num_qubits + 1)
    pp_observable = sparsepauliop_to_pp(observable)
    print(pp_observable)
    parameters = np.random.uniform(-1, 1, qc.num_parameters).tolist()

    max_weight = 13
    min_abs_coeff = 1e-9

    qgt_circuits = QGT_Circuits(qc)

    qgt_pp_circuits = [
        qc_to_pp(c) for c in tqdm(qgt_circuits.qgt_circuits, "transforming circuits")
    ]

    counter = 0
    for pp_circuit, parameter_map in tqdm(qgt_pp_circuits, "computing obs"):
        # print(parameter_map)
        pauli_sum = propagation(pp_circuit, parameter_map, pp_observable, parameters)
        pp_value = pp.overlapwithzero(pauli_sum)
        qc = qgt_circuits.qgt_circuits[counter]
        qc = dag_to_circuit()
        qiskit_parameters = [parameters[i] for i in parameter_map if isinstance(i, int)]
        print(qc)
        qiskit_value = Statevector(
            qc.assign_parameters(qiskit_parameters)
        ).expectation_value(observable)
        print(pp_value, qiskit_value)
        print(qc == qiskit_value)
        counter += 1
