import numpy as np
from qiskit.circuit.library import TwoLocal
from qiskit.quantum_info import SparsePauliOp, Statevector

from qiskit_paulipropagation import pp_estimator

if __name__ == "__main__":
    qc = TwoLocal(10, ["ry", "rz"], ["cx"]).decompose()
    observable = SparsePauliOp.from_sparse_list([("ZZ", [4, 5], 1.4)], qc.num_qubits)
    parameters = np.random.uniform(-1, 1, qc.num_parameters).tolist()

    max_weight = 13
    min_abs_coeff = 1e-9

    print(
        "PP result",
        pp_estimator(
            qc,
            observable,
            parameters,
            max_weight=max_weight,
            min_abs_coeff=min_abs_coeff,
        ),
    )
    print(
        "Qiskit Result",
        Statevector(qc.assign_parameters(parameters)).expectation_value(observable),
    )
