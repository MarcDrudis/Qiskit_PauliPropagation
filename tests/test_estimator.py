"""Tests for generator analysis."""

from unittest import TestCase

import numpy as np
from qiskit.circuit import parameter
from qiskit.circuit.library import EfficientSU2
from qiskit.circuit.random import random_circuit
from qiskit.compiler import transpile
from qiskit.quantum_info import Statevector, state_fidelity
from qiskit.quantum_info.states.statevector import SparsePauliOp

from qiskit_paulipropagation.julia import pp_estimator, supported_gates


class TestHamiltonians(TestCase):
    """Tests if we can get the generators out of any circuit."""

    def test_random(self):
        """Test if the resulting circuit is indeed equivalent."""

        for i in range(25):
            with self.subTest(f"Random:{i}"):
                qc = random_circuit(10, 5, seed=i)
                obs = SparsePauliOp.from_sparse_list(
                    [("X", [5], 2.0), ("ZZ", [4, 7], -1.0)], qc.num_qubits
                )
                qc = transpile(qc, basis_gates=supported_gates)
                pp_value = pp_estimator(qc, obs, params=[])
                qiskit_value = Statevector(qc).expectation_value(obs)
                # print(qc)
                # print(qiskit_value, pp_value)
                np.testing.assert_almost_equal(pp_value, qiskit_value)

    # def test_efficientSU2(self):
    #     """Test if the resulting circuit is indeed equivalent."""
    #     qc = EfficientSU2(10, entanglement="pairwise").decompose()
    #     obs = SparsePauliOp.from_sparse_list(
    #         [("X", [5], 2.0), ("ZZ", [4, 7], -1.0)], qc.num_qubits
    #     )
    #
    #     for i in range(10):
    #         with self.subTest(f"Random:{i}"):
    #             parameter_values = np.random.uniform(-1, 1, qc.num_parameters)
    #             pp_value = pp_estimator(
    #                 qc,
    #                 obs,
    #                 params=parameter_values,
    #                 # compile=True,
    #             )
    #             qiskit_value = Statevector(
    #                 qc.assign_parameters(parameter_values)
    #             ).expectation_value(obs)
    #             np.testing.assert_almost_equal(pp_value, qiskit_value)
