"""Tests for generator analysis."""

import os
import time
from unittest import TestCase

import numpy as np
from qiskit.circuit import parameter
from qiskit.circuit.library import EfficientSU2, TwoLocal, n_local
from qiskit.circuit.random import random_circuit
from qiskit.compiler import transpile
from qiskit.quantum_info import Statevector, state_fidelity
from qiskit.quantum_info.states.statevector import SparsePauliOp
from qiskit_algorithms.gradients.reverse.reverse_qgt import ReverseQGT

import unittest

from qiskit_paulipropagation.julia import compute_qgt, pp_estimator, supported_gates, pp


class TestHamiltonians(TestCase):
    """Tests if we can get the generators out of any circuit."""

    # np.testing.assert_almost_equal(pp_value, qiskit_value)

    def test_A(self):
        qc = TwoLocal(
            5, ["rx", "rz"], ["cx"], entanglement="pairwise", reps=2
        ).decompose()
        use_test_qgt(qc, max_weight=7, truncation_threshold=0.0001)

    def test_B(self):
        qc = n_local(5, ["ry"], ["cx"], entanglement="pairwise", reps=2)
        use_test_qgt(qc)

    def test_C(self):
        qc = TwoLocal(20, ["ry"], ["cx"], entanglement="pairwise", reps=2).decompose()
        rng = np.random.default_rng(0)
        parameters = rng.uniform(-1, 1, qc.num_parameters).tolist()
        pp_qgt = compute_qgt(qc, parameters, max_weight=5, truncation_threshold=1e-3)


def use_test_qgt(qc, **kwargs):
    rng = np.random.default_rng(0)
    parameters = rng.uniform(-1, 1, qc.num_parameters).tolist()
    expected_qgt = ReverseQGT(False).run([qc], [parameters]).result().qgts[0].real
    pp_qgt = compute_qgt(qc, parameters, **kwargs)
    # compare_qgt(pp_qgt, expected_qgt)
    np.testing.assert_almost_equal(expected_qgt, pp_qgt, 4)


def compare_qgt(a, b):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2)
    ax[0].imshow(a)
    ax[1].imshow(b)
    plt.show()


if __name__ == "__main__":
    unittest.main()
