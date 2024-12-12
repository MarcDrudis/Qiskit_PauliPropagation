from juliacall import Main as jl
from juliacall import Pkg as jlPkg
from juliacall import convert
from qiskit.circuit import QuantumCircuit
from qiskit.circuit.parametervector import ParameterVectorElement
from qiskit.converters import circuit_to_dag
from qiskit.quantum_info import SparsePauliOp

print("Activating PP")
jlPkg.activate("PauliPropagation")
jl.seval("using PauliPropagation")
pp = jl.PauliPropagation


def qc_to_pp(qc: QuantumCircuit) -> tuple[list[tuple[str, list[int]]], list[int]]:
    """Returns a list with all the PP gates in julia and a map with the parameter indices."""
    pauli_rotations = {
        "rx": (convert(jl.Symbol, "X"),),
        "ry": (convert(jl.Symbol, "Y"),),
        "rz": (convert(jl.Symbol, "Z"),),
    }
    clifford_gates = {
        "h": convert(jl.Symbol, "H"),
        "x": convert(jl.Symbol, "X"),
        "y": convert(jl.Symbol, "Y"),
        "z": convert(jl.Symbol, "Y"),
        "s": convert(jl.Symbol, "S"),
        "cx": convert(jl.Symbol, "CNOT"),
        "swap": convert(jl.Symbol, "swap"),
    }
    dag = circuit_to_dag(qc)
    op_nodes = list(dag.topological_op_nodes())
    pp_circuit = []
    parameter_map = []
    for node in op_nodes:
        q_indices = tuple(qarg._index + 1 for qarg in node.qargs)
        name = node.op.name
        if name in pauli_rotations:
            pauli_rot = pp.PauliRotation(pauli_rotations[name], q_indices)
            pp_circuit.append(pp.tofastgates(pauli_rot, qc.num_qubits))
            parameter_map.append(node.op.params[0].index)
        elif name in clifford_gates:
            clifford_gate = pp.CliffordGate(clifford_gates[name], q_indices)
            pp_circuit.append(pp.tofastgates(clifford_gate, qc.num_qubits))
        else:
            print(f"We did not find a gate for {node.op.name}. Skipping Gate.")
    return pp_circuit, parameter_map


def pp_estimator(qc: QuantumCircuit, obs: SparsePauliOp, params: list[float], **kwargs):
    """Retruns the expectation value for an initial state of |0>."""
    pauli_sum = pp_propagation(qc, obs, params, **kwargs)
    return pp.overlapwithzero(pauli_sum)


def sparsepauliop_to_pp(op: SparsePauliOp):
    nqubits = op.num_qubits
    pp_paulisum = pp.PauliSum(nqubits)
    for pauli, qubits, coefficient in op.to_sparse_list():
        pauli_symbols = tuple(convert(jl.Symbol, str(p)) for p in pauli)
        pp_qubits = tuple(q + 1 for q in qubits)
        term = pp.PauliString(
            nqubits,
            pauli_symbols,
            pp_qubits,
            coefficient,
        )
        pp.add_b(pp_paulisum, term)
    return pp_paulisum


def pp_propagation(
    qc: QuantumCircuit, obs: SparsePauliOp, params: list[float], **kwargs
):
    pp_circuit, parameter_map = qc_to_pp(qc)
    pp_params = [params[i] for i in parameter_map]
    pp_observable = sparsepauliop_to_pp(obs)
    pauli_sum = pp.propagate(
        pp_circuit,
        pp_observable,
        pp_params,
        **kwargs,
    )
    return pauli_sum
