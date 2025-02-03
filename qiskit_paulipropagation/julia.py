from juliacall import Main as jl
from multipledispatch import dispatch
import numpy as np
from juliacall import Pkg as jlPkg
from juliacall import convert
from qiskit.circuit import QuantumCircuit
from qiskit.compiler import transpile
from qiskit.converters import circuit_to_dag
from qiskit.dagcircuit import DAGCircuit
from qiskit.quantum_info import SparsePauliOp

print("Activating PP")
jlPkg.activate("PauliPropagation")
jl.seval("using PauliPropagation")
pp = jl.PauliPropagation
pp.include("julia_functions/qgt_compute.jl")

# Here is the mapping between the supported qiskit gates and the corresponding PP gates.
pauli_rotations = {
    "rx": (convert(jl.Symbol, "X"),),
    "ry": (convert(jl.Symbol, "Y"),),
    "rz": (convert(jl.Symbol, "Z"),),
    "rxx": (
        convert(jl.Symbol, "X"),
        convert(jl.Symbol, "X"),
    ),
    "ryy": (
        convert(jl.Symbol, "Y"),
        convert(jl.Symbol, "Y"),
    ),
    "rzz": (
        convert(jl.Symbol, "Z"),
        convert(jl.Symbol, "Z"),
    ),
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

supported_gates = list(clifford_gates.keys()) + list(pauli_rotations.keys())


@dispatch(QuantumCircuit, bool)
def qc_to_pp(
    qc: DAGCircuit, gate_position: bool = False
) -> (
    tuple[list[tuple[str, list[int]]], list[int]]
    | tuple[list[tuple[str, list[int]]], list[int], dict[int, int]]
):
    return qc_to_pp(circuit_to_dag(qc), gate_position)


@dispatch(DAGCircuit, bool)
def qc_to_pp(
    qc: DAGCircuit, gate_position: bool = False
) -> (
    tuple[list[tuple[str, list[int]]], list[int]]
    | tuple[list[tuple[str, list[int]]], list[int], dict[int, int]]
):
    """
    Returns a list of Gates which describes the circuit in the PP package. It also returns a mapping
    between the parameter circuit in each library.
    """

    # The circuit must only contain supported gates.
    op_nodes = list(qc.topological_op_nodes())
    pp_circuit = pp.seval("Vector{Gate}")()
    parameter_map = []
    parameter_position = pp.seval("Dict{Int, Int}")()
    for position, node in enumerate(op_nodes):
        q_indices = tuple(qc.find_bit(qarg).index + 1 for qarg in node.qargs)
        name = node.op.name
        if name in pauli_rotations:
            pauli_rot = pp.PauliRotation(pauli_rotations[name], q_indices)
            pp.push_b(pp_circuit, pauli_rot)
            if isinstance(node.op.params[0], float):
                parameter_map.append(node.op.params[0])
            else:
                parameter_map.append(node.op.params[0].index)
                parameter_position[node.op.params[0].index + 1] = position + 1
        elif name in clifford_gates:
            clifford_gate = pp.CliffordGate(clifford_gates[name], q_indices)
            pp.push_b(pp_circuit, clifford_gate)
        else:
            print(f"We did not find a gate for {node.op.name}. Skipping Gate.")
    if gate_position:
        return pp_circuit, parameter_map, parameter_position

    return pp_circuit, parameter_map


def sparsepauliop_to_pp(op: SparsePauliOp):
    """Returns the PP PauliSum representation of the SparsePauliOp."""
    nqubits = op.num_qubits
    pp_paulisum = pp.PauliSum(nqubits)
    for pauli, qubits, coefficient in op.to_sparse_list():
        pauli_symbols = pp.seval("Vector{Symbol}")()
        for p in pauli:
            jl.push_b(pauli_symbols, convert(jl.Symbol, str(p)))
        pp_qubits = pp.seval("Vector{Int}")()
        for q in qubits:
            jl.push_b(pp_qubits, q + 1)
        # Here I am assuming the coefficient will be real. I could be wrong.
        pp.add_b(pp_paulisum, pauli_symbols, pp_qubits, coefficient.real)
    return pp_paulisum


def pp_estimator(
    qc: QuantumCircuit | DAGCircuit,
    obs: SparsePauliOp,
    params: list[float],
    compile: bool = False,
    **kwargs,
) -> float:
    """Retruns the expectation value for an observable and a
    QuantumCircuit starting at |0>."""
    pauli_sum = pp_propagation(qc, obs, params, compile, **kwargs)
    return pp.overlapwithzero(pauli_sum)


def propagation(
    pp_circuit, parameter_map, pp_observable, params: list[float], **kwargs
):
    """
    Propagates a PP circuit. This method can be used to reduce the overhead of converting circuits from
    Qiskit to PP repeatidly.
    """
    pp_params = [params[i] if isinstance(i, int) else i for i in parameter_map]
    pauli_sum = pp.propagate(
        pp_circuit,
        pp_observable,
        pp_params,
        **kwargs,
    )
    return pauli_sum


def pp_propagation(
    qc: QuantumCircuit | DAGCircuit,
    obs: SparsePauliOp,
    params: list[float],
    compile: bool = False,
    **kwargs,
):
    """
    Returns the PP PauliSum of an observable propagated over a quantum circuit.
    """
    if compile:
        qc = transpile(qc, basis_gates=supported_gates)
    pp_circuit, parameter_map = qc_to_pp(qc)
    pp_observable = sparsepauliop_to_pp(obs)
    return propagation(pp_circuit, parameter_map, pp_observable, params, **kwargs)


def compute_qgt(qc: QuantumCircuit, parameters: list[float], **kwargs):
    pp_circuit, parameter_map, parameter_position = qc_to_pp(qc, True)
    pp_params = [parameters[i] if isinstance(i, int) else i for i in parameter_map]
    print(parameter_map)
    qgt = pp.compute_qgt(
        qc.num_parameters,
        pp_circuit,
        pp_params,
        parameter_map,
        parameter_position,
        **kwargs,
    )
    return np.array(qgt) / 4
