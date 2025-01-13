from juliacall import Main as jl
from juliacall import Pkg as jlPkg
from juliacall import convert
from qiskit.circuit import QuantumCircuit
from qiskit.circuit.parametervector import ParameterVectorElement
from qiskit.compiler import transpile
from qiskit.converters import circuit_to_dag, dag_to_circuit
from qiskit.dagcircuit import DAGCircuit
from qiskit.quantum_info import SparsePauliOp

print("Activating PP")
jlPkg.activate("PauliPropagation")
jl.seval("using PauliPropagation")
pp = jl.PauliPropagation

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


def qc_to_pp(
    qc: QuantumCircuit | DAGCircuit,
) -> tuple[list[tuple[str, list[int]]], list[int]]:
    """Returns a list with all the PP gates in julia and a map with the parameter indices."""

    # The circuit must only contain supported gates.
    dag = qc if isinstance(qc, DAGCircuit) else circuit_to_dag(qc)
    num_qubits = qc.num_qubits() if isinstance(qc, DAGCircuit) else qc.num_qubits
    op_nodes = list(dag.topological_op_nodes())
    pp_circuit = []
    parameter_map = []
    for node in op_nodes:
        q_indices = tuple(qarg._index + 1 for qarg in node.qargs)
        name = node.op.name
        if name in pauli_rotations:
            pauli_rot = pp.PauliRotation(pauli_rotations[name], q_indices)
            # pp_circuit.append(pp.tofastgates(pauli_rot, num_qubits))
            pp_circuit.append(pauli_rot)
            if isinstance(node.op.params[0], float):
                parameter_map.append(node.op.params[0])
            else:
                parameter_map.append(node.op.params[0].index)
        elif name in clifford_gates:
            clifford_gate = pp.CliffordGate(clifford_gates[name], q_indices)
            # pp_circuit.append(pp.tofastgates(clifford_gate, num_qubits))
            pp_circuit.append(clifford_gate)
        else:
            print(f"We did not find a gate for {node.op.name}. Skipping Gate.")
    return pp_circuit, parameter_map


def pp_estimator(
    qc: QuantumCircuit | DAGCircuit,
    obs: SparsePauliOp,
    params: list[float],
    compile: bool = False,
    **kwargs,
):
    """Retruns the expectation value for an initial state of |0>."""
    pauli_sum = pp_propagation(qc, obs, params, compile, **kwargs)
    return pp.overlapwithzero(pauli_sum)


def sparsepauliop_to_pp(op: SparsePauliOp):
    nqubits = op.num_qubits
    pp_paulisum = pp.PauliSum(nqubits)
    for pauli, qubits, coefficient in op.to_sparse_list():
        pauli_symbols = pp.seval("Vector{Symbol}")()
        for p in pauli:
            jl.push_b(pauli_symbols,convert(jl.Symbol,str(p)))
        pp_qubits= pp.seval("Vector{Int}")()
        for q in qubits:
            jl.push_b(pp_qubits,q+1)
        #Here I am assuming the coefficient will be real. Good Luck!
        pp.add_b(pp_paulisum,pauli_symbols,pp_qubits,coefficient.real)
    return pp_paulisum


def propagation(
    pp_circuit, parameter_map, pp_observable, params: list[float], **kwargs
):

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
    pp_circuit, parameter_map = qc_to_pp(qc)
    pp_observable = sparsepauliop_to_pp(obs)
    return propagation(pp_circuit, parameter_map, pp_observable, params, **kwargs)
