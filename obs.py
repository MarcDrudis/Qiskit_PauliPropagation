from qiskit.quantum_info import SparsePauliOp
from qiskit_paulipropagation.julia import sparsepauliop_to_pp,pp,jl
from juliacall import convert

obs= SparsePauliOp.from_sparse_list([("Y",[3],1.0),("X",[1],1.0)],4)


print(sparsepauliop_to_pp(obs))
#
pauli_symbols= pp.seval("Vector{Symbol}")()
qubits= pp.seval("Vector{Int}")()

for p in "XXZ":
    pp.push_b(pauli_symbols,convert(jl.Symbol,str(p)))
    pp.push_b(qubits,1)

print(pauli_symbols)
print(qubits)



psum=pp.PauliSum(4)
print(psum)

pp.add_b(psum,pauli_symbols,qubits,0.5)
print(psum)
