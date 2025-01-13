using PauliPropagation


nl = 4
nq = 50

topology = bricklayertopology(nq; periodic=false)
circuit = hardwareefficientcircuit(nq, nl; topology=topology)
# circuit = Gate[PauliRotation([:X], [1]), PauliRotation([:Z], [1]), PauliRotation([:X], [1]), PauliRotation([:X], [2]), PauliRotation([:Z], [2]), PauliRotation([:X], [2]), PauliRotation([:Y, :Y], [1, 2])]
# circuit = Gate[PauliRotation([:Y, :Y], [1, 2])]

# println(circuit)


#It would make no sense the the control overlaps with a non identity term. 
struct CPauli <: StaticGate
  control::Int  # The control qubit 
  pauli::PauliString
end

#The gate appies a Z if the pauli that is controlled anticommutes with the pstr
#and applies said pauli if the control anticommutes with Z.
function PauliPropagation.apply(gate::CPauli, pstr, coefficient; kwargs...)
  control_pauli = getpauli(pstr, gate.control)

  control_commutes = commutes(control_pauli, :Z)
  target_commutes = commutes(pstr, gate.pauli.term)

  effect_pauli = zero(pstr)
  if !target_commutes
    effect_pauli = setpauli(effect_pauli, :Z, gate.control)
  end
  if !control_commutes
    effect_pauli = effect_pauli + gate.pauli.term
  end

  sign, new_pstr = pauliprod(pstr, effect_pauli)

  if !control_commutes & !target_commutes
    return new_pstr, -real(sign) * coefficient
  end

  return new_pstr, real(sign) * coefficient
end

function get_ith_parameter(i::Int)
  return i
end

function add_control_pauli!(nq, circ, i::Int)
  i_param = get_ith_parameter(i)
  cpauli = CPauli(nq + 1, PauliString(nq + 1, circ[i_param].symbols, circ[i_param].qinds))
  insert!(circ, i_param, cpauli)
end


#We assume that i<j
function qgt_element(nq, circ, thetas, i::Int, j::Int)
  new_circ = deepcopy(circ)

  add_control_pauli!(nq, new_circ, i)
  insert!(new_circ, i + 1, CliffordGate(:X, nq + 1))
  add_control_pauli!(nq, new_circ, j + 2)
  insert!(new_circ, j + 3, CliffordGate(:X, nq + 1))

  # for (index, value) in enumerate(new_circ)
  #   println("$index $value")
  # end
  return overlapwithplus(propagate(new_circ, PauliString(nq + 1, :X, nq + 1), thetas, max_weight=5, truncation_threshold=0.001))
end





nparams = countparameters(circuit)
pstr = PauliString(nq, [:I, :Z], [1, 2])



using Random
Random.seed!(42)
thetas = randn(nparams);

for i in 1:nparams
  for j in i+1:nparams
    println(i, ",", j)
    qgt_element(nq, circuit, thetas, i, j)
  end

end

println(qgt_element(nq, circuit, thetas, 7, 14))

# println(overlapwithzero(propagate(circuit, pstr, thetas)))
# println(overlapwithzero(propagate(new_circuit, pstr, thetas)))



