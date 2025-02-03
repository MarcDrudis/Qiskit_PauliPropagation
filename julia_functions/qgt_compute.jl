using PauliPropagation
using Base.Threads


#It would make no sense the the control overlaps with a non identity term. 
struct CPauli <: StaticGate
  control::Int  # The control qubit 
  pauli::PauliString
end

#The gate appies a Z if the pauli that is controlled anticommutes with the pstr
#and applies said pauli if the control anticommutes with Z.
function PauliPropagation.apply(gate::CPauli, pstr, coefficient; kwargs...)

  #Check that we are not controling a pauli that has a non Identity action
  #assert 0==getpauli(gate.pauli.term, gate.control)

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
function qgt_element(nq::Int, circ, thetas, i::Int, j::Int)
  if i >= j
    i, j = j, i
  end
  new_circ = deepcopy(circ)
  add_control_pauli!(nq, new_circ, i)
  insert!(new_circ, i + 1, CliffordGate(:X, nq + 1))
  add_control_pauli!(nq, new_circ, j + 2)
  insert!(new_circ, j + 3, CliffordGate(:X, nq + 1))
  push!(new_circ, CliffordGate(:H, nq + 1))
  insert!(new_circ, 1, CliffordGate(:H, nq + 1))
  obs = PauliString(nq + 1, :Z, nq + 1)
  return overlapwithzero(propagate(new_circ, obs, thetas, max_weight=7, truncation_threshold=0.000001))
end


function compute_qgt(nq::Int, circ, thetas, parameter_map, parameter_position)
  # Initialize an empty matrix
  nparams = countparameters(circ)
  qgt = zeros(Float64, nparams, nparams)
  # Total number of elements
  total_elements = nparams * nparams

  # Parallel computation over the flattened array
  for idx in 1:total_elements
    i = div(idx - 1, nparams) + 1  # Row index
    j = mod(idx - 1, nparams) + 1  # Column index

    # println("Thread ", threadid(), " is handling operation at (", i, ", ", j, ")")
    # println(nq, circ, thetas, i, j)
    qgt[i, j] = qgt_element(nq, circ, thetas, parameter_position[i], parameter_position[j])  # Call the provided function
  end

  return qgt
end
