using PauliPropagation


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
  println(sign)

  if !control_commutes & !target_commutes
    return new_pstr, -real(sign) * coefficient
  end

  return new_pstr, real(sign) * coefficient
end



new_cnot = CPauli(1, PauliString(2, [:X], [2]))

for symbols in [[:I, :I], [:I, :X], [:I, :Y], [:I, :Z], [:X, :I], [:X, :X], [:X, :Y], [:X, :Z], [:Y, :I], [:Y, :X], [:Y, :Y], [:Y, :Z], [:Z, :I], [:Z, :X], [:Z, :Y], [:Z, :Z]]
  println(symbols)
  test_pauli = PauliString(2, symbols, [1, 2]).term
  propagated_custom = PauliPropagation.apply(new_cnot, test_pauli, 1.0)
  propagated_default = PauliPropagation.apply(CliffordGate(:CNOT, (1, 2)), test_pauli, 1.0)
  println(propagated_custom, "  ", propagated_default)
  # println(PauliString(2, propagated[1], [1, 2]), propagated[2])
end
