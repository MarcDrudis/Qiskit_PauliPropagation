
using PauliPropagation

println("hello 2")

nl = 4
nq = 5


cnot = PauliSum(2)
add!(cnot, [:I, :I], [1, 2], 0.5)
add!(cnot, [:Z, :I], [1, 2], 0.5)
add!(cnot, [:I, :X], [1, 2], 0.5)
add!(cnot, [:Z, :X], [1, 2], -0.5)


pstr = PauliString(nq, [:I, :Z], [1, 2])

println(pstr * pstr)
