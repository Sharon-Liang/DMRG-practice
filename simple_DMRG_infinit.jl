# This is a simple DMRG practice for infinit 1D XXZ model
# H = J(Sx Sx + Sy Sy) + Jz Sz Sz
#   =J/2 (S+ S- + S- S+) + Jz Sz Sz

# packages

# Formation of Hamiltonian (basis: eigenstates of Sz)
  # pauli matrics
  sp = [0 1 ; 0 0] ; sz = 1/2 * [1 0 ; 0 -1];

  # coupling constants
  J = 1 ; Jz = 1;

  # construct Hamiltonian
  # two H-matrix are needed : H and  H_enlarge
  # H : original Hamiltonian of A (B) ;  H_enlarge: Hamiltonian of A .. B
  
