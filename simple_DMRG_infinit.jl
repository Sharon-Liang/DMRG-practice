# This is a simple DMRG practice for infinit 1D XXZ model
# H = J(Sx Sx + Sy Sy) + Jz Sz Sz
#   =J/2 (S+ S- + S- S+) + Jz Sz Sz

# packages : copy from example
using Arpack
using LinearAlgebra
using SparseArrays

# Formation of Hamiltonian (basis: eigenstates of Sz)
  # pauli matrics
  sp = [0 1 ; 0 0] ; sz = 1/2 * [1 0 ; 0 -1];

  # coupling constants
  # use const to claim J and Jz are unchanged.
  # In julia, external variables are changable, that may slow the program
  # down, and local variables don't have this problem.
  const J = 1
  const Jz = 1

  # construct Hamiltonian
  # two H-matrix are needed : H and  H_enlarge
  # H : original Hamiltonian of A (B) ;  H_enlarge: Hamiltonian of A .. B
  ##struct Block
  ##  length :: Int  # length of system L(A)
  ##  basis_number :: Int # number of basis (D) after cutoff
  ##  hamintonain :: AbstractMatrix{Float64}
  ##end

  ##struct Block_enlarge
    ##length :: Int  # length of system A .. B
    ##basis_number :: Int # number of basis (D) after cutoff
    ##hamintonain :: AbstractMatrix{Float64}
  ##end

  # A fuction need to constract Hamiltonian
    # single site Hamiltonian
    H1 = zeros(Float64,2,2)

    # two sites Hamiltonian
    function H2(sp1,sz1,sp2,sz2)
      J/2 * (kron(sp1,sp2') + kron(sp1',sp2)) + Jz * kron(sz1,sz2)
    end

    # complete Hamiltonian
    # only two sites are different


    # starting point:two sites ..

    # ground state : eigenstate with smallest eigenvalue
    Block.hamintonain = H2(sp,sz,sp,sz)
