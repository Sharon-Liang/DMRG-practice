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

  # A fuction need to constract Hamiltonian
    # single site Hamiltonian
    H1 = zeros(Float64,2,2)

    # two sites Hamiltonian
    function H2(sp1,sz1,sp2,sz2)
      J/2 * (kron(sp1,sp2') + kron(sp1',sp2)) + Jz * kron(sz1,sz2)
    end

    # Hamiltonain of A .. B
    ## should add a projection operator
    function H3(len, block)
      # dimesion of A part
      d = 2^len
      # unit matrices needed
      unit2 = Matrix{Float64}(I, 2, 2)
      unit = Matrix{Float64}(I, d, d)
      unit_plus = Matrix{Float64}(I, d * 2, d * 2)
      unit_minus = Matrix{Float64}(I, d / 2, d / 2)
      # hamiltonian of A.
      half = kron(block,unit2) + kron(unit_minus, H2(sp,sz,sp,sz))
      # hamiltonian of A..B
      kron(half,unit_plus) + kron(unit_plus, half) + kron( kron(unit, H2(sp,sz,sp,sz)), unit)
    end

  # circulation
  const Lmax = 20     # we count Lmax * 2 sites in total
  const D = 1000      # we only keep maximum 1000 states, 9 points exact

  for L = 1 : Lmax
    # dim : dimension of total Hilbert space
    dim = 2^L

    if L = 1
      # h is hamiltonian for 4 points . . . .
      h = H3(L,H1)
    else
      h = H3(L,half_new)
    end

    # ground state : eigenstate with smallest eigenvalue
    gs = eigen(h).vectors[:,1]  # ground state
    len = sqrt( length(gs) )    # length of new half: L+1

    # reduced density matrix rho
    U = reshape(gs,len.len)
    rho = U * U'

    # projection operator
    if len < D
      P = psi.vectors
    else
      psi = eigen(rho)   # eigenvalues and eigenvectors : small -> large
      P = psi.vectors[:, len_cut - D -1, len_cut]
    end

    # projected Hamiltonian
    half_new = P'* h * P
  end
    
