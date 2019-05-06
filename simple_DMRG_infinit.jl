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

  # Two struct are needed to store A and A.
  # one is input, one is output
  # length is the length of A or A.
  # basis_dim is the dimension of the Hilbert space
  # basis_dim is used to constract unit matrix afterwards
  struct Block
    length :: Int
    basis_dim :: Int
    operator_dict :: Dict{Symbol,AbstractMatrix{Float64}}
  end

  struct EnlargedBlock
    length :: Int
    basis_dim :: Int
    operator_dict :: Dict{Symbol,AbstractMatrix{Float64}}
  end


  # Constract Hamiltonian
    # single site Hamiltonian
    H1 = zeros(Float64,2,2)

    # two sites Hamiltonian
    function H2(sp1,sz1,sp2,sz2)
      J/2 * (kron(sp1,sp2') + kron(sp1',sp2)) + Jz * kron(sz1,sz2)
    end

    # Hamiltonain of A .
    function enlarge(block :: Block)
             # NOTE: input is the kind of block !
             # We have a lot of imputs !
      # dimesion of A part
      d = block.basis_dim
      o = block.operator_dict

      # unit matrices needed
      I2 = Matrix{Float64}(I, 2, 2)
      Iblock = Matrix{Float64}(I, d, d)

      # hamiltonian of A.
      enlarged_operator_dict = Dict{Symbol,AbstractMatrix{Float64}}(
                             :H => kron(o[:H],I2) + kron(Iblock,H1) + H2(o[:conn_sp],o[:conn_sz],sp,sz)
                             :conn_sp => kron(Iblock,o[:conn_sp])
                             :conn_sz => kron(Iblock,o[:conn_sz])
                              )
      half = kron(block,unit2) + kron(unit_minus, H2(sp,sz,sp,sz))
      return EnlargedBlock(block.length + 1,
                            d * 2,
                            enlarged_operator_dict)
    end

    # initialize Hamiltonian : single site .
    H_initial = Dict{Symbol,AbstractMatrix{Float64}}(
               :H => H1
               :conn_sp => sp
               :conn_sz => sz
    )





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
      P = psi.vectors[:, len - D -1, len]
    end

    # projected Hamiltonian
    half_new = P'* h * P
  end
