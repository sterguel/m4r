package slab

import (
	"errors"
	"math"
	"math/big"
	"sync"
)

// hopefully this just gets inlined
// R is always 2^32
const redcmask uint64 = (1 << 32) - 1

func redc(x uint64, n uint64, ninv uint64) uint64 {
	m := ((x & redcmask) * ninv) & redcmask
	p := m * n
	t := (x + p) >> 32
	if t >= n {
		return t - n
	}
	return t
}
func redadd(x uint64, n uint64) uint64 {
	if x >= n {
		return x - n
	}
	return x
}
func pos_mod(x int64, n int64) int64 {
	if x < 0 {
		return (x % n) + n
	}
	return x % n
}
func Mul_modular_montgomery(A *IntMatrix, B *IntMatrix) (*IntMatrix, error) {
	if A.Cols != B.Rows {
		return nil, errors.New("mismatched dimensions")
	}
	maxv, hasnega := mmax(A)
	maxvb, hasnegb := mmax(B)
	maxbitsa := len(maxv.Bits())
	maxbitsb := len(maxvb.Bits())
	maxbitsab := imax(maxbitsa, maxbitsb)
	maxv.Mul(maxv, maxvb)
	maxv.Abs(maxv)
	maxv.Mul(maxv, big.NewInt(int64(A.Cols)))
	if hasnega || hasnegb {
		maxv.Mul(maxv, big.NewInt(2))
	}

	var moduli []uint64
	moduli = append(moduli, 0)
	modmax := big.NewInt(1 << 32)
	modmax.Mul(modmax, modmax)
	var currentprime int64 = 3
	for modmax.Cmp(maxv) == -1 {
		bound := int64(31.0 / math.Log2(float64(currentprime)))
		mod := int64pow(currentprime, bound)
		modmax.Mul(modmax, big.NewInt(mod))
		moduli = append(moduli, uint64(mod))
		currentprime = nextprime(currentprime)
	}
	halfway_bound := big.NewInt(0).Div(modmax, big.NewInt(2))
	nmoduli := len(moduli)
	reduced_mats_A := make([]*uint64Matrix, nmoduli)
	reduced_mats_B := make([]*uint64Matrix, nmoduli)
	final_mats := make([]*IntMatrix, nmoduli)
	for i := 0; i < nmoduli; i++ {
		reduced_mats_A[i] = make_uint_matrix(A.Rows, A.Cols)
		reduced_mats_B[i] = make_uint_matrix(B.Cols, B.Rows)
		final_mats[i] = Make_Int_matrix(A.Rows, B.Cols)
	}
	//special 2^64 case
	m2matAv := reduced_mats_A[0].Vals
	for i, row := range A.Vals {
		for j, num := range row {
			switch num.Sign() {
			case 1:
				m2matAv[i][j] = uint64(num.Bits()[0])
			case -1:
				m2matAv[i][j] = ^(uint64(num.Bits()[0])) + 1
			case 0:
				m2matAv[i][j] = 0
			}

		}
	}
	m2matBv := reduced_mats_B[0].Vals
	for i, row := range B.Vals {
		for j, num := range row {
			switch num.Sign() {
			case 1:
				m2matBv[j][i] = uint64(num.Bits()[0])
			case -1:
				m2matBv[j][i] = ^(uint64(num.Bits()[0])) + 1
			case 0:
				m2matBv[j][i] = 0
			}

		}
	}
	modinvs := make([]uint64, nmoduli)
	for Mi, m := range moduli[1:] {
		u := Mi + 1
		_, _, t := inteeu(1<<32, int64(m))
		//always positive as the sign bit is the leading bit
		t = t & (1<<32 - 1)
		tu := uint64(t)
		modinvs[u] = (1 << 32) - tu
		var base uint64 = (1 << 32) % m

		powers := make([]uint64, maxbitsab)
		base = (base * base) % m
		var st uint64 = base
		powers[0] = base
		//Precompute (2^64)^i mod m
		for i := 1; i < len(powers); i++ {
			st = (st * base) % m
			powers[i] = st
		}

		Ar := reduced_mats_A[u].Vals
		Br := reduced_mats_B[u].Vals
		//Reduce matrices mod m
		for i, row := range A.Vals {
			for j, n := range row {
				bits := n.Bits()
				var r uint64 = 0
				for bi, x := range bits {
					u := redc((uint64(x)%m)*powers[bi], m, modinvs[u])
					r = redadd(r+u, m)
				}
				if n.Sign() == -1 {
					r = m - r
				}
				Ar[i][j] = r
			}
		}
		for i, row := range B.Vals {
			for j, n := range row {
				bits := n.Bits()
				var r uint64 = 0
				for bi, x := range bits {
					u := redc((uint64(x)%m)*powers[bi], m, modinvs[u])
					r = redadd(r+u, m)
				}
				if n.Sign() == -1 {
					r = m - r
				}
				Br[j][i] = r
			}
		}
	}

	//parallel multiplication
	var exgroup sync.WaitGroup
	for ii := 0; ii < A.Rows; ii += blocksize {
		for jj := 0; jj < B.Cols; jj += blocksize {
			erow := imin(A.Rows, ii+blocksize)
			ecol := imin(B.Cols, jj+blocksize)
			exgroup.Add(1)
			go func(C *IntMatrix, A *uint64Matrix, B *uint64Matrix, srow int, erow int, scol int, ecol int) {
				defer exgroup.Done()
				for i := srow; i < erow; i++ {
					for j := scol; j < ecol; j++ {
						var r uint64 = 0
						for k := 0; k < A.Cols; k++ {
							r += A.Vals[i][k] * B.Vals[j][k]
						}
						C.Vals[i][j].SetUint64(r)
					}
				}

			}(final_mats[0], reduced_mats_A[0], reduced_mats_B[0], ii, erow, jj, ecol)
			for u := 1; u < nmoduli; u++ {
				exgroup.Add(1)
				go func(C *IntMatrix, A *uint64Matrix, B *uint64Matrix, srow int, erow int, scol int, ecol int, mix int) {
					defer exgroup.Done()
					mod := moduli[mix]
					for i := srow; i < erow; i++ {
						for j := scol; j < ecol; j++ {

							var r uint64 = 0
							for k := 0; k < A.Cols; k++ {

								m := redc(A.Vals[i][k]*B.Vals[j][k], mod, modinvs[mix])
								r = redadd(m+r, mod)

							}
							C.Vals[i][j].SetUint64(redc(r, mod, modinvs[mix]))
						}
					}
				}(final_mats[u], reduced_mats_A[u], reduced_mats_B[u], ii, erow, jj, ecol, u)
			}
		}
	}
	exgroup.Wait()

	icoll := 1
	//reconstruction
	modlist := make([]*big.Int, nmoduli)
	scratch, rx, sx := big.NewInt(0), big.NewInt(0), big.NewInt(0)
	modlist[0] = big.NewInt(1 << 32)
	modlist[0] = modlist[0].Mul(modlist[0], modlist[0])
	for i, m := range moduli[1:] {
		modlist[i+1] = big.NewInt(int64(m))
	}

	for nmoduli > icoll {
		for i := 0; i < nmoduli-icoll; i += 2 * icoll {
			moda, modb := modlist[i], modlist[i+icoll]
			scratch.GCD(rx, sx, moda, modb)
			sx.Mul(sx, modb)
			rx.Mul(rx, moda)
			moda.Mul(moda, modb)
			for ci, row := range final_mats[i].Vals {
				for cj, t := range row {

					u := final_mats[i+icoll].Vals[ci][cj]
					t.Mul(sx, t)
					u.Mul(rx, u)
					t.Add(t, u)
					t.Mod(t, moda)
				}
			}
		}
		icoll *= 2
	}
	if hasnega || hasnegb {
		for _, row := range final_mats[0].Vals {
			for _, t := range row {
				if t.Cmp(halfway_bound) == 1 {
					t.Sub(t, modmax)
				}
			}
		}
	}

	return final_mats[0], nil
}
