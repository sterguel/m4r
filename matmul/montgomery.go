package matmul

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
	t := (x + m*n) >> 32
	if t >= n {
		return t - n
	}
	return t
}
func redadd(x uint64, n uint64) uint64 {
	if x > n {
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
	maxv := mmax(A)
	maxvb := mmax(B)
	maxbitsa := len(maxv.Bits())
	maxbitsb := len(maxvb.Bits())
	maxbitsab := imax(maxbitsa, maxbitsb)
	maxv.Mul(maxv, maxvb)
	maxv.Abs(maxv)
	maxv.Mul(maxv, big.NewInt(int64(A.Cols)))
	maxv.Mul(maxv, big.NewInt(2))

	var moduli []uint64
	moduli = append(moduli, 0)
	modmax := big.NewInt(1 << 32)
	modmax.Mul(modmax, modmax)
	var currentprime int64 = 3
	for modmax.Cmp(maxv) == -1 {
		bound := int64(64.0/math.Log2(float64(currentprime))) / 2
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
		reduced_mats_B[i] = make_uint_matrix(B.Rows, B.Cols)
		final_mats[i] = make_Int_matrix(A.Rows, B.Cols)
	}
	//special 2^64 case
	m2matAv := reduced_mats_A[0].Vals
	for i := 0; i < A.Rows; i++ {
		for j := 0; j < A.Cols; j++ {
			num := (A.Vals)[i][j]
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
	for i := 0; i < B.Rows; i++ {
		for j := 0; j < B.Cols; j++ {
			num := (B.Vals)[i][j]
			switch num.Sign() {
			case 1:
				m2matBv[i][j] = uint64(num.Bits()[0])
			case -1:
				m2matBv[i][j] = ^(uint64(num.Bits()[0])) + 1
			case 0:
				m2matBv[i][j] = 0
			}

		}
	}
	modinvs := make([]uint64, nmoduli)
	rinvs := make([]uint64, nmoduli)
	for u := 1; u < nmoduli; u++ {
		m := moduli[u]
		_, s, t := inteeu(1<<32, int64(m))
		rinvs[u] = uint64(pos_mod(s, int64(m)))
		//always positive as the sign bit is the leading bit
		t = t & (1<<32 - 1)
		tu := uint64(t)
		modinvs[u] = (1 << 32) - tu
		var base uint64 = (1 << 32) % m
		var st uint64 = base
		powers := make([]uint64, maxbitsab)
		powers[0] = base

		base = (base * base) % m
		//Precompute (2^64)^i mod m

		for i := 1; i < len(powers); i++ {
			st = (st * base) % m
			powers[i] = st
		}
		Ar := reduced_mats_A[u].Vals
		Br := reduced_mats_B[u].Vals
		//Reduce matrices mod m
		for i := 0; i < A.Rows; i++ {
			for j := 0; j < A.Cols; j++ {
				n := A.Vals[i][j]
				bits := n.Bits()
				var r uint64 = 0
				for bi, x := range bits {
					r = redadd((r + ((uint64(x)%m)*powers[bi])%m), m)
				}
				if n.Sign() == -1 {
					r = m - r
				}
				Ar[i][j] = r
			}
		}
		for i := 0; i < B.Rows; i++ {
			for j := 0; j < B.Cols; j++ {
				n := B.Vals[i][j]
				bits := n.Bits()
				var r uint64 = 0
				for bi, x := range bits {
					r = (r + ((uint64(x)%m)*powers[bi])%m) % m
				}
				if n.Sign() == -1 {
					r = m - r
				}
				Br[i][j] = r
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
							r += A.Vals[i][k] * B.Vals[k][j]
						}
						C.Vals[i][j].SetUint64(uint64(r))
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
								m := redc(A.Vals[i][k]*B.Vals[k][j], mod, modinvs[mix])
								r = redadd(m+r, mod)

							}
							C.Vals[i][j].SetUint64((r * rinvs[mix]) % mod)
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
	for i := 1; i < nmoduli; i++ {
		modlist[i] = big.NewInt(int64(moduli[i]))
	}
	for nmoduli > icoll {
		for i := 0; i < nmoduli-icoll; i += 2 * icoll {
			moda, modb := modlist[i], modlist[i+icoll]
			scratch.GCD(rx, sx, moda, modb)
			sx.Mul(sx, modb)
			rx.Mul(rx, moda)
			moda.Mul(moda, modb)
			for ci := 0; ci < A.Rows; ci++ {
				for cj := 0; cj < B.Cols; cj++ {
					t := final_mats[i].Vals[ci][cj]
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
	for ci := 0; ci < A.Rows; ci++ {
		for cj := 0; cj < B.Cols; cj++ {
			t := final_mats[0].Vals[ci][cj]
			if t.Cmp(halfway_bound) == 1 {
				t.Sub(t, modmax)
			}
		}
	}
	return final_mats[0], nil
}
