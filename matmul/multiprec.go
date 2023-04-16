package matmul

import (
	"errors"
	"math"
	"math/big"
	"sync"
)

type IntMatrix struct {
	Rows int
	Cols int
	Vals [][]*big.Int
}
type uint64Matrix struct {
	Rows int
	Cols int
	Vals [][]uint64
}

// Extended Euclidean Algorithm
func inteeu(a int64, b int64) (gcd int64, s int64, t int64) {
	var s1 int64 = 0
	var s2 int64 = 1
	var t1 int64 = 1
	var t2 int64 = 0
	r1 := b
	r2 := a
	for r1 != 0 {
		q := r2 / r1
		r2, r1 = r1, r2-q*r1
		s2, s1 = s1, s2-q*s1
		t2, t1 = t1, t2-q*t1
	}

	return r2, s2, t2
}

func uint64pow(a uint64, b uint) uint64 {
	var res uint64 = 1
	for b > 0 {
		if (b & 1) != 0 {
			res *= a
		}
		b >>= 1
		a *= a
	}
	return res
}
func int64pow(a int64, b int64) int64 {
	var res int64 = 1
	for b > 0 {
		if (b & 1) != 0 {
			res *= a
		}
		b >>= 1
		a *= a
	}
	return res
}
func reduceint(n *big.Int, m uint) uint {
	bits := n.Bits()
	var base uint = (1 << 32) % m
	base = (base * base) % m
	var acc uint = 0
	var mul uint = 1
	for _, x := range bits {
		acc = (acc + ((uint(x)%m)*mul)%m) % m
		mul = (mul * base) % m
	}
	if n.Sign() == -1 {
		return m - acc
	}
	return acc
}
func mmax(M *IntMatrix) *big.Int {
	v := M.Vals
	m := big.NewInt(0)
	for i := 0; i < M.Rows; i++ {
		for j := 0; j < M.Cols; j++ {
			elt := v[i][j]
			if elt.CmpAbs(m) == 1 {
				m.Set(elt)
			}
		}
	}
	return m
}

func make_uint_matrix(r int, c int) *uint64Matrix {
	//ensures matrix is contiguous
	out_matrix := make([][]uint64, r)
	orows := make([]uint64, r*c)
	for i := 0; i < r; i++ {
		out_matrix[i] = orows[i*c : (i+1)*c]
	}
	return &uint64Matrix{r, c, out_matrix}
}
func make_Int_matrix(r int, c int) *IntMatrix {
	//ensures matrix is contiguous
	out_matrix := make([][]*big.Int, r)
	orows := make([]*big.Int, r*c)
	for i := 0; i < r; i++ {
		out_matrix[i] = orows[i*c : (i+1)*c]
		for j := 0; j < c; j++ {
			out_matrix[i][j] = big.NewInt(0)
		}
	}
	return &IntMatrix{r, c, out_matrix}
}

const blocksize int = 64

func Mul_modular(A *IntMatrix, B *IntMatrix) (*IntMatrix, error) {
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
	for u := 1; u < nmoduli; u++ {
		m := moduli[u]
		var base uint64 = (1 << 32) % m
		base = (base * base) % m
		//Precompute (2^64)^i mod m
		powers := make([]uint64, maxbitsab)
		var st uint64 = 1
		powers[0] = 1
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
					r = (r + ((uint64(x)%m)*powers[bi])%m) % m
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
				go func(C *IntMatrix, A *uint64Matrix, B *uint64Matrix, srow int, erow int, scol int, ecol int, mod uint64) {
					defer exgroup.Done()
					for i := srow; i < erow; i++ {
						for j := scol; j < ecol; j++ {
							var r uint64 = 0
							for k := 0; k < A.Cols; k++ {
								r = ((A.Vals[i][k]*B.Vals[k][j])%mod + r) % mod
							}
							C.Vals[i][j].SetUint64(uint64(r))
						}
					}
				}(final_mats[u], reduced_mats_A[u], reduced_mats_B[u], ii, erow, jj, ecol, moduli[u])
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
