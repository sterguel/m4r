package slab

import (
	"errors"
	"math"
	"math/big"
	"sync"
)

type matwork struct {
	cblock [][]*big.Int
	arows  [][]uint64
	bcols  [][]uint64
	mod    uint64
	inv    uint64
	scol   int
	ecol   int
}

func SetWorkerCount(w int) {
	workercount = w
}

var workercount int = 8

func mod_work(w *matwork) {
	a, b, c, mod, inv, e, s := w.arows, w.bcols, w.cblock, w.mod, w.inv, w.ecol, w.scol
	for i, row := range c {
		for j, val := range row[s:e] {
			var r uint64 = 0
			// trick from math/big to speed up bounds checking
			for k := 0; k < len(a[i]) && k < len(b[j]); k++ {
				m := redc(a[i][k]*b[j][k], mod, inv)
				r = redadd(m+r, mod)
			}
			val.SetUint64(redc(r, mod, inv))
		}
	}
}
func mod2_work(w *matwork) {
	a, b, c, e, s := w.arows, w.bcols, w.cblock, w.ecol, w.scol
	for i, row := range c {
		for j, val := range row[s:e] {
			var r uint64 = 0
			for k := 0; k < len(a[i]) && k < len(b[j]); k++ {
				r += a[i][k] * b[j][k]
			}
			val.SetUint64(r)
		}
	}
}
func Mul_modular_montgomery_jobs(A *IntMatrix, B *IntMatrix) (*IntMatrix, error) {
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

	njobs := ceildiv(A.Rows, blocksize) * ceildiv(B.Cols, blocksize) * nmoduli
	jobs := make([]matwork, njobs)

	ct := 0
	for ii := 0; ii < A.Rows; ii += blocksize {
		for jj := 0; jj < B.Cols; jj += blocksize {
			erow := imin(A.Rows, ii+blocksize)
			ecol := imin(B.Cols, jj+blocksize)
			jobs[ct] = matwork{
				final_mats[0].Vals[ii:erow],
				reduced_mats_A[0].Vals[ii:erow],
				reduced_mats_B[0].Vals[jj:ecol],
				0, 0, jj, ecol,
			}
			ct++
		}
	}
	for k := 1; k < len(moduli) && k < len(modinvs); k++ {
		for ii := 0; ii < A.Rows; ii += blocksize {
			for jj := 0; jj < B.Cols; jj += blocksize {
				erow := imin(A.Rows, ii+blocksize)
				ecol := imin(B.Cols, jj+blocksize)
				jobs[ct] = matwork{
					final_mats[k].Vals[ii:erow],
					reduced_mats_A[k].Vals[ii:erow],
					reduced_mats_B[k].Vals[jj:ecol],
					moduli[k], modinvs[k], jj, ecol,
				}
				ct++
			}
		}
	}
	//parallel multiplication
	var exgroup sync.WaitGroup
	jsplit := ceildiv(njobs, workercount)
	for w := 0; w < workercount; w++ {
		if jsplit*w >= njobs {
			break
		}
		exgroup.Add(1)
		go func(tasks []matwork) {
			defer exgroup.Done()
			for _, task := range tasks {
				if task.mod == 0 {
					mod2_work(&task)
				} else {
					mod_work(&task)
				}
			}
		}(jobs[jsplit*w : imin(jsplit*(w+1), njobs)])
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
