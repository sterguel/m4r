package matmul

import "errors"

const BLOCKSIZE = 32

func Make_int64_matrix(rows int, cols int) *[][]int64 {
	C := make([][]int64, rows)
	Crows := make([]int64, rows*cols)
	for i := 0; i < rows; i++ {
		C[i] = Crows[i*cols : (i+1)*cols]
	}
	return &C
}

func Mul_singlethread_nonblocked(Ap *[][]int64, Bp *[][]int64) (*[][]int64, error) {
	A := *Ap
	B := *Bp
	m := len(A)
	k := len(B[0])
	n1, n2 := len(A[0]), len(B)
	if n1 != n2 {
		return nil, errors.New("mismatched dimensions")
	}
	C := make([][]int64, m)
	Crows := make([]int64, m*k)
	for i := 0; i < m; i++ {
		C[i] = Crows[i*k : (i+1)*k]
	}

	for i := 0; i < m; i++ {
		for j := 0; j < k; j++ {
			var s int64 = 0
			for l := 0; l < n1; l++ {
				s += A[i][l] * B[l][j]
			}
			C[i][j] = s
		}
	}
	return &C, nil
}

func Mul_singlethread_blocked(Ap *[][]int64, Bp *[][]int64) (*[][]int64, error) {
	// (MxN) * (NxP)
	A := *Ap
	B := *Bp
	m := len(A)
	p := len(B[0])
	n1, n2 := len(A[0]), len(B)
	if n1 != n2 {
		return nil, errors.New("mismatched dimensions")
	}
	C := make([][]int64, m)
	Crows := make([]int64, m*p)
	for i := 0; i < m; i++ {
		C[i] = Crows[i*p : (i+1)*p]
	}
	rowblocks := m / BLOCKSIZE
	colblocks := p / BLOCKSIZE
	//loop peel
	for ii := 0; ii < rowblocks; ii++ {
		for jj := 0; jj < colblocks; jj++ {
			for i := (ii * BLOCKSIZE); i < (ii+1)*BLOCKSIZE; i++ {
				for j := (jj * BLOCKSIZE); j < (jj+1)*BLOCKSIZE; j++ {
					var s int64 = 0
					for k := 0; k < n1; k++ {
						s += A[i][k] * B[k][j]
					}
					C[i][j] = s
				}
			}
		}
	}
	for jj := 0; jj < colblocks; jj++ {
		for i := rowblocks * BLOCKSIZE; i < m; i++ {
			for j := (jj * BLOCKSIZE); j < (jj+1)*BLOCKSIZE; j++ {

				var s int64 = 0
				for k := 0; k < n1; k++ {
					s += A[i][k] * B[k][j]
				}
				C[i][j] = s
			}
		}
	}
	for ii := 0; ii < rowblocks; ii++ {
		for i := ii * BLOCKSIZE; i < (ii+1)*BLOCKSIZE; i++ {
			for j := colblocks * BLOCKSIZE; j < p; j++ {

				var s int64 = 0
				for k := 0; k < n1; k++ {
					s += A[i][k] * B[k][j]
				}
				C[i][j] = s
			}
		}
	}
	for i := rowblocks * BLOCKSIZE; i < m; i++ {
		for j := colblocks * BLOCKSIZE; j < p; j++ {

			var s int64 = 0
			for k := 0; k < n1; k++ {
				s += A[i][k] * B[k][j]
			}
			C[i][j] = s
		}
	}
	return &C, nil
}
func intmin(a int, b int) int {
	if a < b {
		return a
	}
	return b
}
func Mul_singlethreaded_blocked_RLM(Ap *[][]int64, Bp *[][]int64) (*[][]int64, error) {
	// (MxN) * (NxP)
	A := *Ap
	B := *Bp
	m := len(A)
	p := len(B[0])
	n1, n2 := len(A[0]), len(B)
	if n1 != n2 {
		return nil, errors.New("mismatched dimensions")
	}
	C := make([][]int64, m)
	Crows := make([]int64, m*p)
	for i := 0; i < m; i++ {
		C[i] = Crows[i*p : (i+1)*p]
	}
	for kk := 0; kk < n1; kk += BLOCKSIZE {
		for jj := 0; jj < p; jj += BLOCKSIZE {
			for i := 0; i < m; i++ {
				for k := kk; k < intmin(kk+BLOCKSIZE, n1); k++ {
					r := A[i][k]
					for j := jj; j < intmin(jj+BLOCKSIZE, p); j++ {
						C[i][j] += r * B[k][j]
					}
				}
			}
		}
	}
	return &C, nil
}

// something is wrong with this one
func Mul_singlethreaded_blocked_peeled_RLM(Ap *[][]int64, Bp *[][]int64) (*[][]int64, error) {
	// (MxN) * (NxP)
	A := *Ap
	B := *Bp
	m := len(A)
	p := len(B[0])
	n1, n2 := len(A[0]), len(B)
	if n1 != n2 {
		return nil, errors.New("mismatched dimensions")
	}
	C := make([][]int64, m)
	Crows := make([]int64, m*p)
	for i := 0; i < m; i++ {
		C[i] = Crows[i*p : (i+1)*p]
	}
	for kk := 0; kk < n1-BLOCKSIZE; kk += BLOCKSIZE {
		for jj := 0; jj < p-BLOCKSIZE; jj += BLOCKSIZE {
			for i := 0; i < m; i++ {
				for k := kk; k < kk+BLOCKSIZE; k++ {
					r := A[i][k]
					for j := jj; j < jj+BLOCKSIZE; j++ {
						C[i][j] += r * B[k][j]
					}
				}
			}
		}
	}
	kblocks := (n1 / BLOCKSIZE) * BLOCKSIZE
	for jj := 0; jj < p-BLOCKSIZE; jj += BLOCKSIZE {
		for i := 0; i < m; i++ {
			for k := kblocks; k < n1; k++ {
				r := A[i][k]
				for j := jj; j < jj+BLOCKSIZE; j++ {
					C[i][j] += r * B[k][j]
				}
			}
		}
	}
	jblocks := (p / BLOCKSIZE) * BLOCKSIZE
	for kk := 0; kk < n1-BLOCKSIZE; kk += BLOCKSIZE {
		for i := 0; i < m; i++ {
			for k := kk; k < kk+BLOCKSIZE; k++ {
				r := A[i][k]
				for j := jblocks; j < p; j++ {
					C[i][j] += r * B[k][j]
				}
			}
		}
	}
	for i := 0; i < m; i++ {
		for k := kblocks; k < n1; k++ {
			r := A[i][k]
			for j := jblocks; j < p; j++ {
				C[i][j] += r * B[k][j]
			}
		}
	}
	return &C, nil
}
