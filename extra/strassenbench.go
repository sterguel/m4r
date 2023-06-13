package main

import (
	"fmt"
	"math/rand"
	"time"
)

const BLOCKSIZE = 32

func min[T uint64 | uint | int64 | int](a T, b T) T {
	if a < b {
		return a
	}
	return b
}
func matmul(A [][]int64, B [][]int64) [][]int64 {
	// (MxN) * (NxP)
	n := len(A)
	C := make([][]int64, n)
	Crows := make([]int64, n*n)
	for i := 0; i < n; i++ {
		C[i] = Crows[i*n : (i+1)*n]
	}
	for kk := 0; kk < n; kk += BLOCKSIZE {
		for jj := 0; jj < n; jj += BLOCKSIZE {
			for i := 0; i < n; i++ {
				for k := kk; k < min(kk+BLOCKSIZE, n); k++ {
					r := A[i][k]
					for j := jj; j < min(jj+BLOCKSIZE, n); j++ {
						C[i][j] += r * B[k][j]
					}
				}
			}
		}
	}
	return C
}
func makemat(n int) [][]int64 {
	C := make([][]int64, n)
	Crows := make([]int64, n*n)
	for i := 0; i < n; i++ {
		C[i] = Crows[i*n : (i+1)*n]
	}
	return C
}
func mat_add(C [][]int64, A [][]int64, B [][]int64) {
	for i := 0; i < len(C); i++ {
		for j := 0; j < len(C); j++ {
			C[i][j] = A[i][j] + B[i][j]
		}
	}
}
func mat_sub(C [][]int64, A [][]int64, B [][]int64) {
	for i := 0; i < len(C); i++ {
		for j := 0; j < len(C); j++ {
			C[i][j] = A[i][j] - B[i][j]
		}
	}
}
func mat_subblock(A [][]int64, oi int, oj int, m int) [][]int64 {
	B := makemat(m)
	for i := oi; i < oi+m; i++ {
		for j := oj; j < oj+m; j++ {
			B[i-oi][j-oj] = A[i][j]
		}
	}
	return B
}
func Mul_strassen(A [][]int64, B [][]int64) [][]int64 {
	//assume square
	n2 := len(A)
	n := n2 / 2
	A1, A2, A3, A4 := mat_subblock(A, 0, 0, n), mat_subblock(A, 0, n, n), mat_subblock(A, n, 0, n), mat_subblock(A, n, n, n)
	B1, B2, B3, B4 := mat_subblock(B, 0, 0, n), mat_subblock(B, 0, n, n), mat_subblock(B, n, 0, n), mat_subblock(B, n, n, n)
	M1 := makemat(n)
	mat_add(M1, A1, A4)
	M2 := makemat(n)
	mat_add(M2, B1, B4)
	M1 = matmul(M1, M2)
	mat_add(M2, A3, A4)
	M2 = matmul(M2, B1)
	M3, M4, M5, M6, M7, Mf := makemat(n), makemat(n), makemat(n), makemat(n), makemat(n), makemat(n)
	mat_sub(M3, B2, B4)
	M3 = matmul(A1, M3)
	mat_sub(M4, B3, B1)
	M4 = matmul(A4, M4)
	mat_add(M5, A1, A2)
	M5 = matmul(M5, B4)
	mat_sub(M6, A3, A1)
	mat_add(M7, B1, B2)
	M6 = matmul(M6, M7)
	mat_sub(M7, A2, A4)
	mat_add(Mf, B3, B4)
	M7 = matmul(M7, Mf)
	C := makemat(n2)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			C[i][j] = M1[i][j] + M4[i][j] - M5[i][j] + M7[i][j]
		}
		for j := n; j < n2; j++ {
			jo := j - n
			C[i][j] = M3[i][jo] + M5[i][jo]
		}
	}
	for io := n; io < n2; io++ {
		i := io - n
		for j := 0; j < n; j++ {
			C[io][j] = M2[i][j] + M4[i][j]
		}
		for j := n; j < n2; j++ {
			jo := j - n
			C[io][j] = M1[i][jo] - M2[i][jo] + M3[i][jo] + M6[i][jo]
		}
	}
	return C
}
func testn(n int) (int64, int64) {
	A := makemat(2 * n)
	B := makemat(2 * n)
	for i := 0; i < 2*n; i++ {
		for j := 0; j < 2*n; j++ {
			A[i][j] = rand.Int63n(2_000_000) - 1_000_000
			B[i][j] = rand.Int63n(2_000_000) - 1_000_000
		}
	}
	startA := time.Now().UnixMicro()
	C1 := matmul(A, B)
	endA := time.Now().UnixMicro() - startA
	startB := time.Now().UnixMicro()
	C2 := Mul_strassen(A, B)
	endB := time.Now().UnixMicro() - startB
	for i := 0; i < 2*n; i++ {
		for j := 0; j < 2*n; j++ {
			if C1[i][j] != C2[i][j] {
				fmt.Println("Strassen produced incorrect matrix")
				break
			}
		}
	}
	return endA, endB
}
func makerange(start int, length int, inc int) []int {
	out := make([]int, length)
	for i := range out {
		out[i] = start + i*inc
	}
	return out
}
func main() {
	rand.Seed(time.Now().UnixNano())
	r := makerange(2, 250, 2)
	mn := make([]int64, 250)
	ms := make([]int64, 250)
	for i, v := range r {
		fmt.Println(2 * i)
		mn[i], ms[i] = testn(v)
	}
	fmt.Println(r)
	fmt.Println(mn)
	fmt.Println(ms)
}
