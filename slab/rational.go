package slab

import (
	"math/big"
	"scale/sfax"
)

func NewRatMatrix(rows int, cols int) RatMatrix {
	zero := &sfax.Rat{N: new(big.Rat)}
	out_matrix := make([][]sfax.Rat, rows)
	orows := make([]sfax.Rat, rows*cols)
	for i := 0; i < rows; i++ {
		out_matrix[i] = orows[i*cols : (i+1)*cols]
		for j := 0; j < cols; j++ {
			out_matrix[i][j] = zero.Copy()
		}
	}
	return out_matrix
}

// Set C = A+B. This will never go out of bounds so can be used to set sub-blocks
func (C RatMatrix) Add(A RatMatrix, B RatMatrix) {
	for i := 0; i < len(C) && i < len(A) && i < len(B); i++ {
		arow, brow, crow := A[i], B[i], C[i]
		for j := 0; j < len(arow) && j < len(brow) && j < len(crow); j++ {
			crow[j].Add(arow[j], brow[j])
		}
	}
}

// no aliasing allowed
func lcm(x *big.Int, a *big.Int, b *big.Int) {
	x.GCD(nil, nil, a, b)
	x.Div(a, x)
	x.Mul(x, b)
}
func qclearrows(M RatMatrix) ([][]*big.Int, []*big.Int) {
	rows, cols := len(M), len(M[0])
	out := NewIntMatrix(rows, cols)
	lcms := make([]*big.Int, rows)
	l2 := new(big.Int) //avoid aliasing issue with lcm
	for i, row := range M {
		lc := new(big.Int).Set(row[0].N.Denom()) //necessary to not change M
		for _, val := range row[1:] {
			lcm(l2, lc, val.N.Denom())
			lc.Set(l2)
		}
		lcms[i] = lc
		for j, u := range out[i] {
			u.Div(lc, row[j].N.Denom())
			u.Mul(u, row[j].N.Num())
		}
	}
	return out, lcms
}
func qclearcols(M RatMatrix) ([][]*big.Int, []*big.Int) {
	rows, cols := len(M), len(M[0])
	out := NewIntMatrix(rows, cols)
	lcms := make([]*big.Int, cols)
	l2 := new(big.Int) //avoid aliasing issue with lcm
	for j := 0; j < cols; j++ {
		lc := new(big.Int).Set(M[0][j].N.Denom())
		for i := 1; i < rows; i++ {
			lcm(l2, lc, M[i][j].N.Denom())
			lc.Set(l2)
		}
		lcms[j] = lc
		for i := 0; i < rows; i++ {
			out[i][j].Div(lc, M[i][j].N.Denom())
			out[i][j].Mul(out[i][j], M[i][j].N.Num())
		}
	}
	return out, lcms
}
func (C RatMatrix) Mul_multimod(A RatMatrix, B RatMatrix) {
	Aclear, rlcms := qclearrows(A)
	Bclear, clcms := qclearcols(B)
	Cclear, _ := Mul_modular_montgomery_jobs(&IntMatrix{len(Aclear), len(Aclear[0]), Aclear},
		&IntMatrix{len(Bclear), len(Bclear[0]), Bclear})
	sc := new(big.Int)
	for i, row := range C {
		for j, val := range row {
			sc.Mul(rlcms[i], clcms[j])
			val.N.SetFrac(Cclear.Vals[i][j], sc)
		}
	}
}
func (C RatMatrix) Mul_cleared_normal(A RatMatrix, B RatMatrix) {
	Aclear, rlcms := qclearrows(A)
	Bclear, clcms := qclearcols(B)
	Cclear, _ := Mul_Int_standard(&IntMatrix{len(Aclear), len(Aclear[0]), Aclear},
		&IntMatrix{len(Bclear), len(Bclear[0]), Bclear})
	sc := new(big.Int)
	for i, row := range C {
		for j, val := range row {
			sc.Mul(rlcms[i], clcms[j])
			val.N.SetFrac(Cclear.Vals[i][j], sc)
		}
	}
}
