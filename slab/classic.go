package slab

import (
	"errors"
	"math/big"
)

func Mul_Int_standard(A *IntMatrix, B *IntMatrix) (*IntMatrix, error) {
	if A.Cols != B.Rows {
		return nil, errors.New("mismatched dimensions")
	}
	C := Make_Int_matrix(A.Rows, B.Cols)
	for i, row := range C.Vals {
		for j, v := range row {
			for k := 0; k < A.Cols; k++ {
				p := big.NewInt(0)
				p.Mul(A.Vals[i][k], B.Vals[k][j])
				v.Add(v, p)
			}
		}
	}
	return C, nil
}

// Set C = AB with C preinitialised and not equal to A or B.
func (C Matrix[T]) Mul_singlethreaded(A Matrix[T], B Matrix[T]) {
	zero := C[0][0].Zero()
	s := C[0][0].Zero()
	for i, row := range C {
		for j, val := range row {
			arow := A[i]
			val.Set(zero)
			for k := 0; k < len(B) && k < len(arow); k++ {
				s.Mul(arow[k], B[k][j])
				val.Add(val, s)
			}
		}
	}
}
