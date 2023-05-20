package slab

import (
	"errors"
	"math/big"
)

func Mul_standard(A *IntMatrix, B *IntMatrix) (*IntMatrix, error) {
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
