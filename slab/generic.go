package slab

import (
	"math/big"
	"scale/sfax"
)

func NewMatrix[T sfax.FieldElement[T]](rows int, cols int, e T) Matrix[T] {
	out_matrix := make([][]T, rows)
	orows := make([]T, rows*cols)
	for i := 0; i < rows; i++ {
		out_matrix[i] = orows[i*cols : (i+1)*cols]
		for j := 0; j < cols; j++ {
			out_matrix[i][j] = e.Zero()
		}
	}
	return out_matrix
}
func NewIntMatrix(rows int, cols int) [][]*big.Int {
	out_matrix := make([][]*big.Int, rows)
	orows := make([]*big.Int, rows*cols)
	for i := 0; i < rows; i++ {
		out_matrix[i] = orows[i*cols : (i+1)*cols]
		for j := 0; j < cols; j++ {
			out_matrix[i][j] = new(big.Int)
		}
	}
	return out_matrix
}
