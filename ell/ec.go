package ell

import "scale/sfax"

// Point on a curve of the form y^2 = x^3 + Ax + B
type ECPoint[T sfax.FieldElement[T]] struct {
	A   T
	B   T
	X   T
	Y   T
	Inf bool
}

func (x *ECPoint[T]) Zero() *ECPoint[T] {
	return &ECPoint[T]{
		x.A,
		x.B,
		x.A.Zero(),
		x.B.Zero(),
		true,
	}
}
