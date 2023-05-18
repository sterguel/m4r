package sfax

import "math/big"

type commonElement[T any] interface {
	Copy() T
	String() string
	Equals(T) bool
	Set(T)
}
type MulMonoidElement[T any] interface {
	commonElement[T]
	One() T
	Mul(T, T)
	Times(T) T
}
type MulGroupElement[T any] interface {
	MulMonoidElement[T]
	Inv() T
	InvR(T)
	Div(T) T
	DivR(T, T)
}
type AddGroupElement[T any] interface {
	commonElement[T]
	Neg() T
	NegR(T)
	Plus(T) T
	Minus(T) T
	Add(T, T)
	Sub(T, T)
	Zero() T
}
type RingElement[T any] interface {
	MulMonoidElement[T]
	AddGroupElement[T]
}

type FieldData struct {
	Char   *big.Int // characteristic of the field
	Degree int      // degree of the field over its characteristic field, 0 for infinite
	Order  *big.Int // order of the field, nil if infinite
}
type FieldElement[T any] interface {
	FieldData() FieldData
	AddGroupElement[T]
	MulGroupElement[T]
}

func Power[T MulGroupElement[T]](x T, y int) T {
	if y < 0 {
		return Power(x.Inv(), -y)
	}
	u := x.Copy()
	res := x.One()
	for y > 0 {
		if (y & 1) != 0 {
			res.Mul(res, u)
		}
		y = y >> 1
		u.Mul(u, u)
	}
	return res
}
