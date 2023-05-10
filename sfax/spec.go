package sfax

import "math/big"

type commonElement[T any] interface {
	Copy() T
	Equals(T) bool
}
type MulMonoidElement[T any] interface {
	commonElement[T]
	MulMonoid() MulMonoid[T]
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
	AddGroup() AddGroup[T]
	Neg() T
	NegR(T)
	Plus(T) T
	Minus(T) T
	Add(T, T)
	Sub(T, T)
}
type RingElement[T any] interface {
	MulMonoidElement[T]
	AddGroupElement[T]
}
type MulMonoid[T any] interface {
	One() T
}
type AddGroup[T any] interface {
	Zero() T
}
type Ring[T any] interface {
	MulMonoid[T]
	AddGroup[T]
}
type Field[T any] interface {
	Ring[T]
	Char() *big.Int
}
type FieldElement[T any] interface {
	Field() Field[T]
	AddGroupElement[T]
	MulGroupElement[T]
}

func Power[T MulGroupElement[T]](x T, y int) T {
	if y < 0 {
		return Power(x.Inv(), -y)
	}
	u := x.Copy()
	res := x.MulMonoid().One()
	for y > 0 {
		if (y & 1) != 0 {
			res.Mul(res, u)
		}
		y = y >> 1
		u.Mul(u, u)
	}
	return res
}
