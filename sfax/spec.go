package sfax

import "math/big"

// Call the variable being acted on x, and the arguments y,z
type commonElement[T any] interface {
	Copy() T //Creates a new copy of x
	String() string
	Equals(T) bool
	Set(T) //Set the value of x in-place to y
}
type MulMonoidElement[T any] interface {
	commonElement[T]
	One() T    //Return a new copy of the multiplicative identity
	Mul(T, T)  //Set x to y * z
	Times(T) T //Return z = x * y
}
type MulGroupElement[T any] interface {
	MulMonoidElement[T]
	Inv() T    //Return the multiplicative inverse of x
	InvR(T)    //Set x to the multiplicative inverse of y
	Div(T) T   //Return z = x / y
	DivR(T, T) //Set x = z / y
}
type AddGroupElement[T any] interface {
	commonElement[T]
	Neg() T    //Return -x
	NegR(T)    //Set x = -y
	Plus(T) T  //Return z = x + y
	Minus(T) T //Return z = x - y
	Add(T, T)  //Set x to y+z
	Sub(T, T)  //Set x to y - z
	Zero() T   //Return a new copy of the additive identity
}
type RingElement[T any] interface {
	MulMonoidElement[T]
	AddGroupElement[T]
}

type FieldData struct {
	Char   *big.Int // Characteristic of the field
	Degree int      // Degree of the field over its characteristic field, 0 for infinite
	Order  *big.Int // Order of the field, nil if infinite
}
type FieldElement[T any] interface {
	FieldData() FieldData //Return the corresponding FieldData struct
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
