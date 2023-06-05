package sfax

import "math/big"

// Wrapper for big.Rat to satisfy the FieldElement structure
// The underlying *big.Rat value can be accessed and modified directly.
type Rat struct {
	N *big.Rat
}

func (x Rat) FieldData() FieldData {
	return FieldData{
		big.NewInt(0),
		1,
		nil,
	}
}
func (x Rat) Copy() Rat {
	return Rat{
		new(big.Rat).Set(x.N),
	}
}
func (x Rat) String() string {
	if x.N.IsInt() {
		return x.N.Num().String()
	}
	return x.N.String()

}
func (x Rat) Equals(y Rat) bool {
	return x.N.Cmp(y.N) == 0
}
func (x Rat) Set(y Rat) {
	x.N.Set(y.N)
}
func (x Rat) One() Rat {
	return Rat{
		big.NewRat(1, 1),
	}
}
func (x Rat) Zero() Rat {
	return Rat{
		big.NewRat(0, 1),
	}
}
func (x Rat) Add(a Rat, b Rat) {
	x.N.Add(a.N, b.N)
}
func (x Rat) Plus(y Rat) Rat {
	return Rat{
		new(big.Rat).Add(x.N, y.N),
	}
}
func (x Rat) Sub(a Rat, b Rat) {
	x.N.Sub(a.N, b.N)
}
func (x Rat) Minus(y Rat) Rat {
	return Rat{
		new(big.Rat).Sub(x.N, y.N),
	}
}
func (x Rat) Mul(a Rat, b Rat) {
	x.N.Mul(a.N, b.N)
}
func (x Rat) Times(y Rat) Rat {
	return Rat{
		new(big.Rat).Mul(x.N, y.N),
	}
}
func (x Rat) DivR(a Rat, b Rat) {
	x.N.Quo(a.N, b.N)
}
func (x Rat) Div(y Rat) Rat {
	return Rat{
		new(big.Rat).Quo(x.N, y.N),
	}
}
func (x Rat) Inv() Rat {
	return Rat{
		new(big.Rat).Inv(x.N),
	}
}
func (x Rat) InvR(a Rat) {
	x.N.Inv(a.N)
}
func (x Rat) Neg() Rat {
	return Rat{
		new(big.Rat).Neg(x.N),
	}
}
func (x Rat) NegR(a Rat) {
	x.N.Neg(a.N)
}
