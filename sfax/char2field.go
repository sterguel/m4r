package sfax

import (
	"math/big"
	"math/bits"
	"strconv"
)

// Special case for char 2 because we can pack bits
// This can be used for non-irreducible polys as well as long as div is never called
// Maximum degree of 32
type Char2Field struct {
	Poly uint64
	Data FieldData
}

type Char2FieldElement struct {
	val   uint64
	Field *Char2Field
}
type el2 = Char2FieldElement

func c2mul(x uint64, y uint64) uint64 {
	var m, n, acc uint64
	if x > y {
		m = x
		n = y
	} else {
		n = x
		m = y
	}
	for n != 0 {
		if n&1 != 0 {
			acc ^= m
		}
		n >>= 1
		m <<= 1
	}
	return acc
}
func c2red(x uint64, p uint64, deg int) uint64 {
	var cmp uint64 = 1 << deg
	for x > cmp {
		diff := bits.Len64(x) - deg - 1
		x ^= (p << diff)
	}
	return x
}
func c2qr(x uint64, p uint64, deg int) (uint64, uint64) {
	var cmp uint64 = 1 << deg
	var q uint64
	for x > cmp {
		diff := bits.Len64(x) - deg - 1
		q ^= (1 << diff)
		x ^= (p << diff)
	}
	return q, x
}
func c2eeu(x uint64, y uint64) uint64 {
	r2, r := x, y
	var s2, s uint64 = 1, 0
	for r != 0 {
		q, rem := c2qr(r2, r, bits.Len64(r))
		r2, r = r, rem
		s2, s = s, s2^c2mul(q, s)
	}
	return s2
}
func NewChar2Field(poly uint64) *Char2Field {
	deg := bits.Len64(poly) - 1
	order := big.NewInt(1)
	order.Lsh(order, uint(deg))
	data := FieldData{
		big.NewInt(2),
		deg,
		order,
	}
	return &Char2Field{
		poly,
		data,
	}
}
func (F *Char2Field) Element(v uint64) *el2 {
	return &el2{
		v,
		F,
	}
}
func (x *el2) Copy() *el2 {
	return &el2{
		x.val,
		x.Field,
	}
}
func (x *el2) String() string {
	return strconv.FormatUint(x.val, 10)
}
func (x *el2) Equals(y *el2) bool {
	return x.val == y.val && x.Field == y.Field
}
func (x *el2) Set(y *el2) {
	x.val = y.val
}
func (x *el2) Zero() *el2 {
	return &el2{
		0,
		x.Field,
	}
}
func (x *el2) One() *el2 {
	return &el2{
		1,
		x.Field,
	}
}
func (x *el2) Plus(y *el2) *el2 {
	return &el2{
		x.val ^ y.val,
		x.Field,
	}
}
func (x *el2) Minus(y *el2) *el2 {
	return &el2{
		x.val ^ y.val,
		x.Field,
	}
}
func (x *el2) Add(a *el2, b *el2) {
	x.val = a.val ^ b.val
}
func (x *el2) Sub(a *el2, b *el2) {
	x.val = a.val ^ b.val
}
func (x *el2) Neg() *el2 {
	return x.Copy()
}
func (x *el2) NegR(y *el2) {
	x.Set(y)
}
func (x *el2) Times(y *el2) *el2 {
	m := c2mul(x.val, y.val)
	nv := c2red(m, x.Field.Poly, x.Field.Data.Degree)
	return &el2{
		nv,
		x.Field,
	}
}
func (x *el2) Mul(a *el2, b *el2) {
	m := c2mul(a.val, b.val)
	x.val = c2red(m, a.Field.Poly, a.Field.Data.Degree)
}
func (x *el2) Inv() *el2 {
	return &el2{
		c2eeu(x.val, x.Field.Poly),
		x.Field,
	}
}
func (x *el2) InvR(y *el2) {
	x.val = c2eeu(y.val, y.Field.Poly)
}
func (x *el2) Div(y *el2) *el2 {
	m := y.Inv()
	m.Mul(m, x)
	return m
}
func (x *el2) DivR(a *el2, b *el2) {
	x.InvR(b)
	x.Mul(a, x)
}
