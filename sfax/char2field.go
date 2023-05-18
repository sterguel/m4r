package sfax

import (
	"math/bits"
)

// Special case for char 2 because we can pack bits
// This can be used for non-irreducible polys as well as long as div is never called
type Char2Field struct {
	Poly   uint64
	Order  int
	Degree int
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
func NewChar2Field(poly uint64) *Char2Field {
	deg := bits.Len64(poly) - 1
	order := (1 << deg)
	return &Char2Field{
		poly,
		order,
		deg,
	}
}
func (F *Char2Field) Element(v uint64) *el2 {
	return &el2{
		v,
		F,
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
func (x *el2) Times(y *el2) *el2 {
	m := c2mul(x.val, y.val)
	nv := c2red(m, x.Field.Poly, x.Field.Degree)
	return &el2{
		nv,
		x.Field,
	}
}
func (x *el2) Mul(a *el2, b *el2) {
	m := c2mul(a.val, b.val)
	x.val = c2red(m, a.Field.Poly, a.Field.Degree)
}
