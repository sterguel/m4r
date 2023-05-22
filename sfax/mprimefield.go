package sfax

import (
	"math/big"
	"strconv"
)

func inveeu(a int64, p int64) int64 {
	var s1 int64 = 0
	var s2 int64 = 1
	r1 := p
	r2 := a
	for r1 != 0 {
		q := r2 / r1
		r2, r1 = r1, r2-q*r1
		s2, s1 = s1, s2-q*s1
	}
	return s2
}
func redadd(a uint64, n uint64) uint64 {
	if a >= n {
		return a - n
	}
	return a
}

const redcmask uint64 = (1 << 32) - 1

func redc(x uint64, n uint64, ninv uint64) uint64 {
	m := ((x & redcmask) * ninv) & redcmask
	t := (x + m*n) >> 32
	if t >= n {
		return t - n
	}
	return t
}

// Prime field for small p (<2^32)
type SmallPrimeField struct {
	p    uint64
	r3   uint64
	r2   uint64
	r    uint64
	pinv uint64
	data FieldData
}
type SmallPrimeFieldElement struct {
	Field *SmallPrimeField
	val   uint64
}
type pel = SmallPrimeFieldElement

func NewSmallPrimeField(p uint64) *SmallPrimeField {
	d := FieldData{big.NewInt(int64(p)), 1, big.NewInt(int64(p))}
	r := (1 << 32) % p
	r2 := (r * r) % p
	r3 := (r2 * r) % p
	ninv := (1 << 32) - uint64(inveeu(int64(p), 1<<32))
	return &SmallPrimeField{p: p, r2: r2, r3: r3, data: d, r: r, pinv: ninv}
}
func (x *SmallPrimeField) Element(n uint64) *pel {
	return &pel{x, redc(n*x.r2, x.p, x.pinv)}
}
func (x *pel) Set(y *pel) {
	x.val = y.val
}
func (x *pel) One() *pel {
	return &pel{x.Field, x.Field.r}
}
func (x *pel) Zero() *pel {
	return &pel{x.Field, 0}
}

func (x *pel) Plus(y *pel) *pel {
	return &pel{x.Field, redadd(x.val+y.val, x.Field.p)}
}
func (x *pel) Add(a *pel, b *pel) {
	x.val = redadd(a.val+b.val, a.Field.p)
}
func (x *pel) Minus(y *pel) *pel {
	var v uint64
	if x.val >= y.val {
		v = x.val - y.val
	} else {
		v = x.Field.p + x.val - y.val
	}
	return &pel{x.Field, v}
}
func (x *pel) Sub(a *pel, b *pel) {
	if a.val >= b.val {
		x.val = a.val - b.val
	} else {
		x.val = a.Field.p + a.val - b.val
	}
}
func (x *pel) Times(y *pel) *pel {
	return &pel{x.Field, redc(x.val*y.val, x.Field.p, x.Field.pinv)}
}
func (x *pel) Mul(a *pel, b *pel) {

	x.val = redc(a.val*b.val, a.Field.p, a.Field.pinv)
}
func (x *pel) Inv() *pel {
	i := inveeu(int64(x.val), int64(x.Field.p))
	return &pel{x.Field, redc(uint64(i)*x.Field.r3, x.Field.p, x.Field.pinv)}
}
func (x *pel) InvR(a *pel) {
	i := inveeu(int64(a.val), int64(a.Field.p))
	x.val = redc(uint64(i)*a.Field.r3, a.Field.p, a.Field.pinv)
}
func (x *pel) Neg() *pel {
	return &pel{x.Field, x.Field.p - x.val}
}
func (x *pel) NegR(a *pel) {
	x.val = a.Field.p - a.val
}
func (x *pel) Copy() *pel {
	return &pel{x.Field, x.val}
}
func (x *pel) Div(y *pel) *pel {
	k := y.Inv()
	k.Mul(k, x)
	return k
}
func (x *pel) DivR(a *pel, b *pel) {
	x.InvR(b)
	x.Mul(x, a)
}
func (x *pel) Val() uint64 {
	return redc(x.val, x.Field.p, x.Field.pinv)
}
func (x *pel) RawVal() uint64 {
	return x.val
}
func (x *pel) Equals(y *pel) bool {
	return (x.val == y.val) && (x.Field == y.Field)
}
func (x *pel) FieldData() FieldData {
	return x.Field.data
}
func (x *pel) String() string {
	return strconv.FormatUint(x.Val(), 10)
}
