package sfax

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
	if a > n {
		return a - n
	}
	return a
}

type SmallPrimeField struct {
	p uint64
}
type SmallPrimeFieldElement struct {
	Field *SmallPrimeField
	val   uint64
}
type pel = SmallPrimeFieldElement

func NewPrimeField(p uint64) *SmallPrimeField {
	return &SmallPrimeField{p}
}
func (x *SmallPrimeField) Element(n uint64) *pel {
	return &pel{x, n}
}

func (x *SmallPrimeField) One() *pel {
	return &pel{x, 1}
}
func (x *SmallPrimeField) Zero() *pel {
	return &pel{x, 0}
}
func (x *SmallPrimeField) Char() int64 {
	return int64(x.p)
}

func (x *pel) Plus(y *pel) *pel {
	return &pel{x.Field, redadd(x.val+y.val, x.Field.p)}
}
func (x *pel) Add(a *pel, b *pel) {
	x.val = redadd(a.val+b.val, a.Field.p)
}
func (x *pel) Minus(y *pel) *pel {
	var v uint64
	if x.val > y.val {
		v = x.val - y.val
	} else {
		v = x.Field.p - x.val - y.val
	}
	return &pel{x.Field, v}
}
func (x *pel) Sub(a *pel, b *pel) {
	if a.val > b.val {
		x.val = a.val - b.val
	} else {
		x.val = a.Field.p - a.val - b.val
	}
}
func (x *pel) Times(y *pel) *pel {
	return &pel{x.Field, (x.val * y.val) % x.Field.p}
}
func (x *pel) Mul(a *pel, b *pel) {
	x.val = (a.val * b.val) % a.Field.p
}
func (x *pel) Inv() *pel {
	i := inveeu(int64(x.val), int64(x.Field.p))
	return &pel{x.Field, uint64(i) % x.Field.p}
}
func (x *pel) InvR(a *pel) {
	x.val = uint64(inveeu(int64(a.val), int64(a.Field.p)))
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
	return x.val
}
func (x *pel) Equals(y *pel) bool {
	return (x.val == y.val) && (x.Field == y.Field)
}
func (x *pel) MulMonoid() *SmallPrimeField {
	return x.Field
}
func (x *pel) AddGroup() *SmallPrimeField {
	return x.Field
}
