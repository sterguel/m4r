package sfax

type RatFunc[T FieldElement[T]] struct {
	num *Polynomial[T]
	den *Polynomial[T]
}

func (x RatFunc[T]) FieldData() FieldData {
	return FieldData{
		x.num.ring.bone.FieldData().Char,
		0,
		nil,
	}
}
func (x RatFunc[T]) mcanon() {
	ld := len(x.den.val)
	fel := x.den.val[ld-1]
	if !fel.Equals(x.den.ring.bone) {
		for _, v := range x.num.val {
			v.DivR(v, fel)
		}
		for _, v := range x.den.val {
			v.DivR(v, fel)
		}
	}
}
func (x RatFunc[T]) gcanon() {
	if len(x.num.val) == 0 {
		x.den.val[0].Set(x.den.ring.bone)
		x.den.val = x.den.val[:1]
	}
	gcd := x.num.Zero()
	gcd.GCD(x.num, x.den)
	sc := x.num.Zero()
	if len(gcd.val) != 1 {
		sc.QR(x.num, x.num, gcd)
		sc.QR(x.den, x.den, gcd)
	}
}
func (x RatFunc[T]) canon() {
	x.gcanon()
	x.mcanon()
}

func (x RatFunc[T]) Set(y RatFunc[T]) {
	x.num.Set(y.num)
	x.den.Set(y.den)
}
func (x RatFunc[T]) SetFrac(num *Polynomial[T], den *Polynomial[T]) {
	x.num.Set(num)
	x.den.Set(den)
	x.canon()
}
func (x RatFunc[T]) Copy() RatFunc[T] {
	return RatFunc[T]{
		x.num.Copy(),
		x.den.Copy(),
	}
}
func (x RatFunc[T]) Equals(y RatFunc[T]) bool {
	return x.num.Equals(y.num) && x.den.Equals(y.den)
}
func (x RatFunc[T]) String() string {
	if len(x.den.val) > 1 {
		return x.num.String() + "/" + x.den.String()
	}
	return x.num.String()
}
func (x RatFunc[T]) Zero() RatFunc[T] {
	return RatFunc[T]{x.num.Zero(), x.num.One()}
}
func (x RatFunc[T]) One() RatFunc[T] {
	return RatFunc[T]{x.num.One(), x.num.One()}
}

func (a RatFunc[T]) Plus(b RatFunc[T]) RatFunc[T] {
	cn := a.num.Times(b.den)
	r := b.num.Times(a.den)
	cn.Add(cn, r)
	r.Mul(a.den, b.den)
	out := RatFunc[T]{cn, r}
	out.canon()
	return out
}
func (a RatFunc[T]) Minus(b RatFunc[T]) RatFunc[T] {
	cn := a.num.Times(b.den)
	r := b.num.Times(a.den)
	cn.Sub(cn, r)
	r.Mul(a.den, b.den)
	out := RatFunc[T]{cn, r}
	out.canon()
	return out
}
func (a RatFunc[T]) Times(b RatFunc[T]) RatFunc[T] {
	n := a.num.Times(b.num)
	d := a.den.Times(b.den)
	out := RatFunc[T]{n, d}
	out.canon()
	return out
}

func (a RatFunc[T]) Neg() RatFunc[T] {
	return RatFunc[T]{
		a.num.Neg(),
		a.den.Copy(),
	}
}
func (a RatFunc[T]) NegR(b RatFunc[T]) {
	a.num.NegR(b.num)
	a.den.Set(b.den)
}
func (x RatFunc[T]) Mul(a RatFunc[T], b RatFunc[T]) {
	x.num.Mul(a.num, b.num)
	x.den.Mul(a.den, b.den)
	x.canon()
}
func (a RatFunc[T]) Inv() RatFunc[T] {
	out := RatFunc[T]{a.den, a.num}
	out.mcanon()
	return out
}
func (a RatFunc[T]) InvR(b RatFunc[T]) {
	if &a.den.val[0] == &b.den.val[0] {
		a.Set(b.Inv())
		return
	}
	a.num.Set(b.den)
	a.den.Set(b.num)
	a.mcanon()
}
func (a RatFunc[T]) Div(b RatFunc[T]) RatFunc[T] {
	q := b.Inv()
	q.Mul(q, a)
	return q
}
func (x RatFunc[T]) DivR(a RatFunc[T], b RatFunc[T]) {
	if &x.den.val[0] == &a.den.val[0] {
		x.Set(a.Div(b))
		return
	}
	x.InvR(b)
	x.Mul(x, a)
}
func (x RatFunc[T]) Add(a RatFunc[T], b RatFunc[T]) {
	if &x.den.val[0] == &a.den.val[0] || &x.den.val[0] == &b.den.val[0] {
		x.Set(a.Plus(b))
		return
	}
	x.den.Mul(b.num, a.den)
	x.num.Mul(a.num, b.den)
	x.num.Add(x.num, x.den)
	x.den.Mul(a.den, b.den)
	x.canon()
}
func (x RatFunc[T]) Sub(a RatFunc[T], b RatFunc[T]) {
	if &x.den.val[0] == &a.den.val[0] || &x.den.val[0] == &b.den.val[0] {
		x.Set(a.Plus(b))
		return
	}
	x.den.Mul(b.num, a.den)
	x.num.Mul(a.num, b.den)
	x.num.Sub(x.num, x.den)
	x.den.Mul(a.den, b.den)
	x.canon()
}
