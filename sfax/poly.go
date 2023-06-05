package sfax

import (
	"strconv"
)

var default_poly_cap = 7 // degree 6
func SetDefaultPolyCapacity(n int) {
	default_poly_cap = n
}

type PolynomialRing[T FieldElement[T]] struct {
	bzero  T
	bone   T
	symbol string
}

type Polynomial[T FieldElement[T]] struct {
	ring *PolynomialRing[T]
	val  []T
}

func NewPolynomialRing[T FieldElement[T]](v string, zero T, one T) *PolynomialRing[T] {
	return &PolynomialRing[T]{
		zero,
		one,
		v,
	}
}

func (x *Polynomial[T]) canon() {
	ar := x.val
	ex := 0
	for i := len(ar) - 1; i >= 0; i-- {
		if !ar[i].Equals(x.ring.bzero) {
			ex = i + 1
			break
		}
	}
	x.val = ar[:ex]
}

func (R *PolynomialRing[T]) Element(x []T) *Polynomial[T] {
	u := &Polynomial[T]{
		R,
		x,
	}
	u.canon()
	return u
}

// math/big style resizing
// Do not let different polynomials share the same underlying array!
func (x *Polynomial[T]) resize(n int) {
	// reuse memory
	if n <= cap(x.val) {
		t := x.val[:n]
		x.val = t
		return
	}
	ln := max(default_poly_cap, n)
	r := make([]T, n, ln)
	rub := r[:ln]
	copy(r, x.val)
	for i := len(x.val); i < ln; i++ {
		rub[i] = x.ring.bzero.Copy()
	}
	x.val = r
}

func (R *PolynomialRing[T]) Zero() *Polynomial[T] {
	t := make([]T, 0)
	return &Polynomial[T]{R, t}
}
func (R *PolynomialRing[T]) One() *Polynomial[T] {
	t := make([]T, 1)
	t[0] = R.bone.Copy()
	return &Polynomial[T]{R, t}
}
func (x *Polynomial[T]) Zero() *Polynomial[T] {
	return x.ring.Zero()
}
func (x *Polynomial[T]) One() *Polynomial[T] {
	return x.ring.One()
}
func (x *Polynomial[T]) Degree() int {
	return len(x.val) - 1
}
func (x *Polynomial[T]) Equals(y *Polynomial[T]) bool {
	if len(x.val) != len(y.val) {
		return false
	}
	for i := 0; i < len(x.val) && i < len(y.val); i++ {
		if !x.val[i].Equals(y.val[i]) {
			return false
		}
	}
	return true
}
func (x *Polynomial[T]) Copy() *Polynomial[T] {
	v := make([]T, len(x.val))
	for i := range v {
		v[i] = x.val[i].Copy()
	}
	return &Polynomial[T]{x.ring, v}
}
func (x *Polynomial[T]) Set(y *Polynomial[T]) {
	x.resize(len(y.val))
	for i := 0; i < len(x.val) && i < len(y.val); i++ {
		x.val[i].Set(y.val[i])
	}
}
func (x *Polynomial[T]) String() string {
	a := ""
	for i, v := range x.val {
		if v.Equals(x.ring.bzero) {
			continue
		}
		a += v.String()
		if i != 0 {
			a += x.ring.symbol
			if i != 1 {
				a += "^" + strconv.Itoa(i)
			}
		}
		if i != len(x.val)-1 {
			a += "+"
		}
	}
	if a == "" {
		return "0"
	}
	return a
}
func (x *Polynomial[T]) Coefs() []T {
	return x.val
}

func (x *Polynomial[T]) Add(a *Polynomial[T], b *Polynomial[T]) {
	av, bv := a.val, b.val //avoid issues if x == a or b
	al, bl := len(av), len(bv)
	if al < bl {
		x.Add(b, a)
		return
	}
	x.resize(al)
	for i, u := range b.val {
		x.val[i].Add(a.val[i], u)
	}
	for i := bl; i < al; i++ {
		x.val[i].Set(a.val[i])
	}
	x.canon()
}
func (x *Polynomial[T]) Plus(y *Polynomial[T]) *Polynomial[T] {
	if len(x.val) < len(y.val) {
		return y.Plus(x)
	}
	l := len(x.val)
	z := make([]T, l)
	for i, u := range y.val {
		z[i] = u.Plus(x.val[i])
	}
	for i := len(y.val); i < l; i++ {
		z[i] = x.val[i].Copy()
	}
	o := &Polynomial[T]{x.ring, z}
	o.canon()
	return o
}
func (x *Polynomial[T]) Sub(a *Polynomial[T], b *Polynomial[T]) {
	av, bv := a.val, b.val
	al, bl := len(av), len(bv)
	if al < bl {
		x.resize(bl)
		for i, val := range av {
			x.val[i].Sub(val, bv[i])
		}
		for i := al; i < bl; i++ {
			x.val[i].NegR(bv[i])
		}
	} else {
		x.resize(al)
		for i, val := range bv {
			x.val[i].Sub(av[i], val)
		}
		for i := bl; i < al; i++ {
			x.val[i].Set(av[i])
		}
	}
	x.canon()
}
func (x *Polynomial[T]) Minus(y *Polynomial[T]) *Polynomial[T] {
	xl, yl := len(x.val), len(y.val)
	var z []T
	if xl < yl {
		z = make([]T, yl)
		for i, val := range x.val {
			z[i] = val.Minus(y.val[i])
		}
		for i := xl; i < yl; i++ {
			z[i] = y.val[i].Neg()
		}
		return &Polynomial[T]{x.ring, z}
	} else {
		z = make([]T, xl)
		for i, val := range y.val {
			z[i] = x.val[i].Minus(val)
		}
		for i := yl; i < xl; i++ {
			z[i] = x.val[i].Copy()
		}
	}
	o := &Polynomial[T]{x.ring, z}
	o.canon()
	return o
}
func (x *Polynomial[T]) Neg() *Polynomial[T] {
	z := make([]T, len(x.val))
	for i, v := range x.val {
		z[i] = v.Neg()
	}
	//canonicalisation unnecessary as x != 0 => -x != 0
	return &Polynomial[T]{x.ring, z}
}
func (x *Polynomial[T]) NegR(y *Polynomial[T]) {
	x.resize(len(y.val))
	for i, v := range y.val {
		x.val[i].NegR(v)
	}
}
func (x *Polynomial[T]) Mul(a *Polynomial[T], b *Polynomial[T]) {
	av, bv := a.val, b.val
	al, bl := len(av), len(bv)
	if al == 0 || bl == 0 {
		x.resize(0)
		return
	}
	x.resize(al + bl - 1)
	if &x.val[0] == &b.val[0] {
		if &x.val[0] == &a.val[0] {
			av = x.Copy().val[:al]
			bv = av
		} else {
			bv = x.Copy().val[:bl]
		}
	} else if &x.val[0] == &a.val[0] {
		av = x.Copy().val[:al]
	}
	s := a.ring.bzero.Copy()
	for _, v := range x.val {
		v.Set(a.ring.bzero)
	}
	for i, u := range av {
		for j, v := range bv {
			s.Mul(u, v)
			x.val[i+j].Add(x.val[i+j], s)
		}
	}
}

// Scalar multiplication (in-place)
func (x *Polynomial[T]) SMul(y T) {
	if y.Equals(x.ring.bzero) {
		x.val = x.val[:0]
		return
	}
	for _, v := range x.val {
		v.Mul(v, y)
	}
}
func (x *Polynomial[T]) Times(y *Polynomial[T]) *Polynomial[T] {
	av, bv := x.val, y.val
	al, bl := len(av), len(bv)
	if al == 0 || bl == 0 {
		return x.Zero()
	}
	z := make([]T, al+bl-1)
	for i := range z {
		z[i] = x.ring.bzero.Copy()
	}
	s := x.ring.bzero.Copy()
	for i, xl := range av {
		for j, yl := range bv {
			s.Mul(xl, yl)
			z[i+j].Add(z[i+j], s)
		}
	}
	//canon unnecessary as fields have no zero divisors
	return &Polynomial[T]{x.ring, z}
}

// Scalar multiplication (out-of-place)
func (x *Polynomial[T]) STimes(y T) *Polynomial[T] {
	z := x.Copy()
	z.SMul(y)
	return z
}

// Reduce a mod n and store the result in x. Does not change a or n.
func (x *Polynomial[T]) Mod(a *Polynomial[T], n *Polynomial[T]) {
	x.Set(a)
	xv := x.val
	nv := n.val
	lnv := len(n.val)
	fn := nv[lnv-1]
	sc := a.ring.bzero.Copy()
	cc := a.ring.bzero.Copy()
	for len(xv) >= lnv {
		d := len(xv)
		sc.DivR(xv[d-1], fn)
		xv = xv[:d-1] //chop off the last element as we don't need to compute a_n - a_n
		for i := 0; i < lnv-1; i++ {
			cc.Mul(nv[lnv-1-i], sc)
			xv[d-2-i].Sub(xv[d-2-i], cc)
		}
		xv = strim(xv, a.ring.bzero)
	}
	x.val = xv
}

// Computes a = qn + r such that deg(r) < deg(n)
func (r *Polynomial[T]) QR(q *Polynomial[T], a *Polynomial[T], n *Polynomial[T]) {
	r.Set(a)
	xv := r.val
	nv := n.val
	lnv := len(n.val)
	degdiff := len(xv) - lnv
	if degdiff < 0 {
		q.val = q.val[:0]
		return
	}
	q.resize(degdiff + 1)
	for _, v := range q.val {
		v.Set(a.ring.bzero)
	}
	fn := nv[lnv-1] //leading term
	cc := a.ring.bzero.Copy()
	for len(xv) >= lnv {
		d := len(xv)
		q.val[d-lnv].DivR(xv[d-1], fn)
		xv = xv[:d-1] //chop off the last element as we don't need to compute a_n - a_n
		for i := 0; i < lnv-1; i++ {
			cc.Mul(nv[lnv-1-i], q.val[d-lnv])
			xv[d-2-i].Sub(xv[d-2-i], cc)
		}
		xv = strim(xv, a.ring.bzero)
	}
	r.val = xv
}

// Destructive GCD - both a and b are destroyed and the GCD is returned.
func (a *Polynomial[T]) DGCD(b *Polynomial[T]) *Polynomial[T] {
	for len(b.val) != 0 {
		a.Mod(a, b)
		b, a = a, b
	}
	return a
}

// Extended Euclidean algorithm for polynomials - computes d, s, t such that as + bt = d.
func (d *Polynomial[T]) XGCD(s *Polynomial[T], t *Polynomial[T], a *Polynomial[T], b *Polynomial[T]) {

}
