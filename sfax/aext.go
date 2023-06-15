package sfax

import (
	"math/big"
	"strconv"
)

// in-place reduction of a polynomial into low by the monic polynomial p
func aextred[T FieldElement[T]](low []T, high []T, p []T) {
	zero := low[0].Zero()
	scratch := zero.Copy()
	for d := len(high) - 1; d >= 0; d-- {
		lt := high[d] //leading term
		if lt.Equals(zero) {
			continue
		}
		pct := len(p) - 2 // (we don't care about doing the leading subtraction)
		ct := d - 1
		for (ct >= 0) && (pct >= 0) {
			scratch.Mul(lt, p[pct])
			high[ct].Sub(high[ct], scratch)
			pct -= 1
			ct -= 1
		}
		ct = len(p) - 2
		for (ct >= 0) && (pct >= 0) {
			scratch.Mul(lt, p[pct])
			low[ct].Sub(low[ct], scratch)
			pct -= 1
			ct -= 1
		}
	}
}

// in-place QR of polynomials
func aextqr[T FieldElement[T]](a []T, p []T) []T {
	la, lp := len(a), len(p)
	degdiff := la - lp
	if degdiff < 0 {
		return []T{}
	}
	q := make([]T, degdiff+1)
	for i := range q {
		dd := degdiff - i
		q[dd] = a[la-1-i].Div(p[lp-1])
		scratch := p[0].Zero()
		for j := dd; j < lp+dd; j++ {
			scratch.Mul(p[j-dd], q[dd])
			a[j].Sub(a[j], scratch)
		}
	}
	return q
}
func aexteeu[T FieldElement[T]](a []T, p []T, S []T) {
	//computes s, such that sA + tP = 1
	//destroys a, p
	z := p[0].Zero()
	one := p[0].One()
	r2, r := strim(a, z), p
	s2 := make([]T, len(S))
	s2[0] = z.One()
	s := S
	for i := range s2[1:] {
		s2[i+1] = z.Zero()
	}
	for len(r) != 0 {
		q := aextqr(r2, r)
		r2 = strim(r2, z)
		r, r2 = r2, r
		scratch := z.Zero()
		for i, v := range q {
			if v.Equals(z) {
				continue
			}
			for j, u := range s {
				if i+j >= len(s2) {
					break
				}
				if u.Equals(z) {
					continue
				}
				scratch.Mul(v, u)
				s2[i+j].Sub(s2[i+j], scratch)
			}
		}
		s, s2 = s2, s
	}
	if !r2[0].Equals(one) {
		for _, v := range s2 {
			v.DivR(v, r2[0])
		}
	}
	if &s2[0] != &S[0] {
		copy(S, s2)
	}

}

type AlgExtField[T FieldElement[T]] struct {
	poly   []T
	symbol string
	data   FieldData
	bzero  T
	bone   T
	Deg    int // this is NOT THE SAME as the total degree!!
}

type AlgExtElement[T FieldElement[T]] struct {
	field *AlgExtField[T]
	val   []T
}

/*
Create a new finite extension of any field given an irreducible polynomial P
and optionally a symbol s (used for conversion to string).
Represented as F[s]/(P(s)).
*/
func NewAlgebraicExtension[T FieldElement[T]](poly []T, symbol string) *AlgExtField[T] {
	//trim ending 0
	if len(poly) == 0 {
		return nil
	}
	if symbol == "" {
		symbol = "x"
	}
	zero, one := poly[0].Zero(), poly[0].One()
	basedata := poly[0].FieldData()
	// strip trailing 0s
	P := strim(poly, zero)
	deg := len(P) - 1
	// convert to monic polynomial
	if !P[deg].Equals(one) {
		i := P[deg].Inv()
		for _, v := range P {
			v.Mul(v, i)
		}
	}

	var order *big.Int = nil
	if basedata.Order != nil {
		order = new(big.Int).Exp(basedata.Order, big.NewInt(int64(deg)), nil)
	}
	var tdeg = basedata.Degree * deg
	fdata := FieldData{
		basedata.Char,
		tdeg,
		order,
	}
	return &AlgExtField[T]{
		P,
		symbol,
		fdata,
		zero,
		one,
		deg,
	}
}
func (x *AlgExtField[T]) Element(r []T) AlgExtElement[T] {
	nc := r
	if len(r) < x.Deg {
		if cap(r) < x.Deg {
			nc = make([]T, x.Deg)
			copy(nc, r)
		} else {
			nc = r[:x.Deg]
		}
	}
	return AlgExtElement[T]{
		x,
		nc,
	}
}
func (x *AlgExtField[T]) One() AlgExtElement[T] {
	o := make([]T, x.Deg)
	o[0] = x.bone.Copy()
	for i := 1; i < len(o); i++ {
		o[i] = x.bzero.Copy()
	}
	return AlgExtElement[T]{
		x,
		o,
	}
}
func (x *AlgExtField[T]) Zero() AlgExtElement[T] {
	o := make([]T, x.Deg)
	for i := range o {
		o[i] = x.bzero.Copy()
	}
	return AlgExtElement[T]{
		x,
		o,
	}
}
func (x AlgExtElement[T]) FieldData() FieldData {
	return x.field.data
}
func (x AlgExtElement[T]) Set(y AlgExtElement[T]) {
	for i := 0; i < len(x.val) && i < len(y.val); i++ {
		x.val[i].Set(y.val[i])
	}
}
func (x AlgExtElement[T]) One() AlgExtElement[T] {
	return x.field.One()
}
func (x AlgExtElement[T]) Zero() AlgExtElement[T] {
	return x.field.Zero()
}
func (x AlgExtElement[T]) Val() []T {
	return x.val
}
func (x AlgExtElement[T]) Copy() AlgExtElement[T] {
	nx := make([]T, x.field.Deg)
	for i, v := range x.val {
		nx[i] = v.Copy()
	}
	return AlgExtElement[T]{x.field, nx}
}
func (x AlgExtElement[T]) Equals(y AlgExtElement[T]) bool {
	for i, v := range x.val {
		if !v.Equals(y.val[i]) {
			return false
		}
	}
	return x.field == y.field
}
func (x AlgExtElement[T]) String() string {
	a := ""
	for i, v := range x.val {
		if v.Equals(x.field.bzero) {
			continue
		}
		a += v.String()
		if i != 0 {
			a += x.field.symbol
			if i != 1 {
				a += "^" + strconv.Itoa(i)
			}
		}
		a += "+"
	}

	if a == "" {
		return "0"
	}
	if a[len(a)-1] == '+' {
		return string(a[:len(a)-1])
	}
	return a
}
func (x AlgExtElement[T]) Plus(y AlgExtElement[T]) AlgExtElement[T] {
	o := make([]T, x.field.Deg)
	for i, xv := range x.val {
		o[i] = xv.Plus(y.val[i])
	}
	return AlgExtElement[T]{
		x.field,
		o,
	}
}
func (x AlgExtElement[T]) Add(a AlgExtElement[T], b AlgExtElement[T]) {
	for i, u := range x.val {
		u.Add(a.val[i], b.val[i])
	}
}
func (x AlgExtElement[T]) Minus(y AlgExtElement[T]) AlgExtElement[T] {
	o := make([]T, x.field.Deg)
	for i, xv := range x.val {
		o[i] = xv.Minus(y.val[i])
	}
	return AlgExtElement[T]{
		x.field,
		o,
	}
}
func (x AlgExtElement[T]) Sub(a AlgExtElement[T], b AlgExtElement[T]) {
	for i, u := range x.val {
		u.Sub(a.val[i], b.val[i])
	}
}
func (x AlgExtElement[T]) Neg() AlgExtElement[T] {
	nx := make([]T, x.field.Deg)
	for i, v := range x.val {
		nx[i] = v.Neg()
	}
	return AlgExtElement[T]{x.field, nx}
}
func (x AlgExtElement[T]) NegR(y AlgExtElement[T]) {
	for i, v := range x.val {
		v.NegR(y.val[i])
	}
}
func (x AlgExtElement[T]) Times(y AlgExtElement[T]) AlgExtElement[T] {
	ui := make([]T, 2*x.field.Deg-1)
	o := ui[:x.field.Deg]
	high := ui[x.field.Deg:]
	scratch := x.field.bzero.Copy()
	for i := range o {
		acc := x.field.bzero.Copy()
		for k := 0; k < i+1; k++ {
			scratch.Mul(x.val[i-k], y.val[k])
			acc.Add(acc, scratch)
		}
		o[i] = acc
	}

	for i := range high {
		acc := x.field.bzero.Copy()
		for k := 0; k < (x.field.Deg - i - 1); k++ {
			scratch.Mul(x.val[x.field.Deg-k-1], y.val[i+k+1])
			acc.Add(acc, scratch)
		}
		high[i] = acc
	}
	aextred(o, high, x.field.poly)
	return AlgExtElement[T]{x.field, o}
}
func (x AlgExtElement[T]) Mul(a AlgExtElement[T], b AlgExtElement[T]) {
	if &x.val[0] == &b.val[0] || &x.val[0] == &a.val[0] {
		x.Set(a.Times(b))
		return
	}
	high := make([]T, a.field.Deg-1)

	scratch := a.field.bzero.Copy()
	for i := range x.val {
		x.val[i].Set(x.field.bzero)
		for k := 0; k < i+1; k++ {
			scratch.Mul(a.val[i-k], b.val[k])
			x.val[i].Add(x.val[i], scratch)
		}
	}

	for i := range high {
		acc := a.field.bzero.Copy()
		for k := 0; k < (a.field.Deg - i - 1); k++ {
			scratch.Mul(a.val[a.field.Deg-k-1], b.val[i+k+1])
			acc.Add(acc, scratch)
		}
		high[i] = acc
	}
	aextred(x.val, high, a.field.poly)
}
func (x AlgExtElement[T]) Inv() AlgExtElement[T] {
	z := x.field.bzero.Copy()
	p := make([]T, x.field.Deg+1)
	u := make([]T, x.field.Deg)
	for i, v := range x.field.poly {
		p[i] = v.Copy()
	}
	for i := range u {
		u[i] = z.Zero()
	}
	y := x.Copy().val
	aexteeu(y, p, u)
	return AlgExtElement[T]{
		x.field,
		u,
	}

}
func (x AlgExtElement[T]) InvR(y AlgExtElement[T]) {
	p := make([]T, x.field.Deg+1)
	for i, v := range x.field.poly {
		p[i] = v.Copy()
	}
	yc := y.Copy().val
	for _, v := range x.val {
		v.Set(x.field.bzero)
	}
	aexteeu(yc, p, x.val)
}
func (x AlgExtElement[T]) Div(y AlgExtElement[T]) AlgExtElement[T] {
	z := y.Inv()
	z.Mul(z, x)
	return z
}
func (x AlgExtElement[T]) DivR(a AlgExtElement[T], b AlgExtElement[T]) {
	if &x.val[0] == &a.val[0] {
		x.Set(a.Div(b))
		return
	}
	x.InvR(b)
	x.Mul(x, a)
}
