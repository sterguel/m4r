package ell

import (
	"scale/sfax"
)

// Point on a curve of the form y^2 = x^3 + Ax + B
type ECPoint[T sfax.FieldElement[T]] struct {
	A   T
	B   T
	X   T
	Y   T
	Inf bool
}

func (x *ECPoint[T]) Discriminant() T {
	l := x.A.Times(x.A)
	l.Mul(l, x.A)
	l = sfax.ZMul(l, 4)
	r := x.B.Times(x.B)
	r = sfax.ZMul(r, 27)
	l.Add(l, r)
	return l
}
func (p *ECPoint[T]) OnCurve() bool {
	c := p.X.Times(p.X)
	c.Mul(c, p.X)
	c.Add(c, p.X.Times(p.A))
	c.Add(c, p.B)
	d := p.Y.Times(p.Y)
	return d.Equals(c)
}
func (x *ECPoint[T]) Copy() *ECPoint[T] {
	return &ECPoint[T]{x.A, x.B, x.X.Copy(), x.Y.Copy(), x.Inf}
}
func (x *ECPoint[T]) String() string {
	if x.Inf {
		return "Point at Infinity"
	}
	return "(" + x.X.String() + "," + x.Y.String() + ")"
}
func (x *ECPoint[T]) Set(y *ECPoint[T]) {
	x.Inf = y.Inf
	if x.Inf {
		return
	}
	x.X.Set(y.X)
	x.Y.Set(y.Y)
}
func (x *ECPoint[T]) Equals(y *ECPoint[T]) bool {
	return x.Inf == y.Inf && x.X.Equals(y.X) && x.Y.Equals(y.Y)
}
func (x *ECPoint[T]) Zero() *ECPoint[T] {
	return &ECPoint[T]{x.A, x.B, x.A.Zero(), x.B.Zero(), true}
}
func (x *ECPoint[T]) Neg() *ECPoint[T] {
	return &ECPoint[T]{x.A, x.B, x.X, x.Y.Neg(), x.Inf}
}
func (x *ECPoint[T]) NegR(y *ECPoint[T]) {
	x.Inf = y.Inf
	if y.Inf {
		return
	}
	x.X.Set(y.X)
	x.Y.NegR(y.Y)
}
func (p *ECPoint[T]) Plus(q *ECPoint[T]) *ECPoint[T] {
	if p.Inf {
		return q.Copy()
	}
	if q.Inf {
		return p.Copy()
	}
	x3, y3, m, c := p.A.Zero(), p.A.Zero(), p.A.Zero(), p.A.Zero()
	if p.X.Equals(q.X) {
		if p.Y.Equals(q.Y) {
			c.Mul(p.X, p.X)
			m.Add(c, c)
			m.Add(m, c)
			m.Add(m, p.A)
			y3.Add(p.Y, p.Y)
			m.DivR(m, y3)
		} else {
			return p.Zero()
		}
	} else {
		m.Sub(p.Y, q.Y)
		y3.Sub(p.X, q.X)
		m.DivR(m, y3)
	}
	c.Mul(p.X, m)
	c.Sub(p.Y, c)
	x3.Mul(m, m)
	x3.Sub(x3, p.X)
	x3.Sub(x3, q.X)
	y3.Mul(m, x3)
	y3.Add(y3, c)
	y3.NegR(y3)
	return &ECPoint[T]{p.A, p.B, x3, y3, false}
}

func (z *ECPoint[T]) Add(p *ECPoint[T], q *ECPoint[T]) {
	if p.Inf {
		z.Set(q)
		return
	}
	if q.Inf {
		z.Set(p)
		return
	}
	s, m, c := p.A.Zero(), p.A.Zero(), p.A.Zero()
	if p.X.Equals(q.X) {
		if p.Y.Equals(q.Y) {
			c.Mul(p.X, p.X)
			m.Add(c, c)
			m.Add(m, c)
			m.Add(m, p.A)
			s.Add(p.Y, p.Y)
			m.DivR(m, s)
		} else {
			z.Inf = true
			return
		}
	} else {
		m.Sub(p.Y, q.Y)
		s.Sub(p.X, q.X)
		m.DivR(m, s)
	}
	c.Mul(p.X, m)
	c.Sub(p.Y, c)
	s.Mul(m, m)
	s.Sub(s, p.X)
	z.X.Sub(s, q.X)
	z.Y.Mul(m, z.X)
	z.Y.Add(z.Y, c)
	z.Y.NegR(z.Y)
	z.Inf = false
}

func (p *ECPoint[T]) Minus(q *ECPoint[T]) *ECPoint[T] {
	o := q.Neg()
	o.Add(o, p)
	return o
}
func (z *ECPoint[T]) Sub(p *ECPoint[T], q *ECPoint[T]) {
	if &z == &p {
		z.Set(p.Minus(q))
		return
	}
	z.NegR(q)
	z.Add(z, p)
}
