package sfax

import (
	"math/big"
)

var bigzero *big.Int = big.NewInt(0)
var bigtwo *big.Int = big.NewInt(2)
var bigone *big.Int = big.NewInt(1)

func IsSquareFin[T FieldElement[T]](x T) bool {
	rp := big.NewInt(1)
	rp.Sub(x.FieldData().Order, rp)
	rp.Div(rp, big.NewInt(2))
	return BigPower(x, rp).Equals(x.One())
}
func SqrtFin[T FieldElement[T]](x T, not_square T) T {
	d := x.FieldData()
	if d.Char.Cmp(bigzero) == 0 || d.Degree == 0 {
		panic("Field is infinite.")
	}
	if d.Char.Cmp(bigtwo) == 0 {
		y := x.Copy()
		for i := 0; i < d.Degree-1; i++ {
			y.Mul(y, y)
		}
		return y
	}
	zero := x.Zero()
	one := x.One()
	if zero.Equals(x) || !IsSquareFin(x) {
		return zero
	}
	Q := big.NewInt(1)
	Q.Sub(d.Order, Q)
	M := 0
	for Q.Bits()[0]&1 == 0 {
		M++
		Q.Rsh(Q, 1)
	}

	c := BigPower(not_square, Q)
	t := BigPower(x, Q)
	b := x.Zero()
	Q.Add(Q, bigone)
	Q.Rsh(Q, 1)
	R := BigPower(x, Q)
	for {
		if t.Equals(one) {
			return R
		}
		i := 0
		b.Set(t)
		for i = 1; i < M; i++ {
			b.Mul(b, b)
			if b.Equals(one) {
				break
			}
		}
		nsq := M - i - 1
		b.Set(c)
		for i := 0; i < nsq; i++ {
			b.Mul(b, b)
		}

		M = i
		c.Mul(b, b)
		t.Mul(t, c)
		R.Mul(R, b)
	}
}
