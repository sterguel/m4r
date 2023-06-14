package sfax

import (
	"math/big"
	"math/bits"
)

func Power[T MulGroupElement[T]](x T, y int) T {
	if y < 0 {
		return Power(x.Inv(), -y)
	}
	u := x.Copy()
	res := x.One()
	for y > 0 {
		if (y & 1) != 0 {
			res.Mul(res, u)
		}
		y >>= 1
		u.Mul(u, u)
	}
	return res
}
func BigPower[T MulGroupElement[T]](x T, y *big.Int) T {
	if y.Cmp(bigzero) == -1 {
		return BigPower(x.Inv(), new(big.Int).Neg(y))
	}
	b := y.Bits()
	wc := len(b)
	res := x.One()
	if wc == 0 {
		return res
	}
	u := x.Copy()
	for _, v := range b[:wc-1] {
		k := uint(v)
		for i := 0; i < bits.UintSize; i++ {
			if (k & 1) != 0 {
				res.Mul(res, u)
			}
			k >>= 1
			u.Mul(u, u)
		}
	}
	k := uint(b[wc-1])
	for k > 0 {
		if (k & 1) != 0 {
			res.Mul(res, u)
		}
		k >>= 1
		u.Mul(u, u)
	}
	return res
}
func PosPower[T MulMonoidElement[T]](x T, y int) T {
	u := x.Copy()
	res := x.One()
	for y > 0 {
		if (y & 1) != 0 {
			res.Mul(res, u)
		}
		y >>= 1
		u.Mul(u, u)
	}
	return res
}

func ZMul[T AddGroupElement[T]](x T, y int) T {
	if y < 0 {
		return ZMul(x.Neg(), -y)
	}
	u := x.Copy()
	res := x.Zero()
	for y > 0 {
		if (y & 1) != 0 {
			res.Add(res, u)
		}
		y >>= 1
		u.Add(u, u)
	}
	return res
}
