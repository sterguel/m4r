package sfax

import "math/big"

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
		y = y >> 1
		u.Mul(u, u)
	}
	return res
}
func BigPower[T MulGroupElement[T]](x T, y *big.Int) T {
	if y.Cmp(bigzero) == -1 {
		return BigPower(x.Inv(), new(big.Int).Neg(y))
	}
	u := x.Copy()
	res := x.One()
	for _, v := range y.Bits() {
		k := uint(v)
		for k > 0 {
			if (k & 1) != 0 {
				res.Mul(res, u)
			}
			k = k >> 1
			u.Mul(u, u)
		}
	}

	return res
}
