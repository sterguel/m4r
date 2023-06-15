package main

import (
	"fmt"
	"math/big"
	"scale/ell"
	"scale/sfax"
)

type ft = sfax.AlgExtElement[*sfax.SmallPrimeFieldElement]
type qi = sfax.AlgExtElement[sfax.Rat]

var Fp *sfax.SmallPrimeField = sfax.NewSmallPrimeField(833191)
var polyF []*sfax.SmallPrimeFieldElement = []*sfax.SmallPrimeFieldElement{
	Fp.Element(833185), Fp.Element(833175), Fp.Element(1),
}
var F *sfax.AlgExtField[*sfax.SmallPrimeFieldElement] = sfax.NewAlgebraicExtension(polyF, "a")
var polyQi []sfax.Rat = []sfax.Rat{
	{N: big.NewRat(1, 1)}, {N: big.NewRat(0, 1)}, {N: big.NewRat(1, 1)},
}
var QQi *sfax.AlgExtField[sfax.Rat] = sfax.NewAlgebraicExtension(polyQi, "i")

func pos_mod(a int, p int) int {
	n := a % p
	if n < 0 {
		return n + p
	}
	return n
}
func Qi(a int64, b int64, c int64, d int64) qi {
	r1 := sfax.Rat{N: big.NewRat(a, b)}
	r2 := sfax.Rat{N: big.NewRat(c, d)}
	return QQi.Element([]sfax.Rat{r1, r2})
}
func Fel(a int, b int) ft {
	au := uint64(pos_mod(a, 833191))
	bu := uint64(pos_mod(b, 833191))
	return F.Element([]*sfax.SmallPrimeFieldElement{Fp.Element(au), Fp.Element(bu)})
}

func F_demo() {
	A := Fel(420, 69)
	B := Fel(61016, 60)
	P := &ell.ECPoint[ft]{A, B, Fel(253579, 0), Fel(388595, 0), false}
	fmt.Printf("P: %s\n", P.String())
	fmt.Printf("Curve Discriminant: %s\n", P.Discriminant().String())
	fmt.Printf("Is P on the curve? %t\n", P.OnCurve())
	fmt.Printf("123456789P: %s\n", sfax.ZMul(P, 131414).String())
	fmt.Printf("216P: %s\n\n", sfax.ZMul(P, 216).String())
	S := &ell.ECPoint[ft]{A, B, Fel(251571, 821099), Fel(286300, 34063), false}
	fmt.Printf("S: %s\n", S.String())
	fmt.Printf("Is S on the curve? %t\n", S.OnCurve())
	fmt.Printf("108S: %s\n\n", sfax.ZMul(S, 108).String())
	R := P.Plus(S)
	fmt.Printf("R = P+S: %s\n", R.String())
	fmt.Printf("216R: %s\n\n", sfax.ZMul(R, 216).String())

}
func Q_demo() {
	A, B := sfax.Rat{big.NewRat(-756, 1)}, sfax.Rat{big.NewRat(4320, 1)}
	P := &ell.ECPoint[sfax.Rat]{A, B, sfax.Rat{big.NewRat(-3, 1)}, sfax.Rat{big.NewRat(81, 1)}, false}
	fmt.Printf("Curve Discriminant: %s\n", P.Discriminant().String())
	fmt.Printf("P: %s\n", P.String())
	fmt.Printf("Is P on the curve? %t\n", P.OnCurve())
	fmt.Printf("4P: %s\n", sfax.ZMul(P, 4))
	T := &ell.ECPoint[sfax.Rat]{A, B, sfax.Rat{big.NewRat(24, 1)}, sfax.Rat{big.NewRat(0, 1)}, false}
	fmt.Printf("T: %s\n", T.String())
	fmt.Printf("Is T on the curve? %t\n", T.OnCurve())
	fmt.Printf("2T: %s\n", sfax.ZMul(T, 2))
	Q := P.Plus(T)
	fmt.Printf("Q = P+T : %s \n", Q.String())
	fmt.Printf("4Q: %s \n\n", sfax.ZMul(Q, 4))
	//fmt.Printf("400P: %s\n\n", sfax.ZMul(P, 400))
}
func Qi_demo() {
	A, B := Qi(-3, 1, -8, 1), Qi(-34, 1, -32, 1)
	P := &ell.ECPoint[qi]{A, B, Qi(7, 1, 2, 1), Qi(16, 1, 6, 1), false}
	fmt.Printf("Curve Discriminant: %s\n", P.Discriminant().String())
	fmt.Printf("P: %s\n", P.String())
	fmt.Printf("Is P on the curve? %t\n", P.OnCurve())
	fmt.Printf("4P: %s\n\n", sfax.ZMul(P, 4))
}

func main() {
	fmt.Println("F = F_833191(a) where a^2 = 16a + 6, y^2 = x^3 + (420+69a)x + (61016+60a)\n")
	F_demo()
	fmt.Println("F = Q, y^2 = x^3 - 756x + 4320\n")
	Q_demo()
	fmt.Println("F = Q[i], y^2 = x^3 - (3+8i)x -(34 + 32i)\n")
	Qi_demo()

}
