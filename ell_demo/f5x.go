package main

import (
	"fmt"
	"scale/ell"
	"scale/sfax"
)

var F5 *sfax.SmallPrimeField = sfax.NewSmallPrimeField(5)
var F5x = sfax.NewPolynomialRing[*sfax.SmallPrimeFieldElement]("z", F5.Zero(), F5.One())

type fx = sfax.RatFunc[*sfax.SmallPrimeFieldElement]
type pel = sfax.SmallPrimeFieldElement

func F5x_demo() {
	A, B := F5x.FracElement([]*pel{F5.Element(3), F5.Element(0), F5.Element(1)}, []*pel{F5.One()}),
		F5x.FracElement([]*pel{F5.Element(0), F5.Element(4), F5.Element(0), F5.Element(3), F5.Element(3), F5.Element(2)}, []*pel{F5.One()})
	P := &ell.ECPoint[fx]{A, B, F5x.FracElement([]*pel{F5.Element(1), F5.Element(1), F5.Element(1)}, []*pel{F5.One()}),
		F5x.FracElement([]*pel{F5.Element(3), F5.Zero(), F5.Zero(), F5.One()}, []*pel{F5.One()}), false,
	}
	fmt.Printf("Curve Discriminant: %s\n", P.Discriminant().String())
	fmt.Printf("P: %s\n", P.String())
	fmt.Printf("Is P on the curve? %t\n", P.OnCurve())
	fmt.Printf("4P: %s\n", sfax.ZMul(P, 2))
	fmt.Println(sfax.ZMul(P, 2).OnCurve())
}
