package slab

import (
	"math/big"
	"math/rand"
	"scale/sfax"
	"testing"
	"time"
)

func t_fillRationalMatrix(A RatMatrix, wl int) {
	rand.Seed(time.Now().UnixNano())
	for _, row := range A {
		for _, val := range row {
			bt := make([]big.Word, wl)
			bb := make([]big.Word, wl)
			for k := range bt {
				bb[k] = big.Word(rand.Uint64())
				bt[k] = big.Word(rand.Uint64())
			}
			if rand.Intn(2) == 0 {
				val.NegR(val)
			}
			top, bottom := new(big.Int).SetBits(bt), new(big.Int).SetBits(bb)
			val.N.SetFrac(top, bottom)
		}
	}
}
func TestRatMul(t *testing.T) {
	p, q, r := 40, 30, 50
	A := NewRatMatrix(p, q)
	B := NewRatMatrix(q, r)
	s := 3
	t_fillRationalMatrix(A, s)
	t_fillRationalMatrix(B, s)
	Ag, Bg := Matrix[sfax.Rat](A), Matrix[sfax.Rat](B)
	Cvout := Matrix[sfax.Rat](NewRatMatrix(p, r))
	Cvout.Mul_singlethreaded(Ag, Bg)
	Ctout := NewRatMatrix(p, r)
	Ctout.Mul_multimod(A, B)
	for i, row := range Ctout {
		for j, ent := range row {
			if !ent.Equals(Cvout[i][j]) {
				t.Error("Incorrect value!!")
				return
			}
		}
	}

}
