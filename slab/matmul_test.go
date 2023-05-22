package slab

import (
	"math/big"
	"math/rand"
	"testing"
	"time"
)

func fillMatrix(A *IntMatrix, wl int) {
	rand.Seed(time.Now().UnixNano())
	for _, row := range A.Vals {
		for _, val := range row {
			b := make([]big.Word, wl)
			for k := range b {
				b[k] = big.Word(rand.Uint64())
			}
			if rand.Intn(2) == 0 {
				val.Neg(val)
			}
			val.SetBits(b)
		}
	}
}

func TestMul(t *testing.T) {
	A := Make_Int_matrix(400, 300)
	B := Make_Int_matrix(300, 450)
	s := 4
	fillMatrix(A, s)
	fillMatrix(B, s)
	Cv, _ := Mul_standard(A, B)
	Ct, _ := Mul_modular_montgomery(A, B)
	for i, row := range Cv.Vals {
		for j, ent := range row {
			if ent.Cmp(Ct.Vals[i][j]) != 0 {
				t.Error("Incorrect value!!")
				return
			}
		}
	}

}

func TestJobsMul(t *testing.T) {
	A := Make_Int_matrix(400, 300)
	B := Make_Int_matrix(300, 450)
	s := 4
	fillMatrix(A, s)
	fillMatrix(B, s)
	Cv, _ := Mul_modular_montgomery(A, B)
	Ct, _ := Mul_modular_montgomery_jobs(A, B)
	for i, row := range Cv.Vals {
		for j, ent := range row {
			if ent.Cmp(Ct.Vals[i][j]) != 0 {
				t.Error("Incorrect value!!")
				return
			}
		}
	}
}
