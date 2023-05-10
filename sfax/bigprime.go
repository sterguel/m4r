package sfax

import "math/big"

type PrimeField struct {
	p *big.Int
}
type PrimeFieldElement struct {
	Field *PrimeField
	val   *big.Int
}
type bpel = PrimeFieldElement
