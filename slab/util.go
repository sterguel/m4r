package slab

import "math/big"

func imax(a int, b int) int {
	if a > b {
		return a
	}
	return b
}
func imin(a int, b int) int {
	if a < b {
		return a
	}
	return b
}
func ceildiv(a int, b int) int {
	return 1 + ((a - 1) / b)
}
func reduceint(n *big.Int, m uint) uint {
	bits := n.Bits()
	var base uint = (1 << 32) % m
	base = (base * base) % m
	var acc uint = 0
	var mul uint = 1
	for _, x := range bits {
		acc = (acc + ((uint(x)%m)*mul)%m) % m
		mul = (mul * base) % m
	}
	if n.Sign() == -1 {
		return m - acc
	}
	return acc
}
