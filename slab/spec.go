package slab

import "matmul/sfax"

type Matrix[T sfax.FieldElement[T]] [][]T

type RatMatrix [][]sfax.Rat
