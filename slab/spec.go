package slab

import "scale/sfax"

type Matrix[T sfax.FieldElement[T]] [][]T

type RatMatrix [][]sfax.Rat
