package sfax

type PolynomialRing[T FieldElement[T]] struct {
	bzero T
	bone  T
}
type Polynomial[T FieldElement[T]] struct {
	ring *PolynomialRing[T]
	val  []T
}

func NewPolynomialRing[T FieldElement[T]](zero T, one T) *PolynomialRing[T] {
	return &PolynomialRing[T]{
		zero,
		one,
	}
}
func (x *Polynomial[T]) Degree() int {
	return len(x.val) - 1
}

func (x *Polynomial[T]) canon() {
	ar := x.val
	zero := x.ring.bzero
	ex := 0
	for i := len(ar) - 1; i >= 0; i-- {
		if !ar[i].Equals(zero) {
			ex = i + 1
			break
		}
	}
	x.val = ar[:ex]
}
