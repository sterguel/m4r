package sfax

type PolynomialRing[R Ring[R]] struct {
	BaseRing *R
}
type Polynomial[T RingElement[T]] []T

func (x *Polynomial[T]) Degree() int {
	return len(*x) - 1
}

func (x *Polynomial[T]) canon() {
	ar := *x
	w := len(*x)
	zero := ar[0].AddGroup().Zero()
	ex := 0
	for i := w - 1; i >= 0; i-- {
		if !ar[i].Equals(zero) {
			ex = i + 1
			break
		}
	}
	*x = ar[:ex]
}
