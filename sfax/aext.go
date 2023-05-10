package sfax

type AlgExtField[T FieldElement[T]] struct {
	poly   []T
	Degree int
}

type AlgExtElement[T FieldElement[T]] struct {
	val   []T
	field *AlgExtField[T]
}

func (x *AlgExtField[T]) Element(r []T) *AlgExtElement[T] {
	if len(r) > x.Degree {
		return nil
	}
	nc := r
	if len(r) < x.Degree {
		if cap(r) < x.Degree {
			nc = make([]T, x.Degree)
			copy(nc, r)
		} else {
			nc = r[:x.Degree]
		}
	}
	return &AlgExtElement[T]{
		nc,
		x,
	}

}
