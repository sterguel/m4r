package sfax

func min[T uint64 | uint | int64 | int](a T, b T) T {
	if a < b {
		return a
	}
	return b
}
func max[T uint64 | uint | int64 | int](a T, b T) T {
	if a > b {
		return a
	}
	return b
}
func strim[T FieldElement[T]](P []T, z T) []T {
	l := len(P)
	m := 0
	for i := l - 1; i >= 0; i-- {
		if !P[i].Equals(z) {
			m = i + 1
			break
		}
	}
	return P[:m]
}
