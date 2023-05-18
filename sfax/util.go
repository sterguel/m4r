package sfax

func min[T uint64 | uint | int64 | int](a T, b T) T {
	if a < b {
		return a
	}
	return b
}
