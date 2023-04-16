package matmul

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
