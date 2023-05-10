package sfax

type FiniteFieldExtension struct {
	p    uint64
	poly []uint64
}
type FiniteFieldExtensionElement struct {
	val   []uint64
	Field *FiniteFieldExtension
}

type ffx = FiniteFieldExtension
type fxl = FiniteFieldExtensionElement
