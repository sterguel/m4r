package sfax

import (
	"math/big"
	"math/bits"
)

const wordsize = bits.UintSize

var wordsizeBig *big.Int = big.NewInt(wordsize)

// In place REDC for big ints
func bigredc(X *big.Int, p *big.Int, pinv *big.Int, wl int) {
	t := new(big.Int).Set(X)
	wx := t.Bits()
	t.SetBits(wx[:min(wl, len(wx))]) //mod R reduction
	t.Mul(t, pinv)
	wx = t.Bits()
	t.SetBits(wx[:min(wl, len(wx))])
	t.Mul(t, p)
	X.Add(X, t)
	wx = X.Bits()
	X.SetBits(wx[wl:])
	if X.CmpAbs(p) != -1 {
		X.Sub(X, p)
	}
}

/*
Prime fields for big p (>2^32).
Do not use for smaller p as the Montgomery reduction step is much less efficient!
*/
type PrimeField struct {
	p    *big.Int
	data FieldData
	wl   int
	r    *big.Int
	r2   *big.Int
	r3   *big.Int
	pinv *big.Int
	rinv *big.Int
}
type PrimeFieldElement struct {
	Field *PrimeField
	val   *big.Int
}
type bpel = PrimeFieldElement

func NewPrimeField(p *big.Int) *PrimeField {
	d := FieldData{p, 1, p}
	wordlength := len(p.Bits())
	wlb := big.NewInt(int64(wordlength))
	wlb.Mul(wlb, wordsizeBig)
	r := big.NewInt(2)
	r.Exp(r, wlb, nil)
	rmodn := new(big.Int).Mod(r, p)
	r2 := new(big.Int).Mul(rmodn, rmodn)
	r2.Mod(r2, p)
	r3 := new(big.Int).Mul(rmodn, r2)
	r3.Mod(r3, p)
	rinv := new(big.Int)
	pinv := new(big.Int)
	new(big.Int).GCD(rinv, pinv, r, p)
	rinv.Mod(rinv, p)
	pinv.Mod(pinv, r)
	pinv.Sub(r, pinv)
	return &PrimeField{
		p,
		d,
		wordlength,
		rmodn,
		r2,
		r3,
		pinv,
		rinv,
	}
}

func (F *PrimeField) Element(n *big.Int) bpel {
	v := new(big.Int)
	v.Mul(F.r2, n)
	bigredc(v, F.p, F.pinv, F.wl)
	return bpel{
		F,
		v,
	}
}
func (F *PrimeField) Zero() bpel {
	return bpel{
		F,
		big.NewInt(0),
	}
}
func (F *PrimeField) One() bpel {
	return bpel{
		F,
		new(big.Int).Set(F.r),
	}
}
func (x bpel) Set(y bpel) {
	x.val.Set(y.val)
}
func (x bpel) FieldData() FieldData {
	return x.Field.data
}
func (x bpel) One() bpel {
	return bpel{
		x.Field,
		new(big.Int).Set(x.Field.r),
	}
}
func (x bpel) Zero() bpel {
	return bpel{
		x.Field,
		big.NewInt(0),
	}
}
func (x bpel) Copy() bpel {
	return bpel{
		x.Field,
		new(big.Int).Set(x.val),
	}
}
func (x bpel) Val() *big.Int {
	a := new(big.Int).Set(x.val)
	bigredc(a, x.Field.p, x.Field.pinv, x.Field.wl)
	return a
}
func (x bpel) RawVal() *big.Int {
	return x.val
}
func (x bpel) Equals(y bpel) bool {
	//cmpabs to avoid sign checking
	return (x.Field == y.Field) && (x.val.CmpAbs(y.val) == 0)
}
func (x bpel) String() string {
	return x.Val().String()
}
func (x bpel) Neg() bpel {
	return bpel{
		x.Field,
		new(big.Int).Sub(x.Field.p, x.val),
	}
}
func (x bpel) NegR(y bpel) {
	x.val.Sub(y.Field.p, y.val)
}
func (x bpel) Add(a bpel, b bpel) {
	x.val.Add(a.val, b.val)
	if x.val.CmpAbs(a.Field.p) != -1 {
		x.val.Sub(x.val, a.Field.p)
	}
}
func (a bpel) Plus(b bpel) bpel {
	t := new(big.Int).Add(a.val, b.val)
	if t.CmpAbs(a.Field.p) != -1 {
		t.Sub(t, a.Field.p)
	}
	return bpel{
		a.Field,
		t,
	}
}
func (x bpel) Sub(a bpel, b bpel) {
	if a.val.CmpAbs(b.val) != 1 {
		x.val.Sub(a.val, b.val)
	} else {
		x.val.Sub(a.Field.p, b.val)
		x.val.Add(x.val, a.val)
	}
}
func (a bpel) Minus(b bpel) bpel {
	t := new(big.Int)
	if a.val.CmpAbs(b.val) != -1 {
		t.Sub(a.val, b.val)
	} else {
		t.Sub(a.Field.p, b.val)
		t.Add(t, a.val)
	}
	return bpel{
		a.Field,
		t,
	}
}
func (x bpel) Mul(a bpel, b bpel) {
	x.val.Mul(a.val, b.val)
	bigredc(x.val, a.Field.p, a.Field.pinv, a.Field.wl)
}
func (a bpel) Times(b bpel) bpel {
	t := new(big.Int)
	t.Mul(a.val, b.val)
	bigredc(t, a.Field.p, a.Field.pinv, a.Field.wl)
	return bpel{a.Field, t}
}
func (x bpel) Inv() bpel {
	o := new(big.Int)
	new(big.Int).GCD(o, nil, x.val, x.Field.p)
	o.Mod(o, x.Field.p)
	o.Mul(o, x.Field.r3)
	bigredc(o, x.Field.p, x.Field.pinv, x.Field.wl)
	return bpel{
		x.Field,
		o,
	}
}
func (x bpel) InvR(y bpel) {
	new(big.Int).GCD(x.val, nil, y.val, y.Field.p)
	x.val.Mod(x.val, y.Field.p)
	x.val.Mul(x.val, y.Field.r3)
	bigredc(x.val, y.Field.p, y.Field.pinv, y.Field.wl)
}
func (a bpel) Div(b bpel) bpel {
	o := b.Inv()
	o.Mul(o, a)
	return o
}
func (x bpel) DivR(a bpel, b bpel) {
	x.InvR(b)
	x.Mul(x, a)
}
