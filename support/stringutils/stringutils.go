package stringutils

// ReverseString returns s with its characters in reverse order. It operates on
// runes rather than bytes, so multi-byte UTF-8 characters are reversed as whole
// code points and remain valid. Note that grapheme clusters made of multiple
// combining code points are not kept together.
func ReverseString(s string) string {
	runes := []rune(s)
	for i, j := 0, len(runes)-1; i < j; i, j = i+1, j-1 {
		runes[i], runes[j] = runes[j], runes[i]
	}
	return string(runes)
}
