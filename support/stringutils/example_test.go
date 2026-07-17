package stringutils_test

import (
	"fmt"

	"github.com/compgenlab/cghts/support/stringutils"
)

func ExampleReverseString() {
	fmt.Println(stringutils.ReverseString("hello"))
	// Output: olleh
}
