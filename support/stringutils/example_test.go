package stringutils_test

import (
	"fmt"

	"github.com/compgenlab/hts/support/stringutils"
)

func ExampleReverseString() {
	fmt.Println(stringutils.ReverseString("hello"))
	// Output: olleh
}
