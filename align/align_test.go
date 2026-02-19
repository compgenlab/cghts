package align

import "testing"

func TestCigarCondense(t *testing.T) {
	tests := []struct {
		name  string
		input string
		want  string
	}{
		{name: "empty", input: "", want: ""},
		{name: "all_match", input: "MMMM", want: "4M"},
		{name: "mixed_with_del", input: "MMDMM", want: "2M1D2M"},
		{name: "mixed_with_ins_del", input: "IIMMMMMDMM", want: "2I5M1D2M"},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := CigarCondense(tt.input)
			if got != tt.want {
				t.Fatalf("CigarCondense(%q) = %q, want %q", tt.input, got, tt.want)
			}
		})
	}
}
