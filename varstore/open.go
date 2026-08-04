package varstore

import (
	"context"
	"errors"
	"fmt"
	"strings"

	"github.com/compgenlab/cghts/iosource"
)

// Store kinds, as accepted by OpenStore and returned by StoreKind.
const (
	KindVcf     = "vcf"
	KindParquet = "parquet"
)

// vcfSuffixes are the spellings that name a VCF rather than a store.
var vcfSuffixes = []string{".vcf", ".vcf.gz", ".vcf.bgz", ".bcf"}

// ErrUnknownStoreKind reports that a locator matched no backend.
//
// It is a distinct error so a caller can tell "I could not recognize this" from
// "the transport failed" and give the right advice. Suggesting a --store flag
// helps with the first and misleads on the second, where the flag would not
// have changed anything and the real problem is an unlinked transport or an
// object that is not there.
var ErrUnknownStoreKind = errors.New("unrecognized store")

// StoreKind reports which backend a locator names.
//
// The rule is deliberately shallow: a VCF is recognized by its suffix, and
// anything else is a store if a manifest is found at it. That works unchanged
// for a remote locator, which is the point -- the previous version of this logic
// lived in the CLI and reached for os.Stat twice, so a URL could never resolve
// to anything but a suffix match.
//
// An unusable transport is reported as such rather than as an unrecognized
// format. "gs://bucket/cohort" is not a mystery file type; it is a scheme
// nobody linked a transport for, and saying so is the difference between a
// user importing a package and a user rechecking their path.
func StoreKind(ctx context.Context, locator string) (string, error) {
	scheme := iosource.Scheme(locator)
	if !iosource.HasTransport(scheme) {
		return "", fmt.Errorf("no transport registered for %q (registered: %s)",
			scheme+"://", strings.Join(append([]string{"http", "https"}, iosource.Schemes()...), ", "))
	}

	lower := strings.ToLower(locator)
	for _, sfx := range vcfSuffixes {
		if strings.HasSuffix(lower, sfx) {
			return KindVcf, nil
		}
	}
	if strings.HasSuffix(lower, ".parquet") || strings.HasSuffix(lower, ManifestFile) {
		return KindParquet, nil
	}
	if _, err := ReadManifestContext(ctx, locator); err == nil {
		return KindParquet, nil
	}
	return "", fmt.Errorf("%w: cannot tell what kind of store %q is, it has no VCF "+
		"suffix and no %s was found in it", ErrUnknownStoreKind, locator, ManifestFile)
}

// OpenStore opens whichever backend the locator names.
//
// kind forces one ("vcf" or "parquet"); an empty kind infers via StoreKind.
func OpenStore(ctx context.Context, locator, kind string) (Store, error) {
	switch strings.ToLower(kind) {
	case KindVcf:
		return OpenVcfContext(ctx, locator)
	case KindParquet:
		return OpenParquetContext(ctx, locator)
	case "":
	default:
		return nil, fmt.Errorf("unknown store kind %q (use %s or %s)", kind, KindVcf, KindParquet)
	}

	inferred, err := StoreKind(ctx, locator)
	if err != nil {
		return nil, err
	}
	return OpenStore(ctx, locator, inferred)
}
