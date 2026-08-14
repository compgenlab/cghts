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
	// No manifest, but a calls member still identifies this as a store -- an
	// unfinished one. Saying so is what makes the failure diagnosable: opening
	// it then reports "no manifest, re-convert it or inspect it", which is
	// actionable, where calling it unrecognizable sends the reader off to check
	// their path. Without this the common case is circular, since the manifest
	// is both what identifies a store and the thing it is missing.
	if memberExists(ctx, CallsPath(locator)) {
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
	case KindSet:
		return OpenSet(ctx, locator)
	case "":
	default:
		return nil, fmt.Errorf("unknown store kind %q (use %s or %s)", kind, KindVcf, KindParquet)
	}

	// A SET IS CHECKED FIRST, and by its own marker rather than by a field
	// inside a shared one.
	//
	// The alternative -- one manifest.json.gz carrying a kind -- is tidier and
	// fails worse. Every reader that already exists, this package's older
	// versions included, treats the presence of manifest.json.gz as "this is a
	// store"; a set using that name would be opened as one, parse into a
	// manifest describing no members, and fail somewhere further in with an
	// error about a missing calls file. A distinct marker means an old reader
	// says "not a store", which is true and is the direction a failure should
	// point.
	if IsSet(ctx, locator) {
		return OpenSet(ctx, locator)
	}

	inferred, err := StoreKind(ctx, locator)
	if err != nil {
		return nil, err
	}
	return OpenStore(ctx, locator, inferred)
}

// KindSet names a varset: several stores, disjoint by chromosome, read as one.
const KindSet = "varset"
