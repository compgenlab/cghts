# hts is a library (no binary output); this Makefile drives testing and checks.

GO ?= go
PKGS ?= ./...

# Match the project convention for a writable build cache.
GOCACHE ?= /tmp/go-build-cache
export GOCACHE

# hts has no local module dependencies, so it always builds standalone;
# ignore any ambient go.work from parent directories (e.g. a dev workspace).
export GOWORK := off

.PHONY: test test-race cover vet fmt fmt-check tidy build check doc clean

# Run the full test suite.
test:
	$(GO) test $(PKGS)

# Run tests with the race detector.
test-race:
	$(GO) test -race $(PKGS)

# Generate a coverage profile and print the total.
cover:
	$(GO) test -coverprofile=coverage.out $(PKGS)
	$(GO) tool cover -func=coverage.out | tail -1

# Static analysis.
vet:
	$(GO) vet $(PKGS)

# Format all sources in place.
fmt:
	gofmt -w .

# Fail if any file is not gofmt-clean (for CI / pre-commit).
fmt-check:
	@out=$$(gofmt -l .); if [ -n "$$out" ]; then echo "gofmt needed:"; echo "$$out"; exit 1; fi

# Prune and verify module requirements.
tidy:
	$(GO) mod tidy

# Compile every package (catches build breaks without producing a binary).
build:
	$(GO) build $(PKGS)

# One-shot gate: compile, vet, format check, and test.
check: build vet fmt-check test

# Preview the pkg.go.dev documentation locally (requires network for the tool).
doc:
	$(GO) run golang.org/x/pkgsite/cmd/pkgsite@latest -open .

clean:
	rm -f coverage.out
	$(GO) clean -testcache
