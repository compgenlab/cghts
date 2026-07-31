# hts is a library (no binary output); this Makefile drives testing and checks.

GO ?= go
PKGS ?= ./...

# Branch that releases are cut from; overridable so release-check can be
# exercised from a topic branch.
RELEASE_BRANCH ?= main

# Match the project convention for a writable build cache.
GOCACHE ?= /tmp/go-build-cache
export GOCACHE

# hts has no local module dependencies, so it always builds standalone;
# ignore any ambient go.work from parent directories (e.g. a dev workspace).
export GOWORK := off

.PHONY: test test-race cover vet fmt fmt-check tidy build check doc clean release-check

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

# Gate to run BEFORE cutting a tag. It does not tag anything.
#
# A tag is effectively irreversible -- the module proxy and any consumer can
# resolve it the moment it is pushed -- so the state being tagged has to be the
# state that was actually tested. The failure this exists to prevent is subtler
# than a broken build: work that was written and verified but left uncommitted,
# then excluded by `git commit --amend` (which takes only what is staged) or
# destroyed by `git checkout <file>`. Tests pass, the tag ships without the fix.
# That happened, and a released tag contained a known bug as a result.
#
# Usage:
#   make release-check              # gate main
#   make release-check RELEASE_BRANCH=topic   # gate a topic branch instead
#   make release-check VERSION=v1.2.3   # also refuse if that tag already exists
release-check:
	@fail=0; 	if [ -n "$$(git status --porcelain)" ]; then 		echo "REFUSING: working tree is not clean."; 		echo "  Uncommitted changes mean the tested state is not the committed state."; 		git status --short | sed 's/^/    /'; 		fail=1; 	fi; 	branch=$$(git rev-parse --abbrev-ref HEAD); 	if [ "$$branch" != "main" ]; then 		echo "REFUSING: on branch '$$branch', not main. Tags are cut from main."; 		fail=1; 	fi; 	git fetch -q origin "$$branch" 2>/dev/null || true; 	if [ -n "$$(git rev-parse --verify --quiet origin/$$branch)" ]; then 		if [ "$$(git rev-parse HEAD)" != "$$(git rev-parse origin/$$branch)" ]; then 			echo "REFUSING: HEAD differs from origin/$$branch."; 			echo "  Unpushed or stale commits mean the tag would not match what others see."; 			git --no-pager log --oneline origin/$$branch..HEAD | sed 's/^/    ahead:  /'; 			git --no-pager log --oneline HEAD..origin/$$branch | sed 's/^/    behind: /'; 			fail=1; 		fi; 	fi; 	if [ -n "$(VERSION)" ]; then 		if git rev-parse --verify --quiet "refs/tags/$(VERSION)" >/dev/null; then 			echo "REFUSING: tag $(VERSION) already exists locally."; 			echo "  Re-cutting a published tag changes what consumers already resolved."; 			fail=1; 		fi; 		if git ls-remote --exit-code --tags origin "refs/tags/$(VERSION)" >/dev/null 2>&1; then 			echo "REFUSING: tag $(VERSION) already exists on origin."; 			fail=1; 		fi; 	fi; 	if [ $$fail -ne 0 ]; then echo; echo "release-check FAILED -- do not tag."; exit 1; fi; 	echo "state: clean, on $(RELEASE_BRANCH), in sync with origin$(if $(VERSION), and $(VERSION) is unused,)"
	@$(MAKE) --no-print-directory check
	@echo
	@echo "release-check PASSED. Safe to tag$(if $(VERSION), $(VERSION),)."

# Preview the pkg.go.dev documentation locally (requires network for the tool).
doc:
	$(GO) run golang.org/x/pkgsite/cmd/pkgsite@latest -open .

clean:
	rm -f coverage.out
	$(GO) clean -testcache
