#!/bin/sh
# Gate to run before cutting a release tag. Tags nothing itself.
#
# A tag is effectively irreversible: the module proxy and any consumer can resolve
# it the moment it is pushed. So the state being tagged has to be the state that
# was actually tested -- and the gap between those two is subtler than a broken
# build.
#
# The failure this exists to prevent: work written and verified, then left
# uncommitted, and then either excluded by `git commit --amend` (which takes only
# what is staged) or destroyed by `git checkout <file>`. The test suite still
# passes, because it ran against the working tree. The tag ships without the fix.
# That has happened here, and a released tag contained a known bug as a result.
#
# Every check reports rather than exiting early, so one run tells you everything
# that is wrong.
#
# Usage: release-check.sh [release-branch] [tag]

set -u

branch_expected="${1:-main}"
version="${2:-}"
fail=0

say() { printf '%s\n' "$*"; }

# 1. A dirty tree means the tested state is not the committed state. This is the
#    check that would have caught every incident so far.
if [ -n "$(git status --porcelain)" ]; then
	say "REFUSING: working tree is not clean."
	say "  Uncommitted changes mean what you tested is not what you would tag."
	git status --short | sed 's/^/    /'
	fail=1
fi

# 2. Releases come from one branch, so a tag cannot accidentally capture a topic
#    branch that was never reviewed.
branch="$(git rev-parse --abbrev-ref HEAD)"
if [ "$branch" != "$branch_expected" ]; then
	say "REFUSING: on branch '$branch', not '$branch_expected'."
	say "  Tags are cut from '$branch_expected'; pass RELEASE_BRANCH to override."
	fail=1
fi

# 3. Local and remote must agree, or the tag names a commit nobody else has (or
#    omits one they do).
git fetch -q origin "$branch" 2>/dev/null || true
if git rev-parse --verify --quiet "origin/$branch" >/dev/null 2>&1; then
	if [ "$(git rev-parse HEAD)" != "$(git rev-parse "origin/$branch")" ]; then
		say "REFUSING: HEAD differs from origin/$branch."
		say "  A tag here would not match what others can fetch."
		git --no-pager log --oneline "origin/$branch..HEAD" 2>/dev/null | sed 's/^/    unpushed: /'
		git --no-pager log --oneline "HEAD..origin/$branch" 2>/dev/null | sed 's/^/    missing:  /'
		fail=1
	fi
else
	say "NOTE: origin/$branch not found; skipping the sync check."
fi

# 4. Never re-cut an existing tag. Deleting and re-pushing one changes what
#    consumers have already resolved, and has caused version churn here before.
if [ -n "$version" ]; then
	if git rev-parse --verify --quiet "refs/tags/$version" >/dev/null 2>&1; then
		say "REFUSING: tag $version already exists locally."
		fail=1
	fi
	if git ls-remote --exit-code --tags origin "refs/tags/$version" >/dev/null 2>&1; then
		say "REFUSING: tag $version already exists on origin."
		say "  Re-cutting it would change what consumers already resolved."
		fail=1
	fi
fi

if [ "$fail" -ne 0 ]; then
	say ""
	say "release-check FAILED -- do not tag."
	exit 1
fi

if [ -n "$version" ]; then
	say "state: clean, on $branch_expected, in sync with origin; $version is unused."
else
	say "state: clean, on $branch_expected, in sync with origin."
fi
