# Development workflow (branches, parallel sessions, merging)

How several development lines work on this repository at the same time: Xuehan Gao on a local machine and one
or more Claude Code cloud sessions. Every line works on its own branch and merges into `main`. Read this before
the first commit of a session.

## Repository

- GitHub `gao-xh/ZULF_Model` (`gao-xh/zulf_model` is the same repository; names are case-insensitive).
- `main` is the default branch and always holds a working state: tests of the changed area pass and the docs are
  current. Clone it, branch from it, merge back into it.
- Commit author: Xuehan Gao <gao.xh@berkeley.edu> only, no co-author lines, session links or "generated with"
  footers (AGENTS.md). In a new clone or cloud container:

```bash
git config user.name "Xuehan Gao"
git config user.email "gao.xh@berkeley.edu"
```

## One branch per development line

| Line | Branch |
|---|---|
| Xuehan, local machine | `xuehan/<topic>`, e.g. `xuehan/fit-scripts` |
| Claude cloud session | the branch the session is assigned (`claude/...`); the session may push only there and, when the user asks, to `main` |

Two lines never commit to the same branch at the same time. Sharing a branch makes every push depend on pulling
the other side first, and a forgotten pull is rejected or produces needless merges.

## Start of a session

```bash
git fetch origin
git checkout -B <your-branch> origin/main       # fresh work: start from the latest main
# or, to continue earlier work on the branch:  git checkout <your-branch> && git merge origin/main
pip install -e ".[dev]"
python scripts/check_ascii.py
python -m unittest tests.test_processing tests.test_j_tuner     # quick check, about 1.5 min
```

Then claim the work in `docs/PLAN.md`: mark the item `[~]` with "in progress: <branch>", commit and push this
line first, so the other lines see it after their next `git fetch` and pick a different item.

## During the work

- Commit small, self-contained steps; push the branch after each one. A cloud container is reclaimed when the
  session idles, and everything not pushed is lost (results under `runs/` are lost in any case; record the
  numbers and commands in the analysis log, as AGENTS.md requires).
- Before each push run the checks of the changed area: `python scripts/check_ascii.py` and the matching test
  modules (`tests/test_<area>.py`). The full suite takes several hours; run it before larger merges.
- Bring the other lines' work in regularly: `git fetch origin && git merge origin/main`.
- Split the work by area to keep conflicts rare, e.g. one line in `zulf_core/physics`, the other in `scripts/`.

## Merging into main

1. `git fetch origin && git merge origin/main` on the branch; resolve conflicts; run the checks again.
2. Bring `main` to the branch:
   - fast-forward (no review): `git push origin HEAD:main`. It succeeds only when `main` is an ancestor of the
     branch, i.e. after step 1; if it is rejected, `main` moved in between: repeat step 1.
   - with review: open a pull request from the branch into `main` and merge it on GitHub.
3. Update `docs/PLAN.md` (item done, measured result, commit) in the same merge.

Never force-push `main` and never rewrite history on a branch someone else uses (no rebase, amend or
`push --force` there).

## Conflicts

- Code: keep both changes where they are independent; when both sides changed the same logic, decide which
  behaviour is meant, and run the tests of both changes.
- Append-only records (`docs/PLAN.md`, `docs/ANALYSIS_LOG.md`, `docs/DECISIONS.md`, `docs/analysis/README.md`):
  the conflict is almost always two new entries in the same place. Keep both, in date order. Decision numbers
  (D47, D48, ...) can collide: renumber the later one and fix its references.

## Limits of a cloud session (observed 2026-10-06)

- Pushing to its own branch works; pushing to `main` worked when the user asked for it.
- Repository settings (default branch) and branch deletion are refused by the session's GitHub proxy
  (HTTP 403, or the connection is dropped). The owner does these on github.com: Settings -> General -> Default
  branch; Branches page -> delete.
- The session does not open pull requests unless asked.

## Record: branch setup of 2026-10-06

- `main` created at e9a1580 (identical to the former working branch `claude/eager-franklin-tp94gf`).
- Default branch switched to `main` and `claude/eager-franklin-tp94gf` deleted, both by the owner on github.com.
- Author rule added to AGENTS.md (9a09e53) and fast-forwarded into `main`.
