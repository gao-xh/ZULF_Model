From: eager-franklin-t
To: all
Re: 2026-10-06T0751Z_eager-franklin-t.md
Subject: merged into main at 97cfd78: AGENTS.md note on the cloud stop hook (author rule unchanged)

main fast-forwarded from 8ab5c84 to 97cfd78 (one commit, docs only), on the user's request.

- AGENTS.md: the author rule stays (Xuehan Gao <gao.xh@berkeley.edu> only). In Claude Code cloud
  containers a stop hook reports these commits as "Unverified" and asks to reset the author with
  --amend / rebase. That is expected: do not change the author, and never amend, rebase or force-push
  pushed commits. This is the owner's decision of 2026-10-06.
- docs/DEVELOPMENT.md: the same note under "Limits of a cloud session"; the record now says
  claude/eager-franklin-tp94gf was recreated from main at 8ab5c84.

Nothing else changes for the other lines; merge origin/main at your next sync.
