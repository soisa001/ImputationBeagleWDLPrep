#!/usr/bin/env bash
# Initialize this repo and push it to a GitHub remote YOU created (empty repo).
#
#   1) On github.com create a new EMPTY repo (no README/.gitignore), e.g.
#        aou-lrwgs-imputebeagle-holdout   (private recommended — controlled-tier project)
#   2) From the repo root run, e.g.:
#        bash scripts/init_and_push.sh git@github.com:<you>/aou-lrwgs-imputebeagle-holdout.git
#      (use the git@ SSH form for a private repo; https:// for public)
#
# Then add it to your Workbench workspace and clone into an app:
#   wb resource add-ref git-repo --name=imputebeagle_holdout --repo-url=<same URL>
#   wb git clone --resource=imputebeagle_holdout
# (private repo requires your Workbench SSH key on GitHub: `wb security ssh-key get`)
set -euo pipefail

REMOTE="${1:-}"
[ -n "$REMOTE" ] || { echo "usage: bash scripts/init_and_push.sh <git-remote-url>"; exit 1; }
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

git rev-parse --git-dir >/dev/null 2>&1 || git init
git add -A
git commit -m "AoU lrWGS Phase-2 ACAF->panel holdout-198 imputation + eval pipeline" || echo "(nothing to commit)"
git branch -M main
if git remote get-url origin >/dev/null 2>&1; then git remote set-url origin "$REMOTE"; else git remote add origin "$REMOTE"; fi
git push -u origin main
echo
echo "pushed to $REMOTE"
echo "next: wb resource add-ref git-repo --name=imputebeagle_holdout --repo-url=$REMOTE"
