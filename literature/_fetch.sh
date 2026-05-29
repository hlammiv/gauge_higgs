#!/usr/bin/env bash
# Download the curated, arXiv-verified literature set into literature/.
# Polite: 3s between requests, descriptive filenames, skip already-present.
set -u
cd "$(dirname "$0")"
UA="Mozilla/5.0 (gauge-higgs-hmc litreview; mailto:hanklammiv@gmail.com)"
python3 - <<'PY' > _fetch_jobs.txt
import json
for x in json.load(open("_download_list.json")):
    print(f"{x['arxiv']}\t{x['filename']}")
PY
ok=0; fail=0
while IFS=$'\t' read -r id fname; do
  [ -z "$id" ] && continue
  if [ -s "$fname" ]; then echo "SKIP  $fname (exists)"; ok=$((ok+1)); continue; fi
  url="https://arxiv.org/pdf/${id}"
  echo "GET   $url -> $fname"
  if curl -sL --max-time 120 -A "$UA" -o "$fname" "$url" \
     && [ "$(file -b --mime-type "$fname" 2>/dev/null)" = "application/pdf" ]; then
    echo "  OK  $(du -h "$fname" | cut -f1)  $fname"; ok=$((ok+1))
  else
    echo "  FAIL $id (got $(file -b --mime-type "$fname" 2>/dev/null), removing)"; rm -f "$fname"; fail=$((fail+1))
  fi
  sleep 3
done < _fetch_jobs.txt
echo "=== DONE: $ok ok, $fail failed ==="
ls -la *.pdf 2>/dev/null | wc -l
