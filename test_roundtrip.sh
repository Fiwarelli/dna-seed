#!/bin/bash
# Roundtrip verification for dna_seed.
# Compresses + decompresses a few representative inputs and checks byte-exact match.
set -e

DNA=./dna_seed
TMP=$(mktemp -d)
trap "rm -rf $TMP" EXIT

fail() { echo "FAIL: $*" >&2; exit 1; }
ok()   { echo "ok:   $*"; }

check() {
  local label="$1" src="$2"
  $DNA compress   "$src"         "$TMP/out.dna"     > /dev/null
  $DNA decompress "$TMP/out.dna" "$TMP/restored"    > /dev/null
  cmp -s "$src" "$TMP/restored" || fail "$label: roundtrip mismatch"
  ok "$label ($(stat -c%s "$src") bytes → $(stat -c%s "$TMP/out.dna") bytes)"
}

# 1. Tiny input
printf 'hello world\n' > "$TMP/tiny.txt"
check "tiny"  "$TMP/tiny.txt"

# 2. Empty-ish input (single byte)
printf 'a' > "$TMP/single.txt"
check "single-byte" "$TMP/single.txt"

# 3. Highly repetitive
yes "repeat " 2>/dev/null | head -c 65536 > "$TMP/repeat.txt"
check "repetitive-64k" "$TMP/repeat.txt"

# 4. Mixed text (this source file)
check "source-self" "dna_seed.c"

# 5. Random-ish binary (hard case)
head -c 32768 /dev/urandom > "$TMP/random.bin"
check "random-32k" "$TMP/random.bin"

# 6. Larger text (if present)
if [ -s /usr/share/dict/words ]; then
  check "words-dict" /usr/share/dict/words
fi

echo
echo "all roundtrips passed."
