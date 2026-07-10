#!/usr/bin/env bash
set -u

DIR="/Volumes/Extreme SSD/关键work/Data/XeniumofGC"
PIDFILE="$DIR/aria2.pid"
LOG="$DIR/aria2-download.log"
SVS="$DIR/BS06-9313-8_Tumor__2.svs"
ZIP="$DIR/BS06-9313-8_Tumor.zip"

SVS_BYTES=1239136585
ZIP_BYTES=29795534461
SVS_MD5="a0c1eff7a6310fe7f7e046daafa9fe4e"
ZIP_MD5="941740d6436fcef3f23ca0fbd6c2547c"

size_of() {
  local path="$1"
  if [[ -f "$path" ]]; then
    stat -f "%z" "$path"
  else
    echo 0
  fi
}

human_size() {
  local bytes="$1"
  awk -v b="$bytes" 'BEGIN {
    split("B KiB MiB GiB TiB", u, " ");
    i=1;
    while (b >= 1024 && i < 5) { b/=1024; i++ }
    printf "%.2f %s", b, u[i]
  }'
}

report_file() {
  local label="$1" path="$2" expected="$3"
  local actual pct aria2_marker
  actual="$(size_of "$path")"
  pct="$(awk -v a="$actual" -v e="$expected" 'BEGIN { if (e>0) printf "%.2f", (a/e)*100; else print "0.00" }')"
  aria2_marker=""
  [[ -f "$path.aria2" ]] && aria2_marker="  (.aria2 present: still incomplete or resumable)"
  printf "%s: %s / %s (%s%%)%s\n" "$label" "$(human_size "$actual")" "$(human_size "$expected")" "$pct" "$aria2_marker"
}

echo "Target: $DIR"
echo

if [[ -f "$PIDFILE" ]]; then
  pid="$(cat "$PIDFILE" 2>/dev/null || true)"
  if [[ -n "${pid:-}" ]] && ps -p "$pid" >/dev/null 2>&1; then
    echo "aria2c: running (PID $pid)"
    ps -o pid,ppid,stat,etime,command -p "$pid"
  else
    echo "aria2c: not running (PID file exists: ${pid:-empty})"
  fi
else
  echo "aria2c: no PID file found"
fi

echo
report_file "SVS" "$SVS" "$SVS_BYTES"
report_file "ZIP" "$ZIP" "$ZIP_BYTES"

echo
echo "Recent aria2 log:"
if [[ -f "$LOG" ]]; then
  tail -n 30 "$LOG"
else
  echo "No log found at $LOG"
fi

echo
svs_size="$(size_of "$SVS")"
zip_size="$(size_of "$ZIP")"
if [[ "$svs_size" == "$SVS_BYTES" && "$zip_size" == "$ZIP_BYTES" && ! -f "$SVS.aria2" && ! -f "$ZIP.aria2" ]]; then
  echo "Both files appear complete by size. Running md5 verification..."
  svs_actual="$(md5 -q "$SVS")"
  zip_actual="$(md5 -q "$ZIP")"

  if [[ "$svs_actual" == "$SVS_MD5" ]]; then
    echo "SVS md5: OK"
  else
    echo "SVS md5: FAIL expected=$SVS_MD5 actual=$svs_actual"
  fi

  if [[ "$zip_actual" == "$ZIP_MD5" ]]; then
    echo "ZIP md5: OK"
  else
    echo "ZIP md5: FAIL expected=$ZIP_MD5 actual=$zip_actual"
  fi
else
  echo "Not complete yet, or .aria2 resume file is still present. Re-run later."
fi
