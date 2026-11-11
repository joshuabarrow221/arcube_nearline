# Stage A LRS flash extraction — run notes

Date: 2026-11-11 (local run ~11:39–12:00)
Command form (one URL per call, sequential, via WSL):

```
wsl.exe --cd "<repo>" -- .venv-purity/bin/python -I beam_position_scan/code/beam_position_scan/scripts/extract_lrs_flashes.py --out beam_position_scan/flashes --block-size 32768 <URL>
```

URL prefix: `https://portal.nersc.gov/project/dune/data/2x2/reflows/v11/flow/beam/`
All files: association = `ref`, gap_ns = 80.0. Output: `<name>.flashes.npz` in this folder.

**Result: 19/19 ok, 0 FAIL, no `.failed.txt` files.**

| # | File | Status | n_triggers | n_beam_triggers | n_beam_hits | n_flashes | seconds |
|---|------|--------|-----------:|----------------:|------------:|----------:|--------:|
| 1 | july8_2024/.../packet-0050017-2024_07_08_23_54_32_CDT | ok | 7843 | 453 | 15788 | 1049 | 57.2 |
| 2 | july8_2024/.../packet-0050017-2024_07_09_00_04_33_CDT | ok | 7839 | 454 | 16196 | 1131 | 62.0 |
| 3 | july8_2024/.../packet-0050017-2024_07_09_00_14_34_CDT | ok | 8004 | 454 | 15254 | 1002 | 59.9 |
| 4 | july8_2024/.../packet-0050017-2024_07_09_00_24_35_CDT | ok | 8043 | 454 | 14772 | 1038 | 53.6 |
| 5 | july8_2024/.../packet-0050017-2024_07_09_00_34_36_CDT | ok | 8043 | 455 | 16573 | 1088 | 59.4 |
| 6 | july8_2024/.../packet-0050017-2024_07_09_00_44_37_CDT | ok | 8012 | 458 | 16706 | 1082 | 54.2 |
| 7 | july8_2024/.../packet-0050017-2024_07_09_00_54_38_CDT | ok | 7988 | 440 | 15443 | 1006 | 58.0 |
| 8 | july8_2024/.../packet-0050017-2024_07_09_01_04_39_CDT | ok | 8102 | 454 | 16569 | 1065 | 56.6 |
| 9 | july8_2024/.../packet-0050017-2024_07_09_01_14_40_CDT | ok ⚠ | 16053 | 910 | 30098 | 2008 | 114.5 |
| 10 | july8_2024/.../packet-0050017-2024_07_09_01_24_41_CDT | ok | 8086 | 457 | 15997 | 1084 | 54.7 |
| 11 | july8_2024/.../packet-0050017-2024_07_09_01_34_42_CDT | ok | 7978 | 457 | 16076 | 1032 | 66.1 |
| 12 | july8_2024/.../packet-0050017-2024_07_09_01_44_44_CDT | ok | 7954 | 452 | 16888 | 1123 | 53.7 |
| 13 | july8_2024/.../packet-0050017-2024_07_09_01_54_45_CDT | ok | 8010 | 453 | 15751 | 1030 | 69.8 |
| 14 | july10_2024/.../packet-0050018-2024_07_12_00_00_19_CDT | ok ⚠ | 12156 | 682 | 23406 | 1594 | 78.9 |
| 15 | july10_2024/.../packet-0050018-2024_07_12_00_10_20_CDT | ok | 8074 | 457 | 15789 | 1103 | 51.0 |
| 16 | july10_2024/.../packet-0050018-2024_07_12_00_20_21_CDT | ok | 8031 | 456 | 15892 | 1045 | 59.8 |
| 17 | july10_2024/.../packet-0050018-2024_07_12_00_30_22_CDT | ok | 8056 | 457 | 15595 | 1087 | 56.5 |
| 18 | july10_2024/.../packet-0050018-2024_07_12_00_40_23_CDT | ok | 8105 | 460 | 15932 | 1048 | 58.0 |
| 19 | july10_2024/.../packet-0050018-2024_07_12_00_50_24_CDT | ok | 7972 | 457 | 16673 | 1088 | 54.4 |

Totals: run 50017 (13 files) = 6351 beam triggers, 14738 flashes, 819.7 s;
run 50018 (6 files) = 2969 beam triggers, 6965 flashes, 358.6 s;
all = 9320 beam triggers, 21703 flashes, 1178.3 s (~19.6 min).

## Anomalies to check (not investigated — out of scope for this run)

- **#9 (01_14_40)**: ~2.0× the triggers / beam triggers / flashes of its neighbours, even though
  the filename timestamps are still 10 min apart. Possible duplicated trigger blocks in the
  reflow file, or a longer actual time span. Check for duplicate trigger timestamps before using
  it for per-spill normalisation.
- **#14 (07_12_00_00_19)**: ~1.5× the usual trigger count (12156 vs ~8000; 682 beam triggers vs
  ~455). Same caveat.
- Beam triggers per file are otherwise stable at 440–460, flashes/beam trigger ≈ 2.2–2.5.
