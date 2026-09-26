# Older native analysis results

Moved from `sample_bench/data/analysis/` on 2026-09-26 to keep the active
analysis directory focused on native_5 and experimental datasets.

| Directory | Synthetic dataset |
| --- | --- |
| `271_native_2_wxray` | `native_2` |
| `271_native_2_wxray_ref` | `native_2` |
| `272_correct_w_wxray` | `native_2` |
| `272_correct_w_wxray_ref` | `native_2` |
| `273_native_3_wxray_ref` | `native_3` |
| `274_native_3_corr_w_ref` | `native_3` |
| `276_native_4_ref` | `native_4` |

All files were moved intact, with SHA-256 checks confirming identical contents.
Original parameter CSVs, input structures, and reflection datasets remain in
place. To restore a dataset, move its directory back into
`sample_bench/data/analysis/` after checking that the destination does not exist.
The historical `copy_pdb_files.py` helper now points to this archive.
