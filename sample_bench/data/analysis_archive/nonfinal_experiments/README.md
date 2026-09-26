# Analysis results outside the current final set

Moved from `sample_bench/data/analysis/` on 2026-09-26 at the user's request
to archive analyses outside the final dataset collection.

The retained collection is defined by the current temperature/volume plots,
the manuscript R-factor source manifest, and the native_5 known-weight control.
These archived results include scientifically distinct comparisons and alternate
refinement pipelines; archive placement does not imply invalid or duplicate data.

| Directory | Reason for archiving |
| --- | --- |
| `279_exp_phenix_ref` | Earlier two-state experimental tests |
| `279_exp_ref` | Earlier two-state experimental tests |
| `280_exp_all_2` | Earlier two-state condition-pair sweep |
| `280_exp_all_2_phenix_ref` | Earlier two-state condition-pair sweep |
| `280_exp_all_2_ref` | Earlier two-state condition-pair sweep |
| `281_exp_4_test` | Earlier four-state tests |
| `281_exp_4_test_phenix_ref` | Earlier four-state tests |
| `282_test_w_phenix_ref` | Earlier state-count and restraint-weight tests |
| `283_2_cond_phenix_ref` | Alternate refinement pipeline; not a current manuscript source |
| `283_2_cond_phenix_ref_2` | Alternate refinement pipeline; not a current manuscript source |
| `283_2_cond_ref_phenix_ref_copy` | Alternate refinement results; not a current manuscript source |
| `284_2_cond_4_state_phenix_ref_2` | Four-state comparison outside the current final tables |
| `285_2_state_3_cond_phenix_ref` | Alternate refinement pipeline; not a current manuscript source |
| `286_3_state_3_cond_phenix_ref` | Three-state / three-condition comparison outside the current final tables |
| `287_3_state_2_cond_phenix_ref` | Three-state / two-condition comparison outside the current final tables |

The current working copies, including existing local modifications, were moved
intact. SHA-256 checks verified every file's contents after the move. All retained
analysis files were also verified unchanged. Input parameter files and raw run
outputs were not moved.

To restore a dataset, move its directory back to `sample_bench/data/analysis/`
after checking that the destination does not exist. Current plot and manuscript
source paths remain valid.
