# Reproduction of the scoliosis candidate-nomination study

## Run offline

Python 3.12 and R are required (tested with R 4.4.2). R uses only standard packages.

```sh
python -m pip install -r requirements.txt
python reproduce.py --output /path/to/a/fresh_output
```

The output directory must be empty and outside this repository. No network request occurs by default. Readable outputs include:

- `inputs/variants_121.csv`: complete 121-record source inventory and disposition.
- `stability_summary_102.csv`: all 102 scores aggregated from 106 model records, with primary eligibility and both EvoEF references.
- `inputs/panel_35.csv` and `inputs/assemblies_70.csv`: complete panel/entry inventories, selection and provenance.
- `mechanics/*.csv`: recomputed group effects, geometry/redox calibration, constant-WT baseline, Cobb/LOOCV and 20,000 identity-blocked case-bootstrap results, residual/distribution/influence checks and HC3 intervals.
- Threshold/chance-overlap, sequence-correspondence, coverage and DynaMut2 CSVs; a diagnostic PDF; `constructed_structures/*.pdb` (eight files); an R log and `summary.json`.

Generated outputs do not need to be committed to GitHub.

## Optional public-source downloads

```sh
python reproduce.py --output /path/to/another_fresh_output --fetch-public assemblies
python reproduce.py --output /path/to/another_fresh_output --fetch-public all
```

`assemblies` retrieves 70 assembly-1 mmCIF files. `all` retrieves those plus four metadata/mapping files per entry and the original animal workbook (351 objects). URLs and SHA256 values are embedded in the data. A changed/unavailable object stops retrieval instead of being accepted as the old snapshot. These current downloads are not authenticated historical atom-level scoring inputs. Exact historical versions of the two public AlphaFold WT inputs were not recorded; accession lookup URLs and actual input hashes are retained, without asserting that today's API returns identical coordinates or automatically replacing those models.

## Files and analytical scope

`reproduce.py` is the entry point; implementation is separated by analytical responsibility:

| Code file | Responsibility |
|---|---|
| `reproduce.py` | CLI, output-directory checks and workflow orchestration |
| `io_utils.py` | Compressed inputs, checksum/schema checks and CSV output |
| `variant_analysis.py` | Frozen stability aggregation, eligibility and exact overlap tests |
| `structure_analysis.py` | Assembly/coverage summaries, PDB reconstruction and optional downloads |
| `sequence_annotations.py` | Human/fish sequence alignment and saved DynaMut2 response checks |
| `mechanics.R` | Mechanical calibration, bootstrap/LOOCV, diagnostics and HC3 calculations |

The separate **base R** file retains the original RNG and seeds, including seed 261028 for the 20,000 identity-blocked case bootstraps. No R packages require installation. Biopython 1.88 is the sole additional Python dependency (`requirements.txt`). The two compressed data files are byte-for-byte unchanged from the verified compact package.

`data.json.gz` is compressed JSON with named schemas: the necessary processed animal observations, complete variant/model/panel/entry records, frozen interfaces and coverage, sequences, dated candidate annotations, four unmodified completed DynaMut2 responses and public-source/software provenance. Every table is exported to readable CSV. Both original EvoEF reference conventions, NMR model means and exact chance-overlap tests are recomputed. All 32 sequence comparisons are rerun; live ClinVar/gnomAD/AlphaMissense queries are not.

`coordinates.json.gz` contains shared atom bases and local differences. The script reconstructs all eight constructed WT/mutant PDB files and requires each original SHA256 to match. The actual completed DynaMut2 responses are identity-checked and sign-harmonized; jobs are not resubmitted. Its negative-destabilization convention is preserved in the raw response.

This package **does not rerun historical ThermoMPNN/EvoEF neural/energy inference, atom-level solvent-accessible surfaces, live database annotation, original Mol* pixels or all manuscript figure layouts**. It reaggregates frozen prediction/interface inputs. Third-party repository revisions and available weight/executable hashes are retained. Authors keep the full provenance archive separately. The portable diagnostic PDF represents the same model checks, without claiming identical manuscript typesetting.

## Sources and limits

Animal measurements: Pumputis et al. (2025), [10.1038/s41467-025-63742-2](https://doi.org/10.1038/s41467-025-63742-2), workbook `41467_2025_63742_MOESM13_ESM.xlsx`. Variant sources: Haller et al. ([10.1093/hmg/ddv463](https://doi.org/10.1093/hmg/ddv463)), Jiang et al. ([10.1136/jmedgenet-2019-106411](https://doi.org/10.1136/jmedgenet-2019-106411)), Charng et al. ([10.1136/jmg-2024-110586](https://doi.org/10.1136/jmg-2024-110586)) and Yu et al. ([10.7554/eLife.89762](https://doi.org/10.7554/eLife.89762)). Record-specific source/table labels and PDB/SIFTS provenance remain attached. Original resources retain their own attribution/terms; no blanket relicensing is implied.

Both katnb1 control groups are heterozygotes; `WT` is an internal processed-table reference code. Row IDs are not a cross-sheet animal registry. Segment-to-animal assignments are missing. Internal resampling intervals do not establish external accuracy. Geometry loses to the constant-WT baseline in both withheld groups; original two-tool overlap equals its simple chance expectation. No human-variant-to-fish mechanism, clinical pathogenicity or public repository publication is established by running this local package.
