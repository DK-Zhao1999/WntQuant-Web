# WntQuant Web

**WntQuant Web is a ready-to-use, browser-based tool for quantifying direction-aware Wnt/β-catenin pathway activity from any human gene expression matrix.**

No R installation, no package dependencies, and no command line are required. Users upload a gene expression matrix, run the calculation, and download a plain-text activity score file. All computation happens locally in the browser, so expression data never leaves the user's device.

## Live tool

**Access the tool here:**

<https://dk-zhao1999.github.io/WntQuant-Web/>

Source repository:

<https://github.com/DK-Zhao1999/WntQuant-Web>

## What WntQuant Web does

WntQuant Web maps a curated collection of 63 direction-aware Wnt gene sets onto an uploaded expression matrix and produces a signed, sample-level Wnt activity score.

The gene-set collection includes:

- 36 `WPAGS` activation gene sets (`pos` direction)
- 27 `WPIGS` inhibition gene sets (`neg` direction)

For every sample, the tool reports:

```text
activity_score = mean(activation score) - mean(inhibition score)
```

The default scoring method is `relative_ssGSEA`, matching the directional scoring approach used by the WntQuant-Paper pipeline. If a matrix contains only one sample, WntQuant Web automatically switches to an `absolute` within-sample z-score method because relative ssGSEA requires at least two samples for cross-sample standardization.

## Why use this tool

- **Immediate use:** open the website and run an analysis without installing software.
- **Direction-aware:** activation and inhibition are modeled separately before being combined into one signed score.
- **Privacy-preserving:** all calculations run locally; no expression matrix is uploaded.
- **Standard input:** accepts common `CSV`, `TSV`, and `TXT` expression matrices.
- **Reproducible output:** returns a tab-delimited text file that can be used directly in downstream analyses.
- **Self-contained example:** includes a one-click example dataset for demonstration.

## Quick start

1. Open <https://dk-zhao1999.github.io/WntQuant-Web/>.
2. Click **Load example data** to see a completed run immediately.
3. Or drag and drop your own expression matrix into the upload area.
4. Click **Calculate Wnt activity**.
5. Review the table and click **Download TXT**.

## Input format

The input matrix must be arranged as genes in rows and samples in columns.

```text
Gene        Sample_1    Sample_2    Sample_3
AXIN2       8.12        7.95        8.31
CTNNB1      9.03        8.77        9.20
ZNRF3       6.41        6.22        6.58
...
```

Supported details:

- File extensions: `.txt`, `.tsv`, `.csv`
- Delimiters: tab, comma, semicolon, or whitespace
- Header cell: `Gene`, `Symbol`, `ID`, `Ensembl`, `feature`, or empty
- Gene identifiers: human gene symbols are recommended
- Duplicate gene symbols are averaged automatically
- Rows with missing values are skipped and reported

## Output format

The downloaded `WntQuant_activity_scores.txt` file has four tab-separated columns:

| Column | Description |
| --- | --- |
| `Sample` | Sample name from the input matrix |
| `activity_score` | Final signed Wnt activity score |
| `pos_score` | Mean standardized activation score |
| `neg_score` | Mean standardized inhibition score |

## Method summary

The ssGSEA implementation follows the GSVA `ssgsea` kernel and uses:

- average tie ranks
- `alpha = 0.25`
- global min/max normalization
- cross-sample standardization of each gene-set score
- final score = mean activation score − mean inhibition score

The numerical output has been validated against the R implementation using `GSVA::ssgseaParam`.

## Run locally

Because the site is a static web application, it can be opened directly from a local file:

```bash
open index.html
```

Alternatively, serve the directory with a lightweight HTTP server:

```bash
python3 -m http.server 4173
```

Then visit <http://127.0.0.1:4173/>.

## Repository structure

```text
.
├── index.html                    # Main tool interface
├── styles.css                    # Responsive styling
├── app.js                        # Parsing, scoring, and download logic
├── genesets.js                   # Embedded 63 directional Wnt gene sets
├── example_expression_matrix.tsv # Ready-to-use example input
├── README.md
├── LICENSE
└── .github/workflows/deploy.yml  # GitHub Pages deployment
```

## Deployment

The repository is configured to deploy automatically to GitHub Pages with GitHub Actions. Any push to `main` triggers a fresh deployment.

To deploy your own copy:

1. Fork or clone this repository.
2. Push the code to a GitHub repository.
3. In **Settings → Pages**, set **Source** to **GitHub Actions**.
4. The included workflow will publish the site automatically.

## Related projects

- [WntQuant](https://github.com/FangZY-Lab/WntQuant)
- [WntQuant-Paper](https://github.com/FangZY-Lab/WntQuant-Paper)

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE).
