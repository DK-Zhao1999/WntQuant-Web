# WntQuant Web

WntQuant Web is a client-side tool for quantifying direction-aware Wnt/β-catenin pathway activity from any human gene expression matrix. Upload a CSV, TSV, or TXT file and download the resulting activity scores as a TXT file. All calculations run locally in the browser; no expression data is uploaded.

## Live site

The live site is:

<https://www.WntQuant.com>

Source repository: <https://github.com/DK-Zhao1999/WntQuant-Web>

GitHub Pages is configured to deploy the `main` branch automatically with the included workflow.

## What it does

WntQuant Web maps 63 direction-aware Wnt gene sets to your expression matrix:

- 36 `WPAGS` activation gene sets (`pos` direction)
- 27 `WPIGS` inhibition gene sets (`neg` direction)

For each sample it returns:

```text
Sample    activity_score    pos_score    neg_score
```

The default method is `relative_ssGSEA`, matching the scoring approach used by the WntAct pipeline. For a single-sample matrix, the tool automatically switches to an `absolute` within-sample z-score method because `relative_ssGSEA` requires at least two samples for cross-sample standardization.

## Try it locally

Open `index.html` directly in a browser, or serve the directory with a simple HTTP server:

```bash
python3 -m http.server 4173
```

Then visit <http://127.0.0.1:4173/>.

For a quick demo, click **Load example data**, or download `example_expression_matrix.tsv`.

## Input format

- Rows are genes; columns are samples.
- The first cell in the header row may be `Gene`, `Symbol`, `ID`, or empty.
- Supported extensions: `.txt`, `.tsv`, `.csv`.
- Delimiters are detected automatically (tab, comma, semicolon, or whitespace).
- Gene identifiers should be human gene symbols. Duplicate symbols are averaged.

## Method notes

The ssGSEA implementation follows the GSVA `ssgsea` kernel with:

- average tie ranks,
- `alpha = 0.25`,
- global min/max normalization,
- cross-sample standardization of each gene-set score,
- final score = mean activation score − mean inhibition score.

The numerical output was validated against the R implementation using `GSVA::ssgseaParam`.

## Related projects

- [WntQuant](https://github.com/FangZY-Lab/WntQuant)
- [WntAct](https://github.com/FangZY-Lab/WntAct)

## Deploying to GitHub Pages with the custom domain

1. Push this repository to GitHub.
2. In **Settings → Pages**, set the source to **GitHub Actions** (the included workflow deploys the site).
3. In **Settings → Pages → Custom domain**, enter `www.WntQuant.com`.
4. At your DNS provider, add a `CNAME` record:

```text
www  ->  <your-username>.github.io
```

GitHub Pages will then serve the site at `https://www.WntQuant.com`.

## License

MIT. See [LICENSE](LICENSE).
