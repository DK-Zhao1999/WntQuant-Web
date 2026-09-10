(function () {
  "use strict";

  const ALPHA = 0.25;
  const GENE_SETS = window.WNT_GENE_SETS || [];

  const $ = (selector) => document.querySelector(selector);

  const dropzone = $("#dropzone");
  const fileInput = $("#file-input");
  const exampleBtn = $("#example-btn");
  const downloadExampleBtn = $("#download-example-btn");
  const resetBtn = $("#reset-btn");
  const computeBtn = $("#compute-btn");
  const downloadBtn = $("#download-btn");
  const parsePanel = $("#parse-panel");
  const warningBox = $("#warning-box");
  const resultsSection = $("#results");
  const resultsBody = $("#results-body");
  const fileName = $("#file-name");
  const fileSize = $("#file-size");
  const summaryGenes = $("#summary-genes");
  const summarySamples = $("#summary-samples");
  const summarySets = $("#summary-sets");
  const parseNote = $("#parse-note");
  const methodBadge = $("#method-badge");
  const resultNote = $("#result-note");

  const state = {
    parsed: null,
    result: null,
    file: null,
  };

  function formatNumber(value, digits) {
    if (!Number.isFinite(value)) return "—";
    return Number(value).toFixed(digits);
  }

  function formatFileSize(bytes) {
    if (bytes < 1024) return bytes + " B";
    if (bytes < 1024 * 1024) return (bytes / 1024).toFixed(1) + " KB";
    return (bytes / (1024 * 1024)).toFixed(1) + " MB";
  }

  function detectDelimiter(line) {
    const tabs = (line.match(/\t/g) || []).length;
    const commas = (line.match(/,/g) || []).length;
    const semicolons = (line.match(/;/g) || []).length;
    const best = Math.max(tabs, commas, semicolons);
    if (best === 0) return "space";
    if (best === tabs) return "\t";
    if (best === commas) return ",";
    return ";";
  }

  function splitLine(line, delimiter) {
    if (delimiter === "space") return line.trim().split(/\s+/);
    return line.split(delimiter);
  }

  function cleanField(value) {
    return value.trim().replace(/^"|"$/g, "");
  }

  function isNumeric(value) {
    const text = value.trim();
    return text.length > 0 && Number.isFinite(Number(text));
  }

  function uniqueSampleNames(names) {
    const seen = new Set();
    return names.map((name, index) => {
      const base = name && name.trim() ? name.trim() : "Sample_" + (index + 1);
      let candidate = base;
      let suffix = 2;
      while (seen.has(candidate)) {
        candidate = base + "_" + suffix;
        suffix += 1;
      }
      seen.add(candidate);
      return candidate;
    });
  }

  function parseMatrix(text) {
    const lines = text
      .replace(/^\uFEFF/, "")
      .replace(/\r\n/g, "\n")
      .replace(/\r/g, "\n")
      .split("\n")
      .map((line) => line.replace(/\s+$/, ""))
      .filter((line) => line.trim().length > 0);

    if (lines.length === 0) {
      throw new Error("The file is empty. Please check the file and try again.");
    }

    const delimiter = detectDelimiter(lines[0]);
    const firstRow = splitLine(lines[0], delimiter).map(cleanField);

    if (firstRow.length < 2) {
      throw new Error("No sample columns were detected. Make sure the first row contains sample names.");
    }

    const afterFirst = firstRow.slice(1);
    const firstCellLooksLikeHeader = /^(gene|genes|symbol|id|ids|feature|probe|ensembl|sample|samples|rownames)$/i.test(
      firstRow[0],
    );
    const anyNonNumeric = afterFirst.some((cell) => !isNumeric(cell));
    const hasHeader = firstCellLooksLikeHeader || anyNonNumeric;

    const rawSampleNames = hasHeader
      ? afterFirst
      : afterFirst.map((_, index) => "Sample_" + (index + 1));
    const sampleNames = uniqueSampleNames(rawSampleNames);
    const dataStart = hasHeader ? 1 : 0;
    const sampleCount = sampleNames.length;

    const accumulator = new Map();
    let totalRows = 0;
    let duplicateRows = 0;

    for (let lineIndex = dataStart; lineIndex < lines.length; lineIndex += 1) {
      const fields = splitLine(lines[lineIndex], delimiter).map(cleanField);
      if (fields.length < 2) continue;

      const symbol = fields[0];
      if (!symbol) continue;

      totalRows += 1;
      const key = symbol.toUpperCase();

      if (!accumulator.has(key)) {
        accumulator.set(key, {
          symbol,
          sums: new Array(sampleCount).fill(0),
          counts: new Array(sampleCount).fill(0),
        });
      } else {
        duplicateRows += 1;
      }

      const record = accumulator.get(key);
      for (let i = 0; i < sampleCount; i += 1) {
        const value = Number(fields[i + 1]);
        if (Number.isFinite(value)) {
          record.sums[i] += value;
          record.counts[i] += 1;
        }
      }
    }

    if (accumulator.size === 0) {
      throw new Error("No gene expression values could be parsed from the file.");
    }

    const genes = [];
    const matrix = [];
    let droppedMissing = 0;

    accumulator.forEach((record) => {
      const values = record.sums.map((sum, i) =>
        record.counts[i] > 0 ? sum / record.counts[i] : NaN,
      );
      if (values.some((value) => !Number.isFinite(value))) {
        droppedMissing += 1;
        return;
      }
      genes.push(record.symbol);
      matrix.push(values);
    });

    if (genes.length === 0) {
      throw new Error("The matrix contains no complete numeric rows. Please check the data format.");
    }

    return {
      genes,
      samples: sampleNames,
      matrix,
      totalRows,
      duplicateRows,
      droppedMissing,
    };
  }

  function matchGeneSets(parsed) {
    const geneIndex = new Map();
    parsed.genes.forEach((gene, index) => {
      geneIndex.set(gene.toUpperCase(), index);
    });

    const sets = GENE_SETS.map((set) => ({
      name: set.name,
      direction: set.direction,
      indices: set.genes
        .map((gene) => geneIndex.get(gene.toUpperCase()))
        .filter((index) => index !== undefined),
    }));

    const overlappingSets = sets.filter((set) => set.indices.length > 0);
    const matchedGenes = new Set(overlappingSets.flatMap((set) => set.indices));

    return {
      sets,
      overlappingSets: overlappingSets.length,
      matchedGenes: matchedGenes.size,
    };
  }

  function rankAverage(values) {
    const n = values.length;
    const order = Array.from({ length: n }, (_, i) => i);
    order.sort((a, b) => values[a] - values[b] || a - b);

    const ranks = new Array(n);
    let i = 0;
    while (i < n) {
      let j = i + 1;
      while (j < n && values[order[j]] === values[order[i]]) j += 1;
      const averageRank = (i + j - 1) / 2 + 1;
      for (let k = i; k < j; k += 1) ranks[order[k]] = averageRank;
      i = j;
    }
    return ranks;
  }

  function computeSSGSEA(matrix, genes, sets) {
    const geneCount = genes.length;
    const sampleCount = matrix[0].length;
    const rawScores = sets.map(() => new Array(sampleCount).fill(NaN));

    for (let sample = 0; sample < sampleCount; sample += 1) {
      const values = genes.map((_, gene) => matrix[gene][sample]);
      const ranks = rankAverage(values);
      const descendingOrder = genes
        .map((_, i) => i)
        .sort((a, b) => values[b] - values[a] || a - b);
      const positionOfGene = new Array(geneCount);
      descendingOrder.forEach((gene, position) => {
        positionOfGene[gene] = position;
      });

      sets.forEach((set, setIndex) => {
        let sumRaWeight = 0;
        let sumRa = 0;
        let sumWeight = 0;
        let overlap = 0;

        set.indices.forEach((gene) => {
          const rankValue = ranks[gene];
          const position = positionOfGene[gene];
          const weight = geneCount - position;
          const rankPower = Math.pow(rankValue, ALPHA);
          sumRaWeight += rankPower * weight;
          sumRa += rankPower;
          sumWeight += weight;
          overlap += 1;
        });

        if (overlap === 0 || overlap === geneCount) {
          rawScores[setIndex][sample] = NaN;
          return;
        }

        const inSetCdf = sumRaWeight / sumRa;
        const outSetCdf =
          (geneCount * (geneCount + 1)) / 2 - sumWeight;
        rawScores[setIndex][sample] =
          inSetCdf - outSetCdf / (geneCount - overlap);
      });
    }

    const finiteScores = [];
    rawScores.forEach((row) => {
      row.forEach((value) => {
        if (Number.isFinite(value)) finiteScores.push(value);
      });
    });

    if (finiteScores.length === 0) {
      throw new Error("No Wnt gene sets matched the input genes. Please make sure the input uses gene symbols.");
    }

    let min = Infinity;
    let max = -Infinity;
    finiteScores.forEach((value) => {
      if (value < min) min = value;
      if (value > max) max = value;
    });

    const range = max - min;
    return rawScores.map((row) =>
      row.map((value) => {
        if (!Number.isFinite(value)) return NaN;
        return range > 1e-12 ? value / range : 0;
      }),
    );
  }

  function computeAbsolute(matrix, genes, sets) {
    const geneCount = genes.length;
    const sampleCount = matrix[0].length;
    const scores = sets.map(() => new Array(sampleCount).fill(NaN));

    for (let sample = 0; sample < sampleCount; sample += 1) {
      let mean = 0;
      for (let gene = 0; gene < geneCount; gene += 1) {
        mean += matrix[gene][sample];
      }
      mean /= geneCount;

      let sumSquares = 0;
      for (let gene = 0; gene < geneCount; gene += 1) {
        const diff = matrix[gene][sample] - mean;
        sumSquares += diff * diff;
      }
      const sd = Math.sqrt(sumSquares / (geneCount - 1));

      sets.forEach((set, setIndex) => {
        if (set.indices.length === 0 || !Number.isFinite(sd) || sd < 1e-12) {
          scores[setIndex][sample] = NaN;
          return;
        }
        let sum = 0;
        set.indices.forEach((gene) => {
          sum += (matrix[gene][sample] - mean) / sd;
        });
        scores[setIndex][sample] = sum / set.indices.length;
      });
    }

    return scores;
  }

  function buildFinalScores(scoresBySet, sets, samples, standardize) {
    const sampleCount = samples.length;
    const keptSets = sets
      .map((set, index) => ({ ...set, sourceIndex: index }))
      .filter((set) =>
        scoresBySet[set.sourceIndex].some((value) => Number.isFinite(value)),
      );

    if (keptSets.length === 0) {
      throw new Error("No computable Wnt gene sets were found.");
    }

    const z = Array.from({ length: sampleCount }, () =>
      new Array(keptSets.length).fill(0),
    );

    if (standardize) {
      keptSets.forEach((set, column) => {
        const values = scoresBySet[set.sourceIndex];
        let mean = 0;
        values.forEach((value) => {
          mean += value;
        });
        mean /= sampleCount;

        let sumSquares = 0;
        values.forEach((value) => {
          const diff = value - mean;
          sumSquares += diff * diff;
        });
        const sd = Math.sqrt(sumSquares / (sampleCount - 1));

        for (let sample = 0; sample < sampleCount; sample += 1) {
          z[sample][column] =
            sd > 1e-12 ? (values[sample] - mean) / sd : 0;
        }
      });
    } else {
      keptSets.forEach((set, column) => {
        for (let sample = 0; sample < sampleCount; sample += 1) {
          z[sample][column] = scoresBySet[set.sourceIndex][sample];
        }
      });
    }

    const positiveColumns = [];
    const negativeColumns = [];
    keptSets.forEach((set, column) => {
      if (set.direction === "pos") positiveColumns.push(column);
      else negativeColumns.push(column);
    });

    if (positiveColumns.length === 0 || negativeColumns.length === 0) {
      throw new Error("Gene-set direction information is incomplete, so activity scores cannot be calculated.");
    }

    const meanOfColumns = (columns, sample) => {
      let sum = 0;
      columns.forEach((column) => {
        sum += z[sample][column];
      });
      return sum / columns.length;
    };

    return samples.map((sample, index) => {
      const positive = meanOfColumns(positiveColumns, index);
      const negative = meanOfColumns(negativeColumns, index);
      return {
        sample,
        activity: positive - negative,
        positive,
        negative,
      };
    });
  }

  function computeWntScores(parsed) {
    const matched = matchGeneSets(parsed);
    if (matched.overlappingSets === 0) {
      throw new Error("No genes overlap with the Wnt gene sets. Please use gene symbols.");
    }

    const method = parsed.samples.length === 1 ? "absolute" : "relative_ssGSEA";
    const scoresBySet =
      method === "absolute"
        ? computeAbsolute(parsed.matrix, parsed.genes, matched.sets)
        : computeSSGSEA(parsed.matrix, parsed.genes, matched.sets);

    const scores = buildFinalScores(
      scoresBySet,
      matched.sets,
      parsed.samples,
      method !== "absolute",
    );

    return {
      method,
      scores,
      overlappingSets: matched.overlappingSets,
      matchedGenes: matched.matchedGenes,
    };
  }

  function buildTxt(result) {
    const lines = ["Sample\tactivity_score\tpos_score\tneg_score"];
    result.scores.forEach((score) => {
      lines.push(
        [
          score.sample,
          Number.isFinite(score.activity)
            ? score.activity.toFixed(8)
            : "NaN",
          Number.isFinite(score.positive) ? score.positive.toFixed(8) : "NaN",
          Number.isFinite(score.negative) ? score.negative.toFixed(8) : "NaN",
        ].join("\t"),
      );
    });
    return lines.join("\n") + "\n";
  }

  function downloadText(filename, text) {
    const blob = new Blob([text], { type: "text/plain;charset=utf-8" });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.download = filename;
    document.body.appendChild(link);
    link.click();
    link.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  }

  async function loadTextAsFile(filename, text) {
    const file = new File([text], filename, { type: "text/plain;charset=utf-8" });
    await handleFile(file);
    if (state.parsed) {
      computeBtn.click();
    }
  }

  function makeExampleMatrix() {
    const setGenes = Array.from(
      new Set(GENE_SETS.flatMap((set) => set.genes)),
    );
    const genes = setGenes.concat(
      Array.from({ length: 120 }, (_, i) => "RANDOM" + (i + 1)),
    );
    const samples = [
      "Sample_A",
      "Sample_B",
      "Sample_C",
      "Sample_D",
      "Sample_E",
      "Sample_F",
    ];
    let seed = 987654321;
    const random = () => {
      seed = (seed * 1664525 + 1013904223) >>> 0;
      return seed / 4294967296;
    };
    const gaussian = () => {
      const u = Math.max(random(), 1e-9);
      const v = Math.max(random(), 1e-9);
      return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
    };

    const rows = [["Gene", ...samples].join("\t")];
    genes.forEach((gene, geneIndex) => {
      const isWnt = geneIndex < setGenes.length;
      const values = samples.map((_, sampleIndex) => {
        const signal = isWnt ? (sampleIndex - 2.5) * 0.42 : 0;
        return (signal + gaussian() * 0.9).toFixed(5);
      });
      rows.push([gene, ...values].join("\t"));
    });
    return rows.join("\n");
  }

  function showWarning(message) {
    warningBox.textContent = message;
    warningBox.classList.remove("hidden");
  }

  function clearWarning() {
    warningBox.textContent = "";
    warningBox.classList.add("hidden");
  }

  function renderParsed(file, parsed, matched) {
    state.parsed = parsed;
    state.file = file;
    state.result = null;

    fileName.textContent = file.name;
    fileSize.textContent = formatFileSize(file.size);
    summaryGenes.textContent = String(parsed.genes.length);
    summarySamples.textContent = String(parsed.samples.length);
    summarySets.textContent = matched.overlappingSets + " / " + GENE_SETS.length;

    const notes = [];
    if (parsed.duplicateRows > 0) {
      notes.push(parsed.duplicateRows + " duplicate gene symbols were averaged");
    }
    if (parsed.droppedMissing > 0) {
      notes.push(parsed.droppedMissing + " genes with missing values were skipped");
    }
    if (matched.matchedGenes > 0) {
      notes.push(matched.matchedGenes + " Wnt genes matched");
    }
    parseNote.textContent = notes.join(" · ") || "File parsed successfully.";

    parsePanel.classList.remove("hidden");
    resultsSection.classList.add("hidden");
    computeBtn.disabled = false;
    downloadBtn.disabled = true;
    clearWarning();
  }

  function renderResults(result) {
    state.result = result;

    const methodLabels = {
      relative_ssGSEA: "Relative ssGSEA",
      absolute: "Absolute z-score",
    };
    methodBadge.textContent = methodLabels[result.method] || result.method;

    const noteParts = [
      result.overlappingSets + " gene sets matched",
      "Direction score = activation score − inhibition score",
    ];
    if (result.method === "relative_ssGSEA") {
      noteParts.push("Gene-set scores standardized across samples");
    } else {
      noteParts.push("Single-sample mode: within-sample absolute z-score");
    }
    resultNote.textContent = noteParts.join(" · ");

    resultsBody.innerHTML = "";
    result.scores.forEach((score) => {
      const row = document.createElement("tr");

      const sampleCell = document.createElement("td");
      sampleCell.textContent = score.sample;

      const activityCell = document.createElement("td");
      activityCell.textContent = formatNumber(score.activity, 6);
      activityCell.className =
        Number.isFinite(score.activity) && score.activity >= 0
          ? "score-pos"
          : Number.isFinite(score.activity)
            ? "score-neg"
            : "";

      const positiveCell = document.createElement("td");
      positiveCell.textContent = formatNumber(score.positive, 6);
      const negativeCell = document.createElement("td");
      negativeCell.textContent = formatNumber(score.negative, 6);

      row.append(sampleCell, activityCell, positiveCell, negativeCell);
      resultsBody.appendChild(row);
    });

    resultsSection.classList.remove("hidden");
    downloadBtn.disabled = false;
    resultsSection.scrollIntoView({ behavior: "smooth", block: "start" });
  }

  async function handleFile(file) {
    clearWarning();
    try {
      const text = await file.text();
      const parsed = parseMatrix(text);
      const matched = matchGeneSets(parsed);
      renderParsed(file, parsed, matched);
    } catch (error) {
      showWarning(error.message || "Could not parse the file. Please check the format.");
      resetToUpload();
    }
  }

  function resetToUpload() {
    state.parsed = null;
    state.result = null;
    state.file = null;
    fileInput.value = "";
    parsePanel.classList.add("hidden");
    resultsSection.classList.add("hidden");
    computeBtn.disabled = true;
    downloadBtn.disabled = true;
  }

  dropzone.addEventListener("click", () => fileInput.click());
  dropzone.addEventListener("keydown", (event) => {
    if (event.key === "Enter" || event.key === " ") {
      event.preventDefault();
      fileInput.click();
    }
  });

  ["dragenter", "dragover"].forEach((eventName) => {
    dropzone.addEventListener(eventName, (event) => {
      event.preventDefault();
      dropzone.classList.add("dragover");
    });
  });

  ["dragleave", "drop"].forEach((eventName) => {
    dropzone.addEventListener(eventName, (event) => {
      event.preventDefault();
      dropzone.classList.remove("dragover");
    });
  });

  dropzone.addEventListener("drop", (event) => {
    const file = event.dataTransfer.files && event.dataTransfer.files[0];
    if (file) handleFile(file);
  });

  fileInput.addEventListener("change", () => {
    const file = fileInput.files && fileInput.files[0];
    if (file) handleFile(file);
  });

  resetBtn.addEventListener("click", () => {
    resetToUpload();
    clearWarning();
  });

  downloadExampleBtn.addEventListener("click", () => {
    downloadText("example_expression_matrix.tsv", makeExampleMatrix());
  });

  exampleBtn.addEventListener("click", () => {
    loadTextAsFile("example_expression_matrix.tsv", makeExampleMatrix());
  });

  computeBtn.addEventListener("click", () => {
    if (!state.parsed) return;
    clearWarning();
    try {
      const result = computeWntScores(state.parsed);
      renderResults(result);
    } catch (error) {
      showWarning(error.message || "Calculation failed. Please check the input data.");
    }
  });

  downloadBtn.addEventListener("click", () => {
    if (!state.result) return;
    downloadText("WntQuant_activity_scores.txt", buildTxt(state.result));
  });

  computeBtn.disabled = true;
  downloadBtn.disabled = true;
})();
