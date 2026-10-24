"""Module for creating latex table from empagliflozin computational models table."""

import re
from pathlib import Path

import pandas as pd

# Characters that need escaping in LaTeX
LATEX_ESCAPE_MAP = [
    ("\\", r"\textbackslash{}"),
    ("&", r"\&"),
    ("%", r"\%"),
    ("$", r"\$"),
    ("#", r"\#"),
    ("_", r"\_"),
    ("{", r"\{"),
    ("}", r"\}"),
    ("~", r"\textasciitilde{}"),
    ("^", r"\textasciicircum{}"),
]

# Unicode characters replaced with LaTeX math equivalents
UNICODE_REPLACE_MAP = [
    ("≈", r"$\approx$"),
    ("→", r"$\rightarrow$"),
]

# a link with a label (`GitHub: <url>`) or a bare link
NAMED_URL_PATTERN = r"(GitHub|Zenodo|PK-DB)\s*:?\s*(https?://[^\s,;]+)"
URL_SPLIT_PATTERN = re.compile(
    r"((?:GitHub|Zenodo|PK-DB)\s*:?\s*https?://[^\s,;]+|https?://[^\s,;]+)",
    re.IGNORECASE,
)

# The table is written line by line: some lines end in a space, which is kept so
# that the generated table does not change.
LATEX_HEADER = "\n".join(
    [
        r"\begin{landscape}",
        r"\begin{table}[H]",
        r"\centering",
        r"\tabcolsep=3pt",
        r"\renewcommand{\arraystretch}{1.3}",
        r"\tiny",
        r"\begin{threeparttable} ",
        "",
        r"\caption{\scriptsize{\textbf{Summary of published computational models for empagliflozin.} ",
        "Overview of published computational models including model type, software/platform, ",
        "reproducibility criteria (open software, open model, open code, open data, open license,",
        "reproducibility, FAIR, long-term storage), resources, clinical data sources, and model scope.}}",
        "",
        r"\begin{tabularx}{\linewidth}{",
        r"    >{\centering\arraybackslash}m{1.5cm}  % Study",
        r"    >{\centering\arraybackslash}m{0.8cm}  % PubMed ID",
        r"    >{\centering\arraybackslash}m{1.5cm}  % Model Type",
        r"    >{\centering\arraybackslash}m{1.3cm}  % Platform/Software",
        r"    >{\centering\arraybackslash}m{0.7cm}  % Open Software",
        r"    >{\centering\arraybackslash}m{0.7cm}  % Open Model",
        r"    >{\centering\arraybackslash}m{0.7cm}  % Open Code",
        r"    >{\centering\arraybackslash}m{0.7cm}  % Open Data",
        r"    >{\centering\arraybackslash}m{0.7cm}  % Open License",
        r"    >{\centering\arraybackslash}m{0.7cm}  % Reproducibility",
        r"    >{\centering\arraybackslash}m{0.7cm}  % FAIR",
        r"    >{\centering\arraybackslash}m{0.7cm}  % Longterm Storage",
        r"    >{\centering\arraybackslash}m{2.3cm}  % Resources",
        r"    >{\centering\arraybackslash}m{1.8cm}  % Studies",
        r"    >{\centering\arraybackslash}m{2.3cm}  % Clinical Data Used",
        r"    >{\centering\arraybackslash}m{4.4cm}  % Scope",
        "}",
        r"\toprule",
    ]
)
LATEX_FOOTER = "\n".join(
    [
        r"\bottomrule",
        r"\end{tabularx}",
        "",
        r"\end{threeparttable} ",
        r"\end{table}",
        r"\end{landscape}",
    ]
)

# Columns that get Yes/No cell colouring, by index
COLOR_COLUMNS = {
    "Open Software": 4,
    "Open Model": 5,
    "Open Code": 6,
    "Open Data": 7,
    "Open License": 8,
    "Reproducibility": 9,
    "FAIR": 10,
    "Longterm Storage": 11,
}


def escape_latex(text: str) -> str:
    """Escape the special characters of LaTeX, unicode symbols become math."""
    for char, replacement in UNICODE_REPLACE_MAP:
        text = text.replace(char, replacement)
    for char, replacement in LATEX_ESCAPE_MAP:
        text = text.replace(char, replacement)
    # Restore math snippets that had their $ escaped
    text = text.replace(r"\$\approx\$", r"$\approx$")
    return text.replace(r"\$\rightarrow\$", r"$\rightarrow$")


def process_cell(text: str) -> str:
    r"""Escape the text of a cell and convert its links to `\href`.

    A link with a label (`GitHub: <url>`) gets the label as text, a bare link
    the link; consecutive links are separated by a comma.
    """
    parts = URL_SPLIT_PATTERN.split(text)
    result = []
    for i, part in enumerate(parts):
        if i % 2 == 0:
            result.append(escape_latex(part))
        else:
            named = re.match(NAMED_URL_PATTERN, part, re.IGNORECASE)
            if named:
                label = named.group(1)
                url = named.group(2).rstrip(".,;:)")
                result.append(rf"\href{{{url}}}{{{label}}}")
            else:
                url = part.rstrip(".,;:)")
                result.append(rf"\href{{{url}}}{{{url}}}")
    joined = "".join(result)
    return re.sub(
        r"(\\href\{[^}]+\}\{[^}]+\})\s+(\\href\{)",
        r"\1, \2",
        joined,
    )


def create_latex_table(df: pd.DataFrame) -> str:
    """Transform DataFrame into LaTeX table."""
    column_names = df.columns.tolist()
    header_cells = []
    for col in column_names:
        if col == "Reproducibility":
            header_cells.append(r"\textbf{Reprodu\-cibility}")
        else:
            header_cells.append(f"\\textbf{{{escape_latex(col)}}}")
    column_headers = " & ".join(header_cells) + r" \\"

    latex_body = f"{column_headers}\n\\midrule\n"

    n_rows = len(df)

    for row_idx, (_, row) in enumerate(df.iterrows()):
        raw_values = list(row.astype(str).values)

        # Replace 'nan' / empty with empty string
        raw_values = ["" if v in ("nan", "NaN") else v for v in raw_values]

        # Process each cell: escape LaTeX + convert URLs
        values = [process_cell(v) for v in raw_values]

        is_last_row = row_idx == n_rows - 1

        # Column 0 (Study): add \cite (skip for last row)
        cite_key = raw_values[0]
        if not is_last_row:
            values[0] = f"{values[0]} \\cite{{{cite_key}}}"

        # Column 1 (PubMed ID): clickable link, skip for last row
        pmid = raw_values[1].strip()
        if pmid and pmid not in ("-", "") and not is_last_row:
            values[1] = rf"\href{{https://pubmed.ncbi.nlm.nih.gov/{pmid}/}}{{{pmid}}}"

        # Apply Yes/No cell colouring
        for col_idx in COLOR_COLUMNS.values():
            if col_idx < len(values):
                val = raw_values[col_idx].strip()
                if val == "Yes":
                    values[col_idx] = r"\cellcolor{green!20}Yes"
                elif val == "No":
                    values[col_idx] = r"\cellcolor{red!20}No"

        if is_last_row:
            latex_body += r"\midrule" + "\n"
        row_str = " & ".join(values) + r" \\"
        latex_body += row_str + "\n"
        latex_body += r"\addlinespace[1pt]" + "\n"

    return LATEX_HEADER + "\n" + latex_body + LATEX_FOOTER


if __name__ == "__main__":
    tsv_path = Path(__file__).parent / "computational_models.tsv"
    latex_path = Path(__file__).parent / "computational_models.tex"

    df = pd.read_csv(tsv_path, sep="\t")
    df = df.fillna("")

    latex_str = create_latex_table(df)

    with open(latex_path, "w") as f_tex:
        f_tex.write(latex_str)

    print(f"LaTeX table written to: {latex_path}")
