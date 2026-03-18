"""Module for creating latex table from empagliflozin computational models table."""

import re
from pathlib import Path
import pandas as pd

# Characters that need escaping in LaTeX
LATEX_ESCAPE_MAP = [
    ('\\', r'\textbackslash{}'),
    ('&', r'\&'),
    ('%', r'\%'),
    ('$', r'\$'),
    ('#', r'\#'),
    ('_', r'\_'),
    ('{', r'\{'),
    ('}', r'\}'),
    ('~', r'\textasciitilde{}'),
    ('^', r'\textasciicircum{}'),
]

# Unicode characters replaced with LaTeX math equivalents
UNICODE_REPLACE_MAP = [
    ('≈', r'$\approx$'),
    ('→', r'$\rightarrow$'),
]


def escape_latex(text: str) -> str:
    for char, replacement in UNICODE_REPLACE_MAP:
        text = text.replace(char, replacement)
    for char, replacement in LATEX_ESCAPE_MAP:
        text = text.replace(char, replacement)
    # Restore math snippets that had their $ escaped
    text = text.replace(r'\$\approx\$', r'$\approx$')
    text = text.replace(r'\$\rightarrow\$', r'$\rightarrow$')
    return text


def urls_to_named_hrefs(text: str) -> str:
    named_pattern = re.compile(
        r'(GitHub|Zenodo|PK-DB)\s*:\s*(https?://[^\s,;]+)',
        re.IGNORECASE,
    )

    def replace_named(match):
        label = match.group(1)
        url = match.group(2).rstrip('.,;:)')
        return rf'\href{{{url}}}{{{label}}}'

    text = named_pattern.sub(replace_named, text)

    # Replace any remaining bare URLs
    bare_pattern = re.compile(r'(https?://[^\s,;]+)')

    def replace_bare(match):
        url = match.group(1).rstrip('.,;:)')
        return rf'\href{{{url}}}{{{url}}}'

    text = bare_pattern.sub(replace_bare, text)
    return text


def process_cell(text: str) -> str:
    splitter = re.compile(
        r'((?:GitHub|Zenodo|PK-DB)\s*:?\s*https?://[^\s,;]+|https?://[^\s,;]+)',
        re.IGNORECASE,
    )
    parts = splitter.split(text)
    result = []
    for i, part in enumerate(parts):
        if i % 2 == 0:
            result.append(escape_latex(part))
        else:
            named = re.match(
                r'(GitHub|Zenodo|PK-DB)\s*:?\s*(https?://[^\s,;]+)',
                part,
                re.IGNORECASE,
            )
            if named:
                label = named.group(1)
                url = named.group(2).rstrip('.,;:)')
                result.append(rf'\href{{{url}}}{{{label}}}')
            else:
                url = part.rstrip('.,;:)')
                result.append(rf'\href{{{url}}}{{{url}}}')
    joined = ''.join(result)
    import re as _re
    joined = _re.sub(
        r'(\\href\{[^}]+\}\{[^}]+\})\s+(\\href\{)',
        r'\1, \2',
        joined,
    )
    return joined


def create_latex_table(df: pd.DataFrame) -> str:
    """Transform DataFrame into LaTeX table."""

    latex_header = r"""
\begin{landscape}
\begin{table}[H]
\centering
\tabcolsep=3pt
\renewcommand{\arraystretch}{1.3}
\tiny
\begin{threeparttable} 

\caption{\scriptsize{\textbf{Summary of published computational models for empagliflozin.} 
Overview of published computational models including model type, software/platform, 
reproducibility criteria (open software, open model, open code, open data, open license,
reproducibility, FAIR, long-term storage), resources, clinical data sources, and model scope.}}

\begin{tabularx}{\linewidth}{
    >{\centering\arraybackslash}m{1.5cm}  % Study
    >{\centering\arraybackslash}m{0.8cm}  % PubMed ID
    >{\centering\arraybackslash}m{1.5cm}  % Model Type
    >{\centering\arraybackslash}m{1.3cm}  % Platform/Software
    >{\centering\arraybackslash}m{0.7cm}  % Open Software
    >{\centering\arraybackslash}m{0.7cm}  % Open Model
    >{\centering\arraybackslash}m{0.7cm}  % Open Code
    >{\centering\arraybackslash}m{0.7cm}  % Open Data
    >{\centering\arraybackslash}m{0.7cm}  % Open License
    >{\centering\arraybackslash}m{0.7cm}  % Reproducibility
    >{\centering\arraybackslash}m{0.7cm}  % FAIR
    >{\centering\arraybackslash}m{0.7cm}  % Longterm Storage
    >{\centering\arraybackslash}m{2.3cm}  % Resources
    >{\centering\arraybackslash}m{1.8cm}  % Studies
    >{\centering\arraybackslash}m{2.3cm}  % Clinical Data Used
    >{\centering\arraybackslash}m{4.4cm}  % Scope
}
\toprule
""".strip('\n')

    column_names = df.columns.tolist()
    header_cells = []
    for col in column_names:
        if col == 'Reproducibility':
            header_cells.append(r'\textbf{Reprodu\-cibility}')
        else:
            header_cells.append(f"\\textbf{{{escape_latex(col)}}}")
    column_headers = " & ".join(header_cells) + r" \\"

    latex_body = f"{column_headers}\n\\midrule\n"

    # Columns that get Yes/No cell colouring
    color_columns = {
        'Open Software': 4,
        'Open Model': 5,
        'Open Code': 6,
        'Open Data': 7,
        'Open License': 8,
        'Reproducibility': 9,
        'FAIR': 10,
        'Longterm Storage': 11,
    }

    n_rows = len(df)

    for row_idx, (_, row) in enumerate(df.iterrows()):
        raw_values = list(row.astype(str).values)

        # Replace 'nan' / empty with empty string
        raw_values = ['' if v in ('nan', 'NaN') else v for v in raw_values]

        # Process each cell: escape LaTeX + convert URLs
        values = [process_cell(v) for v in raw_values]

        is_last_row = (row_idx == n_rows - 1)

        # Column 0 (Study): add \cite (skip for last row)
        cite_key = raw_values[0]
        if not is_last_row:
            values[0] = f"{values[0]} \\cite{{{cite_key}}}"

        # Column 1 (PubMed ID): clickable link — skip for last row
        pmid = raw_values[1].strip()
        if pmid and pmid not in ('-', '') and not is_last_row:
            values[1] = (
                rf"\href{{https://pubmed.ncbi.nlm.nih.gov/{pmid}/}}{{{pmid}}}"
            )

        # Apply Yes/No cell colouring
        for col_name, col_idx in color_columns.items():
            if col_idx < len(values):
                val = raw_values[col_idx].strip()
                if val == 'Yes':
                    values[col_idx] = r'\cellcolor{green!20}Yes'
                elif val == 'No':
                    values[col_idx] = r'\cellcolor{red!20}No'

        if is_last_row:
            latex_body += r"\midrule" + "\n"
        row_str = " & ".join(values) + r" \\"
        latex_body += row_str + "\n"
        latex_body += r"\addlinespace[1pt]" + "\n"

    latex_footer = r"""
\bottomrule
\end{tabularx}

\end{threeparttable} 
\end{table}
\end{landscape}
""".strip('\n')

    full_latex = latex_header + "\n" + latex_body + latex_footer
    return full_latex


if __name__ == "__main__":
    tsv_path = Path(__file__).parent / 'computational_models.tsv'
    latex_path = Path(__file__).parent / 'computational_models.tex'

    df = pd.read_csv(tsv_path, sep="\t")
    df = df.fillna("")

    latex_str = create_latex_table(df)

    with open(latex_path, 'w') as f_tex:
        f_tex.write(latex_str)

    print(f"LaTeX table written to: {latex_path}")