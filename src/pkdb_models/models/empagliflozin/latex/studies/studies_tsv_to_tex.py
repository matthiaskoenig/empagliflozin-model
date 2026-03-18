from pathlib import Path
import pandas as pd


def create_latex_table(df: pd.DataFrame) -> str:
    """Transform .tsv into .tex table."""

    latex_header = r"""\begin{landscape}
\begin{table}[H]
\tiny
\centering
\tabcolsep=3.5pt\relax
\caption{\textbf{Summary of studies for modeling.}
Overview of study identifiers, PubMed IDs, PK-DB IDs, administered substance, route, dosing,
and subject characteristics, including health status (\emph{H}), renal impairment (\emph{RI}),
hepatic impairment (\emph{HI}), fasting status, type 2 diabetes (\emph{T2}),
urinary glucose excretion (\emph{UGE}), fasting plasma glucose (\emph{FPG}),
renal threshold for glucose (\emph{RTG}).
\emph{EMP P} = empagliflozin plasma, \emph{EMP U} = empagliflozin urine,
\emph{EMP B} = empagliflozin bile, \emph{EMP F} = empagliflozin feces}
\label{table:curated_data_overview_empa}
\rowcolors{2}{white}{gray!8}
\begin{tabularx}{\linewidth}{
    >{\raggedright\arraybackslash}m{2.3cm}  % Study
    >{\centering\arraybackslash}m{1.2cm}  % PubMed
    >{\centering\arraybackslash}m{1.4cm}  % PK-DB ID
    >{\centering\arraybackslash}m{2.0cm}  % Substance
    >{\centering\arraybackslash}m{0.5cm}  % Route
    >{\centering\arraybackslash}m{1.3cm}  % Dosing
    >{\centering\arraybackslash}m{1.5cm}  % Dose [mg]
    >{\centering\arraybackslash}m{0.7cm}  % Fast
    c                                      % EMP P
    c                                      % EMP U
    c                                      % EMP B
    c                                      % EMP F
    >{\centering\arraybackslash}m{0.9cm}  % H
    >{\centering\arraybackslash}m{0.9cm}  % RI
    >{\centering\arraybackslash}m{0.9cm}  % HI
    >{\centering\arraybackslash}m{0.9cm}  % T2
    >{\centering\arraybackslash}m{0.9cm}  % UGE
    >{\centering\arraybackslash}m{0.9cm}  % FPG
    >{\centering\arraybackslash}m{0.7cm}  % RTG
}
\arrayrulecolor{black}\toprule"""

    # Rename columns
    df = df.rename(columns={
        "study": "Study",
        "pmid": "PubMed",
        "pkdb": "PK-DB",
        "substance": "Substance",
        "route": "Route",
        "dosing": "Dosing",
        "dose": "Dose [mg]",
        "healthy": "H",
        "renal impairment": "RI",
        "hepatic impairment": "HI",
        "t2dm": "T2",
        "empagliflozin plasma": "EMP P",
        "empagliflozin urine": "EMP U",
        "empagliflozin bile": "EMP B",
        "empagliflozin feces": "EMP F",
        "fasting": "Fast",
        "UGE": "UGE",
        "FPG": "FPG",
        "RTg": "RTG",
    })

    columns_bool = [
        "H", "RI", "HI", "T2",
        "EMP P", "EMP U", "EMP B", "EMP F",
        "UGE", "FPG", "RTG",
    ]
    for col in columns_bool:
        df[col] = df[col].apply(lambda x: r"\checkmark" if str(x).strip().upper() == "TRUE" else "")

    # Select and order columns for output
    output_columns = [
        "Study", "PubMed", "PK-DB", "Substance", "Route", "Dosing", "Dose [mg]", "Fast",
        "EMP P", "EMP U", "EMP B", "EMP F",
        "H", "RI", "HI", "T2",
        "UGE", "FPG", "RTG",
    ]
    df = df[output_columns]

    # Header row
    column_headers = " & ".join([f"\\textbf{{{col}}}" for col in df.columns.tolist()]) + r" \\"
    latex_body = (
        "\n"
        "\\hiderowcolors\n"
        f"{column_headers}\n"
        "\\showrowcolors\n"
        "\\arrayrulecolor{black}\\midrule\n"
    )

    for _, row in df.iterrows():
        values = list(row.astype(str).values)

        # Study -> \cite
        values[0] = f"{values[0]} \\cite{{{values[0]}}}"

        # PubMed -> href
        pmid = values[1].split(".")[0] if "." in values[1] else values[1]
        values[1] = f"\\href{{https://pubmed.ncbi.nlm.nih.gov/{pmid}/}}{{{pmid}}}"

        # PK-DB -> href
        pkdb = values[2]
        values[2] = f"\\href{{https://identifiers.org/pkdb:{pkdb}}}{{{pkdb}}}"

        row_str = " & ".join(
            v.replace("TRUE", r"\checkmark").replace("-", "")
            for v in values
        ) + r" \\"

        # fix for long name
        row_str = row_str.replace("vanderAartvanderBeek2020 ",
                                  r"vanderAart\-vanderBeek2020 ")

        latex_body += row_str + "\n"

    latex_footer = (
        "\\arrayrulecolor{black}\\bottomrule\n"
        "\\end{tabularx}\n"
        "\\end{table}\n"
        "\\end{landscape}"
    )

    return latex_header + latex_body + latex_footer


if __name__ == "__main__":
    tsv_path = Path(__file__).parent / 'studies_for_modeling.tsv'
    latex_path = Path(__file__).parent / 'studies_for_modeling.tex'

    columns = [
        "study",
        "pmid",
        "pkdb",

        "substance",
        "route",
        "dosing",
        "dose",
        "combination",
        "species",
        "fasting",

        "healthy",
        "renal impairment",
        "hepatic impairment",
        "t2dm",

        "empagliflozin plasma",
        "empagliflozin urine",
        "empagliflozin bile",
        "empagliflozin feces",

        "UGE",
        "FPG",
        "glucose plasma",
        "RTG",
    ]

    df = pd.read_csv(
        tsv_path, sep="\t", skiprows=1,
        usecols=columns,
        engine="python",
        dtype={"pmid": str},
    )[columns]

    df = df.fillna("")

    latex_str = create_latex_table(df)

    with open(latex_path, 'w') as f_tex:
        f_tex.write(latex_str)