from pathlib import Path


def require(condition, message):
    if not condition:
        raise SystemExit(message)


def strip_fortran_comments(source):
    cleaned_lines = []
    in_single_quote = False
    in_double_quote = False

    for line in source.splitlines(keepends=True):
        comment_start = None
        index = 0
        while index < len(line):
            character = line[index]
            next_character = line[index + 1] if index + 1 < len(line) else ""
            if character == "'" and not in_double_quote:
                if in_single_quote and next_character == "'":
                    index += 2
                    continue
                in_single_quote = not in_single_quote
            elif character == '"' and not in_single_quote:
                if in_double_quote and next_character == '"':
                    index += 2
                    continue
                in_double_quote = not in_double_quote
            elif character == "!" and not in_single_quote and not in_double_quote:
                comment_start = index
                break
            index += 1

        if comment_start is None:
            cleaned_lines.append(line)
        else:
            if line.endswith("\r\n"):
                line_end = "\r\n"
            elif line.endswith(("\n", "\r")):
                line_end = line[-1]
            else:
                line_end = ""
            cleaned_lines.append(line[:comment_start] + line_end)

    return "".join(cleaned_lines)


def marker_positions(source, marker):
    positions = []
    search_from = 0
    while True:
        position = source.find(marker, search_from)
        if position < 0:
            return positions
        positions.append(position)
        search_from = position + len(marker)


def extract_section(source, start_marker, end_marker, section_name):
    start_matches = marker_positions(source, start_marker)
    require(
        start_matches,
        f"{section_name}: missing start marker {start_marker!r}",
    )
    require(
        len(start_matches) == 1,
        f"{section_name}: ambiguous repeated start marker {start_marker!r}",
    )

    end_matches = marker_positions(source, end_marker)
    require(
        end_matches,
        f"{section_name}: missing end marker {end_marker!r}",
    )
    require(
        len(end_matches) == 1,
        f"{section_name}: ambiguous repeated end marker {end_marker!r}",
    )

    start = start_matches[0] + len(start_marker)
    end = end_matches[0]
    require(
        start <= end,
        f"{section_name}: end marker precedes start marker",
    )
    return source[start:end]


steady_raw = Path("src/mod_SteadyImplicit.f90").read_text()
solve_raw = Path("src/mod_Solve.f90").read_text()
unsteady_raw = Path("src/mod_UnsteadyImplicit.f90").read_text()
diagnostics_raw = Path("src/mod_ImplicitDiagnostics.f90").read_text()
steady = strip_fortran_comments(steady_raw)
solve = strip_fortran_comments(solve_raw)
unsteady = strip_fortran_comments(unsteady_raw)
diagnostics = strip_fortran_comments(diagnostics_raw)

find_body = extract_section(
    steady,
    "subroutine FindAdmissibleMicroAlpha(",
    "end subroutine FindAdmissibleMicroAlpha",
    "FindAdmissibleMicroAlpha",
)
require(
    "public :: SelectGlobalAcceptedStateAdmissibilityMaxima, WriteAcceptedStateAdmissibilityDiagnostics"
    in diagnostics,
    "accepted-state diagnostic writer must be public",
)
require(
    "call FindAdmissibleMicroAlpha(implicit_alpha,ok,.false.)" in steady,
    "ImplicitSteady must disable the compatibility negative-fraction gate explicitly",
)
require(
    "enforce_negative_fraction_active=.true." in find_body,
    "FindAdmissibleMicroAlpha must default enforce_negative_fraction to .true.",
)
require(
    "call FindAdmissibleMicroAlpha(micro_alpha,ok)" in unsteady,
    "ImplicitUnsteady must retain its no-argument admissibility call",
)

output_section = extract_section(
    solve_raw,
    "!Output the residuals",
    "!Determine whether this accepted state is terminal",
    "Nstep_show output dispatch",
)
output_section = strip_fortran_comments(output_section)
require(
    "call WriteAcceptedStateAdmissibilityDiagnostics('ImplicitSteady')" in output_section,
    "ImplicitSteady must call the accepted-state diagnostic writer",
)
require(
    "call WriteAcceptedStateAdmissibilityDiagnostics('Explicit')" in output_section,
    "Explicit must call the accepted-state diagnostic writer",
)
unsteady_output_branch = extract_section(
    output_section,
    "elseif(IsImplicitUnsteady()) then",
    "else\n",
    "ImplicitUnsteady output branch",
)
require(
    "WriteAcceptedStateAdmissibilityDiagnostics" not in unsteady_output_branch,
    "ImplicitUnsteady output branch must not call the accepted-state diagnostic writer",
)

print("PASS: Shakhov admissibility policy source contract")
