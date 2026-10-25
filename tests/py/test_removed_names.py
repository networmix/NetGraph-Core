"""One rule for the names this release removed (the old spellings of removed_names.txt's tables:
the two-column table of public names and the three-column tables of members and C++
declarations), over every kind of text: the public docstrings (the callable behind a staticmethod
or classmethod included), the C++ comments and string literals of include/, src/, bindings/,
tests/cpp/ and dev/perf/ (not the frozen fixture replay, whose text names the frozen data's keys),
the comments, docstrings and string literals of the Python files and stubs (.py, .pyi) of
python/, tests/ and dev/ (not the frozen reference under tests/reference/ and not the fixture
replay's key mapping), the comments of the shell scripts, workflows, Makefile, CMake and TOML
files, and the prose of the Markdown documents (fenced code excluded: the docs tests run it
against the current API; not CHANGELOG.md, which names the removed 0.11 API, nor the frozen
reference's README).

1. A removed name that the API no longer has anywhere and that is not an ordinary word of this
   domain's prose is a stale name wherever it stands as a word.
2. A removed name that is an ordinary word (ENGLISH) or still a name elsewhere in the API is a
   stale name as a word in the public text of its old or new owner: the owner's docstrings, its
   declaration and the comment above it in a public header, the paragraphs of the package
   docstring and of the docs that name it.
3. In a C++ comment, a snake_case name with an underscore, or a name written as a call, that no code
   of the repository declares or uses (C++ and Python identifiers, the API, file names) is a stale
   name: a renamed internal helper or member (removed_names.txt lists only the public ones).
4. No code identifier (C++ outside comments and string literals, not a std:: name; Python names,
   attributes, definitions, parameters and keywords) of src/, include/, bindings/, tests/ (not the
   frozen fixture replay or reference), python/ and dev/ holds a rule-1 name as itself or as a
   snake_case part (`bundle_of`), unless that name is part of a current API name
   (`edge_capacity_units` holds `capacity_units`). Its reviewed rows have the place `code:<path>`
   and the identifier itself as their substring. The parameter shorthands of SHORTHANDS stay as
   ordinary locals.

A word is the name or its plural, capitalized too (a sentence-initial `Placed`), beside a hyphen
(`record-bundle`), also after a dot (`ledger.allocated()`; a one-letter name after
a letter and a dot is "e.g."), also beside a "/" once the file paths are taken out of the line
(`src/dst`). Every remaining use is one row of removed_names_reviewed.txt, "place | word | substring
of the line | reason": it excuses that word on the one line of that place (a docstring's owner such
as `PlacementResult.call_delivered`, or a file path) that holds the substring, and nothing else."""

from __future__ import annotations

import ast
import inspect
import io
import os
import re
import sys
import tokenize
from pathlib import Path

import pytest

import netgraph_core as ngc

ROOT = Path(__file__).resolve().parents[2]
REVIEWED_FILE = Path(__file__).with_name("removed_names_reviewed.txt")
REMOVED_FILE = Path(__file__).with_name("removed_names.txt")

# Planted text for the rule's own controls: {relative path: text read in place of the file's}.
PLANTED = {}


def read(path):
    """A file's text (a planted text in its place when the controls set one)."""
    path = Path(path)
    rel = str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else None
    return PLANTED[rel] if rel in PLANTED else path.read_text()


SKIPPED_DIRS = {
    ".git",
    "build",
    "venv",
    "design",
    "_deps",
    "__pycache__",
    ".ruff_cache",
    ".lake",
}


def tree_files():
    """Every file of the checkout outside build output, environments, Lake's dependencies and
    builds (proofs/.lake), and the design tree."""
    for top, dirs, files in os.walk(ROOT):
        dirs[:] = sorted(
            d for d in dirs if d not in SKIPPED_DIRS and not d.endswith(".noindex")
        )
        for name in sorted(files):
            yield Path(top) / name


# Removed names that are ordinary words of this domain's prose and common in it (each used more
# than 16 times across the texts), with the meaning they keep. A rarer one (allocated, placed,
# volume, ...) is not listed: like any removed name it is stale wherever it stands, and each of its
# ordinary uses is a reviewed line.
ENGLISH = {
    "bottleneck": "a path's smallest score under a HIGHER_* tie rule",
    "capacity": "an edge's capacity",
    "count": "the tunnel count",
    "dropped": "traffic dropped on a link or at a black hole",
    "drops": "losses on links or at a source",
    "entry": "an entry of a program, an array or a CSR list",
    "index": "a tunnel's entry index",
    "loss": "LOSSLESS admits a total 'with no loss'",
    "offsets": "CSR offsets",
    "owner": "the owner of a hold (another policy)",
    "placement": "placing traffic",
    "planned": "T2's planned bandwidth",
    "remaining": "what remains",
    "requested": "a requested amount",
    "via": "the via edge of a DAG entry",
    "weights": "WEIGHTED's member weights",
}


def _split(spelling):
    """(owner, leaf) of an old or new spelling: `Owner.leaf`, `Owner::leaf`, `Owner.method:param`,
    `Owner::method:param`, `function:param`, `Owner(leaf=...)`; a bare name is (None, name)."""
    m = re.match(r"(\w+)\((\w+)=", spelling)
    if m:
        return m.group(1), m.group(2)
    spelling = spelling.split("(")[0]
    m = re.match(r"(\w+)(?:\.|::)(?:\w+:)?(\w+)$", spelling)
    if m:
        return m.group(1), m.group(2)
    m = re.match(r"(\w+):(\w+)$", spelling)
    if m:
        return m.group(1), m.group(2)
    return None, re.split(r"::|\.", spelling)[-1]


def removed():
    """(old owner, leaf, new owners) for every old spelling of removed_names.txt's tables (files
    excluded)."""
    out = []
    for row in REMOVED_FILE.read_text().splitlines():
        cells = [c.strip() for c in row.strip().strip("|").split("|")]
        if (
            not row.startswith("|")
            or len(cells) < 2
            or cells[0] in ("Kind", "file", "Name")
            or set(cells[0]) <= {"-"}
        ):
            continue
        # | Name | Replacement | (public names) or | Kind | Name | Replacement |.
        old_cell, new_cell = (cells[0], cells[1]) if len(cells) == 2 else cells[1:3]
        news = {_split(n)[0] for n in re.findall(r"`([^`]+)`", new_cell)} - {None}
        for old in re.findall(r"`([^`]+)`", old_cell):
            owner, leaf = _split(old)
            if leaf and leaf != "__init__":
                out.append((owner, leaf, news))
    return out


def api_names():
    """Every name of the installed stub (classes, functions, members, parameters), and per owner
    (a class or a module function) its own names."""
    tree = ast.parse(read(Path(ngc.__file__).parent / "__init__.pyi"))
    names, per_owner = set(ngc.__all__), {}
    for node in tree.body:
        if isinstance(node, (ast.ClassDef, ast.FunctionDef)):
            own = per_owner.setdefault(node.name, set())
            for sub in ast.walk(node):
                if isinstance(sub, (ast.FunctionDef, ast.ClassDef)):
                    own.add(sub.name)
                elif isinstance(sub, ast.arg):
                    own.add(sub.arg)
                elif isinstance(sub, ast.AnnAssign) and isinstance(
                    sub.target, ast.Name
                ):
                    own.add(sub.target.id)
            names |= own | {node.name}
    return names, per_owner


# A name inside a path (src/flow_policy.cpp) is not the name: the word pattern excludes "/".

# A file path (with an extension, or a directory ending in "/") is not a name: it is taken out of a
# line before its words are read, so "src/flow_policy.cpp" names nothing and "src/dst" names both.
PATH = re.compile(
    r"(?:[\w.-]+/)+[\w.-]*\.(?:hpp|cpp|h|py|pyi|md|txt|json|sh|toml|cmake)\b|(?:[\w.-]+/)+(?!\w)"
)


def word(name):
    """The words of `name`: itself or its plural, capitalized too, after a dot too (not a
    one-letter name after a letter and a dot: "e.g."), beside a hyphen too (`record-bundle`)."""
    before = r"(?<!\w)" + (r"(?<![A-Za-z]\.)" if len(name) == 1 else "")
    # A capitalized one-letter name is a math symbol (G(z)), not the name.
    forms = {name, name[:1].upper() + name[1:]} if len(name) > 1 else {name}
    alternatives = "|".join(re.escape(f) for f in sorted(forms))
    return re.compile(rf"{before}(?:{alternatives})(?:e?s)?(?!\w)")


def _owners_named(lines, owners):
    return {
        o
        for o in owners
        if any(re.search(rf"(?<![\w.]){o}(?!\w)", line) for line in lines)
    }


# -- the texts: (place, line, owners, public) ------------------------------------------------------


def _doc_lines(obj):
    doc = getattr(obj, "__doc__", None) or ""
    # pybind11 prefixes each overload's text with its signature: names, not prose.
    return [
        line
        for line in doc.splitlines()
        if not re.match(r"^\s*(\d+\. )?\w+\(.*\) -> ", line)
    ]


def docstring_texts(owners=()):
    # The package docstring's paragraphs (blank-line separated, or a "- " item) belong to the
    # owners they name, as docs paragraphs do.
    out, paragraph = [], []
    for line in _doc_lines(ngc) + [""]:
        if (not line.strip() or line.lstrip().startswith("- ")) and paragraph:
            named = _owners_named(paragraph, owners)
            out += [("netgraph_core", text, named, True) for text in paragraph]
            paragraph = []
        if line.strip():
            paragraph.append(line)
    for name in ngc.__all__:
        obj = getattr(ngc, name)
        if isinstance(obj, str):
            continue
        out += [(name, line, {name}, True) for line in _doc_lines(obj)]
        if inspect.isclass(obj) and not issubclass(obj, BaseException):
            for member, value in vars(obj).items():
                if member.startswith("_") and member != "__init__":
                    continue
                if isinstance(value, (staticmethod, classmethod)):
                    value = (
                        value.__func__
                    )  # the function's own text, not the descriptor's help
                if isinstance(value, property) or callable(value):
                    out += [
                        (f"{name}.{member}", line, {name}, True)
                        for line in _doc_lines(value)
                    ]
    return out


def _comments_and_code(text, strings=None):
    """Per line its comment text (// and /* */, string literals skipped) and the file's code
    identifiers (outside comments and string literals); per line the string literals' text into
    `strings` when given."""
    comments, code = {}, []
    i, n, line = 0, len(text), 1
    while i < n:
        if text.startswith("//", i):
            j = text.find("\n", i)
            j = n if j < 0 else j
            comments[line] = comments.get(line, "") + text[i + 2 : j]
            i = j
        elif text.startswith("/*", i):
            j = text.find("*/", i + 2)
            j = n if j < 0 else j
            for k, part in enumerate(text[i + 2 : j].split("\n")):
                comments[line + k] = comments.get(line + k, "") + part
            line += text[i:j].count("\n")
            i = j + 2
        elif text[i] in "\"'":  # a string or character literal, to its closing quote
            quote, j = text[i], i + 1
            while j < n and text[j] != quote and text[j] != "\n":
                j += 2 if text[j] == "\\" else 1
            if strings is not None and quote == '"':
                strings.setdefault(line, []).append(text[i + 1 : j])
            i = j + 1
        else:
            if text[i] == "\n":
                line += 1
            code.append(text[i])
            i += 1
    return comments, set(re.findall(r"[A-Za-z_]\w*", "".join(code)))


def _owner_regions(lines, owners):
    """Per line of a public header, the owners whose text it is: a struct or class to its closing
    brace, a function's declaration, each with the comment block above it."""
    out = {}
    for i, line in enumerate(lines):
        m = re.match(
            r"^(\s*)(?:struct|class)\s+(?:\[\[[^\]]*\]\]\s*)?([A-Z]\w*)\b[^;]*$", line
        )
        names = []
        if m:
            end = i
            while end + 1 < len(lines) and not re.match(
                rf"^{m.group(1)}}};", lines[end]
            ):
                end += 1
            names.append((m.group(2), end))
        for owner in owners:
            if owner[:1].islower() and re.search(
                rf"^\s*(?:\[\[nodiscard\]\]\s*)?[\w:<>,\s*&]*\b{owner}\(", line
            ):
                end = i
                while end + 1 < len(lines) and not re.search(r"[;{]\s*$", lines[end]):
                    end += 1
                names.append((owner, end))
        for owner, end in names:
            start = i
            while start > 0 and lines[start - 1].strip().startswith("//"):
                start -= 1
            for k in range(start, end + 1):
                out.setdefault(k + 1, set()).add(owner)
    return out


FIXTURE_REPLAY = "tests/cpp/flow_policy_fixture_tests.cpp"


def _cpp_files():
    for directory in ("include", "src", "bindings", "tests/cpp", "dev/perf"):
        for path in sorted((ROOT / directory).rglob("*")):
            if (
                path.suffix in (".hpp", ".cpp", ".h")
                and str(path.relative_to(ROOT)) != FIXTURE_REPLAY
            ):
                yield path


FIXTURE_KEY_MAPPING = "tests/py/test_flow_policy_fixtures.py"


def _python_texts(path, skip_module_docstring=False, strings=None):
    """(line number, text) of a Python file's comments and docstrings; its other string constants'
    lines into `strings` when given."""
    text = read(path)
    out = [
        (tok.start[0], tok.string[1:])
        for tok in tokenize.generate_tokens(io.StringIO(text).readline)
        if tok.type == tokenize.COMMENT
    ]
    tree = ast.parse(text)
    docstrings = set()
    for node in ast.walk(tree):
        if not isinstance(
            node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)
        ):
            continue
        body = node.body
        if body and isinstance(body[0], ast.Expr):
            docstrings.add(id(body[0].value))
        if isinstance(node, ast.Module) and skip_module_docstring:
            continue
        body = node.body
        if (
            body
            and isinstance(body[0], ast.Expr)
            and isinstance(body[0].value, ast.Constant)
        ):
            if isinstance(body[0].value.value, str):
                first = body[0].lineno
                out += [
                    (first + k, line)
                    for k, line in enumerate(body[0].value.value.splitlines())
                ]
    if strings is not None:
        for node in ast.walk(tree):
            if (
                isinstance(node, ast.Constant)
                and isinstance(node.value, str)
                and id(node) not in docstrings
            ):
                strings += [
                    (node.lineno + k, line)
                    for k, line in enumerate(node.value.splitlines())
                ]
    return out


def _python_files():
    """The Python files and stubs of python/, tests/ and dev/ (not the design tree, not the frozen
    reference, not the fixture replay's key mapping)."""
    for directory in ("python", "tests", "dev"):
        for path in sorted((ROOT / directory).rglob("*")):
            rel = path.relative_to(ROOT)
            if (
                path.suffix in (".py", ".pyi")
                and "design" not in rel.parts
                and not str(rel).startswith("tests/reference/")
                and str(rel) != FIXTURE_KEY_MAPPING
            ):
                yield path


# The generated stub: its docstrings are the module's, read at run time with their owners.
GENERATED_STUB = "python/netgraph_core/__init__.pyi"


def python_texts():
    """The comments, docstrings and string constants of the Python files and stubs (the package's
    own docstring is read at run time; of the generated stub only the comments)."""
    out = []
    for path in _python_files():
        rel = str(path.relative_to(ROOT))
        strings = []
        skip = rel == "python/netgraph_core/__init__.py"
        texts = _python_texts(path, skip, strings)
        if rel == GENERATED_STUB:
            text = read(path)
            texts = [
                (tok.start[0], tok.string[1:])
                for tok in tokenize.generate_tokens(io.StringIO(text).readline)
                if tok.type == tokenize.COMMENT
            ]
            strings = []
        out += [(rel, line, set(), False) for _, line in sorted(texts + strings)]
    return out


# The files whose "#" comments are text: shell scripts, workflows, the Makefile, CMake and TOML.
HASH_COMMENTED = (".sh", ".yml", ".yaml", ".toml", ".cmake")
HASH_NAMES = ("Makefile", "CMakeLists.txt")


def hash_comment_texts():
    """The "#" comments (from a "#" at the line's start or after a blank) of the shell scripts,
    workflows, Makefile, CMake and TOML files outside build/, venv/ and the design tree."""
    out = []
    for path in tree_files():
        rel = path.relative_to(ROOT)
        if path.suffix not in HASH_COMMENTED and path.name not in HASH_NAMES:
            continue
        for line in read(path).splitlines():
            m = re.search(r"(?:^|\s)#(.*)$", line)
            if m:
                out.append((str(rel), m.group(1), set(), False))
    return out


def comment_texts(owners):
    out = []
    for path in _cpp_files():
        rel = str(path.relative_to(ROOT))
        text = read(path)
        comments, _ = _comments_and_code(text)
        public = rel.startswith("include/netgraph/core/") and "/detail/" not in rel
        regions = _owner_regions(text.splitlines(), owners) if public else {}
        out += [
            (rel, comment, regions.get(n, set()), public)
            for n, comment in sorted(comments.items())
        ]
    return out


def cpp_string_texts(docstrings):
    """The string literals of the C++ files, one text each, but a line of a public docstring
    (`docstrings`: read at run time with its owners)."""
    out = []
    for path in _cpp_files():
        strings = {}
        _comments_and_code(read(path), strings)
        rel = str(path.relative_to(ROOT))
        for _, literals in sorted(strings.items()):
            for text in literals:
                if text.replace("\\n", "").strip() not in docstrings:
                    out.append((rel, text, set(), False))
    return out


def _paragraphs(lines):
    """The prose paragraphs of a Markdown file (blank-line separated; fenced code skipped)."""
    fenced, paragraph = False, []
    for line in lines:
        if line.startswith("```"):
            fenced = not fenced
            if paragraph:
                yield paragraph
            paragraph = []
        elif fenced:
            continue
        elif line.strip():
            paragraph.append(line)
        elif paragraph:
            yield paragraph
            paragraph = []
    if paragraph:
        yield paragraph


# The documents read: every Markdown file but CHANGELOG.md (it names the removed 0.11 API) and the
# frozen reference's README.
def _markdown_files():
    for path in tree_files():
        rel = path.relative_to(ROOT)
        if path.suffix != ".md":
            continue
        if str(rel) == "CHANGELOG.md" or str(rel).startswith("tests/reference/"):
            continue
        yield path


def markdown_texts(owners):
    out = []
    for path in _markdown_files():
        rel = str(path.relative_to(ROOT))
        for paragraph in _paragraphs(read(path).splitlines()):
            named = _owners_named(paragraph, owners)
            out.extend((rel, line, named, True) for line in paragraph)
    return out


def code_names():
    """Every identifier the repository's code declares or uses: C++ code (outside comments and
    string literals) of the library, bindings, tests and dev/perf; Python names, attributes,
    parameters and identifier-like strings of python/, tests/ and dev/ (the frozen reference
    included); CMake's targets and options; the API; file stems."""
    out = set(api_names()[0])
    for path in _cpp_files():
        out |= _comments_and_code(read(path))[1]
    out |= set(re.findall(r"[A-Za-z_]\w*", (ROOT / "CMakeLists.txt").read_text()))
    for directory in ("python", "tests", "dev"):
        for path in sorted((ROOT / directory).rglob("*.py")):
            for node in ast.walk(ast.parse(path.read_text())):
                if isinstance(node, ast.Name):
                    out.add(node.id)
                elif isinstance(node, ast.Attribute):
                    out.add(node.attr)
                elif isinstance(node, (ast.FunctionDef, ast.ClassDef)):
                    out.add(node.name)
                elif isinstance(node, ast.arg):
                    out.add(node.arg)
                elif (
                    isinstance(node, ast.Constant)
                    and isinstance(node.value, str)
                    and node.value.isidentifier()
                ):
                    out.add(node.value)
    for path in tree_files():
        out.add(path.stem)
    return out


# A snake_case name with an underscore, or a name written as a call (clause 3).
IDENTIFIER = re.compile(
    r"(?<![\w.$-])([a-z][a-z0-9]*(?:_[a-z0-9]+)+(?![\w(])|[a-z][a-z0-9_]*(?=\())"
)


# The removed parameter shorthands kept as ordinary local names (not API): one reason each.
SHORTHANDS = {
    "src": "a local or parameter holding a source node",
    "dst": "a local or parameter holding a destination node",
    "g": "a local or parameter holding a graph",
    "opts": "a local or parameter holding an options value",
    "elig": "a local or parameter holding an Eligibility",
}


def _code_text(text):
    """The text outside comments and string literals (newlines kept)."""
    out, i, n = [], 0, len(text)
    while i < n:
        if text.startswith("//", i):
            j = text.find("\n", i)
            i = n if j < 0 else j
        elif text.startswith("/*", i):
            j = text.find("*/", i + 2)
            i = n if j < 0 else j + 2
        elif text[i] in "\"'":
            quote, j = text[i], i + 1
            while j < n and text[j] != quote and text[j] != "\n":
                j += 2 if text[j] == "\\" else 1
            out.append(" ")
            i = j + 1
        else:
            out.append(text[i])
            i += 1
    return "".join(out)


def _code_identifiers():
    """{place: identifiers} of the code clause 4 reads (place: "code:<path>")."""
    out = {}
    for path in _cpp_files():
        code = _code_text(read(path))
        names = set()
        for m in re.finditer(r"[A-Za-z_]\w*", code):
            if code[max(0, m.start() - 5) : m.start()] != "std::":
                names.add(m.group(0))
        out[f"code:{path.relative_to(ROOT)}"] = names
    for directory in ("python", "tests/py", "dev"):
        for path in sorted((ROOT / directory).rglob("*.py")):
            if "design" in path.relative_to(ROOT).parts:
                continue
            names = set()
            for node in ast.walk(ast.parse(path.read_text())):
                if isinstance(node, (ast.Import, ast.ImportFrom)):
                    continue
                for attr in ("id", "attr", "name", "arg"):
                    value = getattr(node, attr, None)
                    if isinstance(value, str):
                        names.add(value)
            out[f"code:{path.relative_to(ROOT)}"] = names
    return out


def code_findings(dead, api):
    """Clause 4: (place, name, identifier) for every code identifier holding a rule-1 name."""

    def holds(name, ident):
        return (
            ident == name
            or re.search(rf"(?:^|_){re.escape(name)}(?:_|$)", ident) is not None
        )

    names = {
        n for n in dead if n not in SHORTHANDS and not any(holds(n, a) for a in api)
    }
    return [
        (place, name, ident)
        for place, idents in _code_identifiers().items()
        for ident in idents
        for name in names
        if holds(name, ident)
    ]


def reviewed():
    """(place, word, substring, reason) per row of the reviewed file."""
    rows = []
    for row in REVIEWED_FILE.read_text().splitlines():
        if row.strip() and not row.startswith("#"):
            place, name, sub, reason = (x.strip() for x in row.split(" | ", 3))
            rows.append((place, name, sub, reason))
    return rows


def find_stale():
    """Every (place, name, line) the rule reports, before the reviewed lines."""
    names, per_owner = api_names()
    rows = removed()
    dead = {leaf for _, leaf, _ in rows if leaf not in names and leaf not in ENGLISH}
    owned = {}  # leaf -> owners whose public text must not use it
    for owner, leaf, news in rows:
        if leaf in dead:
            continue
        for o in ({owner} | news) - {None}:
            if leaf not in per_owner.get(o, set()):
                owned.setdefault(leaf, set()).add(o)
    owners = sorted({o for os in owned.values() for o in os})
    comments = comment_texts(owners)
    docstrings = docstring_texts(owners)
    texts = (
        docstrings
        + comments
        + cpp_string_texts({line.strip() for _, line, _, _ in docstrings})
        + python_texts()
        + hash_comment_texts()
        + markdown_texts(owners)
    )
    patterns = {name: word(name) for name in dead | set(owned)}
    found = []
    known = code_names()
    for place, line, _, _ in comments:
        for m in IDENTIFIER.finditer(PATH.sub(" ", line)):
            name = m.group(1)
            if (
                name not in known
                and name.rstrip("_") not in known
                and name + "_" not in known
            ):
                found.append((place, name, line.strip()))
    for place, line, text_owners, public in texts:
        words = PATH.sub(" ", line)
        for name in dead:
            if patterns[name].search(words):
                found.append((place, name, line.strip()))
        if public:
            for name, name_owners in owned.items():
                if text_owners & name_owners and patterns[name].search(words):
                    found.append((place, name, line.strip()))
    found += code_findings(dead, names)
    return sorted(set(found))


def check(found, rows):
    """The findings no row excuses, and the rows that do not excuse exactly one finding: a row
    excuses its word on the one line of its place that holds its substring (that is the line
    itself when a line of the place is exactly the substring: a one-word string literal; a code
    row's substring is always the identifier itself)."""
    whole = {(p, n, line) for p, n, line in found}
    unreviewed, matched = [], {i: set() for i in range(len(rows))}
    for place, name, line in found:
        match = [
            i
            for i, (p, n, sub, _) in enumerate(rows)
            if p == place
            and n == name
            and (
                sub == line
                if place.startswith("code:") or (p, n, sub) in whole
                else sub in line
            )
        ]
        for i in match:
            matched[i].add(line)
        if not match:
            unreviewed.append(f"{place}: [{name}] {line}")
    bad_rows = [
        f"{rows[i][0]} | {rows[i][1]} | {rows[i][2]} ({len(lines)} lines)"
        for i, lines in matched.items()
        if len(lines) != 1
    ]
    return unreviewed, bad_rows


def test_removed_names_are_gone_from_every_text():
    found = find_stale()
    assert len(found) > 0  # the English uses: the rule reads the texts it is meant to
    unreviewed, bad_rows = check(found, reviewed())
    assert not unreviewed, "\n".join(unreviewed)
    assert not bad_rows, (
        "reviewed rows that do not excuse exactly one line:\n" + "\n".join(bad_rows)
    )


def test_a_reviewed_row_excuses_its_word_only():
    # A row names its line and its word: another removed name on that line is still found.
    rows = reviewed()
    place, name, sub, _ = rows[0]
    line = next(
        line for p, n, line in find_stale() if p == place and n == name and sub in line
    )
    other = "elig" if name != "elig" else "opts"
    unreviewed, _ = check(
        [(place, name, line), (place, other, line + f" {other}")], rows
    )
    assert unreviewed == [f"{place}: [{other}] {line} {other}"]


def test_the_rule_reads_static_method_texts():
    # A staticmethod's own text (from_csr's), not the descriptor's generic help.
    texts = {(place, line.strip()) for place, line, _, _ in docstring_texts()}
    assert any(
        place == "PredecessorDAG.from_csr" and "parent_offsets" in line
        for place, line in texts
    )
    assert not any("It can be called either on the class" in line for _, line in texts)


# One planted stale line per kind of text, (file, appended text with {name} and the tag {tag}, the
# name it plants): each must be found on its own line. `DagCacheStats` and `Tunnels` come only from the
# two-column table of public names.
PLANTS = [
    (
        "include/netgraph/core/flow_policy.hpp",
        "// {name} {tag}",
        "FlowPlacement",
    ),
    ("src/graph.cpp", "// {name} {tag}", "DagCacheStats"),
    ("src/flow_policy.cpp", "// {name} {tag}", "Tunnels"),
    (
        "src/forwarding.cpp",
        'const char* kPlanted = "{name} {tag}";',
        "FlowPlacement",
    ),
    ("tests/py/helpers.py", "# {name} {tag}", "FlowPlacement"),
    (
        "tests/py/helpers.py",
        'def planted():\n    """{name} {tag}"""',
        "FlowPlacement",
    ),
    ("tests/py/helpers.py", 'PLANTED = "{name} {tag}"', "FlowPlacement"),
    ("dev/stubs/overrides.pyi", "# {name} {tag}", "FlowPlacement"),
    ("dev/run-checks.sh", "# {name} {tag}", "FlowPlacement"),
    (".github/workflows/tests.yml", "# {name} {tag}", "FlowPlacement"),
    ("Makefile", "# {name} {tag}", "FlowPlacement"),
    ("CMakeLists.txt", "# {name} {tag}", "FlowPlacement"),
    ("pyproject.toml", "# {name} {tag}", "FlowPlacement"),
    ("README.md", "{name} {tag}.", "FlowPlacement"),
    ("CONTRIBUTING.md", "{name} {tag}.", "FlowPlacement"),
]


def test_a_planted_name_is_found_in_every_kind_of_text(monkeypatch):
    if not all((ROOT / rel).is_file() for rel, _, _ in PLANTS):
        pytest.skip("needs a git checkout: the sdist leaves out .github/ and dev/")
    planted = {}
    for i, (rel, text, name) in enumerate(PLANTS):
        planted[rel] = (
            planted.get(rel, read(ROOT / rel))
            + "\n"
            + text.format(name=name, tag=f"plant{i}x")
            + "\n"
        )
    monkeypatch.setattr(sys.modules[__name__], "PLANTED", planted)
    unreviewed, _ = check(find_stale(), reviewed())
    for i, (rel, _, name) in enumerate(PLANTS):
        assert any(
            u.startswith(f"{rel}: [{name}] ") and f"plant{i}x" in u for u in unreviewed
        ), PLANTS[i]
