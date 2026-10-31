"""The documented examples print what the documents show (CONTRIBUTING, "Behaviour
Documentation"): every Python example of docs/semantics.md, run in order in one namespace as the
document says, against the ``text`` block that follows it (an example without one prints
nothing), and the Python blocks of README's quick start against the output written in the comment
of each print."""

from __future__ import annotations

import contextlib
import io
import re
import traceback
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
FENCE = re.compile(r"^```(\w*)\n(.*?)^```$", re.S | re.M)


def blocks(text: str) -> list[tuple[int, str, str, int, int]]:
    """(line, language, body, start, end) of every fenced block, in document order."""
    found = list(FENCE.finditer(text))
    # A different Markdown closing fence may end a block before the parser's closing line.
    assert not any(
        re.search(r"^ {0,3}`{3,}[ \t]*$", m.group(2), re.M) for m in found
    ), "unsupported closing fence inside a parsed block"
    # A fence marker outside the parsed blocks uses an unsupported form or has no closing
    # fence. Other markers inside a block are contents an example may print literally.
    assert all(
        any(m.start() <= fence.start() < m.end() for m in found)
        for fence in re.finditer(r"^[ \t]*(?:`{3,}|~{3,}).*$", text, re.M)
    ), "only unindented triple-backtick fences with a bare language are supported"
    return [
        (text.count("\n", 0, m.start()) + 1, m.group(1), m.group(2), m.start(), m.end())
        for m in found
    ]


def run(code: str, name: str, namespace: dict) -> str:
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        exec(compile(code, name, "exec"), namespace)
    return out.getvalue()


def check_examples(text: str, name: str) -> tuple[int, int]:
    """Runs every Python example of a document in order in one namespace and compares what each
    prints with the ``text`` block directly after it (an example without one prints nothing).
    Returns (the examples, the documented outputs compared); a structural defect or a
    mismatch fails an assertion."""
    found = blocks(text)
    # Only Python examples and their outputs, each output directly after an example: a fence of
    # another language (a typo such as ```py) would otherwise leave an example unchecked.
    assert {lang for _, lang, *_ in found} <= {"python", "text"}
    for i, (line, lang, _, start, _) in enumerate(found):
        if lang == "text":
            prev = found[i - 1] if i else None
            assert prev is not None and prev[1] == "python", (
                f"line {line}: output of no example"
            )
            assert text[prev[4] : start].strip() == "", (
                f"line {line}: text between example and output"
            )
    examples = [b for b in found if b[1] == "python"]
    assert examples

    namespace: dict = {}
    mismatches = []
    outputs = 0  # documented outputs compared with what their example printed
    for i, (line, lang, code, _, end) in enumerate(found):
        if lang != "python":
            continue
        nxt = found[i + 1] if i + 1 < len(found) else None
        expected = (
            nxt[2]
            if nxt and nxt[1] == "text" and text[end : nxt[3]].strip() == ""
            else ""
        )
        got = run(code, f"{name}:{line}", namespace)
        outputs += expected != ""
        if got != expected:
            mismatches.append(
                f"line {line}:\n--- documented\n{expected}--- printed\n{got}"
            )
    assert not mismatches, "\n".join(mismatches)
    assert outputs == sum(lang == "text" for _, lang, *_ in found) > 0
    return len(examples), outputs


def check_quick_start(text: str) -> int:
    """Runs the Python blocks of README's quick start (its heading up to the next ``## ``
    heading; its other blocks are ``text`` pictures) in order in one namespace, and compares
    every line each block prints with the comment of its print. Returns the prints compared."""
    heading = text.index("\n## Quick Start\n")
    end = text.find("\n## ", heading + 1)
    section = [
        b for b in blocks(text) if heading < b[3] < (len(text) if end < 0 else end)
    ]
    langs = {lang for _, lang, *_ in section}
    assert "python" in langs and langs <= {"python", "text"}, "the quick start's blocks"
    namespace: dict = {}
    compared = 0
    for line, lang, code, *_ in section:
        if lang != "python":
            continue
        prints = [ln for ln in code.splitlines() if ln.startswith("print(")]
        expected = [re.fullmatch(r"print\(.*\)\s+# (.*)", ln) for ln in prints]
        assert prints and all(expected), (
            f"README.md:{line}: a print without its output comment"
        )
        got = run(code, f"README.md:{line}", namespace).splitlines()
        assert got == [m.group(1) for m in expected if m]
        compared += len(prints)
    return compared


def test_semantics_examples_print_their_documented_output():
    text = (ROOT / "docs" / "semantics.md").read_text(encoding="utf-8")
    check_examples(text, "docs/semantics.md")


def test_readme_quick_start_prints_its_commented_output():
    text = (ROOT / "README.md").read_text(encoding="utf-8")
    assert check_quick_start(text) > 0


GOOD = "```python\nx = 2\nprint(x)\n```\n\n```text\n2\n```\n\nProse.\n\n```python\ny = x\n```\n"


def test_the_checks_reject_a_broken_document():
    # The checkers pass a well-formed document and fail each defect they exist for, each at its
    # own clause (the failing source line), so a deleted clause is noticed even when a later
    # clause would also fail.
    assert check_examples(GOOD, "good") == (2, 1)
    literal = "```python\nprint('~~~python\\n  ```python')\n```\n\n```text\n~~~python\n  ```python\n```\n"
    assert check_examples(literal, "literal markers") == (1, 1)
    count = 'assert outputs == sum(lang == "text" for _, lang, *_ in found) > 0'
    mismatch = 'assert not mismatches, "\\n".join(mismatches)'
    broken = {
        "a changed output": (GOOD.replace("```text\n2\n", "```text\n3\n"), mismatch),
        "an unseen indented example": (
            GOOD + "\n  ```python\n  print(999)\n  ```\n",
            "assert all(",
        ),
        "an unseen tilde example": (
            GOOD + "\n~~~python\nprint(999)\n~~~\n",
            "assert all(",
        ),
        "an unseen longer fence": (
            GOOD + "\n````python\nprint(999)\n````\n",
            "assert all(",
        ),
        "an unseen fence with extra information": (
            GOOD + "\n```python title\nprint(999)\n```\n",
            "assert all(",
        ),
        "an unseen unclosed fence": (
            GOOD + "\n```python\nprint(999)\n",
            "assert all(",
        ),
        "another fence language": (
            GOOD.replace("```python\ny", "```py\ny"),
            'assert {lang for _, lang, *_ in found} <= {"python", "text"}',
        ),
        "an output of no example": (
            "```text\n2\n```\n\n" + GOOD,
            'assert prev is not None and prev[1] == "python", (',
        ),
        "prose between an example and its output": (
            GOOD + "\n```python\nz = 1\n```\nProse.\n```text\n1\n```\n",
            'assert text[prev[4] : start].strip() == "", (',
        ),
        "an example that prints without an output": (
            GOOD.replace("y = x", "print(x)"),
            mismatch,
        ),
        "no example": ("Prose only.\n", "assert examples"),
        "an empty output block": ("```python\nx = 1\n```\n\n```text\n```\n", count),
        "no documented output": ("```python\nx = 1\n```\n", count),
    }
    for close in ("  ```", "````", "``` "):
        printed = "2\n" + close + "\n~~~python\nprint(999)\n~~~"
        broken[f"an unseen example after {close!r}"] = (
            GOOD.replace("print(x)", f"print({printed!r})").replace(
                "```text\n2\n```", "```text\n" + printed + "\n```"
            ),
            "assert not any(",
        )
    for what, (text, clause) in broken.items():
        with pytest.raises(AssertionError) as info:
            check_examples(text, what)
        assert traceback.extract_tb(info.tb)[-1].line == clause, what
    second = "```python\nx = 3\nprint(x + 1)  # 4\n```\n"
    quick = (
        "# Title\n\n## Quick Start\n\n```text\nA-B\n```\n\n```python\nprint(1 + 1)  # 2\n```\n"
        f"\n### More\n\n{second}\n## Next\n\n```bash\nmake\n```\n"
    )
    assert check_quick_start(quick) == 2
    commented = "assert prints and all(expected), ("
    wrong = "assert got == [m.group(1) for m in expected if m]"
    langs = 'assert "python" in langs and langs <= {"python", "text"}, "the quick start\'s blocks"'
    for text, clause in (
        (quick.replace("# 2", "# 3"), wrong),  # a wrong output comment
        (quick.replace("# 4", "# 5"), wrong),  # a wrong output comment in a later block
        (quick.replace("  # 2", ""), commented),  # a print without its comment
        (quick.replace("print(1 + 1)  # 2", "x = 1"), commented),  # no print at all
        # The quick start has no Python block; a later block is in another language.
        (quick.replace("```python", "```bash"), langs),
        (quick.replace("```python\nx", "```py\nx"), langs),
    ):
        with pytest.raises(AssertionError) as info:
            check_quick_start(text)
        assert traceback.extract_tb(info.tb)[-1].line == clause, text
