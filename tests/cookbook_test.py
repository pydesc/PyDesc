import traceback
from pathlib import Path

import pytest

TEST_PATH = Path(__file__).parent
DOC_PATH = TEST_PATH / ".." / "doc"
COOKBOOK_PATH = DOC_PATH / "cookbook.md"


def read_cookbook_lines():
    with open(COOKBOOK_PATH) as cb_handler:
        cb_lines = [line.rstrip() for line in cb_handler.readlines()]
    cb_lines = [line for line in cb_lines if bool(line)]
    return cb_lines


def make_section_tree(lines):
    root = {"text": []}
    while lines:
        branch = {"text": []}
        root[lines[0]] = branch
        lvl = 1
        lines = fill_section_branch(branch, lvl, lines[1:])
    return root


def fill_section_branch(branch_dct, lvl, lines):
    n = 0
    while lines:
        line = lines[n]
        n += 1
        if not line.startswith("#"):
            branch_dct["text"].append(line)
            continue
        clvl = line.count("#")
        if clvl <= lvl:
            return lines[n - 1 :]
        sub_branch_dct = {"text": []}
        branch_dct[line] = sub_branch_dct
        lines = fill_section_branch(sub_branch_dct, clvl, lines[n:])
        n = 0
    return []


def format_key(key):
    key = key.replace("#", "")
    key = key.strip()
    key = key.split()
    key = [item for item in key if bool(item)]
    key = "_".join(key)
    return key


def iter_tree(key, tree):
    yield format_key(key), tree["text"]
    for sub_key in tree:
        if sub_key == "text":
            continue
        nested_key = "_".join((key, sub_key))
        yield from iter_tree(nested_key, tree[sub_key])


def extract_example(lines):
    examples = []
    marks = []
    add = False
    for line in lines:
        if line.startswith("```python"):
            examples.append("")
            add = True
            if "web" in line:
                marks.append(True)
            else:
                marks.append(False)
            continue
        elif line.startswith("```"):
            add = False
        if add:
            examples[-1] += line
            examples[-1] += "\n"
    return zip(examples, marks)


def get_examples():
    lines = read_cookbook_lines()
    tree = make_section_tree(lines)
    for section, lines in iter_tree("Cookbook", tree):
        examples = extract_example(lines)
        examples = [example for example in examples if bool(example[0])]
        for n, (example, mark) in enumerate(examples, 1):
            test_name = f"{section}_{n}"
            yield test_name, ExecutableTest(example, mark)


class CookbookTestFail(Exception):
    pass


class ExecutableTest:
    def __init__(self, code, web_mark=False):
        self.code = code
        self.web = web_mark

    def is_web(self):
        return self.web

    def execute(self):
        try:
            exec(self.code)
        except Exception:
            e_msg = traceback.format_exc()
            msg = f"Failed test with code:\n{self.code}\n" f"Error:\n{e_msg}"
            raise CookbookTestFail(msg)


@pytest.mark.parametrize("name,test_code", get_examples())
def test_cookbook(request, name, test_code):
    markings = request.config.getoption("keyword")
    markings += request.config.getoption("markexpr")
    if "not web" in markings:
        if test_code.is_web():
            pytest.xfail("Skipping web test.")
    test_code.execute()
