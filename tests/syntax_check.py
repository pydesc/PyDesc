# Copyright 2017 Pawel Daniluk
#
# This file is part of PyDesc.
#
# PyDesc is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# PyDesc is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with PyDesc.  If not, see <http://www.gnu.org/licenses/>.

"""
Syntax checks and some helper functions.

Helper module for designing unittests.
"""

import tabnanny
import pep8
import unittest
import random
import string
import operator
import warnings
import sys
import contextlib
import re
import functools
import inspect
import os

from collections import defaultdict

import types

module = None

has_test = defaultdict(list)
no_test = defaultdict(list)


try:
    import bitten.util.testrunner
    from bitten.util import xmlio
    from unittest import TextTestRunner

    def run(self, test):
        """
        Monkey patch improving the way bitten extracts test and fixture names from result of TestCase.__str__. It fixes
        the case when test name contains a dot as in failure in module load (unittest.loader.LoadTestsFailure) (not a test
        really, but should be reported any way.

        Works with Bitten 0.7dev-r1026.
        """
        result = TextTestRunner.run(self, test)
        if not self.xml_stream:
            return result

        root = xmlio.Element('unittest-results')
        for testcase, filename, timetaken, stdout, stderr in result.tests:
            status = 'success'
            tb = None

            if testcase in [e[0] for e in result.errors]:
                status = 'error'
                tb = [e[1] for e in result.errors if e[0] is testcase][0]
            elif testcase in [f[0] for f in result.failures]:
                status = 'failure'
                tb = [f[1] for f in result.failures if f[0] is testcase][0]

            name = str(testcase)
            fixture = None
            description = testcase.shortDescription() or ''
            if description.startswith('doctest of '):
                name = 'doctest'
                fixture = description[11:]
                description = None
            else:
                match = re.match(r'(\S+)\s+\(([\w.]+)\)', name)  # this regexp is changed from '(\w+)\s+\(([\w.]+)\)'
                if match:
                    name = match.group(1)
                    fixture = match.group(2)

            test_elem = xmlio.Element('test', file=filename, name=name,
                                      fixture=fixture, status=status,
                                      duration=timetaken)
            if description:
                test_elem.append(xmlio.Element('description')[description])
            if stdout:
                test_elem.append(xmlio.Element('stdout')[stdout])
            if stderr:
                test_elem.append(xmlio.Element('stdout')[stderr])
            if tb:
                test_elem.append(xmlio.Element('traceback')[tb])
            root.append(test_elem)

        root.write(self.xml_stream, newlines=True)
        return result

    bitten.util.testrunner.XMLTestRunner.run = run


except ImportError:
    pass


@contextlib.contextmanager
def capture():
    """ Context manager capturing stdout and stderr. Returns a tuple with two streams (out and err respectively). """
    from cStringIO import StringIO
    oldout, olderr = sys.stdout, sys.stderr
    try:
        out = [StringIO(), StringIO()]
        sys.stdout, sys.stderr = out
        yield out
    finally:
        sys.stdout, sys.stderr = oldout, olderr


def randomcoord(n=3):
    """ Return a vector of n random reals from [-5,5]. """
    return [random.uniform(-5, 5) for dummy in range(n)]


def randomword(length):
    """ Return a random lowercase string of given lehgth. """
    return ''.join(random.choice(string.lowercase) for i in range(length))


def dot(v1, v2):
    """ Return a dot product of vectors. """
    return sum(map(operator.mul, v1, v2))


class expando(object):

    """ Empty class. """

    pass


def test_deprec(testcase, function, name):
    """Assert whether function raises a DeprecationWarning."""
    msg = "%s should be deprecated." % name
    with warnings.catch_warnings(record=True) as w:
        function()
        testcase.assertEqual(len(w), 1, msg)
        testcase.assertTrue(issubclass(w[-1].category, DeprecationWarning), msg)


def warning_message(testcase, wlist):
    """Fail a testcase displaying warnings from wlist."""
    if len(wlist) > 0:
        msg = '\n'.join([str(w.category) + ': ' + str(w) for w in wlist[:10]])
        testcase.fail("Warnings occured (%d total, %d shown):\n%s" % (len(wlist), len(wlist[:10]), msg))


def testing(item, mod=None):
    """ A decorator to register items tested by the method. Returns unmodified function. """
    if mod is None:
        mod = module
    if mod is None:
        raise Exception("Module not given.")

    if isinstance(item, types.MethodType):
        item = item.im_func

    def wrapper(f):
        """ Register test. """
        test(item, mod)
        return f

    return wrapper


def _modify_testcase_name(name, *args):
    """ Append args to a unittest name (returned by TestCase.__str__) in a way which doesn't disturb bitten unittest XML reporting.

        args are converted to strings (classes are converted to a fully qualified names) and added to test name.
        All characters not matching '\\w' are converted to underscores.

        If name cannot be split into test name and fixture, args are appended to name in parenthesesized comma separated list.
    """
    cargs = []

    for a in args:
        if isinstance(a, type):
            ca = '%s.%s' % (a.__module__, a.__name__)
        else:
            ca = str(a)
        cargs.append(re.sub(r'\W', '_', ca))

    match = re.match(r'(\w+)\s+\(([\w.]+)\)', name)
    if match:
        name = '%s_%s' % (match.group(1), '_'.join(cargs))
        fixture = match.group(2)

        return '%s (%s)' % (name, fixture)

    return name[:-1] + '(%s))' % (', '.join(cargs))


def test_name_append(*args):
    """ A decorator for unittest.TestCase. Mofifies __str__ method to append args to its return."""
    def wrapper(cls):
        def str_wrapper(old_str):
            def new_str(self):
                return _modify_testcase_name(old_str(self), *args)
            return new_str

        cls.__str__ = str_wrapper(cls.__str__)
        return cls
    return wrapper


def _regtest(d, item, mod=None):
    """ Register item in dict d. Return none. """
    if mod is None:
        mod = module
    if mod is None:
        raise Exception("Module not given.")

    if isinstance(item, types.MethodType):
        item = item.im_func
    d[mod].append(item)

test = functools.partial(_regtest, has_test)
notest = functools.partial(_regtest, no_test)


def module_syntax(mod=None):  # pylint: disable=R0912
    """
    Create and return a TestCase for syntax checking in module.


    Also checks for environment variables FAST and SKIP_SLOW and
    sets respective attributes in tested module. This is an ugly hack
    to avoid modifying all existing modules.
    """
    if mod is None:
        mod = module
    if mod is None:
        raise Exception("Module not given.")

    filename = mod.__file__
    if filename.endswith('.pyc'):
        filename = filename.rstrip('c')

    skiplist = ['__str__', '__repr__']

    calling_mod = inspect.getmodule(inspect.currentframe().f_back)

    calling_mod.fast = from_env('fast')
    calling_mod.skip_slow = from_env('skip_slow')

    calling_mod = calling_mod.__name__


    def traverse(test_func, el, path=None, skip=None):
        """ Traverse el object and run test_func function for each suitable item.

        Returns list of elements for which test_func fails.

        """
        if path is None:
            path = []

        if skip is None:
            skip = skiplist

        res = []
        private = lambda name: name.startswith('_') and not (name.startswith('__') and name.endswith('__'))

        for name in el.__dict__:
            if not private(name):
                item = el.__dict__[name]

                if item and not name in skip and getattr(item, '__module__', mod.__name__) == mod.__name__:
                    p1 = path + [name]
                    if not test_func(item):
                        res.append('.'.join(p1))
                    if isinstance(item, types.ClassType) or isinstance(item, types.TypeType) and not item in no_test[mod]:
                        res.extend(traverse(test_func, item, p1, skip))
        return res

    class TestSyntax(unittest.TestCase):

        """ Various syntax tests."""

        def __str__(self):
            return unittest.TestCase.__str__(self).replace(__name__, calling_mod)

        def test_tabnanny(self):
            """ Tabnanny check """
            self.longMessage = False
            try:
                tabnanny.check(filename)
            except tabnanny.NannyNag:
                self.fail("File %s fails tabnanny.check." % filename)

        def test_pep8_conformance(self):
            """ PEP8 conformance """
            pep8style = pep8.StyleGuide(quiet=True, reporter=pep8.StandardReport, ignore=['E501'])
            with capture() as (out, dummy_err):
                result = pep8style.check_files([filename])
                msg = '\n'.join(out.getvalue().split('\n')[:10])
                self.assertEqual(result.total_errors, 0, "Found %d code style errors (and warnings) (first 10 shown):\n%s" % (result.total_errors, msg))

        def test_pylint(self):
            """ PyLint conformance """

            class WritableObject(object):
                "dummy output stream for pylint"
                def __init__(self):
                    self.content = []

                def write(self, st):
                    "dummy write"
                    self.content.append(st)

                def read(self):
                    "dummy read"
                    wholetext = ''.join(self.content)
                    return re.findall(r"^[^I]:[\s\d]+,[\s\d]+:.*\n", wholetext, re.MULTILINE)

            from pylint import lint
            from pylint.reporters.text import TextReporter
            ARGS = ["-d I,C0301,R0901,R0902,R0903,R0904,R0913,R0915,W0141,W0142,W0232,W0613"]
            pylint_output = WritableObject()
            with capture() as (dummy_out, dummy_err):
                lint.Run([filename] + ARGS, reporter=TextReporter(pylint_output), exit=False)
            msg = pylint_output.read()
            self.assertEqual(len(msg), 0, "Pylint detected %d errors (and warnings) (first 10 shown):\n%s" % (len(msg), ''.join(msg[:10])))

        def test_docs(self):
            """ Docstrings (internal procedure) """

            def test_func(item):
                """ Check if item as a nonempty docstring. """
                return not item.__doc__ in [None, ""] or isinstance(item, types.ModuleType)

            res = traverse(test_func, mod, [mod.__name__])
            self.assertEqual(res, [], "Following elements lack a docstring (%d total): %s" % (len(res), str(sorted(res))))

        def test_has_test(self):
            """ Unittests for all elements """
            testables = [types.ClassType, types.MethodType, types.FunctionType, types.GetSetDescriptorType, types.MemberDescriptorType, property, types.TypeType]
            testables = [types.ClassType, types.MethodType, types.FunctionType, property, staticmethod, classmethod, types.TypeType]

            def test_func(item):
                """ Check if item has a registered unittest if applicable. """
                for t in testables:
                    if isinstance(item, t):
                        if isinstance(item, staticmethod) or isinstance(item, classmethod):
                            item = item.__func__

                        if item in has_test[mod] or item in no_test[mod]:
                            return True

                        if getattr(item, '__isabstractmethod__', False):
                            return True

                        if hasattr(item, 'im_class'):
                            return item.im_class in no_test[mod]

                        return False

                return True

            res = traverse(test_func, mod, [mod.__name__])
            self.assertEqual(res, [], "Following elements lack a test (%d total): %s" % (len(res), str(sorted(res))))

    return TestSyntax


def variants_tested(tests_performed):
    """ Return a testcase, which asserts that all keys in tests_performed dict are True. """

    calling_mod = inspect.getmodule(inspect.currentframe().f_back).__name__

    class VariantsTested(unittest.TestCase):
        """ Testcase checking if all tests were performed (in case it depends on random or outside data. """

        def __str__(self):
            return unittest.TestCase.__str__(self).replace(__name__, calling_mod)

        def test_variants(self):
            """ All variants tested """
            notperformed = [x[0] for x in tests_performed.items() if not x[1]]

            self.assertEqual(notperformed, [], "Did not test the following: %s" % str(notperformed))

    return VariantsTested


def rip_argv(opt):
    """ Check if opt exists in sys.argv and remove it. Returns True if opt was there. """
    if opt in sys.argv:
        sys.argv.remove(opt)
        return True
    return False


def from_env(opt):
    """ Check if the exists an environment variable named opt.upper(), and return it. """

    try:
        return os.environ[opt.upper()]
    except KeyError:
        return None
