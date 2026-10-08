"""Fast packaging/entry-point checks; this does not replay either numerical proof."""
import copy
import json
import os
from pathlib import Path
import runpy
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]

class PackagingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix='kg-packaging-')
        cls.root = Path(cls.temporary.name)/'relocated'
        cls.root.mkdir()
        for line in (ROOT/'MANIFEST.sha256').read_text().splitlines():
            _, name = line.split(None, 1)
            dest = cls.root/name
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(ROOT/name, dest)
        shutil.copyfile(ROOT/'MANIFEST.sha256', cls.root/'MANIFEST.sha256')

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def run_wrapper(self, *arguments, optimized=False):
        env = os.environ.copy()
        env['PYTHONDONTWRITEBYTECODE'] = '1'
        command = [sys.executable]
        if optimized:
            command.append('-O')
        command += [str(self.root/'verify.py'), *arguments]
        return subprocess.run(command, cwd=self.temporary.name, env=env,
                              capture_output=True, text=True)

    def assert_rejected(self, result, message):
        self.assertNotEqual(result.returncode, 0, result.stdout+result.stderr)
        self.assertIn(message, result.stderr)
        self.assertFalse(list(self.root.rglob('ACCEPTED.json')))

    def test_relocated_data_only(self):
        result = self.run_wrapper('--data-only')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('EXACT_ROTATION_TAIL_CHECKS=70', result.stdout)
        self.assertIn('NO_ANALYTIC_NUMERICAL_REPLAY=1', result.stdout)
        self.assertFalse(list(self.root.rglob('ACCEPTED.json')))

    def test_tampered_input(self):
        target = self.root/'data/lower_rows.json'
        original = target.read_bytes()
        try:
            target.write_bytes(original+b'\n')
            self.assert_rejected(self.run_wrapper('--data-only'), 'File hash mismatch')
        finally:
            target.write_bytes(original)

    def test_missing_source(self):
        target = self.root/'upper/source/certification/certify_upper.py'
        original = target.read_bytes()
        try:
            target.unlink()
            self.assert_rejected(self.run_wrapper('--data-only'), 'Missing/invalid file')
        finally:
            target.write_bytes(original)

    def test_unlisted_source(self):
        target = self.root/'upper/source/unlisted.py'
        try:
            target.write_text('raise RuntimeError("unlisted")\n')
            self.assert_rejected(self.run_wrapper('--data-only'), 'Source/input inventory differs')
        finally:
            target.unlink()

    def test_optimized_python(self):
        self.assert_rejected(self.run_wrapper('--data-only', optimized=True),
                             'Run ordinary Python without -O')

    def test_output_inside_inputs(self):
        self.assert_rejected(self.run_wrapper('--lower','--output',str(self.root/'data/runs')),
                             'Output must be outside source/input directories')

    def test_all_70_corrupted_tail_weights(self):
        checker = runpy.run_path(str(self.root/'lower/check-tail.py'))['check']
        data = json.loads((self.root/'data/lower_rows.json').read_text())
        for row in range(35):
            for rotation in range(2):
                invalid = copy.deepcopy(data)
                invalid['rows'][row]['angular'][rotation]['weights'][0] = '0'
                with self.subTest(row=row, rotation=rotation):
                    with self.assertRaisesRegex(ValueError, 'Dominance failed'):
                        checker(invalid)

if __name__ == '__main__':
    unittest.main(verbosity=2)
