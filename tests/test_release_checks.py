"""Publication guards: failures here must stop release-side mutations."""
import copy
import importlib.util
import json
import re
from pathlib import Path
import tempfile
import tomllib
import unittest

ROOT = Path(__file__).resolve().parents[1]

def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'scripts' / (name + '.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

checks = load('check_release')
packages = load('package_release')

class ReleasePolicyTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        for name in ('.zenodo.json', 'CITATION.cff', 'Cargo.toml', 'LICENSE', 'koth_ff/LICENSE'):
            (self.root / name).parent.mkdir(parents=True, exist_ok=True)
            (self.root / name).write_bytes((ROOT / name).read_bytes())
        version = tomllib.loads((ROOT / 'Cargo.toml').read_text())['workspace']['package']['version']
        self.event = {'action': 'published', 'repository': {'full_name': 'pgarrett-scripps/koth', 'private': False},
                      'release': {'tag_name': 'v' + version, 'draft': False}}

    def test_valid_metadata_and_release(self):
        checks.check_metadata(self.root)
        checks.check_event(self.root, self.event)

    def test_rejects_wrong_tag_private_repo_draft_and_nonrelease(self):
        for section, key, value in [('release', 'tag_name', 'v999.0.0'),
                                    ('release', 'draft', True),
                                    ('repository', 'private', True),
                                    ('repository', 'full_name', 'someone/fork')]:
            event = copy.deepcopy(self.event)
            event[section][key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                checks.check_event(self.root, event)
        self.event['action'] = 'edited'
        with self.assertRaises(ValueError):
            checks.check_event(self.root, self.event)

    def test_rejects_version_date_and_nested_funding(self):
        path = self.root / '.zenodo.json'
        original = json.loads(path.read_text())
        for mutation in [{'version': '0.0.1'}, {'publication_date': '2020-01-01'},
                         {'grants': []}, {'extra': {'funding': 'not allowed'}}]:
            path.write_text(json.dumps(original | mutation))
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                checks.check_metadata(self.root)
        path.write_text(json.dumps(original))
        cff = self.root / 'CITATION.cff'
        text = cff.read_text()
        for mutated in [re.sub(r'(?m)^version: .*$', 'version: "0.0.1"', text),
                        re.sub(r'(?m)^date-released: .*\n', '', text),
                        text + '\nfunding: not allowed\n']:
            cff.write_text(mutated)
            with self.subTest(mutated=mutated[-40:]), self.assertRaises(ValueError):
                checks.check_metadata(self.root)

    def test_rejects_creator_drift(self):
        path = self.root / '.zenodo.json'
        metadata = json.loads(path.read_text())
        metadata['creators'][0]['name'] = 'Someone, Else'
        path.write_text(json.dumps(metadata))
        with self.assertRaises(ValueError):
            checks.check_metadata(self.root)

    def test_archives_include_binaries_metadata_and_documentation(self):
        import tarfile
        import zipfile
        from unittest.mock import patch
        for name in ('README.md', 'CHANGELOG.md', 'CONTRIBUTING.md', 'RELEASE.md',
                     'example_config.toml', 'example_config_align.toml', 'docs/README.md'):
            path = self.root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text('test content')
        output = self.root / 'dist'
        with patch.object(packages, 'ROOT', self.root):
            for target in packages.TARGETS:
                extension = '.exe' if 'windows' in target else ''
                for binary in ('koth_ff', 'koth_align'):
                    path = self.root / 'target' / target / 'release' / (binary + extension)
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.write_bytes(b'test executable')
                    path.chmod(0o755)
                packages.package(target, output)
                archive = output / packages.archive_name(target)
                if extension:
                    with zipfile.ZipFile(archive) as stream:
                        names = stream.namelist()
                else:
                    with tarfile.open(archive) as stream:
                        names = stream.getnames()
                        binary = next(item for item in stream.getmembers() if item.name.endswith('/koth_ff'))
                        self.assertTrue(binary.mode & 0o111)
                for name in ('koth_ff' + extension, 'koth_align' + extension, 'LICENSE',
                             'CITATION.cff', '.zenodo.json', 'docs/README.md'):
                    self.assertTrue(any(n.endswith('/' + name) for n in names), name)
            packages.verify_dist(output)

    def test_incomplete_or_corrupt_downloads_cannot_publish(self):
        import hashlib
        for target in packages.TARGETS:
            name = packages.archive_name(target)
            (self.root / name).write_bytes(b'test archive')
            (self.root / (name + '.sha256')).write_text(hashlib.sha256(b'test archive').hexdigest() + '  ' + name + '\n')
        import shutil
        shutil.rmtree(self.root / 'koth_ff')
        for name in ('.zenodo.json', 'CITATION.cff', 'Cargo.toml', 'LICENSE'):
            (self.root / name).unlink()
        packages.verify_dist(self.root)
        name = packages.archive_name(packages.TARGETS[0])
        (self.root / name).write_bytes(b'corrupt archive')
        with self.assertRaises(ValueError):
            packages.verify_dist(self.root)
        (self.root / name).unlink()
        with self.assertRaises(ValueError):
            packages.verify_dist(self.root)

if __name__ == '__main__':
    unittest.main()
